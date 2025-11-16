#!/bin/bash
set -euo pipefail

# ========================================
# Generate Enhanced GitHub Pages
# ========================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

cd "${ROOT}"

DB_PATH="catalog/db/ocp_results.sqlite"
PAGES_DIR="docs/pages"
RESULTS_DIR="results"

mkdir -p "${PAGES_DIR}/runs"
mkdir -p "${PAGES_DIR}/assets"

TIMESTAMP=$(date '+%Y-%m-%d %H:%M:%S')

echo "========================================="
echo "Generating Enhanced GitHub Pages"
echo "========================================="
echo "Database : ${DB_PATH}"
echo "Output   : ${PAGES_DIR}"
echo "Timestamp: ${TIMESTAMP}"
echo ""

# ========================================
# Main Index Page - Results Table
# ========================================
echo "[1/3] Generating main index page..."

cat > "${PAGES_DIR}/index.md" <<EOF
---
title: OCP Results Dashboard
layout: default
---

# 🧬 OCP Results Dashboard

**Last Updated:** ${TIMESTAMP}

## 📊 Latest Results

EOF

# Get all results with details
cat >> "${PAGES_DIR}/index.md" <<'HEADER'
| ZINC ID | SMILES | 3D | Vina Score (kcal/mol) | GROMACS Status | MD Status | Details |
|---------|--------|----|-----------------------|----------------|-----------|---------|
HEADER

sqlite3 "${DB_PATH}" -separator '|' <<'SQL' | while IFS='|' read -r zinc_id smiles has_3d affinity gmx_status md_status run_id; do
SELECT 
    l.zinc_id,
    COALESCE(substr(l.smiles, 1, 30), 'N/A'),
    CASE WHEN l.has_3d = 1 THEN '✓' ELSE '✗' END,
    COALESCE(printf('%.4f', vr.affinity), 'N/A'),
    COALESCE(vr.gromacs_prep_status, 'pending'),
    COALESCE(vr.md_status, 'pending'),
    r.id
FROM vina_results vr
JOIN ligands l ON vr.ligand_id = l.id
JOIN runs r ON vr.run_id = r.id
WHERE vr.mode_rank = 1
ORDER BY r.started_at DESC, vr.affinity ASC
LIMIT 50;
SQL
    if [[ "${smiles}" == *"..."* ]] || [[ ${#smiles} -gt 30 ]]; then
        smiles="${smiles:0:27}..."
    fi
    echo "| ${zinc_id} | \`${smiles}\` | ${has_3d} | ${affinity} | ${gmx_status} | ${md_status} | [View](runs/run_${run_id}_${zinc_id}.html) |" >> "${PAGES_DIR}/index.md"
done

# Replace timestamp placeholder
sed -i "s|\${TIMESTAMP}|${TIMESTAMP}|g" "${PAGES_DIR}/index.md"

cat >> "${PAGES_DIR}/index.md" <<'EOF'

## 📈 Database Statistics

| Metric | Value |
|--------|-------|
EOF

sqlite3 "${DB_PATH}" -separator '|' <<'SQL' | while IFS='|' read -r metric value; do
SELECT 
    'Total Runs',
    COUNT(*)
FROM runs
UNION ALL
SELECT 
    'Total Ligands',
    COUNT(*)
FROM ligands
UNION ALL
SELECT 
    'Ligands with 3D',
    COUNT(*)
FROM ligands WHERE has_3d = 1
UNION ALL
SELECT 
    'Vina Results',
    COUNT(*)
FROM vina_results WHERE mode_rank = 1
UNION ALL
SELECT 
    'GROMACS Completed',
    COUNT(*)
FROM vina_results WHERE gromacs_prep_status = 'completed'
UNION ALL
SELECT 
    'MD Completed',
    COUNT(*)
FROM vina_results WHERE md_status = 'completed';
SQL
    echo "| ${metric} | ${value} |" >> "${PAGES_DIR}/index.md"
done

cat >> "${PAGES_DIR}/index.md" <<'EOF'

---

## 🔬 Recent Runs

EOF

sqlite3 "${DB_PATH}" <<'SQL' >> "${PAGES_DIR}/index.md"
SELECT 
    '- **RUN ' || r.id || '** (' || r.run_uuid || ')  ' || char(10) ||
    '  Started: ' || r.started_at || '  ' || char(10) ||
    '  Target: ' || t.name || ' (' || t.code || ')  ' || char(10) ||
    '  Results: ' || COUNT(vr.ligand_id) || ' ligands'
FROM runs r
JOIN targets t ON r.target_id = t.id
LEFT JOIN vina_results vr ON r.id = vr.run_id AND vr.mode_rank = 1
GROUP BY r.id
ORDER BY r.started_at DESC
LIMIT 10;
SQL

echo "[✓] Main index page generated"

# ========================================
# Individual Result Detail Pages
# ========================================
echo "[2/3] Generating individual result pages..."

DETAIL_COUNT=0

# Get all unique RUN_ID + LIGAND_ID combinations
sqlite3 "${DB_PATH}" "
SELECT DISTINCT r.id, l.zinc_id
FROM vina_results vr
JOIN ligands l ON vr.ligand_id = l.id
JOIN runs r ON vr.run_id = r.id
WHERE vr.mode_rank = 1
ORDER BY r.id DESC, vr.affinity ASC;
" | while IFS='|' read -r RUN_ID ZINC_ID; do
    
    DETAIL_FILE="${PAGES_DIR}/runs/run_${RUN_ID}_${ZINC_ID}.md"
    
    # Get detailed info
    DETAILS=$(sqlite3 "${DB_PATH}" <<SQL
SELECT 
    l.zinc_id,
    l.smiles,
    vr.affinity,
    vr.pose_file,
    vr.gromacs_prep_status,
    vr.gromacs_prep_dir,
    vr.md_status,
    vr.md_output_dir,
    r.run_uuid,
    r.started_at,
    t.code,
    t.name,
    t.pdb_id,
    t.receptor_pdbqt
FROM vina_results vr
JOIN ligands l ON vr.ligand_id = l.id
JOIN runs r ON vr.run_id = r.id
JOIN targets t ON r.target_id = t.id
WHERE r.id = ${RUN_ID}
  AND l.zinc_id = '${ZINC_ID}'
  AND vr.mode_rank = 1;
SQL
)
    
    # Parse details
    SMILES=$(echo "${DETAILS}" | cut -d'|' -f2)
    AFFINITY=$(echo "${DETAILS}" | cut -d'|' -f3)
    POSE_FILE=$(echo "${DETAILS}" | cut -d'|' -f4)
    GMX_STATUS=$(echo "${DETAILS}" | cut -d'|' -f5)
    GMX_DIR=$(echo "${DETAILS}" | cut -d'|' -f6)
    MD_STATUS=$(echo "${DETAILS}" | cut -d'|' -f7)
    MD_DIR=$(echo "${DETAILS}" | cut -d'|' -f8)
    RUN_UUID=$(echo "${DETAILS}" | cut -d'|' -f9)
    STARTED_AT=$(echo "${DETAILS}" | cut -d'|' -f10)
    TARGET_CODE=$(echo "${DETAILS}" | cut -d'|' -f11)
    TARGET_NAME=$(echo "${DETAILS}" | cut -d'|' -f12)
    PDB_ID=$(echo "${DETAILS}" | cut -d'|' -f13)
    RECEPTOR=$(echo "${DETAILS}" | cut -d'|' -f14)
    
    # Create detail page
    cat > "${DETAIL_FILE}" <<DETAIL
---
title: ${ZINC_ID} - RUN ${RUN_ID}
layout: default
---

# 🔬 Detailed Results: ${ZINC_ID}

**Generated:** ${TIMESTAMP}  
**RUN ID:** ${RUN_ID}  
**RUN UUID:** ${RUN_UUID}  
**Started:** ${STARTED_AT}

---

## 📋 Summary

| Property | Value |
|----------|-------|
| **ZINC ID** | ${ZINC_ID} |
| **SMILES** | \`${SMILES}\` |
| **Vina Affinity** | **${AFFINITY} kcal/mol** |
| **GROMACS Status** | ${GMX_STATUS:-pending} |
| **MD Status** | ${MD_STATUS:-pending} |

---

## 🎯 Target Information

| Property | Value |
|----------|-------|
| **Target Name** | ${TARGET_NAME} |
| **Target Code** | ${TARGET_CODE} |
| **PDB ID** | [${PDB_ID}](https://www.rcsb.org/structure/${PDB_ID}) |
| **Receptor File** | \`${RECEPTOR}\` |

---

## 📊 Docking Results

### Vina Pose
- **Pose File:** \`${POSE_FILE}\`
- **Binding Affinity:** ${AFFINITY} kcal/mol

DETAIL

    # Add Vina log if exists
    if [[ -n "${POSE_FILE}" ]]; then
        LOG_FILE="${POSE_FILE%.pdbqt}.log"
        LOG_FILE="${LOG_FILE/poses/logs}"
        
        if [[ -f "${ROOT}/${LOG_FILE}" ]]; then
            cat >> "${DETAIL_FILE}" <<DETAIL

### Docking Log

\`\`\`
$(tail -30 "${ROOT}/${LOG_FILE}")
\`\`\`
DETAIL
        fi
    fi
    
    # Add GROMACS results if exists
    if [[ "${GMX_STATUS}" == "completed" ]] && [[ -n "${GMX_DIR}" ]]; then
        cat >> "${DETAIL_FILE}" <<DETAIL

---

## 🧪 GROMACS MD Preparation

**Status:** ✅ Completed  
**Output Directory:** \`${GMX_DIR}\`

### Generated Files

DETAIL
        
        if [[ -d "${ROOT}/${GMX_DIR}" ]]; then
            cat >> "${DETAIL_FILE}" <<DETAIL

\`\`\`
$(ls -lh "${ROOT}/${GMX_DIR}"/*.{gro,top,itp,pdb} 2>/dev/null | awk '{print $9, $5}' || echo "No files found")
\`\`\`
DETAIL
        fi
    fi
    
    # Add MD results if exists
    if [[ "${MD_STATUS}" == "completed" ]] && [[ -n "${MD_DIR}" ]]; then
        MD_JSON="${ROOT}/${MD_DIR}/md_summary.json"
        
        if [[ -f "${MD_JSON}" ]]; then
            cat >> "${DETAIL_FILE}" <<DETAIL

---

## 🔬 MD Simulation Results

**Status:** ✅ Completed

### Summary

\`\`\`json
$(cat "${MD_JSON}")
\`\`\`

DETAIL
            
            # Extract energy if available
            if command -v jq &> /dev/null; then
                FINAL_ENERGY=$(jq -r '.final_potential_energy // "N/A"' "${MD_JSON}")
                MD_STAGE=$(jq -r '.md_stage // "N/A"' "${MD_JSON}")
                
                cat >> "${DETAIL_FILE}" <<DETAIL

### Key Metrics

| Metric | Value |
|--------|-------|
| **MD Stage** | ${MD_STAGE} |
| **Final Energy** | ${FINAL_ENERGY} kJ/mol |
| **Vina Affinity** | ${AFFINITY} kcal/mol |

DETAIL
            fi
        fi
        
        # Check for analysis plots
        if [[ -d "${ROOT}/${MD_DIR}" ]]; then
            PLOTS=$(find "${ROOT}/${MD_DIR}" -name "*.png" -o -name "*.svg" 2>/dev/null | head -5)
            
            if [[ -n "${PLOTS}" ]]; then
                cat >> "${DETAIL_FILE}" <<DETAIL

### Analysis Plots

DETAIL
                echo "${PLOTS}" | while read -r plot; do
                    REL_PLOT=$(realpath --relative-to="${ROOT}/${PAGES_DIR}" "${plot}")
                    PLOT_NAME=$(basename "${plot}")
                    echo "![${PLOT_NAME}](../../${REL_PLOT})" >> "${DETAIL_FILE}"
                    echo "" >> "${DETAIL_FILE}"
                done
            fi
        fi
    fi
    
    # Add back navigation
    cat >> "${DETAIL_FILE}" <<DETAIL

---

[← Back to Dashboard](../index.html)
DETAIL
    
    ((DETAIL_COUNT++))
done

echo "[✓] Generated ${DETAIL_COUNT} detail pages"

# ========================================
# Summary Statistics Page
# ========================================
echo "[3/3] Generating statistics page..."

cat > "${PAGES_DIR}/statistics.md" <<EOF
---
title: Statistics
layout: default
---

# 📈 OCP Statistics

**Generated:** ${TIMESTAMP}

## Pipeline Performance

EOF

sqlite3 "${DB_PATH}" <<'SQL' >> "${PAGES_DIR}/statistics.md"
.mode markdown

SELECT 
    t.name AS Target,
    COUNT(DISTINCT r.id) AS Runs,
    COUNT(DISTINCT vr.ligand_id) AS Ligands,
    printf('%.4f', AVG(vr.affinity)) AS 'Avg Vina Score',
    printf('%.4f', MIN(vr.affinity)) AS 'Best Score',
    SUM(CASE WHEN vr.gromacs_prep_status = 'completed' THEN 1 ELSE 0 END) AS 'GROMACS Done',
    SUM(CASE WHEN vr.md_status = 'completed' THEN 1 ELSE 0 END) AS 'MD Done'
FROM targets t
LEFT JOIN runs r ON t.id = r.target_id
LEFT JOIN vina_results vr ON r.id = vr.run_id AND vr.mode_rank = 1
GROUP BY t.id
ORDER BY COUNT(DISTINCT r.id) DESC;
SQL

cat >> "${PAGES_DIR}/statistics.md" <<'HEADER'

## Success Rates

| Stage | Total | Success | Rate |
|-------|-------|---------|------|
HEADER

sqlite3 "${DB_PATH}" -separator '|' <<'SQL' | while IFS='|' read -r stage total success rate; do
SELECT 
    'Vina Docking',
    COUNT(*),
    SUM(CASE WHEN affinity IS NOT NULL THEN 1 ELSE 0 END),
    printf('%.1f%%', 100.0 * SUM(CASE WHEN affinity IS NOT NULL THEN 1 ELSE 0 END) / COUNT(*))
FROM vina_results WHERE mode_rank = 1
UNION ALL
SELECT 
    'GROMACS Prep',
    COUNT(*),
    SUM(CASE WHEN gromacs_prep_status = 'completed' THEN 1 ELSE 0 END),
    printf('%.1f%%', 100.0 * SUM(CASE WHEN gromacs_prep_status = 'completed' THEN 1 ELSE 0 END) / COUNT(*))
FROM vina_results WHERE mode_rank = 1
UNION ALL
SELECT 
    'MD Simulation',
    COUNT(*),
    SUM(CASE WHEN md_status = 'completed' THEN 1 ELSE 0 END),
    printf('%.1f%%', 100.0 * SUM(CASE WHEN md_status = 'completed' THEN 1 ELSE 0 END) / COUNT(*))
FROM vina_results WHERE mode_rank = 1;
SQL
    echo "| ${stage} | ${total} | ${success} | ${rate} |" >> "${PAGES_DIR}/statistics.md"
done

cat >> "${PAGES_DIR}/statistics.md" <<'EOF'

---

[← Back to Dashboard](index.html)
EOF

echo "[✓] Statistics page generated"

# ========================================
# Summary
# ========================================
echo ""
echo "========================================="
echo "Pages Generation Complete!"
echo "========================================="
echo "Main index    : ${PAGES_DIR}/index.md"
echo "Detail pages  : ${DETAIL_COUNT} pages in ${PAGES_DIR}/runs/"
echo "Statistics    : ${PAGES_DIR}/statistics.md"
echo "Timestamp     : ${TIMESTAMP}"
echo ""
echo "Files generated:"
find "${PAGES_DIR}" -name "*.md" -type f | wc -l
echo "========================================="
