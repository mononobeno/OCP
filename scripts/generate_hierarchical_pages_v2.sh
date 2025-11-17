#!/bin/bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DB="${ROOT}/catalog/db/ocp_results.sqlite"
PAGES_DIR="${ROOT}/docs/pages"

echo "Generating hierarchical pages..."

mkdir -p "${PAGES_DIR}/targets" "${PAGES_DIR}/compounds"

# トップページ
cat > "${PAGES_DIR}/index.md" <<'EOF'
---
layout: default
title: OCP Pipeline
---

# Open Compound Pipeline

## Pipeline Architecture

```
ZINC20 → RDKit 3D → Vina Docking → GROMACS Prep → MD Simulation → Analysis
```

## Kubernetes Jobs

| Job | Status | Description |
|-----|--------|-------------|
| vina-prep | ✅ | Vina前処理 |
| vina-runner | ✅ | GPU Docking |
| md-production | ✅ | GPU MD |
| md-analysis | ✅ | 軌跡解析 |

## Targets

| PDB | Name | Vina Results | MD Done | Link |
|-----|------|--------------|---------|------|
EOF

sqlite3 -separator '|' "${DB}" "SELECT t.id, t.pdb_id, t.name, COUNT(v.rowid), SUM(CASE WHEN v.md_status='completed' THEN 1 ELSE 0 END) FROM targets t LEFT JOIN runs r ON r.target_id=t.id LEFT JOIN vina_results v ON v.run_id=r.id GROUP BY t.id;" | while IFS='|' read -r tid pdb name vcnt mcnt; do
    echo "| ${pdb} | ${name} | ${vcnt} | ${mcnt} | [View](./targets/target_${tid}.html) |" >> "${PAGES_DIR}/index.md"
done

echo "" >> "${PAGES_DIR}/index.md"
echo "_Updated: $(date '+%Y-%m-%d %H:%M:%S')_" >> "${PAGES_DIR}/index.md"

# ターゲットページ
sqlite3 -list -separator '|' "${DB}" "SELECT DISTINCT t.id, t.pdb_id, t.name, t.description FROM targets t;" | while IFS='|' read -r tid pdb name desc; do
    cat > "${PAGES_DIR}/targets/target_${tid}.md" <<TEOF
---
layout: default
title: ${name}
---

# ${name} (${pdb})

[← Back](../index.html)

## Protein Info

- **PDB ID**: ${pdb}
- **Description**: ${desc}
- **Structure**: [RCSB](https://www.rcsb.org/structure/${pdb})

## Compounds

| ZINC ID | SMILES | Vina (kcal/mol) | MD | Link |
|---------|--------|-----------------|----|----|
TEOF

    sqlite3 -separator '|' "${DB}" "SELECT v.rowid, l.zinc_id, l.smiles, v.affinity_kcal, v.md_status FROM vina_results v JOIN ligands l ON v.ligand_id=l.id JOIN runs r ON v.run_id=r.id WHERE r.target_id=${tid} ORDER BY v.affinity_kcal LIMIT 100;" | while IFS='|' read -r vrid zid smi aff mds; do
        smi_short="${smi:0:30}"
        [[ ${#smi} -gt 30 ]] && smi_short="${smi_short}..."
        md_badge="⏳"
        [[ "${mds}" == "completed" ]] && md_badge="✅"
        echo "| ${zid} | \`${smi_short}\` | ${aff} | ${md_badge} | [View](../compounds/compound_${vrid}.html) |" >> "${PAGES_DIR}/targets/target_${tid}.md"
    done
done

# 化合物詳細ページ
sqlite3 -separator '|' "${DB}" "SELECT v.rowid, l.zinc_id, l.smiles, v.affinity_kcal, v.md_output_dir, v.md_rmsd_avg, v.md_simulation_time_ns, v.md_performance_nsday, t.pdb_id, t.name FROM vina_results v JOIN ligands l ON v.ligand_id=l.id JOIN runs r ON v.run_id=r.id JOIN targets t ON r.target_id=t.id WHERE v.md_status='completed';" | while IFS='|' read -r vrid zid smi aff mddir rmsd simtime perf pdb tname; do
    cat > "${PAGES_DIR}/compounds/compound_${vrid}.md" <<CEOF
---
layout: default
title: ${zid}
---

# ${zid} - MD Analysis

[← Back to ${tname}](../targets/target_$(sqlite3 "${DB}" "SELECT r.target_id FROM vina_results v JOIN runs r ON v.run_id=r.id WHERE v.rowid=${vrid};").html)

## Compound Info

- **ZINC ID**: ${zid}
- **SMILES**: \`${smi}\`
- **Target**: ${tname} (${pdb})
- **Vina Score**: ${aff} kcal/mol

## MD Results

- **Status**: ✅ Completed
- **Time**: ${simtime} ns
- **RMSD**: ${rmsd} nm
- **Performance**: ${perf} ns/day

---

## RMSD Analysis

<div id="rmsd-plot"></div>

<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
<script>
CEOF

    # RMSD データ
    if [[ -f "${mddir}/rmsd_backbone.xvg" ]]; then
        echo "const rmsdData = {x: [" >> "${PAGES_DIR}/compounds/compound_${vrid}.md"
        awk '/^[^@#]/ {printf "%s,", $1}' "${mddir}/rmsd_backbone.xvg" >> "${PAGES_DIR}/compounds/compound_${vrid}.md"
        echo "], y: [" >> "${PAGES_DIR}/compounds/compound_${vrid}.md"
        awk '/^[^@#]/ {printf "%s,", $2}' "${mddir}/rmsd_backbone.xvg" >> "${PAGES_DIR}/compounds/compound_${vrid}.md"
        echo "]};" >> "${PAGES_DIR}/compounds/compound_${vrid}.md"
        echo "Plotly.newPlot('rmsd-plot', [{x: rmsdData.x, y: rmsdData.y, type: 'scatter', mode: 'lines', name: 'RMSD'}], {title: 'Backbone RMSD', xaxis: {title: 'Time (ns)'}, yaxis: {title: 'RMSD (nm)'}});" >> "${PAGES_DIR}/compounds/compound_${vrid}.md"
    fi

    echo "</script>" >> "${PAGES_DIR}/compounds/compound_${vrid}.md"

    # H-bond
    cat >> "${PAGES_DIR}/compounds/compound_${vrid}.md" <<'HBEOF'

## Hydrogen Bonds

<div id="hbond-plot"></div>

<script>
HBEOF

    if [[ -f "${mddir}/hbond.xvg" ]]; then
        echo "const hbondData = {x: [" >> "${PAGES_DIR}/compounds/compound_${vrid}.md"
        awk '/^[^@#]/ {printf "%s,", $1}' "${mddir}/hbond.xvg" >> "${PAGES_DIR}/compounds/compound_${vrid}.md"
        echo "], y: [" >> "${PAGES_DIR}/compounds/compound_${vrid}.md"
        awk '/^[^@#]/ {printf "%s,", $2}' "${mddir}/hbond.xvg" >> "${PAGES_DIR}/compounds/compound_${vrid}.md"
        echo "]};" >> "${PAGES_DIR}/compounds/compound_${vrid}.md"
        echo "Plotly.newPlot('hbond-plot', [{x: hbondData.x, y: hbondData.y, type: 'scatter', mode: 'lines', name: 'H-bonds', line: {color: 'orange'}}], {title: 'Hydrogen Bonds', xaxis: {title: 'Time (ns)'}, yaxis: {title: 'Count'}});" >> "${PAGES_DIR}/compounds/compound_${vrid}.md"
    fi

    echo "</script>" >> "${PAGES_DIR}/compounds/compound_${vrid}.md"

    # 平衡化チェック
    cat >> "${PAGES_DIR}/compounds/compound_${vrid}.md" <<'EQEOF'

---

## Equilibration Check

### Energy Minimization

<div id="em-plot"></div>

### NVT Temperature

<div id="nvt-temp-plot"></div>

### NPT Pressure

<div id="npt-pressure-plot"></div>

<script>
// EM/NVT/NPT plots will be added here
</script>

_Generated: $(date)_
EQEOF

done

echo "✅ Pages generated: ${PAGES_DIR}"
