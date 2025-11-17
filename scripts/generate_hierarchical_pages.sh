#!/bin/bash
# 新しいページ構成でGitHub Pagesを生成
# 階層: トップ → ターゲットタンパク質 → 化合物詳細

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
DB="${ROOT}/catalog/db/ocp_results.sqlite"
DOCS_DIR="${ROOT}/docs"
PAGES_DIR="${DOCS_DIR}/pages"

echo "==========================================="
echo "GitHub Pages生成 (階層構造版)"
echo "==========================================="

# ディレクトリ作成
mkdir -p "${PAGES_DIR}/targets"
mkdir -p "${PAGES_DIR}/compounds"

# ===========================
# 1. トップページ生成
# ===========================
cat > "${PAGES_DIR}/index.md" <<'TOP_HEADER'
---
layout: default
title: OCP - Drug Discovery Pipeline
---

# Open Compound Pipeline (OCP)

GPU加速による大規模ドラッグディスカバリーパイプライン

## Pipeline Architecture

```
┌─────────────────┐
│  ZINC20 Library │
│  (ML Subset)    │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│  RDKit 3D Gen   │
│  (ETKDG)        │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│  AutoDock Vina  │
│  (GPU)          │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│  GROMACS Prep   │
│  (ACPYPE)       │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│  MD Simulation  │
│  (GPU, GROMACS) │
└────────┬────────┘
         │
         ▼
┌─────────────────┐
│  Analysis &     │
│  Visualization  │
└─────────────────┘
```

## Kubernetes Job Status

| Job Type | Status | Description |
|----------|--------|-------------|
| `ligand-selector` | ✅ Active | 3D未生成リガンド選択 |
| `vina-prep` | ✅ Active | Vina前処理バッチ |
| `vina-runner` | ✅ Active | GPU Vinaドッキング |
| `gromacs-prep` | 🔧 Dev | GROMACS topology生成 |
| `md-equilibration` | 🔧 Dev | NVT/NPT平衡化 |
| `md-production` | ✅ Validated | 本番MD (GPU) |
| `md-analysis` | ✅ Validated | 軌跡解析・ページ生成 |

## Target Proteins

TOP_HEADER

# ターゲットタンパク質一覧取得
TARGETS=$(sqlite3 -separator '|' "${DB}" <<'SQL'
SELECT DISTINCT 
    t.id,
    t.pdb_id,
    t.name,
    t.description,
    COUNT(DISTINCT r.id) as run_count,
    COUNT(DISTINCT v.rowid) as result_count,
    SUM(CASE WHEN v.md_status = 'completed' THEN 1 ELSE 0 END) as md_completed
FROM targets t
LEFT JOIN runs r ON r.target_id = t.id
LEFT JOIN vina_results v ON v.run_id = r.id
GROUP BY t.id
ORDER BY t.id;
SQL
)

if [[ -n "${TARGETS}" ]]; then
    echo "" >> "${PAGES_DIR}/index.md"
    echo "| PDB ID | Target Name | Runs | Vina Results | MD Completed | Details |" >> "${PAGES_DIR}/index.md"
    echo "|--------|-------------|------|--------------|--------------|---------|" >> "${PAGES_DIR}/index.md"
    
    echo "${TARGETS}" | while IFS='|' read -r target_id pdb_id name desc run_count result_count md_completed; do
        echo "| ${pdb_id} | ${name} | ${run_count} | ${result_count} | ${md_completed} | [View](./targets/target_${target_id}.html) |" >> "${PAGES_DIR}/index.md"
    done
else
    echo "⚠️  No targets found" >> "${PAGES_DIR}/index.md"
fi

cat >> "${PAGES_DIR}/index.md" <<'TOP_FOOTER'

## System Statistics

```sql
-- Database: catalog/db/ocp_results.sqlite
```

| Metric | Count |
|--------|-------|
TOP_FOOTER

# 統計情報追加
sqlite3 -separator '|' "${DB}" <<'SQL' | while IFS='|' read -r metric value; do
    echo "| ${metric} | ${value} |" >> "${PAGES_DIR}/index.md"
done
SELECT 'Total Ligands', COUNT(*) FROM ligands
UNION ALL SELECT '3D Generated', COUNT(*) FROM ligands WHERE has_3d = 1
UNION ALL SELECT 'Vina Results', COUNT(*) FROM vina_results
UNION ALL SELECT 'MD Completed', COUNT(*) FROM vina_results WHERE md_status = 'completed';
SQL

echo "" >> "${PAGES_DIR}/index.md"
echo "_Last updated: $(date '+%Y-%m-%d %H:%M:%S')_" >> "${PAGES_DIR}/index.md"

echo "✅ トップページ生成: ${PAGES_DIR}/index.md"

# ===========================
# 2. ターゲットタンパク質詳細ページ生成
# ===========================
echo ""
echo "Generating target pages..."

echo "${TARGETS}" | while IFS='|' read -r target_id pdb_id name desc run_count result_count md_completed; do
    TARGET_PAGE="${PAGES_DIR}/targets/target_${target_id}.md"
    
    cat > "${TARGET_PAGE}" <<TARGET_HEADER
---
layout: default
title: ${name} (${pdb_id})
---

# Target: ${name}

[← Back to Top](../index.html)

## Protein Information

- **PDB ID**: ${pdb_id}
- **Name**: ${name}
- **Description**: ${desc}
- **3D Structure**: [RCSB PDB](https://www.rcsb.org/structure/${pdb_id})

## Calculation Summary

- **Total Runs**: ${run_count}
- **Vina Results**: ${result_count}
- **MD Completed**: ${md_completed}

## Compound Results

| ZINC ID | SMILES | Vina Score (kcal/mol) | MD Status | Detail |
|---------|--------|-----------------------|-----------|--------|
TARGET_HEADER

    # このターゲットの化合物結果を取得
    sqlite3 -separator '|' "${DB}" <<SQL | while IFS='|' read -r vina_rowid zinc_id smiles affinity md_status; do
        smiles_short="\${smiles:0:40}"
        [[ \${#smiles} -gt 40 ]] && smiles_short="\${smiles_short}..."
        
        md_badge="⏳ Pending"
        [[ "\${md_status}" == "completed" ]] && md_badge="✅ Done"
        [[ "\${md_status}" == "running" ]] && md_badge="🔄 Running"
        
        echo "| \${zinc_id} | \\\`\${smiles_short}\\\` | \${affinity} | \${md_badge} | [View](../compounds/compound_\${vina_rowid}.html) |" >> "${TARGET_PAGE}"
    done
SELECT 
    v.rowid,
    l.zinc_id,
    l.smiles,
    v.affinity_kcal,
    v.md_status
FROM vina_results v
JOIN ligands l ON v.ligand_id = l.id
JOIN runs r ON v.run_id = r.id
WHERE r.target_id = ${target_id}
ORDER BY v.affinity_kcal ASC
LIMIT 100;
SQL

    echo "" >> "${TARGET_PAGE}"
    echo "_Generated: $(date '+%Y-%m-%d %H:%M:%S')_" >> "${TARGET_PAGE}"
    
    echo "  ✓ Target ${pdb_id} (${result_count} compounds)"
done

echo "✅ ターゲットページ生成完了"

# ===========================
# 3. 化合物詳細ページ生成
# ===========================
echo ""
echo "Generating compound detail pages with analysis..."

COMPOUNDS=$(sqlite3 -separator '|' "${DB}" <<'SQL'
SELECT 
    v.rowid,
    l.zinc_id,
    l.smiles,
    v.affinity_kcal,
    v.md_output_dir,
    v.md_status,
    v.md_rmsd_avg,
    v.md_simulation_time_ns,
    v.md_performance_nsday,
    t.pdb_id,
    t.name as target_name,
    r.id as run_id
FROM vina_results v
JOIN ligands l ON v.ligand_id = l.id
JOIN runs r ON v.run_id = r.id
JOIN targets t ON r.target_id = t.id
WHERE v.md_status = 'completed'
ORDER BY v.rowid;
SQL
)

if [[ -n "${COMPOUNDS}" ]]; then
    echo "${COMPOUNDS}" | while IFS='|' read -r vina_rowid zinc_id smiles affinity md_dir md_status rmsd sim_time perf pdb_id target_name run_id; do
        COMPOUND_PAGE="${PAGES_DIR}/compounds/compound_${vina_rowid}.md"
        
        # 解析データ読み込み
        ANALYSIS_JSON="${md_dir}/analysis.json"
        HAS_ANALYSIS=false
        [[ -f "${ANALYSIS_JSON}" ]] && HAS_ANALYSIS=true
        
        cat > "${COMPOUND_PAGE}" <<COMPOUND_HEADER
---
layout: default
title: ${zinc_id} - MD Analysis
---

# Compound: ${zinc_id}

[← Back to ${target_name}](../targets/target_$(sqlite3 "${DB}" "SELECT target_id FROM runs WHERE id = ${run_id};").html)

## Compound Information

- **ZINC ID**: ${zinc_id}
- **SMILES**: \`${smiles}\`
- **Target**: ${target_name} (${pdb_id})
- **Vina Affinity**: **${affinity} kcal/mol**

## MD Simulation Results

- **Status**: ✅ Completed
- **Simulation Time**: ${sim_time} ns
- **RMSD (Backbone Avg)**: **${rmsd} nm**
- **Performance**: ${perf} ns/day
- **Output Directory**: \`${md_dir}\`

---

## Trajectory Analysis

### RMSD (Root Mean Square Deviation)

Backbone RMSD over simulation time shows protein stability.

\`\`\`plotly
{
  "data": [
    {
      "type": "scatter",
      "mode": "lines",
      "name": "Backbone RMSD",
      "x": [],
      "y": [],
      "line": {"color": "rgb(31, 119, 180)", "width": 2}
    }
  ],
  "layout": {
    "title": "RMSD vs Time",
    "xaxis": {"title": "Time (ns)"},
    "yaxis": {"title": "RMSD (nm)"},
    "hovermode": "closest"
  }
}
\`\`\`

COMPOUND_HEADER

        # XVGデータをJSONに変換してページに埋め込み
        if [[ -f "${md_dir}/rmsd_backbone.xvg" ]]; then
            echo "<script>" >> "${COMPOUND_PAGE}"
            echo "// RMSD data from XVG" >> "${COMPOUND_PAGE}"
            echo "const rmsdData = [" >> "${COMPOUND_PAGE}"
            awk '/^[^@#]/ {printf "[%s, %s],\n", $1, $2}' "${md_dir}/rmsd_backbone.xvg" >> "${COMPOUND_PAGE}"
            echo "];" >> "${COMPOUND_PAGE}"
            echo "// Update plot with actual data" >> "${COMPOUND_PAGE}"
            echo "if (typeof Plotly !== 'undefined') {" >> "${COMPOUND_PAGE}"
            echo "  const x = rmsdData.map(d => d[0]);" >> "${COMPOUND_PAGE}"
            echo "  const y = rmsdData.map(d => d[1]);" >> "${COMPOUND_PAGE}"
            echo "  Plotly.newPlot('rmsd-plot', [{x: x, y: y, type: 'scatter', mode: 'lines', name: 'Backbone RMSD', line: {color: 'rgb(31, 119, 180)', width: 2}}], {title: 'RMSD vs Time', xaxis: {title: 'Time (ns)'}, yaxis: {title: 'RMSD (nm)'}});" >> "${COMPOUND_PAGE}"
            echo "}" >> "${COMPOUND_PAGE}"
            echo "</script>" >> "${COMPOUND_PAGE}"
            echo '<div id="rmsd-plot" style="width:100%;height:400px;"></div>' >> "${COMPOUND_PAGE}"
        fi

        cat >> "${COMPOUND_PAGE}" <<'HBOND_SECTION'

### Hydrogen Bonds

Number of hydrogen bonds between protein and ligand over time.

<div id="hbond-plot" style="width:100%;height:400px;"></div>

HBOND_SECTION

        if [[ -f "${md_dir}/hbond.xvg" ]]; then
            echo "<script>" >> "${COMPOUND_PAGE}"
            echo "const hbondData = [" >> "${COMPOUND_PAGE}"
            awk '/^[^@#]/ {printf "[%s, %s],\n", $1, $2}' "${md_dir}/hbond.xvg" >> "${COMPOUND_PAGE}"
            echo "];" >> "${COMPOUND_PAGE}"
            echo "if (typeof Plotly !== 'undefined') {" >> "${COMPOUND_PAGE}"
            echo "  const x = hbondData.map(d => d[0]);" >> "${COMPOUND_PAGE}"
            echo "  const y = hbondData.map(d => d[1]);" >> "${COMPOUND_PAGE}"
            echo "  Plotly.newPlot('hbond-plot', [{x: x, y: y, type: 'scatter', mode: 'lines', name: 'H-bonds', line: {color: 'rgb(255, 127, 14)', width: 2}}], {title: 'Hydrogen Bonds vs Time', xaxis: {title: 'Time (ns)'}, yaxis: {title: 'Number of H-bonds'}});" >> "${COMPOUND_PAGE}"
            echo "}" >> "${COMPOUND_PAGE}"
            echo "</script>" >> "${COMPOUND_PAGE}"
        fi

        cat >> "${COMPOUND_PAGE}" <<'EQUILIBRATION_SECTION'

---

## Equilibration Quality Check

### Energy Minimization

<div id="em-energy-plot" style="width:100%;height:300px;"></div>

### NVT Temperature

<div id="nvt-temp-plot" style="width:100%;height:300px;"></div>

### NPT Pressure & Density

<div id="npt-pressure-plot" style="width:100%;height:300px;"></div>
<div id="npt-density-plot" style="width:100%;height:300px;"></div>

EQUILIBRATION_SECTION

        # EM/NVT/NPTグラフデータ埋め込み（存在する場合）
        if [[ -f "${md_dir}/../em.edr" ]]; then
            echo "<!-- EM data available -->" >> "${COMPOUND_PAGE}"
        fi

        cat >> "${COMPOUND_PAGE}" <<FOOTER

---

## Additional Analysis

\`\`\`bash
# RMSF (residue flexibility)
gmx rmsf -s ${md_dir}/md_*.tpr -f ${md_dir}/md_*.xtc -o rmsf.xvg

# Radius of gyration
gmx gyrate -s ${md_dir}/md_*.tpr -f ${md_dir}/md_*.xtc -o gyrate.xvg

# Distance analysis
gmx distance -s ${md_dir}/md_*.tpr -f ${md_dir}/md_*.xtc -select ...
\`\`\`

_Generated: $(date '+%Y-%m-%d %H:%M:%S')_

<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
FOOTER

        echo "  ✓ Compound ${zinc_id}"
    done
    
    echo "✅ 化合物詳細ページ生成完了 ($(echo "${COMPOUNDS}" | wc -l) compounds)"
else
    echo "⚠️  MD完了化合物なし"
fi

echo ""
echo "==========================================="
echo "✅ GitHub Pages生成完了"
echo "==========================================="
echo "  Top: ${PAGES_DIR}/index.md"
echo "  Targets: ${PAGES_DIR}/targets/"
echo "  Compounds: ${PAGES_DIR}/compounds/"
