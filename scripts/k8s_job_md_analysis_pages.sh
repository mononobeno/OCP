#!/bin/bash
# MD解析結果をGitHub Pagesとして生成

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
DB="${ROOT}/catalog/db/ocp_results.sqlite"

echo "==========================================="
echo "MD解析結果ページ生成"
echo "==========================================="

# 完了したMD結果を取得
RESULTS=$(sqlite3 -separator '|' "${DB}" <<'SQL'
SELECT 
    v.rowid,
    l.zinc_id,
    l.smiles,
    v.affinity_kcal,
    v.md_output_dir,
    v.md_rmsd_avg,
    v.md_simulation_time_ns,
    v.md_performance_nsday
FROM vina_results v
JOIN ligands l ON v.ligand_id = l.id
WHERE v.md_status = 'completed'
ORDER BY v.md_rmsd_avg ASC;
SQL
)

if [[ -z "${RESULTS}" ]]; then
    echo "⚠️  完了したMD結果がありません"
    exit 0
fi

# MD結果テーブルページ生成
PAGES_DIR="${ROOT}/docs/pages/md"
mkdir -p "${PAGES_DIR}"

cat > "${PAGES_DIR}/index.md" <<'HEADER'
---
layout: default
title: MD Simulation Results
---

# MD Simulation Results

分子動力学シミュレーション完了結果の一覧

| ZINC ID | SMILES | Vina (kcal/mol) | RMSD (nm) | Time (ns) | Performance (ns/day) | Detail |
|---------|--------|-----------------|-----------|-----------|---------------------|--------|
HEADER

echo "${RESULTS}" | while IFS='|' read -r rowid zinc_id smiles affinity md_dir rmsd sim_time perf; do
    # SMILESを30文字に切り詰め
    smiles_short="${smiles:0:30}"
    [[ ${#smiles} -gt 30 ]] && smiles_short="${smiles_short}..."
    
    # テーブル行追加
    echo "| ${zinc_id} | \`${smiles_short}\` | ${affinity} | ${rmsd} | ${sim_time} | ${perf} | [詳細](./result_${rowid}.html) |" >> "${PAGES_DIR}/index.md"
    
    # 個別結果ページ生成
    cat > "${PAGES_DIR}/result_${rowid}.md" <<DETAIL_EOF
---
layout: default
title: MD Result - ${zinc_id}
---

# MD Simulation Result: ${zinc_id}

[← Back to MD Results](./index.html)

## 基本情報

- **ZINC ID**: ${zinc_id}
- **SMILES**: \`${smiles}\`
- **Vina Affinity**: ${affinity} kcal/mol

## MD Simulation

- **Status**: ✅ Completed
- **Output Directory**: \`${md_dir}\`
- **RMSD平均**: **${rmsd} nm**
- **シミュレーション時間**: ${sim_time} ns
- **計算性能**: ${perf} ns/day

## 軌跡解析

\`\`\`bash
# RMSD plot
gmx rms -s ${md_dir}/md_*.tpr -f ${md_dir}/md_*.xtc -o rmsd.xvg

# RMSF plot
gmx rmsf -s ${md_dir}/md_*.tpr -f ${md_dir}/md_*.xtc -o rmsf.xvg

# Radius of gyration
gmx gyrate -s ${md_dir}/md_*.tpr -f ${md_dir}/md_*.xtc -o gyrate.xvg
\`\`\`

DETAIL_EOF

done

echo ""
echo "✅ MD解析ページ生成完了"
echo "   Pages: ${PAGES_DIR}/index.md"
echo "   Results: $(echo "${RESULTS}" | wc -l) 件"
