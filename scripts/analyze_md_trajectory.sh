#!/bin/bash
# MD軌跡から解析データを抽出してJSON形式で出力

set -euo pipefail

MD_DIR="${1:-}"
OUTPUT_JSON="${2:-}"

if [[ -z "${MD_DIR}" || -z "${OUTPUT_JSON}" ]]; then
    echo "Usage: $0 <md_output_dir> <output_json>"
    exit 1
fi

cd "${MD_DIR}"

# 必要なファイルを確認
TPR_FILE=$(ls md_*.tpr 2>/dev/null | head -1 || ls *.tpr 2>/dev/null | head -1 || echo "")
XTC_FILE=$(ls md_*.xtc 2>/dev/null | head -1 || ls *.xtc 2>/dev/null | head -1 || echo "")
EDR_FILE=$(ls md_*.edr 2>/dev/null | head -1 || ls *.edr 2>/dev/null | head -1 || echo "")

if [[ -z "${TPR_FILE}" || -z "${XTC_FILE}" ]]; then
    echo "Error: TPR or XTC file not found in ${MD_DIR}"
    exit 1
fi

echo "Analyzing MD trajectory..."

# 1. RMSD計算 (Backbone)
if [[ ! -f rmsd_backbone.xvg ]]; then
    gmx rms -s "${TPR_FILE}" -f "${XTC_FILE}" -o rmsd_backbone.xvg -tu ns <<EOF >/dev/null 2>&1
4
4
EOF
fi

# 2. RMSD計算 (Ligand) - 存在する場合のみ
LIGAND_RMSD=""
if gmx make_ndx -f "${TPR_FILE}" -o index_temp.ndx <<EOF >/dev/null 2>&1
q
EOF
then
    # Ligand グループがあるか確認（簡易版）
    if grep -q "Other" index_temp.ndx 2>/dev/null; then
        LIGAND_RMSD="available"
    fi
    rm -f index_temp.ndx
fi

# 3. 水素結合解析 (Protein-Ligand or Protein-Protein)
HBOND_DATA=""
if [[ ! -f hbond.xvg ]] && gmx hbond -s "${TPR_FILE}" -f "${XTC_FILE}" -num hbond.xvg -tu ns <<EOF >/dev/null 2>&1
1
1
EOF
then
    HBOND_DATA="available"
fi

# 4. エネルギー解析 (EM, NVT, NPT)
ENERGY_EM=""
ENERGY_NVT=""
ENERGY_NPT=""

# EM解析
if [[ -f ../em.edr ]]; then
    gmx energy -f ../em.edr -o energy_em.xvg <<EOF >/dev/null 2>&1
Potential
EOF
    ENERGY_EM="available"
fi

# NVT解析
if [[ -f ../nvt.edr ]]; then
    gmx energy -f ../nvt.edr -o energy_nvt_temp.xvg -o energy_nvt_potential.xvg <<EOF >/dev/null 2>&1
Temperature
Potential
EOF
    ENERGY_NVT="available"
fi

# NPT解析
if [[ -f ../npt.edr ]]; then
    gmx energy -f ../npt.edr -o energy_npt_pressure.xvg -o energy_npt_density.xvg <<EOF >/dev/null 2>&1
Pressure
Density
EOF
    ENERGY_NPT="available"
fi

# 5. Radius of gyration
if [[ ! -f gyrate.xvg ]]; then
    gmx gyrate -s "${TPR_FILE}" -f "${XTC_FILE}" -o gyrate.xvg -tu ns <<EOF >/dev/null 2>&1
1
EOF
fi

# 6. RMSF計算
if [[ ! -f rmsf.xvg ]]; then
    gmx rmsf -s "${TPR_FILE}" -f "${XTC_FILE}" -o rmsf.xvg -res <<EOF >/dev/null 2>&1
4
EOF
fi

# JSONデータ生成
cat > "${OUTPUT_JSON}" <<JSON_EOF
{
  "md_dir": "${MD_DIR}",
  "tpr_file": "${TPR_FILE}",
  "xtc_file": "${XTC_FILE}",
  "edr_file": "${EDR_FILE}",
  "analysis": {
    "rmsd_backbone": "rmsd_backbone.xvg",
    "hbond": $([ -n "${HBOND_DATA}" ] && echo '"hbond.xvg"' || echo 'null'),
    "gyrate": "gyrate.xvg",
    "rmsf": "rmsf.xvg",
    "energy_em": $([ -n "${ENERGY_EM}" ] && echo '"energy_em.xvg"' || echo 'null'),
    "energy_nvt_temp": $([ -n "${ENERGY_NVT}" ] && echo '"energy_nvt_temp.xvg"' || echo 'null'),
    "energy_nvt_potential": $([ -n "${ENERGY_NVT}" ] && echo '"energy_nvt_potential.xvg"' || echo 'null'),
    "energy_npt_pressure": $([ -n "${ENERGY_NPT}" ] && echo '"energy_npt_pressure.xvg"' || echo 'null'),
    "energy_npt_density": $([ -n "${ENERGY_NPT}" ] && echo '"energy_npt_density.xvg"' || echo 'null')
  }
}
JSON_EOF

echo "✅ Analysis complete: ${OUTPUT_JSON}"
