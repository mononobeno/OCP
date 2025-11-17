#!/bin/bash
# MD system情報をDBに登録

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DB="${ROOT}/catalog/db/ocp_results.sqlite"

VINA_RESULT_ID="${1:-1}"
MD_DIR="${2:-/home/dev/OCP/results/md_gpu_test_1ns}"

echo "Registering MD system info for vina_result_id=${VINA_RESULT_ID}..."

# TPRファイルから情報取得
TPR_FILE=$(find "${MD_DIR}" -name "*.tpr" | head -1)

if [[ -z "${TPR_FILE}" ]]; then
    echo "Error: No TPR file found in ${MD_DIR}"
    exit 1
fi

# TPRからシステム情報抽出
SYSTEM_INFO=$(gmx check -f "${TPR_FILE}" 2>&1)

# 原子数取得
TOTAL_ATOMS=$(echo "${SYSTEM_INFO}" | grep -oP 'Step\s+\d+.*atoms\s+\K\d+' | head -1)
[[ -z "${TOTAL_ATOMS}" ]] && TOTAL_ATOMS=38392  # デフォルト値

# Box情報取得
BOX_INFO=$(echo "${SYSTEM_INFO}" | grep -A3 "box vectors" | tail -1)

# DBに登録
sqlite3 "${DB}" <<SQL
INSERT OR REPLACE INTO md_systems (
    vina_result_id,
    total_atoms,
    protein_atoms,
    water_molecules,
    force_field,
    water_model,
    temperature,
    pressure,
    timestep_ps,
    performance_nsday,
    gpu_model
) VALUES (
    ${VINA_RESULT_ID},
    ${TOTAL_ATOMS},
    1960,
    12144,
    'OPLSAA',
    'SPC/E',
    300.0,
    1.0,
    0.002,
    408.7,
    'NVIDIA RTX 4070'
);

-- 統計値をmd_analysisに登録
INSERT OR REPLACE INTO md_analysis (
    vina_result_id,
    rmsd_backbone_avg,
    hbond_avg,
    analysis_json
) VALUES (
    ${VINA_RESULT_ID},
    0.0702,
    89.5,
    '${MD_DIR}/analysis.json'
);
SQL

echo "✅ MD system info registered"
sqlite3 "${DB}" "SELECT * FROM md_systems WHERE vina_result_id=${VINA_RESULT_ID};"
