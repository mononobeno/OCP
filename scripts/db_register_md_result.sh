#!/bin/bash
# MD結果をDBに登録

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
DB="${ROOT}/catalog/db/ocp_results.sqlite"

VINA_RESULT_ID="${1:-}"
MD_OUTPUT_DIR="${2:-}"
RMSD_AVG="${3:-}"
MD_SIMULATION_TIME="${4:-}"
MD_PERFORMANCE="${5:-}"

if [[ -z "${VINA_RESULT_ID}" || -z "${MD_OUTPUT_DIR}" ]]; then
    echo "Usage: $0 <vina_result_id> <md_output_dir> [rmsd_avg] [simulation_time] [performance]"
    exit 1
fi

# MD結果を登録
sqlite3 "${DB}" <<SQL
UPDATE vina_results
SET 
    md_status = 'completed',
    md_output_dir = '${MD_OUTPUT_DIR}',
    md_rmsd_avg = ${RMSD_AVG:-NULL},
    md_simulation_time_ns = ${MD_SIMULATION_TIME:-NULL},
    md_performance_nsday = ${MD_PERFORMANCE:-NULL}
WHERE rowid = ${VINA_RESULT_ID};
SQL

echo "✅ MD結果登録完了: vina_results.id=${VINA_RESULT_ID}"
echo "   MD出力: ${MD_OUTPUT_DIR}"
[[ -n "${RMSD_AVG}" ]] && echo "   RMSD平均: ${RMSD_AVG} nm"
[[ -n "${MD_SIMULATION_TIME}" ]] && echo "   シミュレーション時間: ${MD_SIMULATION_TIME} ns"
[[ -n "${MD_PERFORMANCE}" ]] && echo "   性能: ${MD_PERFORMANCE} ns/day"
