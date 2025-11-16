#!/usr/bin/env bash
#
# k8s_job_ligand_selector.sh
#
# 役割:
#   - ZINC 化合物ライブラリ（RDKit 整形済み & DB 登録済み）から、
#     対象ライブラリ / RUN 用の ligand を選択し、
#     /workspace/ligands_raw にコピーするステージ。
#   - どの ligand を選ぶかは DB を見て決定する。
#
# 共通仕様:
#   - RUN_ID は第1引数、または環境変数 RUN_ID から取得する。
#   - DB_PATH は環境変数 DB_PATH_OVERRIDE で上書き可能。
#   - DB は catalog/db/ocp_results.sqlite を前提とし、
#     ここから対象 ligand / target / run を決める（DB = source of truth）。
#
set -euo pipefail

RUN_ID="${1:-${RUN_ID:-}}"
if [[ -z "${RUN_ID}" ]]; then
  echo "Usage: $0 <RUN_ID>  # or set RUN_ID env" >&2
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="${SCRIPT_DIR%/scripts}"
DB_PATH="${DB_PATH_OVERRIDE:-/home/dev/OCP/catalog/db/ocp_results.sqlite}"

WORKSPACE="${WORKSPACE:-/workspace}"
LIGANDS_RAW_DIR="${WORKSPACE}/ligands_raw"

mkdir -p "${LIGANDS_RAW_DIR}"

if [[ ! -f "${DB_PATH}" ]]; then
  echo "ERROR: DB not found: ${DB_PATH}" >&2
  exit 1
fi

if ! command -v sqlite3 >/dev/null 2>&1; then
  echo "ERROR: sqlite3 command not found." >&2
  exit 1
fi

echo "[INFO] ligand-selector: RUN_ID=${RUN_ID}"
echo "[INFO] DB_PATH       : ${DB_PATH}"
echo "[INFO] LIGANDS_RAW   : ${LIGANDS_RAW_DIR}"

echo "[INFO] Selecting ligands for RUN_ID=${RUN_ID} from DB..."

# 想定スキーマ:
#   runs(id, library_id, target_id, ...)
#   run_ligands(run_id, ligand_id, order_index, ...)
#   ligands(id, zinc_id, source_file, has_3d, conformer_method, ...)
#
# RDKit による整形済み ligand のみを対象:
#   - has_3d = 1
#   - conformer_method = 'rdkit_etkdg'
#
sqlite3 "${DB_PATH}" <<SQL
.headers on
.mode column
SELECT
  rl.order_index,
  l.id AS ligand_id,
  l.zinc_id,
  l.source_file,
  l.has_3d,
  l.conformer_method
FROM run_ligands rl
JOIN ligands l ON l.id = rl.ligand_id
WHERE rl.run_id = ${RUN_ID}
  AND l.has_3d = 1
  AND l.conformer_method = 'rdkit_etkdg'
ORDER BY rl.order_index;
SQL

echo "[INFO] Copying selected ligands into ${LIGANDS_RAW_DIR} (if source_file exists)..."

sqlite3 -csv "${DB_PATH}" <<SQL | while IFS=',' read -r order_index ligand_id zinc_id source_file has_3d conformer_method; do
SELECT
  rl.order_index,
  l.id AS ligand_id,
  l.zinc_id,
  l.source_file,
  l.has_3d,
  l.conformer_method
FROM run_ligands rl
JOIN ligands l ON l.id = rl.ligand_id
WHERE rl.run_id = ${RUN_ID}
  AND l.has_3d = 1
  AND l.conformer_method = 'rdkit_etkdg'
ORDER BY rl.order_index;
SQL
do
  if [[ -z "${source_file}" || "${source_file}" = "NULL" ]]; then
    echo "[WARN] ligand_id=${ligand_id} zinc_id=${zinc_id}: source_file is NULL, skipping." >&2
    continue
  fi
  SRC_PATH="${REPO_ROOT}/${source_file}"
  if [[ ! -f "${SRC_PATH}" ]]; then
    echo "[WARN] ligand_id=${ligand_id} zinc_id=${zinc_id}: source file not found: ${SRC_PATH}, skipping." >&2
    continue
  fi
  BASENAME="$(basename "${SRC_PATH}")"
  DEST_PATH="${LIGANDS_RAW_DIR}/${zinc_id}_${BASENAME}"
  echo "[INFO] copy: ${SRC_PATH} -> ${DEST_PATH}"
  cp "${SRC_PATH}" "${DEST_PATH}"
done

echo "[INFO] ligand-selector finished."
