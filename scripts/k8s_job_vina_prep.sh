#!/usr/bin/env bash
#
# k8s_job_vina_prep.sh
#
# 役割:
#   - Receptor PDB → PDBQT 変換（MGLTools を想定）
#   - RDKit 整形済み ligand PDB/SDF から ligand PDBQT 生成
#   - /workspace/vina_input に receptor.pdbqt + ligands PDBQT を並べる。
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
VINA_INPUT_DIR="${WORKSPACE}/vina_input"
LIGANDS_RAW_DIR="${WORKSPACE}/ligands_raw"

mkdir -p "${VINA_INPUT_DIR}"

if [[ ! -f "${DB_PATH}" ]]; then
  echo "ERROR: DB not found: ${DB_PATH}" >&2
  exit 1
fi

echo "[INFO] vina-prep: RUN_ID=${RUN_ID}"
echo "[INFO] DB_PATH   : ${DB_PATH}"
echo "[INFO] VINA_INPUT: ${VINA_INPUT_DIR}"

# Receptor 情報取得:
#   runs.target_id -> targets.id -> targets.code, targets.receptor_pdb_path などを想定。
sqlite3 "${DB_PATH}" <<SQL
.headers on
.mode column
SELECT
  r.id         AS run_id,
  t.id         AS target_id,
  t.code       AS target_code,
  t.receptor_pdb_path,
  t.receptor_pdbqt_path,
  t.receptor_preparation_method
FROM runs r
JOIN targets t ON t.id = r.target_id
WHERE r.id = ${RUN_ID};
SQL

# 実際の receptor.pdb -> receptor.pdbqt の変換は MGLTools に委ねる。
# ここではパスの存在チェックと雛形コマンドのみ。
RECEPTOR_PDB_PATH="$(sqlite3 "${DB_PATH}" "SELECT t.receptor_pdb_path FROM runs r JOIN targets t ON t.id = r.target_id WHERE r.id = ${RUN_ID} LIMIT 1;")"

if [[ -z "${RECEPTOR_PDB_PATH}" || "${RECEPTOR_PDB_PATH}" = "NULL" ]]; then
  echo "[WARN] receptor_pdb_path is not set in DB; please update targets.receptor_pdb_path." >&2
else
  ABS_RECEPTOR_PDB="${REPO_ROOT}/${RECEPTOR_PDB_PATH}"
  echo "[INFO] receptor PDB: ${ABS_RECEPTOR_PDB}"
  if [[ -f "${ABS_RECEPTOR_PDB}" ]]; then
    echo "[INFO] (NOTE) ここで MGLTools の prepare_receptor4.py を呼び出して receptor.pdbqt を生成する想定です。"
    echo "[INFO] 例:"
    echo "  prepare_receptor4.py -r ${ABS_RECEPTOR_PDB} -o ${VINA_INPUT_DIR}/receptor.pdbqt"
  else
    echo "[WARN] receptor PDB not found at: ${ABS_RECEPTOR_PDB}" >&2
  fi
fi

echo "[INFO] Preparing ligand PDBQT from ${LIGANDS_RAW_DIR} into ${VINA_INPUT_DIR} ..."

if ! command -v obabel >/dev/null 2>&1; then
  echo "[WARN] obabel not found; skipping actual PDBQT conversion. This is a skeleton." >&2
else
  for f in "${LIGANDS_RAW_DIR}"/*; do
    [[ -e "$f" ]] || continue
    base="$(basename "$f")"
    zinc_id="${base%%_*}"
    out_pdbqt="${VINA_INPUT_DIR}/${zinc_id}.pdbqt"
    echo "[INFO] obabel: $f -> ${out_pdbqt}"
    obabel "$f" -O "${out_pdbqt}" >/dev/null 2>&1 || echo "[WARN] obabel failed for $f" >&2
  done
fi

echo "[INFO] vina-prep finished."
