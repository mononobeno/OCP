#!/usr/bin/env bash
#
# k8s_job_vina2gromacs_prep.sh
#
# 役割:
#   - Vina 結果 (vina_results) を DB から読み、
#     top ヒットの pose を GROMACS 用入力に変換する準備ステージ。
#   - 具体的には:
#       - top N (例: 1〜10) の pose PDBQT を取り出す
#       - receptor + ligand 複合体 PDB を生成
#       - ligand の topology (acpype 等) を作る
#   - ここでは DB 駆動 & ファイル配置の骨格のみ実装する。
#
# 共通仕様:
#   - RUN_ID は第1引数、または環境変数 RUN_ID から取得する。
#   - DB_PATH は環境変数 DB_PATH_OVERRIDE で上書き可能。
#   - DB は catalog/db/ocp_results.sqlite を前提とし、
#     ここから対象 ligand / target / run を決める（DB = source of truth）。
#
set -euo pipefail

RUN_ID="${1:-${RUN_ID:-}}"
TOP_N="${TOP_N:-1}"  # デフォルトは top 1 ヒット
if [[ -z "${RUN_ID}" ]]; then
  echo "Usage: $0 <RUN_ID>  # or set RUN_ID env" >&2
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="${SCRIPT_DIR%/scripts}"
DB_PATH="${DB_PATH_OVERRIDE:-/home/dev/OCP/catalog/db/ocp_results.sqlite}"

WORKSPACE="${WORKSPACE:-/workspace}"
VINA_OUTPUT_DIR="${WORKSPACE}/vina_output"
GMX_INPUT_DIR="${WORKSPACE}/gmx_input"

mkdir -p "${GMX_INPUT_DIR}"

if [[ ! -f "${DB_PATH}" ]]; then
  echo "ERROR: DB not found: ${DB_PATH}" >&2
  exit 1
fi

echo "[INFO] vina2gromacs-prep: RUN_ID=${RUN_ID}"
echo "[INFO] DB_PATH          : ${DB_PATH}"
echo "[INFO] VINA_OUTPUT_DIR  : ${VINA_OUTPUT_DIR}"
echo "[INFO] GMX_INPUT_DIR    : ${GMX_INPUT_DIR}"
echo "[INFO] TOP_N            : ${TOP_N}"

echo "[INFO] Top ${TOP_N} hits from vina_results (by affinity_kcal ASC):"
sqlite3 "${DB_PATH}" <<SQL
.headers on
.mode column
SELECT
  vr.id,
  vr.run_id,
  vr.ligand_id,
  l.zinc_id,
  vr.affinity_kcal,
  vr.pose_path
FROM vina_results vr
JOIN ligands l ON l.id = vr.ligand_id
WHERE vr.run_id = ${RUN_ID}
ORDER BY vr.affinity_kcal ASC
LIMIT ${TOP_N};
SQL

# ここで本来は:
#   - pose_path (PDBQT) を PDB に変換 (obabel)
#   - receptor PDB と merge して複合体 PDB を作成
#   - acpype で ligand トポロジー生成
#   を行う。ここではファイルパスを workspace にコピーする骨格のみ。

sqlite3 -csv "${DB_PATH}" <<SQL | while IFS=',' read -r vr_id run_id ligand_id zinc_id affinity_kcal pose_path; do
SELECT
  vr.id,
  vr.run_id,
  vr.ligand_id,
  l.zinc_id,
  vr.affinity_kcal,
  vr.pose_path
FROM vina_results vr
JOIN ligands l ON l.id = vr.ligand_id
WHERE vr.run_id = ${RUN_ID}
ORDER BY vr.affinity_kcal ASC
LIMIT ${TOP_N};
SQL
do
  if [[ -z "${pose_path}" || "${pose_path}" = "NULL" ]]; then
    echo "[WARN] pose_path is NULL for zinc_id=${zinc_id}, skipping." >&2
    continue
  fi
  SRC_PDBQT="${VINA_OUTPUT_DIR}/${pose_path##*/}"
  if [[ ! -f "${SRC_PDBQT}" ]]; then
    # DB に相対パスが入っている想定 (vina_out/poses/...), REPO_ROOT 由来のケースも考慮
    ALT_PDBQT="${REPO_ROOT}/${pose_path}"
    if [[ -f "${ALT_PDBQT}" ]]; then
      SRC_PDBQT="${ALT_PDBQT}"
    else
      echo "[WARN] pose file not found for zinc_id=${zinc_id}: ${SRC_PDBQT} or ${ALT_PDBQT}" >&2
      continue
    fi
  fi
  DEST_PDBQT="${GMX_INPUT_DIR}/${zinc_id}_pose.pdbqt"
  echo "[INFO] copy pose: ${SRC_PDBQT} -> ${DEST_PDBQT}"
  cp "${SRC_PDBQT}" "${DEST_PDBQT}"
done

echo "[INFO] (NOTE) ここで obabel により PDBQT -> PDB 変換、および acpype によるトポロジー生成を行う想定です。"
echo "[INFO] vina2gromacs-prep finished."
