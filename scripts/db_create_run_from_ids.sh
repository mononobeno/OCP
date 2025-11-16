#!/usr/bin/env bash
# scripts/db_create_run_from_ids.sh
# ZINC ID リストから runs を作成し、
# RUN_ID / RUN_UUID を echo するユーティリティ。
# 実際のスキーマ:
#   libraries(id, code, description)
#   targets(id, code, name, pdb_id, receptor_path)
#   ligands(id, zinc_id, library_id, ...)
#   runs(id, run_uuid, started_at, target_id, library_id, notes)
#   vina_results(run_id, ligand_id, mode_rank, affinity_kcal, ...)
# 
# NOTE: run_ligands テーブルは存在しないため、
#       ZINC ID リストを環境変数またはファイルとして保存し、
#       vina-runner がそれを読み込む方式とします。

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="${SCRIPT_DIR%/scripts}"
DEFAULT_DB_PATH="${REPO_ROOT}/catalog/db/ocp_results.sqlite"

usage() {
cat <<USAGE
Usage: $0 <LIBRARY_CODE> <TARGET_CODE> <ZINC_ID_LIST_FILE>

Example:
  $0 zinc20_ml_v1 apoh catalog/libraries/zinc_2d_smi_v1/processed/ligands_zinc_ids.txt

環境変数:
  DB_PATH_OVERRIDE ... SQLite DB のパスを上書き指定できる
  RUN_NOTES ... runs.notes に入れる任意のメモ
  
NOTE: 現在のスキーマには run_ligands テーブルがないため、
      ZINC ID リストは vina-runner がファイルまたは DB から直接読み込みます。
USAGE
exit 1
}

if [[ $# -ne 3 ]]; then
  usage
fi

LIB_CODE="$1"
TARGET_CODE="$2"
ZINC_LIST_FILE="$3"

DB_PATH="${DB_PATH_OVERRIDE:-$DEFAULT_DB_PATH}"

if [[ ! -f "$DB_PATH" ]]; then
  echo "ERROR: DB not found: $DB_PATH" >&2
  exit 1
fi

if [[ ! -f "$ZINC_LIST_FILE" ]]; then
  echo "ERROR: ZINC ID list file not found: $ZINC_LIST_FILE" >&2
  exit 1
fi

if ! command -v sqlite3 >/dev/null 2>&1; then
  echo "ERROR: sqlite3 command not found." >&2
  exit 1
fi

generate_uuid() {
  if command -v uuidgen >/dev/null 2>&1; then
    uuidgen
  else
    python3 - <<'PY'
import uuid
print(uuid.uuid4())
PY
  fi
}

RUN_UUID="$(generate_uuid)"
RUN_NOTES_DEFAULT="Vina run: ${LIB_CODE} × ${TARGET_CODE} ($(wc -l < "$ZINC_LIST_FILE") ligands)"
RUN_NOTES="${RUN_NOTES:-$RUN_NOTES_DEFAULT}"

# ライブラリ / ターゲット ID を取得
LIB_ID="$(
  sqlite3 "$DB_PATH" "SELECT id FROM libraries WHERE code = '$LIB_CODE' LIMIT 1;" || true
)"
TARGET_ID="$(
  sqlite3 "$DB_PATH" "SELECT id FROM targets WHERE code = '$TARGET_CODE' LIMIT 1;" || true
)"

if [[ -z "$LIB_ID" ]]; then
  echo "ERROR: library not found for code=${LIB_CODE}" >&2
  exit 1
fi
if [[ -z "$TARGET_ID" ]]; then
  echo "ERROR: target not found for code=${TARGET_CODE}" >&2
  exit 1
fi

echo "[INFO] DB : $DB_PATH" >&2
echo "[INFO] LIB_CODE : $LIB_CODE (id=${LIB_ID})" >&2
echo "[INFO] TARGET : $TARGET_CODE (id=${TARGET_ID})" >&2
echo "[INFO] ZINC LIST: $ZINC_LIST_FILE" >&2
echo "[INFO] RUN_NOTES: $RUN_NOTES" >&2

# runs に 1 レコード INSERT し、RUN_ID を取得
RUN_ID="$(
  sqlite3 "$DB_PATH" <<SQL
BEGIN;
INSERT INTO runs (run_uuid, started_at, target_id, library_id, notes)
VALUES (
  '$RUN_UUID',
  datetime('now'),
  $TARGET_ID,
  $LIB_ID,
  '$RUN_NOTES'
);
SELECT last_insert_rowid();
COMMIT;
SQL
)"

if [[ -z "$RUN_ID" ]]; then
  echo "ERROR: failed to insert run record." >&2
  exit 1
fi

echo "[INFO] Created run: id=${RUN_ID}, uuid=${RUN_UUID}" >&2

# ZINC ID の存在確認（オプション）
TOTAL_LINES=0
VALID_LIGANDS=0
SKIP_NO_LIGAND=0

while IFS= read -r ZINC_ID; do
  ZINC_ID="${ZINC_ID%%[[:space:]]*}"
  if [[ -z "$ZINC_ID" ]]; then
    continue
  fi
  if [[ "$ZINC_ID" =~ ^# ]]; then
    continue
  fi
  
  TOTAL_LINES=$((TOTAL_LINES + 1))

  LIG_ID="$(
    sqlite3 "$DB_PATH" "SELECT id FROM ligands WHERE zinc_id = '$ZINC_ID' AND library_id = $LIB_ID LIMIT 1;" || true
  )"

  if [[ -n "$LIG_ID" ]]; then
    VALID_LIGANDS=$((VALID_LIGANDS + 1))
  else
    echo "[WARN] ligand not found in DB: zinc_id=${ZINC_ID}" >&2
    SKIP_NO_LIGAND=$((SKIP_NO_LIGAND + 1))
  fi
done < "$ZINC_LIST_FILE"

echo "[INFO] ZINC ID validation: ${VALID_LIGANDS}/${TOTAL_LINES} found in DB (skipped ${SKIP_NO_LIGAND})." >&2

# シェルから eval しやすいように echo
echo "RUN_ID=${RUN_ID}"
echo "RUN_UUID=${RUN_UUID}"
