#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DB="$ROOT_DIR/catalog/db/ocp_results.sqlite"

if ! command -v sqlite3 >/dev/null 2>&1; then
  echo "ERROR: sqlite3 が見つかりません。" >&2
  exit 1
fi

if [ $# -lt 1 ]; then
  echo "Usage: $0 ZINC_ID_LIST_TSV [LIBRARY_CODE]" >&2
  echo "  ZINC_ID_LIST_TSV: 1列目に ZINC ID が入っているテキスト（TSV/改行区切り）" >&2
  echo "  LIBRARY_CODE: 省略時は zinc_2d_smi_v1" >&2
  exit 1
fi

LIST_FILE="$1"
LIB_CODE="${2:-zinc_2d_smi_v1}"

if [ ! -f "$LIST_FILE" ]; then
  echo "ERROR: file not found: $LIST_FILE" >&2
  exit 1
fi

LIB_ID="$(sqlite3 "$DB" "SELECT id FROM libraries WHERE code = '$LIB_CODE';")"
if [ -z "$LIB_ID" ]; then
  echo "ERROR: libraries.code='$LIB_CODE' が DB にありません。" >&2
  exit 1
fi

echo "[INFO] DB            = $DB"
echo "[INFO] LIB_CODE      = $LIB_CODE (id=$LIB_ID)"
echo "[INFO] ZINC ID list  = $LIST_FILE"

COUNT_BEFORE="$(sqlite3 "$DB" "SELECT COUNT(*) FROM ligands WHERE library_id=$LIB_ID;")"
echo "[INFO] ligands (before) = $COUNT_BEFORE"

ADDED=0
while IFS=$'\t' read -r ZINC_ID REST; do
  ZINC_ID_TRIM="$(echo "$ZINC_ID" | tr -d '[:space:]')"
  [ -z "$ZINC_ID_TRIM" ] && continue

  sqlite3 "$DB" << EOSQL
INSERT OR IGNORE INTO ligands (zinc_id, library_id)
VALUES ('$ZINC_ID_TRIM', $LIB_ID);
EOSQL

  ADDED=$((ADDED+1))
done < "$LIST_FILE"

COUNT_AFTER="$(sqlite3 "$DB" "SELECT COUNT(*) FROM ligands WHERE library_id=$LIB_ID;")"

echo "[INFO] processed lines = $ADDED"
echo "[INFO] ligands (after)  = $COUNT_AFTER"
