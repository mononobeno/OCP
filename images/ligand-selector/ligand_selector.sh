#!/usr/bin/env bash
set -euo pipefail

echo "[INFO] ligand-selector starting..."

DB="${DB_PATH:-/zinc-library/db/ocp_results.sqlite}"
WS="${WORKSPACE:-/workspace}"
TARGET_CODE="${TARGET_CODE:-aps_apoh}"
LIB_CODE="${LIBRARY_CODE:-zinc_2d_smi_v1}"
NUM_LIGANDS="${NUM_LIGANDS:-50}"

echo "[INFO] DB          = $DB"
echo "[INFO] WORKSPACE   = $WS"
echo "[INFO] TARGET_CODE = $TARGET_CODE"
echo "[INFO] LIB_CODE    = $LIB_CODE"
echo "[INFO] NUM_LIGANDS = $NUM_LIGANDS"

if [ ! -f "$DB" ]; then
  echo "[ERROR] DB not found: $DB" >&2
  exit 1
fi

mkdir -p "$WS/ligands_raw"

TARGET_ID="$(sqlite3 "$DB" "SELECT id FROM targets   WHERE code = '$TARGET_CODE';")"
LIB_ID="$(sqlite3 "$DB" "SELECT id FROM libraries WHERE code = '$LIB_CODE';")"

if [ -z "$TARGET_ID" ]; then
  echo "[ERROR] targets.code='$TARGET_CODE' not found in DB" >&2
  exit 1
fi
if [ -z "$LIB_ID" ]; then
  echo "[ERROR] libraries.code='$LIB_CODE' not found in DB" >&2
  exit 1
fi

echo "[INFO] TARGET_ID = $TARGET_ID"
echo "[INFO] LIB_ID    = $LIB_ID"

# 新しい run を登録
RUN_UUID="$(cat /proc/sys/kernel/random/uuid)"
NOW_UTC="$(date -u +"%Y-%m-%dT%H:%M:%SZ")"

sqlite3 "$DB" <<EOSQL
INSERT INTO runs (run_uuid, started_at, target_id, library_id, notes)
VALUES ('$RUN_UUID', '$NOW_UTC', $TARGET_ID, $LIB_ID, 'ligand-selector');
EOSQL

RUN_ID="$(sqlite3 "$DB" "SELECT id FROM runs WHERE run_uuid = '$RUN_UUID';")"

echo "[INFO] RUN_ID    = $RUN_ID"
echo "[INFO] RUN_UUID  = $RUN_UUID"

OUT_TSV="$WS/ligands_raw/ligands_selected.tsv"
OUT_IDS="$WS/ligands_raw/ligands_zinc_ids.txt"
META_JSON="$WS/ligands_raw/run_info.json"

echo "[INFO] selecting unused ligands ..."
sqlite3 -separator $'\t' "$DB" "
SELECT l.id, l.zinc_id
FROM ligands l
WHERE l.library_id = $LIB_ID
  AND NOT EXISTS (
    SELECT 1
    FROM runs r
    JOIN vina_results vr ON vr.run_id = r.id
    WHERE r.target_id  = $TARGET_ID
      AND r.library_id = $LIB_ID
      AND vr.ligand_id = l.id
  )
ORDER BY RANDOM()
LIMIT $NUM_LIGANDS;
" > "$OUT_TSV"

cut -f2 "$OUT_TSV" > "$OUT_IDS" || true

cat > "$META_JSON" <<EOF_JSON
{
  "run_id": $RUN_ID,
  "run_uuid": "$RUN_UUID",
  "target_code": "$TARGET_CODE",
  "library_code": "$LIB_CODE",
  "num_ligands": $NUM_LIGANDS,
  "generated_at_utc": "$NOW_UTC",
  "ligands_tsv": "ligands_raw/ligands_selected.tsv",
  "ligands_ids": "ligands_raw/ligands_zinc_ids.txt"
}
EOF_JSON

echo "[INFO] selection done."
echo "[INFO]  - $OUT_TSV"
echo "[INFO]  - $OUT_IDS"
echo "[INFO]  - $META_JSON"

LINES="$(wc -l < "$OUT_TSV" || echo 0)"
echo "[INFO] selected ligands: $LINES"

if [ "$LINES" -eq 0 ]; then
  echo "[WARN] no unused ligands were found. Job will still succeed, but downstream stages may have nothing to do."
fi

echo "[INFO] ligand-selector finished."
