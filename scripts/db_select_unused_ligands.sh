#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DB="$ROOT_DIR/catalog/db/ocp_results.sqlite"

if ! command -v sqlite3 >/dev/null 2>&1; then
  echo "ERROR: sqlite3 が見つかりません。" >&2
  exit 1
fi

TARGET_CODE="${1:-aps_apoh}"
LIB_CODE="${2:-zinc_2d_smi_v1}"
LIMIT="${3:-10}"

TARGET_ID="$(sqlite3 "$DB" "SELECT id FROM targets   WHERE code = '$TARGET_CODE';")"
LIB_ID="$(sqlite3 "$DB" "SELECT id FROM libraries WHERE code = '$LIB_CODE';")"

if [ -z "$TARGET_ID" ]; then
  echo "ERROR: targets.code='$TARGET_CODE' が DB にありません。" >&2
  exit 1
fi
if [ -z "$LIB_ID" ]; then
  echo "ERROR: libraries.code='$LIB_CODE' が DB にありません。" >&2
  exit 1
fi

echo "[INFO] DB          = $DB"
echo "[INFO] TARGET_CODE = $TARGET_CODE (id=$TARGET_ID)"
echo "[INFO] LIB_CODE    = $LIB_CODE (id=$LIB_ID)"
echo "[INFO] LIMIT       = $LIMIT"
echo
echo "# Selecting unused ligands with has_3d=1"
echo "# ligand_id    zinc_id"

sqlite3 "$DB" "
SELECT l.id, l.zinc_id
FROM ligands l
WHERE l.library_id = $LIB_ID
  AND l.has_3d = 1
  AND NOT EXISTS (
    SELECT 1
    FROM runs r
    JOIN vina_results vr ON vr.run_id = r.id
    WHERE r.target_id  = $TARGET_ID
      AND r.library_id = $LIB_ID
      AND vr.ligand_id = l.id
  )
ORDER BY RANDOM()
LIMIT $LIMIT;
"
