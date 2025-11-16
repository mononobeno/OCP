#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DB="$ROOT_DIR/catalog/db/ocp_results.sqlite"
LIB_CODE="${1:-zinc_2d_smi_v1}"
RUN_ID="${2:-}"

LIB_DIR="$ROOT_DIR/catalog/libraries/$LIB_CODE"
OUT_DIR="$LIB_DIR/processed/pdbqt"
SMILES_FILE="$LIB_DIR/raw/library_prepared/smiles_all_01.txt"

echo "[INFO] ROOT_DIR    = $ROOT_DIR"
echo "[INFO] DB          = $DB"
echo "[INFO] LIB_CODE    = $LIB_CODE"
echo "[INFO] RUN_ID      = ${RUN_ID:-<none>}"
echo "[INFO] OUT_DIR     = $OUT_DIR"
echo "[INFO] SMILES_FILE = $SMILES_FILE"

if [ ! -f "$SMILES_FILE" ]; then
  echo "[ERROR] SMILES ファイルが見つかりません: $SMILES_FILE" >&2
  exit 1
fi

if ! command -v obabel >/dev/null 2>&1; then
  echo "[ERROR] obabel が見つかりません。" >&2
  exit 1
fi

LIB_ID="$(sqlite3 "$DB" "SELECT id FROM libraries WHERE code='$LIB_CODE';" || true)"
if [ -z "$LIB_ID" ]; then
  echo "[ERROR] libraries.code='$LIB_CODE' が DB にありません。" >&2
  exit 1
fi

mkdir -p "$OUT_DIR"

# ---- 対象 ZINC ID リストを DB から取得 ----
TMP_LIST="$(mktemp)"

if [ -n "$RUN_ID" ]; then
  echo "[INFO] DB から run_id=$RUN_ID に紐づく ligands を取得"
  sqlite3 "$DB" << EOSQL > "$TMP_LIST"
SELECT DISTINCT l.zinc_id
FROM run_ligands rl
JOIN ligands l ON l.id = rl.ligand_id
WHERE rl.run_id = $RUN_ID
  AND l.library_id = $LIB_ID
  AND l.zinc_id IS NOT NULL;
EOSQL
else
  echo "[INFO] RUN_ID が指定されていないので、library 全体から最新 100 件取得"
  sqlite3 "$DB" << EOSQL > "$TMP_LIST"
SELECT zinc_id
FROM ligands
WHERE library_id = $LIB_ID
  AND zinc_id IS NOT NULL
ORDER BY id DESC
LIMIT 100;
EOSQL
fi

TOTAL="$(wc -l < "$TMP_LIST" || echo 0)"
echo "[INFO] targets to gen3d = $TOTAL"

if [ "$TOTAL" -eq 0 ]; then
  echo "[WARN] 対象リガンドが 0 件です。"
  rm -f "$TMP_LIST"
  exit 0
fi

# ---- 各リガンドに対して 3D 生成 ----
i=0
while read -r ZINC_ID; do
  ZINC_ID_TRIM="$(echo "$ZINC_ID" | tr -d '[:space:]')"
  [ -z "$ZINC_ID_TRIM" ] && continue

  # ZINC_ID から _1 などのサフィックスを除去してベース名取得
  ZINC_BASE="${ZINC_ID_TRIM%%_*}"
  
  # SMILESファイルから該当行をgrep（高速）
  SMILES_LINE="$(grep -m1 "^[^[:space:]]*[[:space:]]${ZINC_BASE}_" "$SMILES_FILE" || true)"
  
  if [ -z "$SMILES_LINE" ]; then
    echo "[WARN] [$i/$TOTAL] SMILES not found for $ZINC_ID_TRIM, skipping"
    continue
  fi
  
  # SMILES部分を抽出（最初のフィールド）
  SMILES="$(echo "$SMILES_LINE" | awk '{print $1}')"
  
  if [ -z "$SMILES" ]; then
    echo "[WARN] [$i/$TOTAL] SMILES empty for $ZINC_ID_TRIM, skipping"
    continue
  fi

  i=$((i+1))
  OUT_PDBQT="$OUT_DIR/${ZINC_ID_TRIM}.pdbqt"

  echo "[INFO] [$i/$TOTAL] $ZINC_ID_TRIM -> 3D PDBQT"

  # 一時 SMILES ファイル
  TMP_SMI="$(mktemp /tmp/ocp_gen3d_XXXXXX.smi)"
  echo "$SMILES" > "$TMP_SMI"

  # ここで --gen3d を ON（1分子ごとに3D化）
  if obabel -ismi "$TMP_SMI" -O "$OUT_PDBQT" --gen3d 2>/dev/null; then
    echo "  [OK] $OUT_PDBQT を生成"
  else
    echo "  [WARN] obabel --gen3d 失敗: $ZINC_ID_TRIM"
    rm -f "$OUT_PDBQT"
  fi

  rm -f "$TMP_SMI"
done < "$TMP_LIST"

rm -f "$TMP_LIST"

echo "[INFO] gen3d_for_selected_ligands.sh finished."
