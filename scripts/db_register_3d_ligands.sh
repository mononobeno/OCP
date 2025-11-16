#!/usr/bin/env bash
set -euo pipefail

# 既存の3D座標付きPDBQTファイルをDBに登録して has_3d=1 にする

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DB="$ROOT_DIR/catalog/db/ocp_results.sqlite"
LIB_CODE="${1:-zinc_2d_smi_v1}"

echo "[INFO] ROOT_DIR = $ROOT_DIR"
echo "[INFO] DB       = $DB"
echo "[INFO] LIB_CODE = $LIB_CODE"

if [ ! -f "$DB" ]; then
  echo "[ERROR] DB not found: $DB" >&2
  exit 1
fi

# library_id取得
LIB_ID="$(sqlite3 "$DB" "SELECT id FROM libraries WHERE code='$LIB_CODE';" || true)"
if [ -z "$LIB_ID" ]; then
  echo "[ERROR] library '$LIB_CODE' not found in DB" >&2
  exit 1
fi

echo "[INFO] LIB_ID = $LIB_ID"

PDBQT_DIR="$ROOT_DIR/catalog/libraries/$LIB_CODE/processed/pdbqt"

if [ ! -d "$PDBQT_DIR" ]; then
  echo "[ERROR] PDBQT dir not found: $PDBQT_DIR" >&2
  exit 1
fi

echo "[INFO] PDBQT_DIR = $PDBQT_DIR"

# 有効なPDBQTファイル（サイズ>1KB）を検索
VALID_PDBQT_FILES="$(find "$PDBQT_DIR" -name "*.pdbqt" -size +1k)"

TOTAL="$(echo "$VALID_PDBQT_FILES" | wc -l || echo 0)"
echo "[INFO] Found $TOTAL valid PDBQT files (>1KB)"

if [ "$TOTAL" -eq 0 ]; then
  echo "[WARN] No valid PDBQT files found"
  exit 0
fi

# 各ファイルに対してhas_3dフラグを立てる
i=0
echo "$VALID_PDBQT_FILES" | while read -r pdbqt_file; do
  i=$((i+1))
  
  BASENAME="$(basename "$pdbqt_file" .pdbqt)"
  
  # ZINC000000000001_1 -> ZINC000000000001 にベース名取得
  ZINC_BASE="${BASENAME%%_*}"
  
  # DBに該当レコードがあるか確認
  EXISTS="$(sqlite3 "$DB" "SELECT COUNT(*) FROM ligands WHERE zinc_id='$BASENAME' AND library_id=$LIB_ID;" || echo "0")"
  
  if [ "$EXISTS" -gt 0 ]; then
    # has_3d=1 に更新
    sqlite3 "$DB" << EOSQL
UPDATE ligands 
SET has_3d = 1, 
    conformer_method = 'obabel_gen3d'
WHERE zinc_id = '$BASENAME' 
  AND library_id = $LIB_ID;
EOSQL
    
    if [ $((i % 50)) -eq 0 ]; then
      echo "[INFO] [$i/$TOTAL] Updated: $BASENAME"
    fi
  else
    # レコードが無い場合はスキップ（または挿入）
    if [ $((i % 100)) -eq 0 ]; then
      echo "[WARN] [$i/$TOTAL] Not in DB: $BASENAME"
    fi
  fi
done

echo ""
echo "[INFO] Update complete!"
echo "[INFO] Check result:"
echo "  sqlite3 $DB \"SELECT COUNT(*) FROM ligands WHERE has_3d=1;\""
echo "  sqlite3 $DB \"SELECT zinc_id, has_3d, conformer_method FROM ligands WHERE has_3d=1 LIMIT 10;\""
