#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# SMILES -> PDBQT Parallel Processor (v2)
# GNU Parallel を使って高速にSMILESをPDBQTに変換し、DBに登録
###############################################################################

if [ "$#" -lt 2 ]; then
  echo "Usage: $0 <smiles_file> <library_code> [max_count] [jobs]"
  echo "  smiles_file  : SMILES file (space-separated: SMILES ZINC_ID)"
  echo "  library_code : library code (e.g. zinc_2d_smi_v1)"
  echo "  max_count    : 処理する最大行数 (default: 0=全件)"
  echo "  jobs         : 並列ジョブ数 (default: CPU cores)"
  exit 1
fi

INPUT="$1"
LIB_CODE="$2"
MAX_COUNT="${3:-0}"
JOBS="${4:-$(nproc)}"

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
DB="$ROOT_DIR/catalog/db/ocp_results.sqlite"

# ライブラリIDを取得
LIB_ID=$(sqlite3 "$DB" "SELECT id FROM libraries WHERE code='$LIB_CODE';" 2>/dev/null || echo "")
if [ -z "$LIB_ID" ]; then
  echo "[ERROR] Library '$LIB_CODE' not found in DB"
  exit 1
fi

OUT_DIR="$ROOT_DIR/catalog/libraries/$LIB_CODE/processed/pdbqt"
mkdir -p "$OUT_DIR"

echo "[INFO] ROOT_DIR   = $ROOT_DIR"
echo "[INFO] INPUT      = $INPUT"
echo "[INFO] LIB_CODE   = $LIB_CODE (id=$LIB_ID)"
echo "[INFO] OUT_DIR    = $OUT_DIR"
echo "[INFO] MAX_COUNT  = $MAX_COUNT"
echo "[INFO] JOBS       = $JOBS (並列処理)"

# 処理する行を一時ファイルに保存
echo "[INFO] Extracting SMILES data..."
TMP_DIR=$(mktemp -d -t ocp_parallel_XXXXXX)
WORK_FILE="$TMP_DIR/smiles_to_process.tsv"

if [ "$MAX_COUNT" -ne 0 ]; then
  cat "$INPUT" | head -n "$MAX_COUNT" | awk '{print $1 "\t" $2}' > "$WORK_FILE"
else
  cat "$INPUT" | awk '{print $1 "\t" $2}' > "$WORK_FILE"
fi

TOTAL_LINES=$(wc -l < "$WORK_FILE")
echo "[INFO] Total entries to process: $TOTAL_LINES"

# GNU Parallelで並列処理（関数をインラインで実行）
echo "[INFO] Using GNU Parallel with $JOBS jobs"

cat "$WORK_FILE" | parallel --colsep '\t' -j "$JOBS" --line-buffer "
  SMILES={1}
  ZINC_ID={2}
  OUT_DIR='$OUT_DIR'
  DB='$DB'
  LIB_ID=$LIB_ID
  
  ZINC_ID_TRIM=\$(echo \"\$ZINC_ID\" | sed 's/_[0-9]*$//')
  OUT_PDBQT=\"\$OUT_DIR/\${ZINC_ID}.pdbqt\"
  
  if [ -f \"\$OUT_PDBQT\" ]; then
    echo \"[SKIP] \$ZINC_ID_TRIM\"
    exit 0
  fi
  
  TMP_SMI=\$(mktemp --suffix=.smi)
  echo \"\$SMILES\" > \"\$TMP_SMI\"
  
  if obabel -ismi \"\$TMP_SMI\" -O \"\$OUT_PDBQT\" >/dev/null 2>&1; then
    sqlite3 \"\$DB\" \"INSERT OR IGNORE INTO ligands (zinc_id, library_id) VALUES ('\$ZINC_ID_TRIM', \$LIB_ID);\" 2>/dev/null || true
    echo \"[OK] \$ZINC_ID_TRIM\"
    rm -f \"\$TMP_SMI\"
    exit 0
  else
    echo \"[FAIL] \$ZINC_ID_TRIM\" >&2
    rm -f \"\$TMP_SMI\"
    exit 1
  fi
"

# クリーンアップ
rm -rf "$TMP_DIR"

# 結果を確認
PDBQT_COUNT=$(ls -1 "$OUT_DIR"/*.pdbqt 2>/dev/null | wc -l)
echo "[INFO] Processing completed"
echo "[INFO] Total processed: $PDBQT_COUNT pdbqt files"
