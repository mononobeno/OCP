#!/usr/bin/env bash
set -euo pipefail

# ZINC SMILES を並列処理で pdbqt に変換
# GNU Parallel または xargs を使用して高速化

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DB="$ROOT_DIR/catalog/db/ocp_results.sqlite"

if ! command -v obabel >/dev/null 2>&1; then
  echo "ERROR: obabel が見つかりません。" >&2
  echo "conda か apt で Open Babel をインストールしてください。" >&2
  exit 1
fi

if ! command -v sqlite3 >/dev/null 2>&1; then
  echo "ERROR: sqlite3 が見つかりません。" >&2
  exit 1
fi

if [ $# -lt 1 ]; then
  echo "Usage: $0 INPUT_SMI[.gz] [LIBRARY_CODE] [MAX_COUNT] [JOBS]" >&2
  echo "  INPUT_SMI    : ZINC SMILES .smi または .smi.gz" >&2
  echo "  LIBRARY_CODE : DB上のライブラリコード (default: zinc_2d_smi_v1)" >&2
  echo "  MAX_COUNT    : 上限件数 (default: 50, 0なら全件)" >&2
  echo "  JOBS         : 並列ジョブ数 (default: CPUコア数)" >&2
  exit 1
fi

INPUT="$1"
LIB_CODE="${2:-zinc_2d_smi_v1}"
MAX_COUNT="${3:-50}"
JOBS="${4:-$(nproc)}"

LIB_DIR="$ROOT_DIR/catalog/libraries/$LIB_CODE"
OUT_DIR="$LIB_DIR/processed/pdbqt"
TMP_DIR="$(mktemp -d /tmp/ocp_parallel_XXXXXX)"

mkdir -p "$OUT_DIR"

LIB_ID="$(sqlite3 "$DB" "SELECT id FROM libraries WHERE code = '$LIB_CODE';")"
if [ -z "$LIB_ID" ]; then
  echo "ERROR: libraries.code='$LIB_CODE' が DB にありません。" >&2
  exit 1
fi

echo "[INFO] ROOT_DIR   = $ROOT_DIR"
echo "[INFO] INPUT      = $INPUT"
echo "[INFO] LIB_CODE   = $LIB_CODE (id=$LIB_ID)"
echo "[INFO] OUT_DIR    = $OUT_DIR"
echo "[INFO] MAX_COUNT  = $MAX_COUNT"
echo "[INFO] JOBS       = $JOBS (並列処理)"

# 入力をストリームとして読む
if [[ "$INPUT" == *.gz ]]; then
  READ_CMD=(zcat "$INPUT")
else
  READ_CMD=(cat "$INPUT")
fi

# 処理する行を一時ファイルに保存
echo "[INFO] Extracting SMILES data..."
WORK_FILE="$TMP_DIR/smiles_to_process.tsv"
if [ "$MAX_COUNT" -ne 0 ]; then
  "${READ_CMD[@]}" | head -n "$MAX_COUNT" | awk '{print $1 "\t" $2}' > "$WORK_FILE"
else
  "${READ_CMD[@]}" | awk '{print $1 "\t" $2}' > "$WORK_FILE"
fi

TOTAL_LINES=$(wc -l < "$WORK_FILE")
echo "[INFO] Total entries to process: $TOTAL_LINES"
echo "[DEBUG] WORK_FILE: $WORK_FILE"

# 処理関数を定義
process_smiles() {
  local SMILES="$1"
  local ZINC_ID="$2"
  local OUT_DIR="$3"
  local DB="$4"
  local LIB_ID="$5"
  
  SMILES_TRIM="$(echo "$SMILES" | tr -d '[:space:]')"
  ZINC_ID_TRIM="$(echo "$ZINC_ID" | tr -d '[:space:]')"
  
  [ -z "$SMILES_TRIM" ] && return 0
  [ -z "$ZINC_ID_TRIM" ] && return 0
  
  OUT_PDBQT="$OUT_DIR/${ZINC_ID_TRIM}.pdbqt"
  
  # 既存ファイルをスキップ
  if [ -f "$OUT_PDBQT" ]; then
    sqlite3 "$DB" "INSERT OR IGNORE INTO ligands (zinc_id, library_id) VALUES ('$ZINC_ID_TRIM', $LIB_ID);" 2>/dev/null || true
    echo "[SKIP] $ZINC_ID_TRIM"
    return 0
  fi
  
  # 一時SMIファイルを作成
  TMP_SMI="$(mktemp /tmp/ocp_smi_XXXXXX.smi)"
  echo "$SMILES_TRIM" > "$TMP_SMI"
  
  # obabelで変換
  if obabel -ismi "$TMP_SMI" -O "$OUT_PDBQT" --gen3d >/dev/null 2>&1; then
    sqlite3 "$DB" "INSERT OR IGNORE INTO ligands (zinc_id, library_id) VALUES ('$ZINC_ID_TRIM', $LIB_ID);" 2>/dev/null || true
    echo "[OK] $ZINC_ID_TRIM"
    rm -f "$TMP_SMI"
    return 0
  else
    echo "[FAIL] $ZINC_ID_TRIM" >&2
    rm -f "$TMP_SMI"
    return 1
  fi
}

export -f process_smiles
export OUT_DIR DB LIB_ID

# GNU Parallel が利用可能かチェック
if command -v parallel >/dev/null 2>&1; then
  echo "[INFO] Using GNU Parallel with $JOBS jobs"
  cat "$WORK_FILE" | parallel --colsep '\t' -j "$JOBS" --line-buffer \
    process_smiles {1} {2} \"\$OUT_DIR\" \"\$DB\" \"\$LIB_ID\"
else
  echo "[INFO] Using xargs with $JOBS jobs (GNU Parallel not found)"
  # xargs での並列処理（少し遅い）
  cat "$WORK_FILE" | while IFS=$'\t' read -r SMILES ZINC_ID REST; do
    echo "$SMILES	$ZINC_ID"
  done | xargs -P "$JOBS" -I {} bash -c 'process_smiles $(echo {} | cut -f1) $(echo {} | cut -f2) "$OUT_DIR" "$DB" "$LIB_ID"'
fi

# クリーンアップ
rm -rf "$TMP_DIR"

# 結果を集計
SUCCESS=$(ls "$OUT_DIR"/*.pdbqt 2>/dev/null | wc -l)
echo
echo "[INFO] ================================"
echo "[INFO] Processing completed"
echo "[INFO] Total processed: $SUCCESS pdbqt files"
echo "[INFO] ================================"
