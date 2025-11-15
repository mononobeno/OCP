#!/usr/bin/env bash
set -euo pipefail

# ZINC SMILES (.smi or .smi.gz) を読み込み、
# 1列目 SMILES / 2列目 ZINC ID を前提に pdbqt を生成して
# catalog/libraries/<LIB_CODE>/processed/pdbqt/ に出力する。
#
# ついでに ligands テーブルにも登録する。

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
  echo "Usage: $0 INPUT_SMI[.gz] [LIBRARY_CODE] [MAX_COUNT]" >&2
  echo "  INPUT_SMI : ZINC-downloader の .smi または .smi.gz" >&2
  echo "  LIBRARY_CODE : DB上のライブラリコード (default: zinc_2d_smi_v1)" >&2
  echo "  MAX_COUNT : 上限件数 (default: 50, 0なら全件)" >&2
  exit 1
fi

INPUT="$1"
LIB_CODE="${2:-zinc_2d_smi_v1}"
MAX_COUNT="${3:-50}"

LIB_DIR="$ROOT_DIR/catalog/libraries/$LIB_CODE"
OUT_DIR="$LIB_DIR/processed/pdbqt"

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

# 入力をストリームとして読む（.gz なら zcat、普通なら cat）
if [[ "$INPUT" == *.gz ]]; then
  READ_CMD=(zcat "$INPUT")
else
  READ_CMD=(cat "$INPUT")
fi

COUNT=0
SUCCESS=0
SKIP_EXIST=0
FAILED=0

TMP_SMI="$(mktemp /tmp/ocp_smi_XXXXXX.smi)"

"${READ_CMD[@]}" | while IFS=$'\t ' read -r SMILES ZINC_ID REST; do
  SMILES_TRIM="$(echo "$SMILES" | tr -d '[:space:]')"
  ZINC_ID_TRIM="$(echo "$ZINC_ID" | tr -d '[:space:]')"

  [ -z "$SMILES_TRIM" ] && continue
  [ -z "$ZINC_ID_TRIM" ] && continue

  COUNT=$((COUNT+1))

  if [ "$MAX_COUNT" -ne 0 ] && [ "$COUNT" -gt "$MAX_COUNT" ]; then
    break
  fi

  OUT_PDBQT="$OUT_DIR/${ZINC_ID_TRIM}.pdbqt"

  if [ -f "$OUT_PDBQT" ]; then
    echo "[SKIP] exists: $OUT_PDBQT"
    SKIP_EXIST=$((SKIP_EXIST+1))
    # 既に ligands テーブルに登録済みかどうかは問わず、一応 INSERT OR IGNORE しておく
    sqlite3 "$DB" << EOSQL
INSERT OR IGNORE INTO ligands (zinc_id, library_id)
VALUES ('$ZINC_ID_TRIM', $LIB_ID);
EOSQL
    continue
  fi

  echo "$SMILES_TRIM" > "$TMP_SMI"

  echo "[INFO] [$COUNT] $ZINC_ID_TRIM -> $OUT_PDBQT"
  if obabel -ismi "$TMP_SMI" -O "$OUT_PDBQT" --gen3d >/dev/null 2>&1; then
    SUCCESS=$((SUCCESS+1))
    sqlite3 "$DB" << EOSQL
INSERT OR IGNORE INTO ligands (zinc_id, library_id)
VALUES ('$ZINC_ID_TRIM', $LIB_ID);
EOSQL
  else
    echo "[WARN] obabel failed for $ZINC_ID_TRIM"
    FAILED=$((FAILED+1))
    rm -f "$OUT_PDBQT"
  fi

done

rm -f "$TMP_SMI" || true

echo "[INFO] processed lines  = $COUNT"
echo "[INFO] pdbqt generated  = $SUCCESS"
echo "[INFO] skipped (exists) = $SKIP_EXIST"
echo "[INFO] failed           = $FAILED"
