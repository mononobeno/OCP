#!/usr/bin/env bash
set -euo pipefail

echo "[INFO] vina-prep starting..."

WS="${WORKSPACE:-/workspace}"
TARGET_CODE="${TARGET_CODE:-aps_apoh}"
LIB_CODE="${LIBRARY_CODE:-zinc_2d_smi_v1}"

ZINC_ROOT="/zinc-library"
TARGET_DIR="$ZINC_ROOT/targets/$TARGET_CODE"
LIB_DIR="$ZINC_ROOT/libraries/$LIB_CODE"

DB="${DB_PATH:-$ZINC_ROOT/db/ocp_results.sqlite}"

echo "[INFO] WORKSPACE   = $WS"
echo "[INFO] TARGET_CODE = $TARGET_CODE"
echo "[INFO] LIB_CODE    = $LIB_CODE"
echo "[INFO] TARGET_DIR  = $TARGET_DIR"
echo "[INFO] LIB_DIR     = $LIB_DIR"
echo "[INFO] DB          = $DB"

LIG_RAW_DIR="$WS/ligands_raw"
LIG_IDS_FILE="$LIG_RAW_DIR/ligands_zinc_ids.txt"
RUN_INFO="$LIG_RAW_DIR/run_info.json"

OUT_DIR="$WS/vina_input"
OUT_LIG_DIR="$OUT_DIR/ligands"

mkdir -p "$OUT_DIR" "$OUT_LIG_DIR"

# --- receptor 処理 ---
REC_PDB_SRC="$TARGET_DIR/receptor.pdb"
REC_PDB_DST="$OUT_DIR/receptor.pdb"
REC_PDBQT_DST="$OUT_DIR/receptor.pdbqt"

if [ ! -f "$REC_PDB_SRC" ]; then
  echo "[ERROR] receptor.pdb not found: $REC_PDB_SRC" >&2
  exit 1
fi

echo "[INFO] copying receptor: $REC_PDB_SRC -> $REC_PDB_DST"
cp "$REC_PDB_SRC" "$REC_PDB_DST"

echo "[INFO] converting receptor to pdbqt..."
# 必要に応じてオプションは後で調整
obabel "$REC_PDB_DST" -O "$REC_PDBQT_DST" >/dev/null 2>&1 || {
  echo "[WARN] obabel receptor conversion failed, receptor.pdbqt may be missing" >&2
}

# --- ligands 処理 ---
if [ ! -f "$LIG_IDS_FILE" ]; then
  echo "[ERROR] ligand IDs file not found: $LIG_IDS_FILE" >&2
  exit 1
fi

echo "[INFO] reading ligand IDs from: $LIG_IDS_FILE"

COUNT_TOTAL=0
COUNT_OK=0
COUNT_MISS=0

while IFS= read -r ZINC_ID; do
  ZINC_ID_TRIM="$(echo "$ZINC_ID" | tr -d '[:space:]')"
  [ -z "$ZINC_ID_TRIM" ] && continue
  COUNT_TOTAL=$((COUNT_TOTAL+1))

  SRC_PDBQT="$LIB_DIR/processed/pdbqt/${ZINC_ID_TRIM}.pdbqt"
  DST_PDBQT="$OUT_LIG_DIR/${ZINC_ID_TRIM}.pdbqt"

  if [ ! -f "$SRC_PDBQT" ]; then
    echo "[WARN] missing preprocessed ligand pdbqt: $SRC_PDBQT"
    COUNT_MISS=$((COUNT_MISS+1))
    continue
  fi

  mkdir -p "$(dirname "$DST_PDBQT")"
  cp "$SRC_PDBQT" "$DST_PDBQT"
  COUNT_OK=$((COUNT_OK+1))
done < "$LIG_IDS_FILE"

echo "[INFO] total ZINC IDs: $COUNT_TOTAL"
echo "[INFO] copied pdbqt:   $COUNT_OK"
echo "[INFO] missing pdbqt:  $COUNT_MISS"

# --- run_info.json に vina_input 情報を追記（あれば） ---
if [ -f "$RUN_INFO" ]; then
  echo "[INFO] updating run_info.json with vina_input info..."
  # jq 無し前提なので、簡易的に別ファイルとして書き足しておく
  cat >> "$RUN_INFO" <<EOF2

# vina-prep:
#   receptor_pdb:    vina_input/receptor.pdb
#   receptor_pdbqt:  vina_input/receptor.pdbqt
#   ligands_dir:     vina_input/ligands/
EOF2
fi

echo "[INFO] vina-prep finished."
