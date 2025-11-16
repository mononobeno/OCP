#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RAW_DIR="$ROOT_DIR/catalog/libraries/zinc_2d_smi_v1/raw"
INPUT_SMI="${1:-$RAW_DIR/all_zinc20_ml_subset.smi}"
LIB_CODE="${2:-zinc_2d_smi_v1}"
MAX_COUNT="${3:-100}"

echo "[INFO] ROOT_DIR  = $ROOT_DIR"
echo "[INFO] INPUT_SMI = $INPUT_SMI"
echo "[INFO] LIB_CODE  = $LIB_CODE"
echo "[INFO] MAX_COUNT = $MAX_COUNT"

if [ ! -f "$INPUT_SMI" ]; then
  echo "[ERROR] INPUT_SMI not found: $INPUT_SMI" >&2
  exit 1
fi

if ! command -v python >/dev/null 2>&1; then
  echo "[ERROR] python command not found. Activate your RDKit env." >&2
  exit 1
fi

# RDKit で 3D 生成 → PDBQT 生成 → DB 更新
python "$ROOT_DIR/tools/rdkit_gen3d_to_pdbqt.py" "$INPUT_SMI" "$LIB_CODE" "$MAX_COUNT"

# 確認
echo "[INFO] Sample of PDBQT files:"
ls -lh "$ROOT_DIR/catalog/libraries/$LIB_CODE/processed/pdbqt" | head

echo "[INFO] ligands count (lib=$LIB_CODE):"
sqlite3 "$ROOT_DIR/catalog/db/ocp_results.sqlite" \
  "SELECT COUNT(*) FROM ligands WHERE library_id = (SELECT id FROM libraries WHERE code='$LIB_CODE');"

echo "[INFO] ligands with has_3d=1:"
sqlite3 "$ROOT_DIR/catalog/db/ocp_results.sqlite" \
  "SELECT COUNT(*) FROM ligands WHERE library_id = (SELECT id FROM libraries WHERE code='$LIB_CODE') AND has_3d=1;"
