#!/usr/bin/env bash
set -euo pipefail

# PDBファイルをPDBQTに変換（簡易版）
# 実運用ではMGLToolsのprepare_receptorを使用すべき

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJ_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

TARGET="${1:-aps_apoh}"
PDB_FILE="$PROJ_ROOT/catalog/targets/$TARGET/receptor.pdb"
PDBQT_FILE="$PROJ_ROOT/catalog/targets/$TARGET/receptor.pdbqt"

if [ ! -f "$PDB_FILE" ]; then
  echo "[ERROR] PDB file not found: $PDB_FILE"
  exit 1
fi

echo "[INFO] Converting $PDB_FILE to PDBQT..."
echo "[WARN] This is a simplified conversion. For production, use MGLTools prepare_receptor4.py"

# PDBQT形式に変換（AutoDock Vina形式）
grep -E "^(ATOM|HETATM)" "$PDB_FILE" | \
  awk '{
    # PDB形式からフィールドを抽出
    record = substr($0, 1, 6)
    atomNum = substr($0, 7, 5)
    atomName = substr($0, 13, 4)
    resName = substr($0, 18, 3)
    chain = substr($0, 22, 1)
    resNum = substr($0, 23, 4)
    x = substr($0, 31, 8)
    y = substr($0, 39, 8)
    z = substr($0, 47, 8)
    occ = substr($0, 55, 6)
    bfac = substr($0, 61, 6)
    element = substr($0, 77, 2)
    
    # Atom名とelementをトリム
    gsub(/^ +| +$/, "", atomName)
    gsub(/^ +| +$/, "", element)
    
    if (element == "") {
      # Elementがない場合、atom名から推測
      firstChar = substr(atomName, 1, 1)
      if (firstChar ~ /[0-9]/) {
        element = substr(atomName, 2, 1)
      } else {
        element = firstChar
      }
    }
    
    # PDBQT形式で出力（charge, atom_type）
    printf("%-6s%5s %-4s %3s %1s%4s    %8s%8s%8s%6s%6s    %+6.3f %1s \n",
           record, atomNum, atomName, resName, chain, resNum,
           x, y, z, occ, bfac, 0.000, element)
  }' > "$PDBQT_FILE"

echo "[INFO] Created: $PDBQT_FILE"
echo "[WARN] This is NOT production-ready. Install MGLTools for proper conversion:"
echo "  prepare_receptor4 -r receptor.pdb -o receptor.pdbqt"
