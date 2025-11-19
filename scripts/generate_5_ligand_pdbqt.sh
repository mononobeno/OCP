#!/bin/bash
set -e

DB="/home/dev/OCP/catalog/db/ocp_results.sqlite"
LIGAND_IDS="1684,1689,1690,1691,1692"

echo "Generating PDBQT for 5 ligands..."
mkdir -p /home/dev/OCP/catalog/libraries/zinc20_ml/pdbqt

LIGANDS=$(sqlite3 "$DB" "SELECT zinc_id, smiles FROM ligands WHERE id IN ($LIGAND_IDS);")

COUNT=0
while IFS='|' read -r ZINC_ID SMILES; do
  COUNT=$((COUNT + 1))
  echo "[$COUNT/5] $ZINC_ID"
  
  docker run --rm \
    -v /home/dev/OCP:/work \
    -w /work \
    pipeline/ligand-selector:local \
    python3 tools/rdkit_gen3d_to_pdbqt.py \
      --smiles "$SMILES" \
      --output "catalog/libraries/zinc20_ml/pdbqt/${ZINC_ID}.pdbqt" \
      --name "$ZINC_ID"
  
  echo "  ✅ ${ZINC_ID}.pdbqt"
done <<< "$LIGANDS"

echo "✅ All 5 ligand PDBQT files generated"
ls -lh /home/dev/OCP/catalog/libraries/zinc20_ml/pdbqt/
