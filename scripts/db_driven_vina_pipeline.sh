#!/bin/bash
set -e

# DB-driven Vina Docking Pipeline
# Reads receptor and ligand info from DB, runs Vina automatically

DB="/home/dev/OCP/catalog/db/ocp_results.sqlite"
RECEPTOR_ID=1  # aps_prothrombin
LIGAND_IDS="1684,1689,1690,1691,1692"  # 5 ligands

echo "=========================================="
echo "DB-Driven Vina Docking Pipeline"
echo "=========================================="
echo "Start: $(TZ=Asia/Tokyo date '+%Y-%m-%d %H:%M:%S JST')"

# Get receptor info from DB
RECEPTOR_INFO=$(sqlite3 "$DB" "SELECT name, pdb_id, pdbqt_path, box_center_x, box_center_y, box_center_z, box_size_x, box_size_y, box_size_z FROM receptors WHERE id=$RECEPTOR_ID;")

IFS='|' read -r REC_NAME REC_PDB REC_PDBQT BOX_X BOX_Y BOX_Z SIZE_X SIZE_Y SIZE_Z <<< "$RECEPTOR_INFO"

echo "Receptor: $REC_NAME ($REC_PDB)"
echo "Box center: $BOX_X, $BOX_Y, $BOX_Z"
echo "Box size: $SIZE_X x $SIZE_Y x $SIZE_Z"
echo ""

# Get ligand info from DB
LIGANDS=$(sqlite3 "$DB" "SELECT id, zinc_id, smiles FROM ligands WHERE id IN ($LIGAND_IDS) AND has_3d=1;")

COUNT=0
while IFS='|' read -r LIG_ID ZINC_ID SMILES; do
  COUNT=$((COUNT + 1))
  echo "----------------------------------------"
  echo "Ligand $COUNT: $ZINC_ID (ID: $LIG_ID)"
  echo "----------------------------------------"
  
  # Check if ligand PDBQT exists
  LIGAND_PDBQT="/home/dev/OCP/catalog/libraries/zinc20_ml/pdbqt/${ZINC_ID}.pdbqt"
  
  if [ ! -f "$LIGAND_PDBQT" ]; then
    echo "⚠️  PDBQT not found, generating from SMILES..."
    mkdir -p /home/dev/OCP/catalog/libraries/zinc20_ml/pdbqt
    
    # Generate 3D structure and convert to PDBQT
    docker run --rm \
      -v /home/dev/OCP:/work \
      -w /work \
      ligand-selector:latest \
      python3 tools/rdkit_gen3d_to_pdbqt.py \
        --smiles "$SMILES" \
        --output "catalog/libraries/zinc20_ml/pdbqt/${ZINC_ID}.pdbqt" \
        --name "$ZINC_ID"
  fi
  
  # Prepare output directory
  OUTPUT_DIR="/home/dev/OCP/results/vina_output/${ZINC_ID}_${REC_NAME}"
  mkdir -p "$OUTPUT_DIR"
  
  # Create Vina config
  cat > "$OUTPUT_DIR/config.txt" << VINA_CONFIG
receptor = $REC_PDBQT
ligand = $LIGAND_PDBQT

center_x = $BOX_X
center_y = $BOX_Y
center_z = $BOX_Z
size_x = $SIZE_X
size_y = $SIZE_Y
size_z = $SIZE_Z

exhaustiveness = 16
num_modes = 9
energy_range = 3
VINA_CONFIG
  
  # Run Vina
  echo "Running Vina docking..."
  docker run --rm \
    -v /home/dev/OCP:/work \
    -w /work/results/vina_output/${ZINC_ID}_${REC_NAME} \
    vina-runner:latest \
    vina --config config.txt --out docked.pdbqt --log vina.log
  
  # Extract best score
  SCORE=$(grep "^   1 " "$OUTPUT_DIR/vina.log" | awk '{print $2}')
  echo "✅ Docking complete. Best score: $SCORE kcal/mol"
  
  # Update DB with docking result
  sqlite3 "$DB" "UPDATE ligands SET docking_score=$SCORE WHERE id=$LIG_ID;"
  
done <<< "$LIGANDS"

echo ""
echo "=========================================="
echo "Completed: $COUNT compounds docked"
echo "End: $(TZ=Asia/Tokyo date '+%Y-%m-%d %H:%M:%S JST')"
echo "=========================================="
