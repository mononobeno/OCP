#!/bin/bash
set -euo pipefail

# ========================================
# GROMACS Preparation Script
# Vina pose → GROMACS-ready complex
# ========================================

: "${RUN_ID:?RUN_ID required}"
: "${LIGAND_ID:?LIGAND_ID required}"
: "${WORKSPACE:=/workspace}"
: "${DB_PATH:=${WORKSPACE}/catalog/db/ocp_results.sqlite}"

cd "${WORKSPACE}"

RESULTS_DIR="results/gmx_input"
VINA_POSES_DIR="results/vina_output/poses"
OUTPUT_DIR="results/gmx_output/${LIGAND_ID}"

mkdir -p "${OUTPUT_DIR}"

echo "========================================="
echo "GROMACS Preparation"
echo "========================================="
echo "RUN_ID      : ${RUN_ID}"
echo "LIGAND_ID   : ${LIGAND_ID}"
echo "Output      : ${OUTPUT_DIR}"
echo "========================================="

# ----------------------------------------
# 1. Get best Vina pose from DB
# ----------------------------------------
echo "[1/6] Querying best Vina pose from DB..."

QUERY="
SELECT 
    vr.pose_file,
    vr.affinity,
    t.code AS target_code,
    t.receptor_pdbqt
FROM vina_results vr
JOIN ligands l ON vr.ligand_id = l.id
JOIN runs r ON vr.run_id = r.id
JOIN targets t ON r.target_id = t.id
WHERE r.id = ${RUN_ID}
  AND l.zinc_id = '${LIGAND_ID}'
ORDER BY vr.affinity ASC
LIMIT 1;
"

RESULT=$(sqlite3 "${DB_PATH}" "${QUERY}")

if [[ -z "${RESULT}" ]]; then
    echo "ERROR: No Vina result found for RUN_ID=${RUN_ID}, LIGAND_ID=${LIGAND_ID}"
    exit 1
fi

POSE_FILE=$(echo "${RESULT}" | cut -d'|' -f1)
AFFINITY=$(echo "${RESULT}" | cut -d'|' -f2)
TARGET_CODE=$(echo "${RESULT}" | cut -d'|' -f3)
RECEPTOR_PDBQT=$(echo "${RESULT}" | cut -d'|' -f4)

echo "  Pose file  : ${POSE_FILE}"
echo "  Affinity   : ${AFFINITY} kcal/mol"
echo "  Target     : ${TARGET_CODE}"
echo "  Receptor   : ${RECEPTOR_PDBQT}"

# ----------------------------------------
# 2. Convert PDBQT pose to PDB
# ----------------------------------------
echo "[2/6] Converting PDBQT pose to PDB..."

LIGAND_PDBQT="${WORKSPACE}/${POSE_FILE}"
LIGAND_PDB="${OUTPUT_DIR}/ligand.pdb"

if [[ ! -f "${LIGAND_PDBQT}" ]]; then
    echo "ERROR: Pose file not found: ${LIGAND_PDBQT}"
    exit 1
fi

# Extract only the first MODEL from Vina output (best pose)
LIGAND_PDBQT_MODEL1="${OUTPUT_DIR}/ligand_model1.pdbqt"
sed -n '/^MODEL 1$/,/^ENDMDL$/p' "${LIGAND_PDBQT}" > "${LIGAND_PDBQT_MODEL1}"

# Convert to PDB with hydrogen atoms
obabel -ipdbqt "${LIGAND_PDBQT_MODEL1}" -opdb -O "${LIGAND_PDB}" -h

echo "  Generated: $(basename ${LIGAND_PDB}) ($(stat -c%s ${LIGAND_PDB}) bytes)"

# ----------------------------------------
# 3. Convert receptor PDBQT to PDB
# ----------------------------------------
echo "[3/6] Converting receptor PDBQT to PDB..."

RECEPTOR_PDBQT_FULL="${WORKSPACE}/${RECEPTOR_PDBQT}"
RECEPTOR_PDB="${OUTPUT_DIR}/receptor.pdb"

if [[ ! -f "${RECEPTOR_PDBQT_FULL}" ]]; then
    echo "ERROR: Receptor PDBQT not found: ${RECEPTOR_PDBQT_FULL}"
    exit 1
fi

obabel -ipdbqt "${RECEPTOR_PDBQT_FULL}" -opdb -O "${RECEPTOR_PDB}"

echo "  Generated: $(basename ${RECEPTOR_PDB}) ($(stat -c%s ${RECEPTOR_PDB}) bytes)"

# ----------------------------------------
# 4. Create complex PDB
# ----------------------------------------
echo "[4/6] Creating complex PDB (receptor + ligand)..."

COMPLEX_PDB="${OUTPUT_DIR}/complex.pdb"

cat "${RECEPTOR_PDB}" "${LIGAND_PDB}" > "${COMPLEX_PDB}"

echo "  Generated: $(basename ${COMPLEX_PDB}) ($(stat -c%s ${COMPLEX_PDB}) bytes)"

# ----------------------------------------
# 5. Generate ligand topology with ACPYPE
# ----------------------------------------
echo "[5/6] Generating ligand topology with ACPYPE..."

# ACPYPE requires MOL2 format for best results
LIGAND_MOL2="${OUTPUT_DIR}/ligand.mol2"
obabel -ipdb "${LIGAND_PDB}" -omol2 -O "${LIGAND_MOL2}" -h

if [[ ! -f "${LIGAND_MOL2}" ]] || [[ ! -s "${LIGAND_MOL2}" ]]; then
    echo "  Warning: MOL2 conversion failed, trying direct PDB conversion..."
    cd "${OUTPUT_DIR}"
    acpype -i "ligand.pdb" -b "LIG" -o gmx -n 0 2>&1 | tee acpype.log || {
        echo "ERROR: ACPYPE failed with PDB input"
        exit 1
    }
else
    echo "  Converted to MOL2: $(basename ${LIGAND_MOL2}) ($(stat -c%s ${LIGAND_MOL2}) bytes)"
    
    # Run ACPYPE from output directory (it creates files in current dir)
    cd "${OUTPUT_DIR}"
    acpype -i "ligand.mol2" -b "LIG" -o gmx -n 0 2>&1 | tee acpype.log || {
        echo "WARNING: ACPYPE failed with MOL2, trying PDB..."
        acpype -i "ligand.pdb" -b "LIG" -o gmx -n 0 2>&1 | tee acpype_fallback.log || {
            echo "ERROR: ACPYPE failed with both MOL2 and PDB"
            exit 1
        }
    }
fi

# ACPYPE output: LIG.acpype/LIG_GMX.gro, LIG_GMX.top, etc.
if [[ -d "LIG.acpype" ]]; then
    mv LIG.acpype/LIG_GMX.gro ligand_gmx.gro
    mv LIG.acpype/LIG_GMX.top ligand_gmx.top
    mv LIG.acpype/LIG_GMX.itp ligand_gmx.itp
    echo "  Generated GROMACS topology files:"
    echo "    - ligand_gmx.gro (coordinates)"
    echo "    - ligand_gmx.top (topology)"
    echo "    - ligand_gmx.itp (parameters)"
else
    echo "ERROR: ACPYPE output directory not found"
    exit 1
fi

# ----------------------------------------
# 6. Process receptor with pdb2gmx
# ----------------------------------------
echo "[6/6] Processing receptor with GROMACS pdb2gmx..."

cd "${OUTPUT_DIR}"

# Use AMBER99SB-ILDN force field (matches GAFF for ligand)
echo "1" | gmx pdb2gmx \
    -f "${RECEPTOR_PDB}" \
    -o receptor_processed.gro \
    -p receptor_topol.top \
    -i receptor_posre.itp \
    -water tip3p \
    -ff amber99sb-ildn \
    -ignh \
    2>&1 | tee pdb2gmx.log || {
        echo "ERROR: pdb2gmx failed on receptor"
        exit 1
    }

echo "  Generated receptor topology:"
echo "    - receptor_processed.gro"
echo "    - receptor_topol.top"
echo "    - receptor_posre.itp"

# ----------------------------------------
# Summary
# ----------------------------------------
echo ""
echo "========================================="
echo "GROMACS Preparation Complete!"
echo "========================================="
echo "Output directory: ${OUTPUT_DIR}"
echo ""
echo "Generated files:"
ls -lh "${OUTPUT_DIR}"/*.{gro,top,itp,pdb} 2>/dev/null || true
echo ""
echo "Next steps:"
echo "  1. Merge receptor and ligand topologies"
echo "  2. Create simulation box (gmx editconf)"
echo "  3. Solvate system (gmx solvate)"
echo "  4. Add ions (gmx genion)"
echo "  5. Energy minimization (gmx mdrun -v -deffnm em)"
echo "========================================="

# Update DB with preparation status
sqlite3 "${DB_PATH}" <<EOF
UPDATE vina_results
SET gromacs_prep_status = 'completed',
    gromacs_prep_dir = '${OUTPUT_DIR}'
WHERE run_id = ${RUN_ID}
  AND ligand_id = (SELECT id FROM ligands WHERE zinc_id = '${LIGAND_ID}');
EOF

echo "[DB] Updated gromacs_prep_status to 'completed'"
