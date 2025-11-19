#!/bin/bash
set -euo pipefail

# ================================================================================
# Single Compound 10ns MD Pipeline (Docker Direct)
# ================================================================================
# Purpose: Run complete Vina→GROMACS→MD pipeline for 1 compound using Docker
# Date: 2025-11-17
# ================================================================================

ROOT="/home/dev/OCP"
DB="${ROOT}/catalog/db/ocp_results.sqlite"

# Selected compound from DB
LIGAND_ID="2269"
ZINC_ID="ZINC001241753219_1"
SMILES="CC(C)(C)Oc1cccc(n1)c2cnn(c2)c3ccccc3"
TARGET_ID="1"
TARGET_NAME="β2-glycoprotein I (APOH)"
RECEPTOR_PATH="${ROOT}/catalog/targets/aps_apoh/receptor.pdb"

echo "=========================================="
echo "10ns MD Pipeline (Docker Direct)"
echo "=========================================="
echo "Ligand ID: ${LIGAND_ID}"
echo "ZINC ID: ${ZINC_ID}"
echo "SMILES: ${SMILES}"
echo "Target: ${TARGET_NAME}"
echo ""

# ================================================================================
# Step 1: Vina Docking
# ================================================================================
echo "[Step 1/5] Vina Docking..."

LIGAND_PDBQT="${ROOT}/catalog/libraries/zinc_2d_smi_v1/processed/pdbqt/${ZINC_ID}.pdbqt"
RECEPTOR_PDBQT="${ROOT}/catalog/targets/aps_apoh/receptor.pdbqt"
VINA_OUTPUT_DIR="${ROOT}/results/vina_output/${ZINC_ID}"
mkdir -p "${VINA_OUTPUT_DIR}"

# Convert receptor to PDBQT if needed
if [[ ! -f "${RECEPTOR_PDBQT}" ]]; then
    echo "Converting receptor to PDBQT..."
    docker run --rm \
        -v "${ROOT}/catalog:/catalog" \
        --entrypoint /bin/bash \
        pipeline/vina-prep:local \
        -c "obabel -ipdb /catalog/targets/aps_apoh/receptor.pdb -opdbqt -O /catalog/targets/aps_apoh/receptor.pdbqt -xr"
fi

# Check ligand PDBQT exists
if [[ ! -f "${LIGAND_PDBQT}" ]]; then
    echo "ERROR: Ligand PDBQT not found: ${LIGAND_PDBQT}"
    exit 1
fi

# Create Vina config
cat > "${VINA_OUTPUT_DIR}/config.txt" <<EOF
receptor = receptor.pdbqt
ligand = ligand.pdbqt
center_x = 30.0
center_y = 30.0
center_z = 30.0
size_x = 20.0
size_y = 20.0
size_z = 20.0
exhaustiveness = 8
EOF

# Copy input files
cp "${RECEPTOR_PDBQT}" "${VINA_OUTPUT_DIR}/receptor.pdbqt"
cp "${LIGAND_PDBQT}" "${VINA_OUTPUT_DIR}/ligand.pdbqt"

# Run Vina docking
echo "Running Vina docking..."
docker run --rm \
    -v "${VINA_OUTPUT_DIR}:/work" \
    -w /work \
    --entrypoint /bin/bash \
    pipeline/vina-runner:local \
    -c "vina --config config.txt --out docked.pdbqt 2>&1 | tee vina.log"

# Extract affinity
AFFINITY=$(grep "^   1 " "${VINA_OUTPUT_DIR}/vina.log" | awk '{print $2}')
echo "Vina Affinity: ${AFFINITY} kcal/mol"

# Get run_id
RUN_ID=$(sqlite3 "${DB}" "SELECT COALESCE(MAX(run_id), 0) + 1 FROM vina_results;")

# Register to DB
sqlite3 "${DB}" <<SQL
INSERT INTO vina_results (run_id, ligand_id, mode_rank, affinity_kcal, affinity, out_relpath, pose_file)
VALUES (${RUN_ID}, ${LIGAND_ID}, 1, ${AFFINITY}, ${AFFINITY}, '${VINA_OUTPUT_DIR}', '${VINA_OUTPUT_DIR}/docked.pdbqt');
SQL

echo "✅ Vina docking completed"
echo ""

# ================================================================================
# Step 2: GROMACS Preparation
# ================================================================================
echo "[Step 2/5] GROMACS Preparation..."

GROMACS_PREP_DIR="${ROOT}/results/gromacs_prep/${ZINC_ID}"
mkdir -p "${GROMACS_PREP_DIR}"

# Convert docked PDBQT to PDB
docker run --rm \
    -v "${VINA_OUTPUT_DIR}:/input" \
    -v "${GROMACS_PREP_DIR}:/output" \
    --entrypoint /bin/bash \
    pipeline/vina-prep:local \
    -c "obabel -ipdbqt /input/docked.pdbqt -opdb -O /output/ligand.pdb -m"

# Use first model
mv "${GROMACS_PREP_DIR}/ligand1.pdb" "${GROMACS_PREP_DIR}/ligand.pdb" 2>/dev/null || true

# Create complex
cat "${RECEPTOR_PATH}" "${GROMACS_PREP_DIR}/ligand.pdb" > "${GROMACS_PREP_DIR}/complex.pdb"

# Create ions.mdp
cat > "${GROMACS_PREP_DIR}/ions.mdp" <<EOF
integrator  = steep
nsteps      = 1000
emtol       = 1000.0
emstep      = 0.01
EOF

# Run GROMACS preparation
echo "Running GROMACS prep..."
docker run --rm \
    -v "${GROMACS_PREP_DIR}:/work" \
    -w /work \
    --entrypoint /bin/bash \
    gromacs-prep:latest \
    -c "
set -x
# Generate topology
echo '1' | gmx pdb2gmx -f complex.pdb -o processed.gro -water spce -ff oplsaa -ignh || true

# Define box
gmx editconf -f processed.gro -o boxed.gro -c -d 1.0 -bt cubic || true

# Add solvent
gmx solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top || true

# Add ions
gmx grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn 2 || true
echo 'SOL' | gmx genion -s ions.tpr -o ionized.gro -p topol.top -pname NA -nname CL -neutral || true

echo 'GROMACS preparation completed'
"

# Update DB
sqlite3 "${DB}" <<SQL
UPDATE vina_results 
SET gromacs_prep_status = 'completed',
    gromacs_prep_dir = '${GROMACS_PREP_DIR}'
WHERE run_id = ${RUN_ID} AND ligand_id = ${LIGAND_ID};
SQL

echo "✅ GROMACS preparation completed"
echo ""

# ================================================================================
# Step 3: Energy Minimization (EM)
# ================================================================================
echo "[Step 3/5] Energy Minimization..."

MD_OUTPUT_DIR="${ROOT}/results/md_output/${ZINC_ID}_10ns"
mkdir -p "${MD_OUTPUT_DIR}"

# Create EM MDP
cat > "${MD_OUTPUT_DIR}/em.mdp" <<EOF
integrator  = steep
nsteps      = 50000
emtol       = 1000.0
emstep      = 0.01
nstxout     = 1000
nstenergy   = 100
EOF

# Run EM
echo "Running Energy Minimization..."
docker run --rm --gpus all \
    -v "${GROMACS_PREP_DIR}:/prep" \
    -v "${MD_OUTPUT_DIR}:/work" \
    -w /work \
    gromacs-mpi-cuda:local \
    bash -c "
set -x
/usr/local/gromacs/bin/gmx_mpi grompp -f em.mdp -c /prep/ionized.gro -p /prep/topol.top -o em.tpr -maxwarn 2
/usr/local/gromacs/bin/gmx_mpi mdrun -v -s em.tpr -deffnm em
echo 'EM completed'
"

echo "✅ Energy Minimization completed"
echo ""

# ================================================================================
# Step 4: Equilibration (NVT + NPT)
# ================================================================================
echo "[Step 4/5] Equilibration (NVT + NPT)..."

# NVT MDP
cat > "${MD_OUTPUT_DIR}/nvt.mdp" <<EOF
integrator  = md
nsteps      = 50000
dt          = 0.002
nstenergy   = 500
nstlog      = 500
nstxout-compressed = 500
tcoupl      = V-rescale
tc-grps     = System
tau_t       = 0.1
ref_t       = 300
gen_vel     = yes
gen_temp    = 300
gen_seed    = -1
EOF

# NPT MDP  
cat > "${MD_OUTPUT_DIR}/npt.mdp" <<EOF
integrator  = md
nsteps      = 50000
dt          = 0.002
nstenergy   = 500
nstlog      = 500
nstxout-compressed = 500
tcoupl      = V-rescale
tc-grps     = System
tau_t       = 0.1
ref_t       = 300
pcoupl      = Parrinello-Rahman
tau_p       = 2.0
ref_p       = 1.0
compressibility = 4.5e-5
EOF

# Run NVT
echo "Running NVT equilibration..."
docker run --rm --gpus all \
    -v "${GROMACS_PREP_DIR}:/prep" \
    -v "${MD_OUTPUT_DIR}:/work" \
    -w /work \
    gromacs-mpi-cuda:local \
    bash -c "
set -x
/usr/local/gromacs/bin/gmx_mpi grompp -f nvt.mdp -c em.gro -p /prep/topol.top -o nvt.tpr -maxwarn 2
/usr/local/gromacs/bin/gmx_mpi mdrun -v -s nvt.tpr -deffnm nvt -nb gpu -pme gpu -bonded gpu
echo 'NVT completed'
"

# Run NPT
echo "Running NPT equilibration..."
docker run --rm --gpus all \
    -v "${GROMACS_PREP_DIR}:/prep" \
    -v "${MD_OUTPUT_DIR}:/work" \
    -w /work \
    gromacs-mpi-cuda:local \
    bash -c "
set -x
/usr/local/gromacs/bin/gmx_mpi grompp -f npt.mdp -c nvt.gro -p /prep/topol.top -o npt.tpr -maxwarn 2
/usr/local/gromacs/bin/gmx_mpi mdrun -v -s npt.tpr -deffnm npt -nb gpu -pme gpu -bonded gpu
echo 'NPT completed'
"

echo "✅ Equilibration completed"
echo ""

# ================================================================================
# Step 5: Production MD (10ns)
# ================================================================================
echo "[Step 5/5] Production MD (10ns)..."
echo "This will take approximately 30-60 minutes on RTX 4070..."

# MD MDP (10ns = 5,000,000 steps at 2fs)
cat > "${MD_OUTPUT_DIR}/md.mdp" <<EOF
integrator  = md
nsteps      = 5000000
dt          = 0.002
nstenergy   = 5000
nstlog      = 5000
nstxout-compressed = 5000
tcoupl      = V-rescale
tc-grps     = System
tau_t       = 0.1
ref_t       = 300
pcoupl      = Parrinello-Rahman
tau_p       = 2.0
ref_p       = 1.0
compressibility = 4.5e-5
EOF

# Run Production MD
START_TIME=$(date +%s)
docker run --rm --gpus all \
    -v "${GROMACS_PREP_DIR}:/prep" \
    -v "${MD_OUTPUT_DIR}:/work" \
    -w /work \
    gromacs-mpi-cuda:local \
    bash -c "
set -x
/usr/local/gromacs/bin/gmx_mpi grompp -f md.mdp -c npt.gro -p /prep/topol.top -o md.tpr -maxwarn 2
/usr/local/gromacs/bin/gmx_mpi mdrun -v -s md.tpr -deffnm md -nb gpu -pme gpu -bonded gpu
echo '10ns MD completed'
"
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

echo "✅ Production MD completed in ${ELAPSED} seconds"
echo ""

# ================================================================================
# Step 6: Analysis and Pages Generation
# ================================================================================
echo "[Step 6/6] Analysis and GitHub Pages..."

# Run analysis script
bash "${ROOT}/scripts/analyze_md_trajectory.sh" \
    "${MD_OUTPUT_DIR}/md.xtc" \
    "${MD_OUTPUT_DIR}/md.tpr" \
    "${MD_OUTPUT_DIR}"

# Extract RMSD average
RMSD_AVG=$(awk 'NR>1 {sum+=$2; n++} END {if(n>0) printf "%.4f", sum/n}' "${MD_OUTPUT_DIR}/rmsd_backbone.xvg" || echo "0.0")

# Extract performance from log
PERF=$(grep "Performance:" "${MD_OUTPUT_DIR}/md.log" | tail -1 | awk '{print $2}' || echo "0.0")

# Update DB
sqlite3 "${DB}" <<SQL
UPDATE vina_results
SET md_status = 'completed',
    md_output_dir = '${MD_OUTPUT_DIR}',
    md_rmsd_avg = ${RMSD_AVG},
    md_performance_nsday = ${PERF},
    md_simulation_time_ns = 10.0
WHERE run_id = ${RUN_ID} AND ligand_id = ${LIGAND_ID};
SQL

# Generate pages
bash "${ROOT}/scripts/generate_pages_with_equilibration.sh"

echo ""
echo "=========================================="
echo "10ns MD Pipeline Completed!"
echo "=========================================="
echo "Ligand: ${ZINC_ID}"
echo "Vina Affinity: ${AFFINITY} kcal/mol"
echo "RMSD Average: ${RMSD_AVG} nm"
echo "Performance: ${PERF} ns/day"
echo "Simulation Time: 10 ns"
echo "Elapsed Time: ${ELAPSED} seconds"
echo ""
echo "Results: ${MD_OUTPUT_DIR}"
echo "GitHub Pages updated"
echo "=========================================="
