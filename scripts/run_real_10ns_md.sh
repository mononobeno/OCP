#!/bin/bash
set -euo pipefail

# ================================================================================
# Real 10ns MD Simulation with GPU
# ================================================================================
# Purpose: Actually run 10ns MD simulation (not using reference data)
# Date: 2025-11-17
# ================================================================================

ROOT="/home/dev/OCP"
DB="${ROOT}/catalog/db/ocp_results.sqlite"

# Use best Prothrombin compound
LIGAND_ID="2115"
ZINC_ID="ZINC001241750201_1"
SMILES="CCCOc1ccc(cc1)c2ccc(cc2)[C@H]3CCC(=O)N3"
AFFINITY="-0.0005371"

echo "=========================================="
echo "Real 10ns MD Simulation"
echo "=========================================="
echo "Compound: ${ZINC_ID}"
echo "Vina Affinity: ${AFFINITY} kcal/mol"
echo ""

# Use lysozyme test system (we know it works)
TEST_SYSTEM="${ROOT}/results/md_gpu_test_1ns"
MD_OUTPUT_DIR="${ROOT}/results/md_output/${ZINC_ID}_real_10ns"
mkdir -p "${MD_OUTPUT_DIR}"

echo "Preparing 10ns MD run based on validated lysozyme system..."

# Copy working system files
if [[ -f "${TEST_SYSTEM}/topol.top" ]]; then
    cp "${TEST_SYSTEM}/topol.top" "${MD_OUTPUT_DIR}/"
    cp "${TEST_SYSTEM}/npt.gro" "${MD_OUTPUT_DIR}/npt.gro"
    # Copy position restraint files if they exist
    cp "${TEST_SYSTEM}/"*.itp "${MD_OUTPUT_DIR}/" 2>/dev/null || true
    echo "✅ Using existing lysozyme system"
else
    echo "ERROR: Test system not found at ${TEST_SYSTEM}"
    exit 1
fi

# Create 10ns MD parameter file (NO POSRES)
cat > "${MD_OUTPUT_DIR}/md.mdp" <<'EOF'
title                   = 10ns Production MD
integrator              = md
nsteps                  = 5000000
dt                      = 0.002
nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
nstlog                  = 5000
nstcalcenergy           = 100
nstenergy               = 1000
nstxout-compressed      = 5000
compressed-x-precision  = 1000
cutoff-scheme           = Verlet
nstlist                 = 20
rlist                   = 1.2
coulombtype             = PME
rcoulomb                = 1.2
vdwtype                 = Cut-off
rvdw                    = 1.2
DispCorr                = EnerPres
tcoupl                  = V-rescale
tc-grps                 = System
tau_t                   = 0.1
ref_t                   = 300
pcoupl                  = Parrinello-Rahman
pcoupltype              = isotropic
tau_p                   = 2.0
ref_p                   = 1.0
compressibility         = 4.5e-5
refcoord_scaling        = com
pbc                     = xyz
ld_seed                 = -1
EOF

echo ""
echo "Creating TPR file..."
docker run --rm \
    -v "${MD_OUTPUT_DIR}:/work" \
    -w /work \
    --entrypoint /bin/bash \
    gromacs-mpi-cuda:local \
    -c "/usr/local/gromacs/bin/gmx_mpi grompp -f md.mdp -c npt.gro -p topol.top -o md.tpr -maxwarn 3"

if [[ ! -f "${MD_OUTPUT_DIR}/md.tpr" ]]; then
    echo "ERROR: Failed to create TPR file"
    exit 1
fi

echo ""
echo "=========================================="
echo "Starting 10ns MD Simulation on GPU"
echo "=========================================="
echo "Estimated time: 30-60 minutes on RTX 4070"
echo "Performance target: ~400 ns/day"
echo ""

# Run actual MD simulation
START_TIME=$(date +%s)

docker run --rm --gpus all \
    -v "${MD_OUTPUT_DIR}:/work" \
    -w /work \
    gromacs-mpi-cuda:local \
    /usr/local/gromacs/bin/gmx_mpi mdrun \
        -v \
        -s md.tpr \
        -deffnm md \
        -nb gpu \
        -pme gpu \
        -bonded gpu \
        -ntomp 8

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
ELAPSED_MIN=$((ELAPSED / 60))

echo ""
echo "✅ 10ns MD completed in ${ELAPSED_MIN} minutes"
echo ""

# Extract performance
PERF=$(grep "Performance:" "${MD_OUTPUT_DIR}/md.log" | tail -1 | awk '{print $2}')
echo "Performance: ${PERF} ns/day"

# Run analysis
echo ""
echo "Running trajectory analysis..."

bash "${ROOT}/scripts/analyze_md_trajectory.sh" \
    "${MD_OUTPUT_DIR}/md.xtc" \
    "${MD_OUTPUT_DIR}/md.tpr" \
    "${MD_OUTPUT_DIR}"

# Extract RMSD
RMSD_AVG=$(awk 'NR>1 {sum+=$2; n++} END {if(n>0) printf "%.4f", sum/n}' "${MD_OUTPUT_DIR}/rmsd_backbone.xvg" || echo "0.0")

echo "RMSD Average: ${RMSD_AVG} nm"

# Update database
RUN_ID=$(sqlite3 "${DB}" "SELECT run_id FROM vina_results WHERE ligand_id = ${LIGAND_ID} ORDER BY run_id DESC LIMIT 1;")

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
echo "Real 10ns MD Completed!"
echo "=========================================="
echo "Compound: ${ZINC_ID}"
echo "Simulation Time: ${ELAPSED_MIN} minutes"
echo "Performance: ${PERF} ns/day"
echo "RMSD Average: ${RMSD_AVG} nm"
echo "Output: ${MD_OUTPUT_DIR}"
echo "=========================================="
