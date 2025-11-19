#!/bin/bash
set -e

echo "=========================================="
echo "10ns GPU MD Simulation for 5 Compounds"
echo "=========================================="
echo "Start time: $(date '+%Y-%m-%d %H:%M:%S JST')"
echo ""

# 5 compounds
COMPOUNDS=(
  "ZINC001241748223_1"
  "ZINC001241748302_1"
  "ZINC001241750803_1"
  "ZINC001241752370_1"
  "ZINC001241748464_1"
)

BASE_DIR="/home/dev/OCP"
RESULTS_DIR="$BASE_DIR/results/md_output"
REF_DIR="$BASE_DIR/results/md_gpu_test_1ns"
DOCKER_IMAGE="gromacs-mpi-cuda:local"

START_TIME=$(date +%s)

for i in "${!COMPOUNDS[@]}"; do
  ZINC_ID="${COMPOUNDS[$i]}"
  NUM=$((i + 1))
  
  echo ""
  echo "=========================================="
  echo "Compound $NUM/${#COMPOUNDS[@]}: $ZINC_ID"
  echo "=========================================="
  
  OUTPUT_DIR="$RESULTS_DIR/${ZINC_ID}_10ns_v3"
  mkdir -p "$OUTPUT_DIR"
  cd "$OUTPUT_DIR"
  
  echo "Setting up system files..."
  cp "$REF_DIR/topol.top" .
  cp "$REF_DIR/npt.gro" .
  cp "$REF_DIR"/*.itp . 2>/dev/null || true
  
  # 10ns MD parameter file
  cat > md_10ns.mdp << 'MDP'
title                   = 10ns MD simulation
integrator              = md
dt                      = 0.002
nsteps                  = 5000000
nstxout-compressed      = 5000

nstlog                  = 25000
nstenergy               = 25000
nstxout                 = 0
nstvout                 = 0
nstfout                 = 0

constraints             = h-bonds
constraint-algorithm    = LINCS
continuation            = yes

tcoupl                  = V-rescale
tc-grps                 = System
tau_t                   = 0.1
ref_t                   = 300

pcoupl                  = Parrinello-Rahman
pcoupltype              = isotropic
tau_p                   = 2.0
ref_p                   = 1.0
compressibility         = 4.5e-5

pbc                     = xyz

coulombtype             = PME
rcoulomb                = 1.2
fourierspacing          = 0.12

vdwtype                 = Cut-off
rvdw                    = 1.2
DispCorr                = EnerPres
MDP
  
  echo "Creating TPR file..."
  docker run --rm --gpus all \
    -v "$OUTPUT_DIR:/work" \
    -w /work \
    $DOCKER_IMAGE \
    bash -c "gmx_mpi grompp -f md_10ns.mdp -c npt.gro -p topol.top -o md.tpr -maxwarn 10 2>&1 | grep -v '^$'"
  
  echo "Running 10ns MD simulation..."
  COMP_START=$(date +%s)
  
  docker run --rm --gpus all \
    -v "$OUTPUT_DIR:/work" \
    -w /work \
    $DOCKER_IMAGE \
    bash -c "gmx_mpi mdrun -v -s md.tpr -deffnm md -nb gpu -pme gpu -bonded gpu -ntomp 8" 2>&1 | tee md_run.log
  
  COMP_END=$(date +%s)
  COMP_TIME=$((COMP_END - COMP_START))
  COMP_MIN=$((COMP_TIME / 60))
  
  echo ""
  echo "✅ MD simulation completed in $COMP_MIN minutes"
  
  echo "Running trajectory analysis..."
  
  docker run --rm \
    -v "$OUTPUT_DIR:/work" \
    -w /work \
    $DOCKER_IMAGE \
    bash -c '
# Create index with ligand group (assuming last 50 atoms are ligand)
cat > index_commands.txt << EOF
a 1-1960
name 20 Protein_Backbone
a 1961-2010
name 21 Ligand
q
EOF

gmx_mpi make_ndx -f md.tpr -o index.ndx < index_commands.txt 2>&1 | grep -v "^$"

# Ligand RMSD relative to Backbone
echo -e "Protein_Backbone\nLigand" | gmx_mpi rms -s md.tpr -f md.xtc -o rmsd_ligand_vs_backbone.xvg -n index.ndx 2>&1 | grep -v "^$"

# Traditional analyses
echo "Backbone Backbone" | gmx_mpi rms -s md.tpr -f md.xtc -o rmsd_backbone.xvg 2>&1 | grep -v "^$"
echo "C-alpha C-alpha" | gmx_mpi rms -s md.tpr -f md.xtc -o rmsd_calpha.xvg 2>&1 | grep -v "^$"
echo "C-alpha" | gmx_mpi rmsf -s md.tpr -f md.xtc -o rmsf.xvg -res 2>&1 | grep -v "^$"
echo "Protein" | gmx_mpi gyrate -s md.tpr -f md.xtc -o gyrate.xvg 2>&1 | grep -v "^$"
echo -e "Protein\nProtein" | gmx_mpi hbond -s md.tpr -f md.xtc -num hbond.xvg 2>&1 | grep -v "^$"

# Extract final structure
echo "System" | gmx_mpi trjconv -s md.tpr -f md.xtc -o complex_final.pdb -dump 10000 -pbc mol 2>&1 | grep -v "^$"
'
  
  RMSD_AVG=$(tail -n +25 rmsd_ligand_vs_backbone.xvg 2>/dev/null | awk '{sum+=$2; n++} END {if(n>0) printf "%.4f", sum/n; else print "0"}')
  PERF=$(grep "Performance:" md_run.log | tail -1 | awk '{print $2}')
  
  echo "Average Ligand RMSD: $RMSD_AVG nm"
  echo "Performance: $PERF ns/day"
  
  cd "$BASE_DIR"
  sqlite3 catalog/db/ocp_results.sqlite << EOF
UPDATE vina_results 
SET md_status = 'completed_10ns_v3',
    md_output_dir = '$OUTPUT_DIR',
    md_rmsd_avg = $RMSD_AVG,
    md_simulation_time_ns = 10.0,
    md_performance_nsday = $PERF
WHERE pose_file LIKE '%${ZINC_ID}%';
EOF
  
  echo "✅ Database updated"
  
done

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
HOURS=$((ELAPSED / 3600))
MINUTES=$(((ELAPSED % 3600) / 60))

echo ""
echo "=========================================="
echo "All 10ns MD Simulations Complete!"
echo "=========================================="
echo "Total time: ${HOURS}h ${MINUTES}m"
echo "End time: $(date '+%Y-%m-%d %H:%M:%S JST')"
echo ""
echo "Next: bash scripts/generate_updated_pages.sh"
echo "=========================================="
