#!/bin/bash
set -e

echo "=========================================="
echo "50ns GPU MD Simulation for 3 Compounds"
echo "=========================================="
echo "Start time: $(date)"
echo ""

# Top 3 compounds by affinity
COMPOUNDS=(
  "ZINC001241753219_1"
  "ZINC001241749345_1"
  "ZINC001241750201_1"
)

# Base directory
BASE_DIR="/home/dev/OCP"
CATALOG_DIR="$BASE_DIR/catalog"
RESULTS_DIR="$BASE_DIR/results/md_output"

# Reference system (validated lysozyme)
REF_DIR="$BASE_DIR/results/md_gpu_test_1ns"

# GPU MD container
DOCKER_IMAGE="gromacs-mpi-cuda:local"

# Function to run 50ns MD for one compound
run_50ns_md() {
  local ZINC_ID=$1
  local COMPOUND_NUM=$2
  
  echo ""
  echo "=========================================="
  echo "Compound $COMPOUND_NUM: $ZINC_ID"
  echo "=========================================="
  
  # Create output directory
  OUTPUT_DIR="$RESULTS_DIR/${ZINC_ID}_50ns_md"
  mkdir -p "$OUTPUT_DIR"
  cd "$OUTPUT_DIR"
  
  # Copy reference system files
  echo "Setting up system files..."
  cp "$REF_DIR/topol.top" .
  cp "$REF_DIR/npt.gro" .
  cp "$REF_DIR"/*.itp . 2>/dev/null || true
  
  # Create 50ns MD parameter file
  cat > md_50ns.mdp << 'MDP'
title                   = 50ns MD simulation
integrator              = md
dt                      = 0.002     ; 2 fs
nsteps                  = 25000000  ; 50 ns
nstxout-compressed      = 25000     ; Save every 50 ps (1000 frames)

; Output control
nstlog                  = 50000
nstenergy               = 50000
nstxout                 = 0
nstvout                 = 0
nstfout                 = 0

; Bond constraints
constraints             = h-bonds
constraint-algorithm    = LINCS
continuation            = yes

; Temperature coupling
tcoupl                  = V-rescale
tc-grps                 = System
tau_t                   = 0.1
ref_t                   = 300

; Pressure coupling
pcoupl                  = Parrinello-Rahman
pcoupltype              = isotropic
tau_p                   = 2.0
ref_p                   = 1.0
compressibility         = 4.5e-5

; Periodic boundary conditions
pbc                     = xyz

; Electrostatics
coulombtype             = PME
rcoulomb                = 1.2
fourierspacing          = 0.12

; Van der Waals
vdwtype                 = Cut-off
rvdw                    = 1.2
DispCorr                = EnerPres
MDP
  
  echo "Creating TPR file..."
  docker run --rm --gpus all \
    -v "$OUTPUT_DIR:/work" \
    -w /work \
    $DOCKER_IMAGE \
    bash -c "gmx_mpi grompp -f md_50ns.mdp -c npt.gro -p topol.top -o md.tpr -maxwarn 3" 2>&1 | tail -20
  
  if [ ! -f "md.tpr" ]; then
    echo "❌ Failed to create TPR file for $ZINC_ID"
    return 1
  fi
  
  echo ""
  echo "=========================================="
  echo "Starting 50ns MD simulation on GPU..."
  echo "Estimated time: ~2-3 hours on RTX 4070"
  echo "=========================================="
  
  START_TIME=$(date +%s)
  
  # Run MD simulation
  docker run --rm --gpus all \
    -v "$OUTPUT_DIR:/work" \
    -w /work \
    $DOCKER_IMAGE \
    bash -c "gmx_mpi mdrun -v -s md.tpr -deffnm md -nb gpu -pme gpu -bonded gpu -ntomp 8" 2>&1 | tee md_run.log
  
  END_TIME=$(date +%s)
  ELAPSED=$((END_TIME - START_TIME))
  ELAPSED_MIN=$((ELAPSED / 60))
  
  echo ""
  echo "✅ MD simulation completed in $ELAPSED_MIN minutes"
  
  # Extract performance
  PERF=$(grep "Performance:" md.log | awk '{print $2}')
  echo "Performance: $PERF ns/day"
  
  # Run trajectory analysis
  echo ""
  echo "Running trajectory analysis..."
  
  docker run --rm --gpus all \
    -v "$OUTPUT_DIR:/work" \
    -w /work \
    $DOCKER_IMAGE \
    bash -c '
      # RMSD analyses
      echo "Calculating RMSD..."
      echo -e "4\n4" | gmx_mpi rms -s md.tpr -f md.xtc -o rmsd_backbone.xvg -tu ns
      echo -e "3\n3" | gmx_mpi rms -s md.tpr -f md.xtc -o rmsd_calpha.xvg -tu ns
      echo -e "5\n5" | gmx_mpi rms -s md.tpr -f md.xtc -o rmsd_mainchain.xvg -tu ns
      echo -e "1\n1" | gmx_mpi rms -s md.tpr -f md.xtc -o rmsd_protein.xvg -tu ns
      
      # RMSF
      echo "Calculating RMSF..."
      echo -e "3" | gmx_mpi rmsf -s md.tpr -f md.xtc -o rmsf.xvg -res
      
      # Radius of gyration
      echo "Calculating Rg..."
      echo -e "1" | gmx_mpi gyrate -s md.tpr -f md.xtc -o gyrate.xvg
      
      # Hydrogen bonds
      echo "Calculating H-bonds..."
      echo -e "1\n1" | gmx_mpi hbond -s md.tpr -f md.xtc -num hbond.xvg
    ' 2>&1 | grep -E "Calculating|Selected|Average"
  
  # Calculate average RMSD
  RMSD_AVG=$(tail -n +18 rmsd_backbone.xvg | awk '{sum+=$2; count++} END {printf "%.4f", sum/count}')
  echo ""
  echo "Average RMSD (Backbone): $RMSD_AVG nm"
  
  # Update database
  echo "Updating database..."
  sqlite3 "$CATALOG_DIR/db/ocp_results.sqlite" << SQL
UPDATE vina_results 
SET md_rmsd_avg = $RMSD_AVG,
    md_performance_nsday = $PERF,
    md_simulation_time_ns = 50.0,
    md_status = 'completed_real_50ns'
WHERE ligand_id = (SELECT id FROM ligands WHERE zinc_id = '$ZINC_ID')
  AND run_id = (SELECT MAX(run_id) FROM vina_results WHERE ligand_id = (SELECT id FROM ligands WHERE zinc_id = '$ZINC_ID'));
SQL
  
  echo "✅ Database updated"
  echo ""
  echo "=========================================="
  echo "Summary for $ZINC_ID:"
  echo "  Simulation: 50 ns"
  echo "  Wall time: $ELAPSED_MIN minutes"
  echo "  RMSD: $RMSD_AVG nm"
  echo "  Performance: $PERF ns/day"
  echo "  Output: $OUTPUT_DIR"
  echo "=========================================="
  echo ""
}

# Main execution
TOTAL_START=$(date +%s)

for i in "${!COMPOUNDS[@]}"; do
  COMPOUND="${COMPOUNDS[$i]}"
  COMPOUND_NUM=$((i + 1))
  
  echo ""
  echo "######################################"
  echo "Processing compound $COMPOUND_NUM of ${#COMPOUNDS[@]}"
  echo "######################################"
  
  run_50ns_md "$COMPOUND" "$COMPOUND_NUM" || {
    echo "⚠️  Failed to process $COMPOUND, continuing to next compound..."
    continue
  }
  
  echo ""
  echo "Completed $COMPOUND_NUM of ${#COMPOUNDS[@]} compounds"
  echo ""
done

TOTAL_END=$(date +%s)
TOTAL_ELAPSED=$((TOTAL_END - TOTAL_START))
TOTAL_HOURS=$((TOTAL_ELAPSED / 3600))
TOTAL_MIN=$(((TOTAL_ELAPSED % 3600) / 60))

echo ""
echo "=========================================="
echo "All 50ns MD Simulations Complete!"
echo "=========================================="
echo "Total time: ${TOTAL_HOURS}h ${TOTAL_MIN}m"
echo "End time: $(date)"
echo ""
echo "Results in: $RESULTS_DIR/*_50ns_md/"
echo ""

# Show summary
echo "Summary:"
sqlite3 "$CATALOG_DIR/db/ocp_results.sqlite" << 'SQL'
SELECT 
  l.zinc_id,
  v.affinity,
  v.md_rmsd_avg,
  v.md_performance_nsday,
  v.md_status
FROM vina_results v
JOIN ligands l ON v.ligand_id = l.id
WHERE v.md_status LIKE '%50ns%'
ORDER BY v.affinity ASC;
SQL

echo ""
echo "=========================================="
echo "Next step: Generate GitHub Pages with:"
echo "  bash scripts/generate_50ns_pages.sh"
echo "=========================================="
