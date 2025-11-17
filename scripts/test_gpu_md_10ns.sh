#!/bin/bash
set -euo pipefail

# ========================================
# 10ns GPU MD Simulation Test
# ========================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

cd "${ROOT}"

TEST_DIR="results/md_gpu_test"
mkdir -p "${TEST_DIR}"

echo "========================================="
echo "10ns GPU MD Simulation Test"
echo "========================================="
echo "Output: ${TEST_DIR}"
echo ""

# Use receptor from catalog as test system
RECEPTOR_PDB="catalog/targets/aps_apoh/receptor.pdb"

if [[ ! -f "${RECEPTOR_PDB}" ]]; then
    echo "ERROR: Receptor PDB not found: ${RECEPTOR_PDB}"
    exit 1
fi

cd "${TEST_DIR}"

echo "[1/7] Preparing system with pdb2gmx..."
echo "1" | gmx pdb2gmx \
    -f "../../${RECEPTOR_PDB}" \
    -o processed.gro \
    -p topol.top \
    -i posre.itp \
    -water tip3p \
    -ff amber99sb-ildn \
    -ignh \
    2>&1 | tail -15

echo ""
echo "[2/7] Creating simulation box..."
gmx editconf \
    -f processed.gro \
    -o boxed.gro \
    -c -d 1.2 -bt cubic \
    2>&1 | tail -8

echo ""
echo "[3/7] Solvating system..."
gmx solvate \
    -cp boxed.gro \
    -cs spc216.gro \
    -o solvated.gro \
    -p topol.top \
    2>&1 | tail -8

echo ""
echo "[4/7] Adding ions..."
cat > ions.mdp <<'MDP'
integrator  = steep
emtol       = 1000.0
emstep      = 0.01
nsteps      = 500
nstlist     = 10
cutoff-scheme = Verlet
ns_type     = grid
coulombtype = PME
rcoulomb    = 1.0
rvdw        = 1.0
pbc         = xyz
MDP

gmx grompp \
    -f ions.mdp \
    -c solvated.gro \
    -p topol.top \
    -o ions.tpr \
    -maxwarn 2 \
    2>&1 | tail -10

echo "SOL" | gmx genion \
    -s ions.tpr \
    -o solvated_ions.gro \
    -p topol.top \
    -pname NA -nname CL \
    -neutral \
    2>&1 | tail -10

echo ""
echo "[5/7] Energy minimization..."
cat > em.mdp <<'MDP'
integrator  = steep
emtol       = 1000.0
emstep      = 0.01
nsteps      = 5000
nstlist     = 10
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb    = 1.0
rvdw        = 1.0
pbc         = xyz
MDP

gmx grompp \
    -f em.mdp \
    -c solvated_ions.gro \
    -p topol.top \
    -o em.tpr \
    -maxwarn 2 \
    2>&1 | tail -10

echo "Running energy minimization on GPU..."
gmx mdrun \
    -v -deffnm em \
    -nb gpu -pme gpu -bonded gpu \
    -nt 8 -ntmpi 1 -ntomp 8 \
    2>&1 | tee em.log | tail -20

echo ""
echo "[6/7] NVT equilibration (100ps)..."
cat > nvt.mdp <<'MDP'
title                   = NVT equilibration
integrator              = md
nsteps                  = 50000     ; 100 ps
dt                      = 0.002     ; 2 fs
nstxout                 = 500
nstvout                 = 500
nstenergy               = 100
nstlog                  = 100
continuation            = no
constraint_algorithm    = lincs
constraints             = h-bonds
cutoff-scheme           = Verlet
ns_type                 = grid
nstlist                 = 10
rcoulomb                = 1.0
rvdw                    = 1.0
coulombtype             = PME
pme_order               = 4
fourierspacing          = 0.16
tcoupl                  = V-rescale
tc-grps                 = Protein Non-Protein
tau_t                   = 0.1     0.1
ref_t                   = 300     300
pcoupl                  = no
pbc                     = xyz
gen_vel                 = yes
gen_temp                = 300
gen_seed                = -1
MDP

gmx grompp \
    -f nvt.mdp \
    -c em.gro \
    -r em.gro \
    -p topol.top \
    -o nvt.tpr \
    -maxwarn 2 \
    2>&1 | tail -10

echo "Running NVT equilibration on GPU..."
gmx mdrun \
    -v -deffnm nvt \
    -nb gpu -pme gpu -bonded gpu \
    -nt 8 -ntmpi 1 -ntomp 8 \
    2>&1 | tee nvt.log | tail -20

echo ""
echo "[7/7] NPT production (10ns)..."
cat > npt_10ns.mdp <<'MDP'
title                   = NPT production 10ns
integrator              = md
nsteps                  = 5000000   ; 10 ns
dt                      = 0.002     ; 2 fs
nstxout                 = 5000      ; save coords every 10 ps
nstvout                 = 5000
nstfout                 = 0
nstenergy               = 500       ; save energy every 1 ps
nstlog                  = 500
nstxout-compressed      = 5000      ; compressed trajectory every 10 ps
compressed-x-grps       = System
continuation            = yes
constraint_algorithm    = lincs
constraints             = h-bonds
cutoff-scheme           = Verlet
ns_type                 = grid
nstlist                 = 10
rcoulomb                = 1.0
rvdw                    = 1.0
coulombtype             = PME
pme_order               = 4
fourierspacing          = 0.16
tcoupl                  = V-rescale
tc-grps                 = Protein Non-Protein
tau_t                   = 0.1     0.1
ref_t                   = 300     300
pcoupl                  = Parrinello-Rahman
pcoupltype              = isotropic
tau_p                   = 2.0
ref_p                   = 1.0
compressibility         = 4.5e-5
pbc                     = xyz
gen_vel                 = no
MDP

gmx grompp \
    -f npt_10ns.mdp \
    -c nvt.gro \
    -t nvt.cpt \
    -p topol.top \
    -o npt_10ns.tpr \
    -maxwarn 2 \
    2>&1 | tail -10

echo ""
echo "========================================="
echo "Starting 10ns NPT production MD on GPU"
echo "========================================="
echo "This will take several minutes..."
echo ""

START_TIME=$(date +%s)

gmx mdrun \
    -v -deffnm npt_10ns \
    -nb gpu -pme gpu -bonded gpu \
    -nt 8 -ntmpi 1 -ntomp 8 \
    2>&1 | tee npt_10ns.log

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
ELAPSED_MIN=$((ELAPSED / 60))
ELAPSED_SEC=$((ELAPSED % 60))

echo ""
echo "========================================="
echo "10ns MD Simulation Complete!"
echo "========================================="
echo "Elapsed time: ${ELAPSED_MIN}m ${ELAPSED_SEC}s"
echo ""

# Performance analysis
if [[ -f "npt_10ns.log" ]]; then
    NS_PER_DAY=$(grep "Performance:" npt_10ns.log | tail -1 | awk '{print $2}')
    echo "Performance: ${NS_PER_DAY} ns/day"
    echo ""
fi

# Extract final energy
if command -v gmx &> /dev/null; then
    echo "Extracting energy..."
    echo "Potential" | gmx energy -f npt_10ns.edr -o energy.xvg 2>&1 | grep "^Potential" || true
    
    echo ""
    echo "Calculating RMSD..."
    echo "Backbone Backbone" | gmx rms \
        -s npt_10ns.tpr \
        -f npt_10ns.xtc \
        -o rmsd.xvg \
        -tu ns \
        2>&1 | tail -10 || true
fi

echo ""
echo "Output files:"
ls -lh *.{gro,edr,xtc,log,xvg} 2>/dev/null | awk '{print $9, $5}'

echo ""
echo "========================================="
echo "GPU MD Test Summary"
echo "========================================="
echo "System: $(grep 'atoms' topol.top 2>/dev/null | head -1 || echo 'N/A')"
echo "Simulation time: 10 ns"
echo "Timestep: 2 fs"
echo "Total steps: 5,000,000"
echo "Wall time: ${ELAPSED_MIN}m ${ELAPSED_SEC}s"
if [[ -n "${NS_PER_DAY:-}" ]]; then
    echo "Performance: ${NS_PER_DAY} ns/day"
fi
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)"
echo "========================================="

# Save summary JSON
cat > md_10ns_summary.json <<JSON
{
  "simulation_type": "NPT_production",
  "duration_ns": 10.0,
  "timestep_fs": 2.0,
  "total_steps": 5000000,
  "wall_time_seconds": ${ELAPSED},
  "performance_ns_per_day": "${NS_PER_DAY:-N/A}",
  "gpu_used": "$(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)",
  "force_field": "amber99sb-ildn",
  "water_model": "tip3p",
  "temperature_K": 300,
  "pressure_bar": 1.0,
  "status": "completed"
}
JSON

echo ""
echo "Summary saved to: md_10ns_summary.json"
echo ""
