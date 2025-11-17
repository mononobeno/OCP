#!/bin/bash
set -euo pipefail

# ========================================
# Simple 10ns GPU MD Test with Lysozyme
# ========================================

ROOT="/home/dev/OCP"
TEST_DIR="${ROOT}/results/md_gpu_test_simple"
mkdir -p "${TEST_DIR}"
cd "${TEST_DIR}"

echo "========================================="
echo "10ns GPU MD Simulation Test (Lysozyme)"
echo "========================================="
echo ""

# Download lysozyme PDB
echo "[1/8] Downloading test system (Lysozyme)..."
wget -q http://files.rcsb.org/download/1AKI.pdb -O 1aki.pdb || {
    echo "Download failed, using local receptor"
    # Use first 100 residues of existing receptor as backup
    head -1000 "${ROOT}/catalog/targets/aps_apoh/receptor.pdb" | grep "^ATOM" | head -800 > 1aki.pdb
    echo "END" >> 1aki.pdb
}

echo "  ✓ System prepared ($(grep -c '^ATOM' 1aki.pdb || echo 0) atoms)"

echo ""
echo "[2/8] System preparation (pdb2gmx)..."
echo "15" | timeout 60 gmx pdb2gmx \
    -f 1aki.pdb \
    -o processed.gro \
    -p topol.top \
    -i posre.itp \
    -water spce \
    -ff oplsaa \
    -ignh \
    2>&1 | tail -10

if [[ ! -f "processed.gro" ]]; then
    echo "ERROR: pdb2gmx failed"
    exit 1
fi

echo "  ✓ System processed"

echo ""
echo "[3/8] Creating box..."
gmx editconf -f processed.gro -o boxed.gro -c -d 1.0 -bt cubic 2>&1 | tail -5
echo "  ✓ Box created"

echo ""
echo "[4/8] Solvating..."
gmx solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top 2>&1 | tail -5
echo "  ✓ Solvated"

echo ""
echo "[5/8] Energy minimization..."
cat > em.mdp <<'MDP'
integrator  = steep
emtol       = 1000.0
nsteps      = 5000
nstlist     = 10
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb    = 1.0
rvdw        = 1.0
MDP

gmx grompp -f em.mdp -c solvated.gro -p topol.top -o em.tpr -maxwarn 2 2>&1 | tail -5
gmx mdrun -v -deffnm em -nb gpu -pme gpu -nt 4 2>&1 | tail -10
echo "  ✓ Energy minimized"

echo ""
echo "[6/8] NVT equilibration (100ps)..."
cat > nvt.mdp <<'MDP'
integrator              = md
nsteps                  = 50000
dt                      = 0.002
nstxout                 = 5000
nstvout                 = 5000
nstenergy               = 500
nstlog                  = 500
continuation            = no
constraint_algorithm    = lincs
constraints             = h-bonds
cutoff-scheme           = Verlet
nstlist                 = 10
rcoulomb                = 1.0
rvdw                    = 1.0
coulombtype             = PME
tcoupl                  = V-rescale
tc-grps                 = System
tau_t                   = 0.1
ref_t                   = 300
pcoupl                  = no
pbc                     = xyz
gen_vel                 = yes
gen_temp                = 300
MDP

gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 2 2>&1 | tail -5
gmx mdrun -v -deffnm nvt -nb gpu -pme gpu -nt 4 2>&1 | tail -10
echo "  ✓ NVT done (100ps)"

echo ""
echo "[7/8] NPT equilibration (100ps)..."
cat > npt.mdp <<'MDP'
integrator              = md
nsteps                  = 50000
dt                      = 0.002
nstxout                 = 5000
nstenergy               = 500
nstlog                  = 500
continuation            = yes
constraint_algorithm    = lincs
constraints             = h-bonds
cutoff-scheme           = Verlet
nstlist                 = 10
rcoulomb                = 1.0
rvdw                    = 1.0
coulombtype             = PME
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
MDP

gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 2 2>&1 | tail -5
gmx mdrun -v -deffnm npt -nb gpu -pme gpu -nt 4 2>&1 | tail -10
echo "  ✓ NPT done (100ps)"

echo ""
echo "[8/8] Production MD (10ns)..."
cat > md_10ns.mdp <<'MDP'
integrator              = md
nsteps                  = 5000000
dt                      = 0.002
nstxout                 = 50000
nstenergy               = 5000
nstlog                  = 5000
nstxout-compressed      = 5000
compressed-x-grps       = System
continuation            = yes
constraint_algorithm    = lincs
constraints             = h-bonds
cutoff-scheme           = Verlet
nstlist                 = 10
rcoulomb                = 1.0
rvdw                    = 1.0
coulombtype             = PME
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
MDP

gmx grompp -f md_10ns.mdp -c npt.gro -t npt.cpt -p topol.top -o md_10ns.tpr -maxwarn 2 2>&1 | tail -5

echo ""
echo "========================================="
echo "Starting 10ns Production MD on GPU"
echo "========================================="
START=$(date +%s)

gmx mdrun -v -deffnm md_10ns -nb gpu -pme gpu -bonded gpu -nt 8 -ntmpi 1 -ntomp 8 2>&1 | tee md_10ns.log

END=$(date +%s)
ELAPSED=$((END - START))

echo ""
echo "========================================="
echo "10ns MD Complete!"
echo "========================================="
echo "Wall time: $((ELAPSED/60))m $((ELAPSED%60))s"

# Extract performance
PERF=$(grep "Performance:" md_10ns.log | tail -1 | awk '{print $2}')
echo "Performance: ${PERF} ns/day"

# Save summary
cat > summary.json <<JSON
{
  "simulation": "10ns_production",
  "system": "$(basename $(pwd))",
  "steps": 5000000,
  "duration_ns": 10.0,
  "wall_time_s": ${ELAPSED},
  "performance_ns_day": "${PERF}",
  "gpu": "$(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)",
  "status": "completed"
}
JSON

echo ""
cat summary.json
echo ""
echo "Files:"
ls -lh *.{gro,xtc,edr,log} 2>/dev/null | awk '{print $9, $5}'
