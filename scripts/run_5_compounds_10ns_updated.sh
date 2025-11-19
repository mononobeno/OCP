#!/bin/bash
#
# 5化合物の10ns MD計算スクリプト(改良版)
# - リガンドRMSD: バックボーンを基準にリガンドの動きを計算
# - 構造キャプチャ: 初期構造とMD終了時の構造を保存
#

set -e

COMPOUNDS=(
  "ZINC001241748223_1"
  "ZINC001241748302_1"
  "ZINC001241750803_1"
  "ZINC001241752370_1"
  "ZINC001241748464_1"
)

TARGET="prothrombin"
RECEPTOR_PDB="/home/dev/OCP/catalog/targets/aps_prothrombin/receptor.pdb"

echo "=========================================="
echo "10ns MD Simulation for 5 Compounds"
echo "=========================================="
echo "Target: Prothrombin (F2)"
echo "Compounds: ${#COMPOUNDS[@]}"
echo "MD Duration: 10 ns"
echo "Start time: $(date '+%Y-%m-%d %H:%M:%S JST')"
echo "=========================================="
echo ""

START_TIME=$(date +%s)

for i in "${!COMPOUNDS[@]}"; do
    COMPOUND="${COMPOUNDS[$i]}"
    NUM=$((i + 1))
    
    echo ""
    echo "=========================================="
    echo "Processing compound $NUM/${#COMPOUNDS[@]}: $COMPOUND"
    echo "=========================================="
    
    # Paths
    VINA_DIR="results/vina_output/${COMPOUND}_${TARGET}"
    PREP_DIR="results/gromacs_prep/${COMPOUND}_${TARGET}"
    MD_DIR="results/md_output/${COMPOUND}_${TARGET}_10ns_v2"
    
    if [ ! -d "$VINA_DIR" ]; then
        echo "⚠️  Vina output not found: $VINA_DIR"
        continue
    fi
    
    mkdir -p "$MD_DIR"
    cd "$MD_DIR"
    
    # Copy docked ligand
    if [ ! -f "ligand.pdbqt" ]; then
        cp "../../../$VINA_DIR/docked.pdbqt" ligand.pdbqt
    fi
    
    # Convert PDBQT to PDB
    if [ ! -f "ligand.pdb" ]; then
        grep -v "^ROOT\|^ENDROOT\|^TORSDOF\|^BRANCH\|^ENDBRANCH" ligand.pdbqt | \
        sed 's/^ATOM/HETATM/' > ligand.pdb
    fi
    
    # Prepare protein
    if [ ! -f "protein.pdb" ]; then
        grep "^ATOM" "$RECEPTOR_PDB" > protein.pdb
    fi
    
    # Create complex (initial structure for capture)
    if [ ! -f "complex_initial.pdb" ]; then
        cat protein.pdb ligand.pdb > complex_initial.pdb
    fi
    
    echo "Running GROMACS preparation..."
    
    # GROMACS topology generation with Docker
    docker run --rm -v "$PWD:/work" -w /work \
        gromacs-mpi-cuda:local bash -c '
        # Process protein
        echo "1" | gmx_mpi pdb2gmx -f protein.pdb -o processed.gro -p topol.top -ff oplsaa -water spce -ignh 2>&1 | grep -v "^$"
        
        # Add ligand coordinates
        grep "^HETATM" ligand.pdb | awk "{printf \"%-6s%5d %-4s %3s %c%4d    %8.3f%8.3f%8.3f\n\", \"ATOM\", NR, \$3, \"LIG\", \"A\", 1, \$6, \$7, \$8}" > ligand.gro
        
        # Combine protein and ligand
        head -n -1 processed.gro > temp.gro
        cat ligand.gro >> temp.gro
        NATOMS=$(grep -c "^" ligand.gro)
        PROTEIN_ATOMS=$(sed -n "2p" processed.gro)
        TOTAL=$((PROTEIN_ATOMS + NATOMS))
        sed -i "2s/.*/$TOTAL/" temp.gro
        tail -n 1 processed.gro >> temp.gro
        mv temp.gro complex.gro
        
        # Add box
        gmx_mpi editconf -f complex.gro -o boxed.gro -c -d 1.0 -bt cubic 2>&1 | grep -v "^$"
        
        # Solvate
        gmx_mpi solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top 2>&1 | grep -v "^$"
        
        # Add ions
        cat > ions.mdp << EOF
integrator  = steep
nsteps      = 500
emtol       = 1000.0
EOF
        
        gmx_mpi grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn 10 2>&1 | grep -v "^$"
        echo "SOL" | gmx_mpi genion -s ions.tpr -o ions.gro -p topol.top -pname NA -nname CL -neutral 2>&1 | grep -v "^$"
        
        # Energy minimization
        cat > em.mdp << EOF
integrator  = steep
nsteps      = 5000
emtol       = 1000.0
nstlist     = 10
cutoff-scheme = Verlet
ns_type     = grid
rlist       = 1.2
coulombtype = PME
rcoulomb    = 1.2
vdwtype     = cutoff
rvdw        = 1.2
pbc         = xyz
EOF
        
        gmx_mpi grompp -f em.mdp -c ions.gro -p topol.top -o em.tpr -maxwarn 10 2>&1 | grep -v "^$"
        gmx_mpi mdrun -v -deffnm em -nb gpu -pme gpu -bonded gpu -ntomp 8 2>&1 | grep -E "step|Potential|finished"
        
        # NVT equilibration
        cat > nvt.mdp << EOF
integrator  = md
nsteps      = 50000
dt          = 0.002
nstxout     = 0
nstvout     = 0
nstlog      = 5000
nstenergy   = 1000
nstxout-compressed = 0
continuation = no
gen_vel     = yes
gen_temp    = 300
gen_seed    = -1
tcoupl      = V-rescale
tc-grps     = System
tau_t       = 0.1
ref_t       = 300
pcoupl      = no
cutoff-scheme = Verlet
nstlist     = 10
ns_type     = grid
rlist       = 1.2
coulombtype = PME
rcoulomb    = 1.2
vdwtype     = cutoff
rvdw        = 1.2
pbc         = xyz
DispCorr    = EnerPres
EOF
        
        gmx_mpi grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 10 2>&1 | grep -v "^$"
        gmx_mpi mdrun -v -deffnm nvt -nb gpu -pme gpu -bonded gpu -ntomp 8 2>&1 | grep -E "step|finished"
        
        # NPT equilibration
        cat > npt.mdp << EOF
integrator  = md
nsteps      = 50000
dt          = 0.002
nstxout     = 0
nstvout     = 0
nstlog      = 5000
nstenergy   = 1000
nstxout-compressed = 0
continuation = yes
gen_vel     = no
tcoupl      = V-rescale
tc-grps     = System
tau_t       = 0.1
ref_t       = 300
pcoupl      = Parrinello-Rahman
pcoupltype  = isotropic
tau_p       = 2.0
ref_p       = 1.0
compressibility = 4.5e-5
refcoord_scaling = com
cutoff-scheme = Verlet
nstlist     = 10
ns_type     = grid
rlist       = 1.2
coulombtype = PME
rcoulomb    = 1.2
vdwtype     = cutoff
rvdw        = 1.2
pbc         = xyz
DispCorr    = EnerPres
EOF
        
        gmx_mpi grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr -maxwarn 10 2>&1 | grep -v "^$"
        gmx_mpi mdrun -v -deffnm npt -nb gpu -pme gpu -bonded gpu -ntomp 8 2>&1 | grep -E "step|finished"
    '
    
    echo "Running 10ns MD simulation..."
    
    # Production MD (10ns)
    docker run --rm --gpus all -v "$PWD:/work" -w /work \
        gromacs-mpi-cuda:local bash -c '
        cat > md_10ns.mdp << EOF
integrator  = md
nsteps      = 5000000
dt          = 0.002
nstxout     = 0
nstvout     = 0
nstlog      = 5000
nstenergy   = 5000
nstxout-compressed = 5000
compressed-x-grps = System
continuation = yes
gen_vel     = no
tcoupl      = V-rescale
tc-grps     = System
tau_t       = 0.1
ref_t       = 300
pcoupl      = Parrinello-Rahman
pcoupltype  = isotropic
tau_p       = 2.0
ref_p       = 1.0
compressibility = 4.5e-5
cutoff-scheme = Verlet
nstlist     = 10
ns_type     = grid
rlist       = 1.2
coulombtype = PME
rcoulomb    = 1.2
vdwtype     = cutoff
rvdw        = 1.2
pbc         = xyz
DispCorr    = EnerPres
EOF
        
        gmx_mpi grompp -f md_10ns.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr -maxwarn 10 2>&1 | grep -v "^$"
        gmx_mpi mdrun -v -s md.tpr -deffnm md -nb gpu -pme gpu -bonded gpu -ntomp 8 2>&1 | tee md.log
    '
    
    # Extract final structure for capture
    docker run --rm -v "$PWD:/work" -w /work \
        gromacs-mpi-cuda:local bash -c '
        echo "System" | gmx_mpi trjconv -s md.tpr -f md.xtc -o complex_final.pdb -dump 10000 -pbc mol 2>&1 | grep -v "^$"
    '
    
    echo "Running trajectory analysis..."
    
    # Enhanced RMSD analysis: Ligand RMSD relative to Backbone
    docker run --rm -v "$PWD:/work" -w /work \
        gromacs-mpi-cuda:local bash -c '
        # Create index groups
        cat > index_commands.txt << EOF
a 1-2000
name 20 Protein_Backbone
a 2001-2050
name 21 Ligand
q
EOF
        
        gmx_mpi make_ndx -f md.tpr -o index.ndx < index_commands.txt 2>&1 | grep -v "^$"
        
        # RMSD: Backbone as reference, Ligand as target
        echo -e "Protein_Backbone\nLigand" | gmx_mpi rms -s md.tpr -f md.xtc -o rmsd_ligand_vs_backbone.xvg -n index.ndx 2>&1 | grep -v "^$"
        
        # Traditional Backbone RMSD
        echo "Backbone Backbone" | gmx_mpi rms -s md.tpr -f md.xtc -o rmsd_backbone.xvg 2>&1 | grep -v "^$"
        
        # C-alpha RMSD
        echo "C-alpha C-alpha" | gmx_mpi rms -s md.tpr -f md.xtc -o rmsd_calpha.xvg 2>&1 | grep -v "^$"
        
        # RMSF
        echo "C-alpha" | gmx_mpi rmsf -s md.tpr -f md.xtc -o rmsf.xvg -res 2>&1 | grep -v "^$"
        
        # Radius of gyration
        echo "Protein" | gmx_mpi gyrate -s md.tpr -f md.xtc -o gyrate.xvg 2>&1 | grep -v "^$"
        
        # Hydrogen bonds
        echo -e "Protein\nProtein" | gmx_mpi hbond -s md.tpr -f md.xtc -num hbond.xvg 2>&1 | grep -v "^$"
    '
    
    # Calculate average RMSD
    RMSD_AVG=$(tail -n +25 rmsd_ligand_vs_backbone.xvg | awk '{sum+=$2; n++} END {if(n>0) printf "%.4f", sum/n; else print "0"}')
    
    # Extract performance
    PERF=$(grep "Performance:" md.log | tail -1 | awk '{print $2}')
    
    echo ""
    echo "✅ MD simulation completed"
    echo "Average Ligand RMSD (vs Backbone): $RMSD_AVG nm"
    echo "Performance: $PERF ns/day"
    
    # Update database
    cd /home/dev/OCP
    sqlite3 catalog/db/ocp_results.sqlite << EOF
UPDATE vina_results 
SET md_status = 'completed_10ns_v2',
    md_output_dir = '$MD_DIR',
    md_rmsd_avg = $RMSD_AVG,
    md_simulation_time_ns = 10.0,
    md_performance_nsday = $PERF
WHERE pose_file LIKE '%${COMPOUND}%';
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
echo "Results in: results/md_output/*_10ns_v2/"
echo ""
echo "=========================================="
echo "Next step: Generate GitHub Pages with:"
echo "  bash scripts/generate_updated_pages.sh"
echo "=========================================="
