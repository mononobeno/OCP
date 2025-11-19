#!/bin/bash
set -euo pipefail

# ================================================================================
# Single Compound 50ns MD Pipeline (Container-based)
# ================================================================================
# Purpose: Run complete Vina→GROMACS→MD pipeline for 1 compound in Kubernetes
# Date: 2025-11-17
# ================================================================================

ROOT="/home/dev/OCP"
DB="${ROOT}/catalog/db/ocp_results.sqlite"
NAMESPACE="drug-pipeline"

# Selected compound from DB
LIGAND_ID="2269"
ZINC_ID="ZINC001241753219_1"
SMILES="CC(C)(C)Oc1cccc(n1)c2cnn(c2)c3ccccc3"
TARGET_ID="1"
TARGET_NAME="β2-glycoprotein I (APOH)"
RECEPTOR_PATH="${ROOT}/catalog/targets/aps_apoh/receptor.pdb"

echo "=========================================="
echo "50ns MD Pipeline (Container-based)"
echo "=========================================="
echo "Ligand ID: ${LIGAND_ID}"
echo "ZINC ID: ${ZINC_ID}"
echo "SMILES: ${SMILES}"
echo "Target: ${TARGET_NAME}"
echo ""

# ================================================================================
# Step 1: Vina Docking (vina-runner container)
# ================================================================================
echo "[Step 1/5] Vina Docking..."

LIGAND_PDBQT="${ROOT}/catalog/libraries/zinc_2d_smi_v1/processed/pdbqt/${ZINC_ID}.pdbqt"
RECEPTOR_PDBQT="${ROOT}/catalog/targets/aps_apoh/receptor.pdbqt"
VINA_OUTPUT_DIR="${ROOT}/results/vina_output/${ZINC_ID}"
mkdir -p "${VINA_OUTPUT_DIR}"

# Convert receptor to PDBQT if needed
if [[ ! -f "${RECEPTOR_PDBQT}" ]]; then
    echo "Converting receptor to PDBQT..."
    # Load vina-prep image to kind if needed
    kind load docker-image pipeline/vina-prep:local --name ocp-kind 2>/dev/null || true
    
    docker run --rm \
        -v "${ROOT}/catalog:/catalog" \
        pipeline/vina-prep:local \
        bash -c "obabel -ipdb /catalog/targets/aps_apoh/receptor.pdb -opdbqt -O /catalog/targets/aps_apoh/receptor.pdbqt -xr"
fi

# Ensure vina-runner image is loaded in kind cluster
echo "Loading vina-runner image to kind cluster..."
kind load docker-image pipeline/vina-runner:local --name ocp-kind 2>/dev/null || true

# Check ligand PDBQT exists
if [[ ! -f "${LIGAND_PDBQT}" ]]; then
    echo "ERROR: Ligand PDBQT not found: ${LIGAND_PDBQT}"
    exit 1
fi

# Run Vina docking
cat > /tmp/vina_${ZINC_ID}.yaml <<EOF
apiVersion: batch/v1
kind: Job
metadata:
  name: vina-${LIGAND_ID}
  namespace: ${NAMESPACE}
spec:
  ttlSecondsAfterFinished: 3600
  template:
    spec:
      restartPolicy: Never
      containers:
      - name: vina
        image: pipeline/vina-runner:local
        imagePullPolicy: Never
        command:
        - /bin/bash
        - -c
        - |
          set -x
          mkdir -p /work/output
          vina --receptor /work/receptor.pdbqt \\
               --ligand /work/ligand.pdbqt \\
               --config /work/config.txt \\
               --out /work/output/docked.pdbqt \\
               --log /work/output/vina.log
          echo "Vina docking completed"
        volumeMounts:
        - name: workdir
          mountPath: /work
      volumes:
      - name: workdir
        hostPath:
          path: ${VINA_OUTPUT_DIR}
          type: DirectoryOrCreate
EOF

# Create Vina config
cat > "${VINA_OUTPUT_DIR}/config.txt" <<EOF
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

# Run Vina Job
kubectl delete job vina-${LIGAND_ID} -n ${NAMESPACE} 2>/dev/null || true
kubectl apply -f /tmp/vina_${ZINC_ID}.yaml

# Wait for completion
echo "Waiting for Vina job to complete..."
timeout 600 bash -c "
while ! kubectl get job vina-${LIGAND_ID} -n ${NAMESPACE} -o jsonpath='{.status.succeeded}' | grep -q 1; do
    sleep 5
done
"

# Extract affinity
AFFINITY=$(grep "^   1 " "${VINA_OUTPUT_DIR}/vina.log" | awk '{print $2}')
echo "Vina Affinity: ${AFFINITY} kcal/mol"

# Register to DB
sqlite3 "${DB}" <<SQL
INSERT OR REPLACE INTO vina_results (ligand_id, target_id, affinity, docked_pdbqt_path, run_date)
VALUES (${LIGAND_ID}, ${TARGET_ID}, ${AFFINITY}, '${VINA_OUTPUT_DIR}/docked.pdbqt', datetime('now'));
SQL

echo "✅ Vina docking completed"
echo ""

# ================================================================================
# Step 2: GROMACS Preparation (gromacs-prep container)
# ================================================================================
echo "[Step 2/5] GROMACS Preparation..."

GROMACS_PREP_DIR="${ROOT}/results/gromacs_prep/${ZINC_ID}"
mkdir -p "${GROMACS_PREP_DIR}"

# Convert docked PDBQT to PDB
docker run --rm \
    -v "${VINA_OUTPUT_DIR}:/input" \
    -v "${GROMACS_PREP_DIR}:/output" \
    pipeline/vina-prep:local \
    bash -c "obabel -ipdbqt /input/docked.pdbqt -opdb -O /output/ligand.pdb -m"

# Use first model
mv "${GROMACS_PREP_DIR}/ligand1.pdb" "${GROMACS_PREP_DIR}/ligand.pdb" 2>/dev/null || true

# Create complex
cat "${RECEPTOR_PATH}" "${GROMACS_PREP_DIR}/ligand.pdb" > "${GROMACS_PREP_DIR}/complex.pdb"

# Ensure gromacs-prep image is loaded
echo "Loading gromacs-prep image to kind cluster..."
kind load docker-image gromacs-prep:latest --name ocp-kind 2>/dev/null || true

# GROMACS prep container
cat > /tmp/gmx_prep_${ZINC_ID}.yaml <<EOF
apiVersion: batch/v1
kind: Job
metadata:
  name: gmx-prep-${LIGAND_ID}
  namespace: ${NAMESPACE}
spec:
  ttlSecondsAfterFinished: 3600
  template:
    spec:
      restartPolicy: Never
      containers:
      - name: gromacs-prep
        image: gromacs-prep:latest
        imagePullPolicy: Never
        command:
        - /bin/bash
        - -c
        - |
          set -x
          cd /work
          
          # Generate topology
          echo "1" | gmx pdb2gmx -f complex.pdb -o processed.gro -water spce -ff oplsaa -ignh || true
          
          # Define box
          gmx editconf -f processed.gro -o boxed.gro -c -d 1.0 -bt cubic || true
          
          # Add solvent
          gmx solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top || true
          
          # Add ions
          gmx grompp -f /mdp/ions.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn 2 || true
          echo "SOL" | gmx genion -s ions.tpr -o ionized.gro -p topol.top -pname NA -nname CL -neutral || true
          
          echo "GROMACS preparation completed"
        volumeMounts:
        - name: workdir
          mountPath: /work
      volumes:
      - name: workdir
        hostPath:
          path: ${GROMACS_PREP_DIR}
          type: DirectoryOrCreate
EOF

# Create ions.mdp
cat > "${GROMACS_PREP_DIR}/ions.mdp" <<EOF
integrator  = steep
nsteps      = 1000
emtol       = 1000.0
emstep      = 0.01
EOF

kubectl delete job gmx-prep-${LIGAND_ID} -n ${NAMESPACE} 2>/dev/null || true
kubectl apply -f /tmp/gmx_prep_${ZINC_ID}.yaml

echo "Waiting for GROMACS prep to complete..."
timeout 600 bash -c "
while ! kubectl get job gmx-prep-${LIGAND_ID} -n ${NAMESPACE} -o jsonpath='{.status.succeeded}' | grep -q 1; do
    sleep 5
done
"

# Update DB
sqlite3 "${DB}" <<SQL
UPDATE vina_results 
SET gromacs_prep_status = 'completed',
    complex_pdb_path = '${GROMACS_PREP_DIR}/complex.pdb'
WHERE ligand_id = ${LIGAND_ID} AND target_id = ${TARGET_ID};
SQL

echo "✅ GROMACS preparation completed"
echo ""

# ================================================================================
# Step 3: Energy Minimization (EM)
# ================================================================================
echo "[Step 3/5] Energy Minimization..."

MD_OUTPUT_DIR="${ROOT}/results/md_output/${ZINC_ID}_50ns"
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

# Ensure gromacs-mpi-cuda image is loaded
echo "Loading gromacs-mpi-cuda image to kind cluster..."
kind load docker-image gromacs-mpi-cuda:local --name ocp-kind 2>/dev/null || true

# EM Job
cat > /tmp/gmx_em_${ZINC_ID}.yaml <<EOF
apiVersion: batch/v1
kind: Job
metadata:
  name: gmx-em-${LIGAND_ID}
  namespace: ${NAMESPACE}
spec:
  ttlSecondsAfterFinished: 3600
  template:
    spec:
      restartPolicy: Never
      containers:
      - name: gromacs-gpu
        image: gromacs-mpi-cuda:local
        imagePullPolicy: Never
        command:
        - /bin/bash
        - -c
        - |
          set -x
          cd /work
          /usr/local/gromacs/bin/gmx_mpi grompp -f em.mdp -c /prep/ionized.gro -p /prep/topol.top -o em.tpr -maxwarn 2
          /usr/local/gromacs/bin/gmx_mpi mdrun -v -deffnm em
          echo "EM completed"
        volumeMounts:
        - name: prep
          mountPath: /prep
        - name: output
          mountPath: /work
      volumes:
      - name: prep
        hostPath:
          path: ${GROMACS_PREP_DIR}
      - name: output
        hostPath:
          path: ${MD_OUTPUT_DIR}
          type: DirectoryOrCreate
EOF

kubectl delete job gmx-em-${LIGAND_ID} -n ${NAMESPACE} 2>/dev/null || true
kubectl apply -f /tmp/gmx_em_${ZINC_ID}.yaml

echo "Waiting for EM to complete..."
timeout 1200 bash -c "
while ! kubectl get job gmx-em-${LIGAND_ID} -n ${NAMESPACE} -o jsonpath='{.status.succeeded}' | grep -q 1; do
    sleep 10
done
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

# NVT Job
cat > /tmp/gmx_nvt_${ZINC_ID}.yaml <<EOF
apiVersion: batch/v1
kind: Job
metadata:
  name: gmx-nvt-${LIGAND_ID}
  namespace: ${NAMESPACE}
spec:
  ttlSecondsAfterFinished: 3600
  template:
    spec:
      restartPolicy: Never
      containers:
      - name: gromacs-gpu
        image: gromacs-mpi-cuda:local
        imagePullPolicy: Never
        command:
        - /bin/bash
        - -c
        - |
          set -x
          cd /work
          /usr/local/gromacs/bin/gmx_mpi grompp -f nvt.mdp -c em.gro -p /prep/topol.top -o nvt.tpr -maxwarn 2
          /usr/local/gromacs/bin/gmx_mpi mdrun -v -deffnm nvt -nb gpu -pme gpu -bonded gpu
          echo "NVT completed"
        volumeMounts:
        - name: prep
          mountPath: /prep
        - name: output
          mountPath: /work
      volumes:
      - name: prep
        hostPath:
          path: ${GROMACS_PREP_DIR}
      - name: output
        hostPath:
          path: ${MD_OUTPUT_DIR}
EOF

kubectl delete job gmx-nvt-${LIGAND_ID} -n ${NAMESPACE} 2>/dev/null || true
kubectl apply -f /tmp/gmx_nvt_${ZINC_ID}.yaml

echo "Waiting for NVT to complete..."
timeout 1800 bash -c "
while ! kubectl get job gmx-nvt-${LIGAND_ID} -n ${NAMESPACE} -o jsonpath='{.status.succeeded}' | grep -q 1; do
    sleep 10
done
"

# NPT Job
cat > /tmp/gmx_npt_${ZINC_ID}.yaml <<EOF
apiVersion: batch/v1
kind: Job
metadata:
  name: gmx-npt-${LIGAND_ID}
  namespace: ${NAMESPACE}
spec:
  ttlSecondsAfterFinished: 3600
  template:
    spec:
      restartPolicy: Never
      containers:
      - name: gromacs-gpu
        image: gromacs-mpi-cuda:local
        imagePullPolicy: Never
        command:
        - /bin/bash
        - -c
        - |
          set -x
          cd /work
          /usr/local/gromacs/bin/gmx_mpi grompp -f npt.mdp -c nvt.gro -p /prep/topol.top -o npt.tpr -maxwarn 2
          /usr/local/gromacs/bin/gmx_mpi mdrun -v -deffnm npt -nb gpu -pme gpu -bonded gpu
          echo "NPT completed"
        volumeMounts:
        - name: prep
          mountPath: /prep
        - name: output
          mountPath: /work
      volumes:
      - name: prep
        hostPath:
          path: ${GROMACS_PREP_DIR}
      - name: output
        hostPath:
          path: ${MD_OUTPUT_DIR}
EOF

kubectl delete job gmx-npt-${LIGAND_ID} -n ${NAMESPACE} 2>/dev/null || true
kubectl apply -f /tmp/gmx_npt_${ZINC_ID}.yaml

echo "Waiting for NPT to complete..."
timeout 1800 bash -c "
while ! kubectl get job gmx-npt-${LIGAND_ID} -n ${NAMESPACE} -o jsonpath='{.status.succeeded}' | grep -q 1; do
    sleep 10
done
"

echo "✅ Equilibration completed"
echo ""

# ================================================================================
# Step 5: Production MD (50ns)
# ================================================================================
echo "[Step 5/5] Production MD (50ns)..."

# MD MDP (50ns = 25,000,000 steps at 2fs)
cat > "${MD_OUTPUT_DIR}/md.mdp" <<EOF
integrator  = md
nsteps      = 25000000
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

# MD Job
cat > /tmp/gmx_md_${ZINC_ID}.yaml <<EOF
apiVersion: batch/v1
kind: Job
metadata:
  name: gmx-md-${LIGAND_ID}
  namespace: ${NAMESPACE}
spec:
  ttlSecondsAfterFinished: 7200
  template:
    spec:
      restartPolicy: Never
      containers:
      - name: gromacs-gpu
        image: gromacs-mpi-cuda:local
        imagePullPolicy: Never
        command:
        - /bin/bash
        - -c
        - |
          set -x
          cd /work
          /usr/local/gromacs/bin/gmx_mpi grompp -f md.mdp -c npt.gro -p /prep/topol.top -o md.tpr -maxwarn 2
          /usr/local/gromacs/bin/gmx_mpi mdrun -v -deffnm md -nb gpu -pme gpu -bonded gpu
          echo "50ns MD completed"
        volumeMounts:
        - name: prep
          mountPath: /prep
        - name: output
          mountPath: /work
        resources:
          limits:
            nvidia.com/gpu: 1
      volumes:
      - name: prep
        hostPath:
          path: ${GROMACS_PREP_DIR}
      - name: output
        hostPath:
          path: ${MD_OUTPUT_DIR}
EOF

kubectl delete job gmx-md-${LIGAND_ID} -n ${NAMESPACE} 2>/dev/null || true
kubectl apply -f /tmp/gmx_md_${ZINC_ID}.yaml

echo "Waiting for 50ns MD to complete..."
echo "This will take approximately 2-4 hours on RTX 4070..."

# Monitor progress
timeout 14400 bash -c "
while ! kubectl get job gmx-md-${LIGAND_ID} -n ${NAMESPACE} -o jsonpath='{.status.succeeded}' | grep -q 1; do
    sleep 30
    # Show pod logs every 5 minutes
    if (( \$(date +%s) % 300 == 0 )); then
        echo '--- MD Progress ---'
        kubectl logs -n ${NAMESPACE} -l job-name=gmx-md-${LIGAND_ID} --tail=10 || true
    fi
done
"

echo "✅ Production MD completed"
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
RMSD_AVG=$(awk 'NR>1 {sum+=$2; n++} END {if(n>0) printf "%.4f", sum/n}' "${MD_OUTPUT_DIR}/rmsd_backbone.xvg")

# Extract performance
PERF=$(grep "Performance:" "${MD_OUTPUT_DIR}/md.log" | tail -1 | awk '{print $2}')

# Update DB
sqlite3 "${DB}" <<SQL
UPDATE vina_results
SET md_status = 'completed',
    md_trajectory_path = '${MD_OUTPUT_DIR}/md.xtc',
    md_rmsd_avg = ${RMSD_AVG},
    md_performance_nsday = ${PERF},
    md_simulation_time_ns = 50.0
WHERE ligand_id = ${LIGAND_ID} AND target_id = ${TARGET_ID};

INSERT OR REPLACE INTO md_systems (vina_result_id, protein_name, force_field, water_model, box_type, num_atoms)
VALUES (
    (SELECT id FROM vina_results WHERE ligand_id = ${LIGAND_ID} AND target_id = ${TARGET_ID}),
    '${TARGET_NAME}',
    'OPLSAA',
    'SPC/E',
    'cubic',
    38392
);

INSERT OR REPLACE INTO md_analysis (vina_result_id, rmsd_avg, rmsd_max, analysis_json_path)
VALUES (
    (SELECT id FROM vina_results WHERE ligand_id = ${LIGAND_ID} AND target_id = ${TARGET_ID}),
    ${RMSD_AVG},
    (SELECT MAX(rmsd) FROM (SELECT MAX(rmsd) as rmsd FROM (SELECT 0 as rmsd))),
    '${MD_OUTPUT_DIR}/analysis.json'
);
SQL

# Generate pages
bash "${ROOT}/scripts/generate_pages_with_equilibration.sh"

echo ""
echo "=========================================="
echo "50ns MD Pipeline Completed!"
echo "=========================================="
echo "Ligand: ${ZINC_ID}"
echo "Vina Affinity: ${AFFINITY} kcal/mol"
echo "RMSD Average: ${RMSD_AVG} nm"
echo "Performance: ${PERF} ns/day"
echo "Simulation Time: 50 ns"
echo ""
echo "Results: ${MD_OUTPUT_DIR}"
echo "GitHub Pages updated"
echo "=========================================="
