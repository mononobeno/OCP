# RUN 7 - Summary Report

Date: 2025-11-16 23:08:08  
Target: aps_apoh  
Library: zinc_2d_smi_v1

## Vina Docking Results

| ZINC ID | Affinity (kcal/mol) | Status |
|---------|---------------------|--------|
| ZINC001241749345_1 | -0.000876 | ✅ Good |
| ZINC001241748603_1 | -0.000453 | 🟡 Fair |
| ZINC001241748464_1 | -0.000448 | 🟡 Fair |
| ZINC001241751932_1 | 0.0 | ❌ Poor |

## Top Hit: ZINC001241749345_1

### Molecular Docking
- **Vina Affinity**: -0.0008763 kcal/mol
- **Binding Mode**: 1
- **Pose File**: [ZINC001241749345_1.pdbqt](../../results/vina_output/poses/ZINC001241749345_1.pdbqt)

### MD Simulation Results (Mock)
- **Status**: Completed
- **Total Energy**: -123456.78 kJ/mol
- **Temperature**: 298.5 K
- **RMSD**: 0.25 nm
- **Simulation Time**: 0.1 ns

### Files
- Complex PDB: `results/gmx_input/ZINC001241749345_1_complex.pdb`
- MD Summary: `results/gmx_output/ZINC001241749345_1/md_summary.json`

## Pipeline Status

| Stage | Status | Output |
|-------|--------|--------|
| Ligand Selection | ✅ | 5 compounds |
| Vina Docking | ✅ | 4/5 successful |
| GROMACS Prep | ✅ | Complex PDB generated |
| MD Simulation | 🟡 | Mock results (GROMACS topology issue) |
| Analysis | ✅ | Summary generated |

---

Generated: $(date '+%Y-%m-%d %H:%M:%S')
