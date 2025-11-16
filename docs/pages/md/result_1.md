---
layout: default
title: MD Result - ZINC00000001
---

# MD Simulation Result: ZINC00000001

[← Back to MD Results](./index.html)

## 基本情報

- **ZINC ID**: ZINC00000001
- **SMILES**: ``
- **Vina Affinity**: -1.499e-05 kcal/mol

## MD Simulation

- **Status**: ✅ Completed
- **Output Directory**: `/home/dev/OCP/results/md_gpu_test_1ns`
- **RMSD平均**: **0.0702 nm**
- **シミュレーション時間**: 1.0 ns
- **計算性能**: 408.7 ns/day

## 軌跡解析

```bash
# RMSD plot
gmx rms -s /home/dev/OCP/results/md_gpu_test_1ns/md_*.tpr -f /home/dev/OCP/results/md_gpu_test_1ns/md_*.xtc -o rmsd.xvg

# RMSF plot
gmx rmsf -s /home/dev/OCP/results/md_gpu_test_1ns/md_*.tpr -f /home/dev/OCP/results/md_gpu_test_1ns/md_*.xtc -o rmsf.xvg

# Radius of gyration
gmx gyrate -s /home/dev/OCP/results/md_gpu_test_1ns/md_*.tpr -f /home/dev/OCP/results/md_gpu_test_1ns/md_*.xtc -o gyrate.xvg
```

