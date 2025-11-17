---
layout: default
title: OCP Pipeline
---

# Open Compound Pipeline

## Pipeline Architecture

```
ZINC20 → RDKit 3D → Vina Docking → GROMACS Prep → MD Simulation → Analysis
```

## Kubernetes Jobs

| Job | Status | Description |
|-----|--------|-------------|
| vina-prep | ✅ | Vina前処理 |
| vina-runner | ✅ | GPU Docking |
| md-production | ✅ | GPU MD |
| md-analysis | ✅ | 軌跡解析 |

## Targets

| PDB | Name | Vina Results | MD Done | Link |
|-----|------|--------------|---------|------|
| 1C1Z | β2-glycoprotein I (APOH) | 16 | 1 | [View](./targets/target_1.html) |
| 6C2W | Prothrombin (F2) | 0 | 0 | [View](./targets/target_2.html) |
| 3FXI | TLR4/MD-2 complex | 0 | 0 | [View](./targets/target_3.html) |

_Updated: 2025-11-17 11:05:52_
