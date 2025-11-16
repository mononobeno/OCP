#!/usr/bin/env bash
cat << 'EOGI'
# === OCP big data / local-only assets ===
catalog/db/ocp_results.sqlite

catalog/libraries/*/raw/
catalog/libraries/*/processed/pdbqt/
catalog/libraries/*/processed/rdkit_3d/

# ZINC / docking heavy archives
catalog/libraries/*/raw/*.tar.gz
catalog/libraries/*/raw/*.tgz
catalog/libraries/*/raw/*.smi
EOGI
