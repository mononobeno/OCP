#!/bin/bash
set -euo pipefail

# ========================================
# DB Migration: Add GROMACS prep fields
# ========================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

DB_PATH="${ROOT}/catalog/db/ocp_results.sqlite"

echo "========================================="
echo "DB Migration: GROMACS Preparation Fields"
echo "========================================="
echo "Database: ${DB_PATH}"
echo ""

# Check current schema
echo "[1/3] Checking current vina_results schema..."
sqlite3 "${DB_PATH}" "PRAGMA table_info(vina_results);" | head -10

# Add new columns if they don't exist
echo ""
echo "[2/3] Adding GROMACS preparation columns..."

sqlite3 "${DB_PATH}" <<'EOF'
-- Add columns for GROMACS preparation tracking
ALTER TABLE vina_results ADD COLUMN pose_file TEXT;
ALTER TABLE vina_results ADD COLUMN affinity REAL;
ALTER TABLE vina_results ADD COLUMN gromacs_prep_status TEXT DEFAULT 'pending';
ALTER TABLE vina_results ADD COLUMN gromacs_prep_dir TEXT;
ALTER TABLE vina_results ADD COLUMN md_status TEXT DEFAULT 'pending';
ALTER TABLE vina_results ADD COLUMN md_output_dir TEXT;
EOF

echo "  ✓ Added columns: pose_file, affinity, gromacs_prep_status, gromacs_prep_dir, md_status, md_output_dir"

# Verify new schema
echo ""
echo "[3/3] Verifying updated schema..."
sqlite3 "${DB_PATH}" "PRAGMA table_info(vina_results);"

echo ""
echo "========================================="
echo "Migration Complete!"
echo "========================================="
echo ""
echo "New fields in vina_results table:"
echo "  - pose_file           : Path to best Vina pose PDBQT"
echo "  - affinity            : Binding affinity (kcal/mol)"
echo "  - gromacs_prep_status : Status of GROMACS preparation (pending/completed/failed)"
echo "  - gromacs_prep_dir    : Directory with GROMACS input files"
echo "  - md_status           : Status of MD simulation (pending/running/completed/failed)"
echo "  - md_output_dir       : Directory with MD results"
echo "========================================="
