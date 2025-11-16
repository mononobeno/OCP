#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DB="$ROOT_DIR/catalog/db/ocp_results.sqlite"

echo "[INFO] DB = $DB"

if [ ! -f "$DB" ]; then
  echo "ERROR: DB not found: $DB" >&2
  exit 1
fi

add_col_if_missing () {
  local table="$1"
  local col="$2"
  local type="$3"

  local exists
  exists="$(sqlite3 "$DB" "SELECT name FROM pragma_table_info('$table') WHERE name = '$col';")"
  if [ -z "$exists" ]; then
    echo "[INFO] adding column $table.$col ($type)"
    sqlite3 "$DB" "ALTER TABLE $table ADD COLUMN $col $type;"
  else
    echo "[INFO] column already exists: $table.$col"
  fi
}

# 既に追加済みでも安全
add_col_if_missing "ligands" "smiles"           "TEXT"
add_col_if_missing "ligands" "source_file"      "TEXT"
add_col_if_missing "ligands" "has_3d"           "INTEGER DEFAULT 0"
add_col_if_missing "ligands" "conformer_method" "TEXT"

echo "[INFO] migration done."
