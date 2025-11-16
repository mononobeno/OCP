#!/usr/bin/env bash
set -eo pipefail

###############################################################################
# update_smiles_in_db.sh - DB内のligandテーブルにSMILES情報を一括更新
###############################################################################

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

SMI_FILE="catalog/libraries/zinc_2d_smi_v1/raw/all_zinc20_ml_subset.smi"
DB_PATH="catalog/db/ocp_results.sqlite"
PDBQT_DIR="catalog/libraries/zinc_2d_smi_v1/processed/pdbqt"

LIB_ID=$(sqlite3 "$DB_PATH" "SELECT id FROM libraries WHERE code='zinc_2d_smi_v1';")

echo "=========================================="
echo "Update SMILES & has_3d in Database"
echo "=========================================="

count=0
updated=0
inserted=0

while IFS=' ' read -r smiles zinc_id; do
    ((count++))
    
    if (( count % 100 == 0 )); then
        echo "[$count] Updated: $updated, Inserted: $inserted"
    fi
    
    # SMILESエスケープ
    smiles_esc="${smiles//\'/\'\'}"
    
    # PDBQTファイルの存在確認
    pdbqt_file="$PDBQT_DIR/${zinc_id}.pdbqt"
    has_3d=0
    conformer_method=""
    
    if [[ -f "$pdbqt_file" ]] && grep -q "ATOM.*[1-9]" "$pdbqt_file" 2>/dev/null; then
        has_3d=1
        conformer_method="obabel_gen3d"
    fi
    
    # DB更新または挿入
    existing=$(sqlite3 "$DB_PATH" "SELECT id FROM ligands WHERE zinc_id='$zinc_id' AND library_id=$LIB_ID;" 2>/dev/null || echo "")
    
    if [[ -n "$existing" ]]; then
        sqlite3 "$DB_PATH" "UPDATE ligands SET smiles='$smiles_esc', has_3d=$has_3d, conformer_method='$conformer_method' WHERE id=$existing;" 2>/dev/null || true
        ((updated++))
    else
        sqlite3 "$DB_PATH" "INSERT INTO ligands (zinc_id, library_id, smiles, has_3d, conformer_method) VALUES ('$zinc_id', $LIB_ID, '$smiles_esc', $has_3d, '$conformer_method');" 2>/dev/null || true
        ((inserted++))
    fi
    
done < "$SMI_FILE"

echo ""
echo "=========================================="
echo "Update Complete!"
echo "=========================================="
echo "Total processed: $count"
echo "Updated:  $updated"
echo "Inserted: $inserted"
echo "=========================================="

# DB統計
sqlite3 "$DB_PATH" <<EOF
.mode column
.headers on
SELECT 
    COUNT(*) as total_ligands,
    SUM(CASE WHEN has_3d=1 THEN 1 ELSE 0 END) as with_3d,
    SUM(CASE WHEN smiles IS NOT NULL AND smiles!='' THEN 1 ELSE 0 END) as with_smiles
FROM ligands WHERE library_id=$LIB_ID;
EOF

echo ""
echo "[OK] Database update complete!"
