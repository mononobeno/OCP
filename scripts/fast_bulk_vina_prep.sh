#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# fast_bulk_vina_prep.sh - 高速一括SMILES→PDBQT変換
###############################################################################

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

BATCH_SIZE="${1:-100}"
SMI_FILE="catalog/libraries/zinc_2d_smi_v1/raw/all_zinc20_ml_subset.smi"
DB_PATH="catalog/db/ocp_results.sqlite"
PDBQT_DIR="catalog/libraries/zinc_2d_smi_v1/processed/pdbqt"

mkdir -p "$PDBQT_DIR"

LIB_ID=$(sqlite3 "$DB_PATH" "SELECT id FROM libraries WHERE code='zinc_2d_smi_v1';")

echo "=========================================="
echo "Fast Bulk SMILES → PDBQT Conversion"
echo "=========================================="
echo "Processing: $SMI_FILE"
echo "Batch size: $BATCH_SIZE"
echo "=========================================="

count=0
success=0
failed=0
skipped=0

while IFS=' ' read -r smiles zinc_id; do
    ((count++))
    
    # 進捗表示
    if (( count % 50 == 0 )); then
        echo "[$count] Processed: $success success, $failed failed, $skipped skipped"
    fi
    
    # PDBQTチェック
    pdbqt_file="$PDBQT_DIR/${zinc_id}.pdbqt"
    if [[ -f "$pdbqt_file" ]] && grep -q "ATOM.*[1-9]" "$pdbqt_file" 2>/dev/null; then
        ((skipped++))
        
        # DBにSMILES登録（エスケープ処理）
        smiles_esc="${smiles//\'/\'\'}"
        sqlite3 "$DB_PATH" "UPDATE ligands SET smiles='$smiles_esc' WHERE zinc_id='$zinc_id' AND library_id=$LIB_ID;" 2>/dev/null || true
        continue
    fi
    
    # 3D生成
    if echo "$smiles" | obabel -ismi -opdbqt --gen3d -h -O "$pdbqt_file" &>/dev/null; then
        if grep -q "ATOM.*[1-9]" "$pdbqt_file" 2>/dev/null; then
            # DB登録
            smiles_esc="${smiles//\'/\'\'}"
            
            # 既存チェック
            existing=$(sqlite3 "$DB_PATH" "SELECT id FROM ligands WHERE zinc_id='$zinc_id' AND library_id=$LIB_ID;" 2>/dev/null || echo "")
            
            if [[ -n "$existing" ]]; then
                sqlite3 "$DB_PATH" "UPDATE ligands SET smiles='$smiles_esc', has_3d=1, conformer_method='obabel_gen3d' WHERE id=$existing;" 2>/dev/null || true
            else
                sqlite3 "$DB_PATH" "INSERT INTO ligands (zinc_id, library_id, smiles, has_3d, conformer_method) VALUES ('$zinc_id', $LIB_ID, '$smiles_esc', 1, 'obabel_gen3d');" 2>/dev/null || true
            fi
            
            ((success++))
        else
            ((failed++))
        fi
    else
        ((failed++))
    fi
    
done < "$SMI_FILE"

echo ""
echo "=========================================="
echo "Conversion Complete!"
echo "=========================================="
echo "Total:   $count"
echo "Success: $success"
echo "Failed:  $failed"
echo "Skipped: $skipped"
echo "=========================================="

# DB統計
sqlite3 "$DB_PATH" <<EOF
.mode column
SELECT 
    COUNT(*) as total,
    SUM(CASE WHEN has_3d=1 THEN 1 ELSE 0 END) as with_3d,
    SUM(CASE WHEN smiles IS NOT NULL AND smiles!='' THEN 1 ELSE 0 END) as with_smiles
FROM ligands WHERE library_id=$LIB_ID;
EOF

echo ""
echo "[OK] Processing complete!"
