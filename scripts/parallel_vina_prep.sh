#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# parallel_vina_prep.sh - 並列SMILES→PDBQT変換
###############################################################################

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

PARALLEL_JOBS="${1:-4}"
SMI_FILE="catalog/libraries/zinc_2d_smi_v1/raw/all_zinc20_ml_subset.smi"
DB_PATH="catalog/db/ocp_results.sqlite"
PDBQT_DIR="catalog/libraries/zinc_2d_smi_v1/processed/pdbqt"
WORK_DIR="/tmp/vina_prep_$$"

mkdir -p "$PDBQT_DIR" "$WORK_DIR"

LIB_ID=$(sqlite3 "$DB_PATH" "SELECT id FROM libraries WHERE code='zinc_2d_smi_v1';")

echo "=========================================="
echo "Parallel SMILES → PDBQT Conversion"
echo "=========================================="
echo "Parallel jobs: $PARALLEL_JOBS"
echo "Input file: $SMI_FILE"
echo "=========================================="

# 処理関数
process_compound() {
    local smiles="$1"
    local zinc_id="$2"
    local pdbqt_file="$PDBQT_DIR/${zinc_id}.pdbqt"
    local smiles_esc="${smiles//\'/\'\'}"
    
    # 既存チェック
    if [[ -f "$pdbqt_file" ]] && grep -q "ATOM.*[1-9]" "$pdbqt_file" 2>/dev/null; then
        sqlite3 "$DB_PATH" "UPDATE ligands SET smiles='$smiles_esc' WHERE zinc_id='$zinc_id' AND library_id=$LIB_ID;" 2>/dev/null || true
        echo "SKIP:$zinc_id"
        return 0
    fi
    
    # 3D生成
    if echo "$smiles" | obabel -ismi -opdbqt --gen3d -h -O "$pdbqt_file" &>/dev/null; then
        if grep -q "ATOM.*[1-9]" "$pdbqt_file" 2>/dev/null; then
            # DB登録
            existing=$(sqlite3 "$DB_PATH" "SELECT id FROM ligands WHERE zinc_id='$zinc_id' AND library_id=$LIB_ID;" 2>/dev/null || echo "")
            
            if [[ -n "$existing" ]]; then
                sqlite3 "$DB_PATH" "UPDATE ligands SET smiles='$smiles_esc', has_3d=1, conformer_method='obabel_gen3d' WHERE id=$existing;" 2>/dev/null || true
            else
                sqlite3 "$DB_PATH" "INSERT INTO ligands (zinc_id, library_id, smiles, has_3d, conformer_method) VALUES ('$zinc_id', $LIB_ID, '$smiles_esc', 1, 'obabel_gen3d');" 2>/dev/null || true
            fi
            
            echo "OK:$zinc_id"
            return 0
        fi
    fi
    
    echo "FAIL:$zinc_id"
    return 1
}

export -f process_compound
export DB_PATH LIB_ID PDBQT_DIR

# GNU parallelがあればそれを使用、なければxargsで並列化
if command -v parallel &>/dev/null; then
    echo "[INFO] Using GNU parallel"
    cat "$SMI_FILE" | parallel --colsep ' ' -j "$PARALLEL_JOBS" --bar process_compound {1} {2} > "$WORK_DIR/results.txt" 2>&1
else
    echo "[INFO] Using xargs (fallback)"
    cat "$SMI_FILE" | xargs -P "$PARALLEL_JOBS" -I {} bash -c 'process_compound $(echo {} | cut -d" " -f1) $(echo {} | cut -d" " -f2)' > "$WORK_DIR/results.txt" 2>&1
fi

# 結果集計
success=$(grep -c "^OK:" "$WORK_DIR/results.txt" 2>/dev/null || echo 0)
failed=$(grep -c "^FAIL:" "$WORK_DIR/results.txt" 2>/dev/null || echo 0)
skipped=$(grep -c "^SKIP:" "$WORK_DIR/results.txt" 2>/dev/null || echo 0)
total=$((success + failed + skipped))

echo ""
echo "=========================================="
echo "Conversion Complete!"
echo "=========================================="
echo "Total:   $total"
echo "Success: $success"
echo "Failed:  $failed"
echo "Skipped: $skipped"
echo "=========================================="

# DB統計
sqlite3 "$DB_PATH" <<EOF
SELECT 
    'Total: ' || COUNT(*) ||
    ' | With 3D: ' || SUM(CASE WHEN has_3d=1 THEN 1 ELSE 0 END) ||
    ' | With SMILES: ' || SUM(CASE WHEN smiles IS NOT NULL AND smiles!='' THEN 1 ELSE 0 END)
FROM ligands WHERE library_id=$LIB_ID;
EOF

# クリーンアップ
rm -rf "$WORK_DIR"

echo ""
echo "[OK] Processing complete!"
