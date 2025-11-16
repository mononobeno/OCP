#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# bulk_register_and_gen3d_vina.sh
#
# SMILESファイルから一括でDB登録 → 3D生成 → PDBQT変換を行う高速スクリプト
#
# 使用方法:
#   bash scripts/bulk_register_and_gen3d_vina.sh [BATCH_SIZE] [START_LINE] [END_LINE]
#
# 例:
#   bash scripts/bulk_register_and_gen3d_vina.sh 100       # 全件を100件ずつ処理
#   bash scripts/bulk_register_and_gen3d_vina.sh 50 1 200  # 1-200行目を50件ずつ処理
###############################################################################

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

# パラメータ
BATCH_SIZE="${1:-100}"
START_LINE="${2:-1}"
END_LINE="${3:-}"

# パス設定
SMI_FILE="catalog/libraries/zinc_2d_smi_v1/raw/all_zinc20_ml_subset.smi"
DB_PATH="catalog/db/ocp_results.sqlite"
PDBQT_DIR="catalog/libraries/zinc_2d_smi_v1/processed/pdbqt"
LIB_CODE="zinc_2d_smi_v1"

# ディレクトリ作成
mkdir -p "$PDBQT_DIR"

# ライブラリIDを取得
LIB_ID=$(sqlite3 "$DB_PATH" "SELECT id FROM libraries WHERE code='$LIB_CODE';")
if [[ -z "$LIB_ID" ]]; then
    echo "[ERROR] Library '$LIB_CODE' not found in DB"
    exit 1
fi

# 総行数を取得
TOTAL_LINES=$(wc -l < "$SMI_FILE")
if [[ -z "$END_LINE" ]]; then
    END_LINE=$TOTAL_LINES
fi

echo "========================================"
echo "Bulk SMILES → DB → 3D → PDBQT Pipeline"
echo "========================================"
echo "SMI File    : $SMI_FILE"
echo "Total Lines : $TOTAL_LINES"
echo "Processing  : Lines $START_LINE - $END_LINE"
echo "Batch Size  : $BATCH_SIZE"
echo "Library     : $LIB_CODE (id=$LIB_ID)"
echo "========================================"

# カウンター
total_processed=0
total_success=0
total_failed=0
total_skipped=0

# バッチ処理
current_line=$START_LINE
while [[ $current_line -le $END_LINE ]]; do
    batch_end=$((current_line + BATCH_SIZE - 1))
    if [[ $batch_end -gt $END_LINE ]]; then
        batch_end=$END_LINE
    fi
    
    echo ""
    echo "[BATCH] Processing lines $current_line - $batch_end..."
    
    # バッチファイル作成
    batch_file="/tmp/batch_${current_line}_${batch_end}.smi"
    sed -n "${current_line},${batch_end}p" "$SMI_FILE" > "$batch_file"
    
    batch_processed=0
    batch_success=0
    batch_failed=0
    batch_skipped=0
    
    # 1行ずつ処理
    while IFS=' ' read -r smiles zinc_id || [[ -n "$smiles" ]]; do
        # 空行スキップ
        if [[ -z "$smiles" ]] || [[ -z "$zinc_id" ]]; then
            continue
        fi
        
        ((batch_processed++))
        
        # DB登録チェック
        existing_id=$(sqlite3 "$DB_PATH" "SELECT id FROM ligands WHERE zinc_id='$zinc_id' AND library_id=$LIB_ID;" 2>/dev/null || echo "")
        
        # SMILESをエスケープ（シングルクォートを二重化）
        smiles_escaped="${smiles//\'/\'\'}"
        
        if [[ -n "$existing_id" ]]; then
            # 既に登録済み - SMILESを更新
            sqlite3 "$DB_PATH" "UPDATE ligands SET smiles='$smiles_escaped' WHERE id=$existing_id;" 2>/dev/null || true
            ligand_db_id=$existing_id
        else
            # 新規登録
            if ! sqlite3 "$DB_PATH" "INSERT INTO ligands (zinc_id, library_id, smiles, has_3d, conformer_method) VALUES ('$zinc_id', $LIB_ID, '$smiles_escaped', 0, '');" 2>/dev/null; then
                echo "  [WARN] Failed to insert $zinc_id"
                ((batch_failed++))
                continue
            fi
            ligand_db_id=$(sqlite3 "$DB_PATH" "SELECT last_insert_rowid();")
        fi
        
        # PDBQT生成チェック
        pdbqt_file="$PDBQT_DIR/${zinc_id}.pdbqt"
        if [[ -f "$pdbqt_file" ]]; then
            # 座標が有効かチェック
            if grep -q "ATOM.*[1-9]" "$pdbqt_file" 2>/dev/null; then
                ((batch_skipped++))
                continue
            fi
        fi
        
        # 3D生成 (Open Babel)
        echo "$smiles" | obabel -ismi -opdbqt --gen3d -h -O "$pdbqt_file" 2>/dev/null 1>/dev/null
        
        if [[ $? -eq 0 ]] && [[ -f "$pdbqt_file" ]]; then
            # 座標検証
            if grep -q "ATOM.*[1-9]" "$pdbqt_file" 2>/dev/null; then
                # DB更新
                sqlite3 "$DB_PATH" "UPDATE ligands SET has_3d=1, conformer_method='obabel_gen3d' WHERE id=$ligand_db_id;"
                ((batch_success++))
            else
                echo "  [WARN] $zinc_id: Generated PDBQT has zero coordinates"
                ((batch_failed++))
            fi
        else
            echo "  [WARN] $zinc_id: 3D generation failed"
            ((batch_failed++))
        fi
        
    done < "$batch_file"
    
    # バッチ結果表示
    echo "  Processed: $batch_processed | Success: $batch_success | Failed: $batch_failed | Skipped: $batch_skipped"
    
    # 累計更新
    total_processed=$((total_processed + batch_processed))
    total_success=$((total_success + batch_success))
    total_failed=$((total_failed + batch_failed))
    total_skipped=$((total_skipped + batch_skipped))
    
    # クリーンアップ
    rm -f "$batch_file"
    
    # 次のバッチ
    current_line=$((batch_end + 1))
done

echo ""
echo "========================================"
echo "Processing Complete!"
echo "========================================"
echo "Total Processed : $total_processed"
echo "Total Success   : $total_success"
echo "Total Failed    : $total_failed"
echo "Total Skipped   : $total_skipped"
echo "========================================"

# DB統計
echo ""
echo "Database Statistics:"
sqlite3 "$DB_PATH" << EOF
SELECT 
    'Total ligands: ' || COUNT(*) 
FROM ligands WHERE library_id=$LIB_ID;

SELECT 
    'With SMILES: ' || COUNT(*) 
FROM ligands WHERE library_id=$LIB_ID AND smiles IS NOT NULL AND smiles != '';

SELECT 
    'With 3D (has_3d=1): ' || COUNT(*) 
FROM ligands WHERE library_id=$LIB_ID AND has_3d=1;

SELECT 
    'Conformer methods: ' || conformer_method || ' = ' || COUNT(*) 
FROM ligands 
WHERE library_id=$LIB_ID AND has_3d=1 
GROUP BY conformer_method;
EOF

echo ""
echo "[OK] Bulk processing complete!"
