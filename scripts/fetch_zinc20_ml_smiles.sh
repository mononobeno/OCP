#!/usr/bin/env bash
set -euo pipefail

# ============================
# ZINC20-ML SMILES ライブラリ取得スクリプト
# - ZINC20-ML_smiles.tar.gz (8.9G) をダウンロード
# - 中の .smi を一部だけ展開して raw/ に配置
# - 将来の pdbqt 前処理や DB 登録の「元データ」として使う
# ============================

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RAW_DIR="$ROOT_DIR/catalog/libraries/zinc_2d_smi_v1/raw"
TAR_URL="https://files.docking.org/zinc20-ML/ZINC20-ML_smiles.tar.gz"
TAR_PATH="$RAW_DIR/ZINC20-ML_smiles.tar.gz"

# 展開する .smi ファイル数（デフォルト 2個だけ）
NUM_FILES="${ZINC20_ML_NUM_FILES:-2}"

echo "[INFO] ROOT_DIR = $ROOT_DIR"
echo "[INFO] RAW_DIR  = $RAW_DIR"
echo "[INFO] will download: $TAR_URL"
echo "[INFO] ZINC20_ML_NUM_FILES = $NUM_FILES"

mkdir -p "$RAW_DIR"

# ---------- 1. tar.gz のダウンロード ----------
if [ -f "$TAR_PATH" ]; then
  echo "[INFO] tar.gz already exists: $TAR_PATH"
else
  echo "[INFO] downloading ZINC20-ML_smiles.tar.gz ..."
  wget -O "$TAR_PATH" "$TAR_URL"
fi

echo "[INFO] downloaded file:"
ls -lh "$TAR_PATH"

# ---------- 2. tarball 内の .smi 一覧を取得 ----------
TMP_LIST="$(mktemp)"
echo "[INFO] listing .smi inside tarball (this may take a bit)..."
tar -tzf "$TAR_PATH" | grep '\.smi$' | sort > "$TMP_LIST"

TOTAL_SMILES_FILES="$(wc -l < "$TMP_LIST" || echo 0)"
if [ "$TOTAL_SMILES_FILES" -eq 0 ]; then
  echo "[ERROR] No .smi files found inside $TAR_PATH"
  rm -f "$TMP_LIST"
  exit 1
fi

echo "[INFO] total .smi files in tarball: $TOTAL_SMILES_FILES"

# ---------- 3. 先頭 NUM_FILES 個だけ展開 ----------
EXTRACT_LIST="$(mktemp)"
head -n "$NUM_FILES" "$TMP_LIST" > "$EXTRACT_LIST"

echo "[INFO] extracting first $NUM_FILES .smi files to $RAW_DIR ..."
while read -r RELPATH; do
  echo "  - $RELPATH"
  tar -xzf "$TAR_PATH" -C "$RAW_DIR" "$RELPATH"
done < "$EXTRACT_LIST"

# 展開された .smi を RAW_DIR 直下に移動（smiles/ を剥がす）
echo "[INFO] flattening directory structure (moving smiles/*.smi to RAW_DIR)..."
find "$RAW_DIR" -type f -name '*.smi' -printf '%P\n' | while read -r P; do
  # 既に RAW_DIR 直下にあるものはそのまま
  case "$P" in
    *.smi)
      SRC="$RAW_DIR/$P"
      BASENAME="$(basename "$P")"
      if [ "$SRC" != "$RAW_DIR/$BASENAME" ]; then
        mv -f "$SRC" "$RAW_DIR/$BASENAME"
      fi
      ;;
  esac
done

echo "[INFO] .smi files now in $RAW_DIR:"
ls -lh "$RAW_DIR"/*.smi | head

# ---------- 4. （オプション）統合ファイルを作る ----------
MERGED="$RAW_DIR/all_zinc20_ml_subset.smi"
echo "[INFO] merging extracted .smi into $MERGED ..."
cat "$RAW_DIR"/*.smi > "$MERGED"
ls -lh "$MERGED"

# ---------- 5. ライブラリ方針レポートを書いておく ----------
POLICY_DIR="$ROOT_DIR/reports/policies"
mkdir -p "$POLICY_DIR"
POLICY_FILE="$POLICY_DIR/$(date +%Y-%m-%d)_library_source_policy.md"

cat > "$POLICY_FILE" << 'EOPOL'
# OCP 化合物ライブラリ方針（ZINC20-ML ベース）

## データソース

- 由来: ZINC20-ML ライブラリ（ZINC20 database から派生した Deep Docking 用 SMILES セット）
- 公式配布場所: https://files.docking.org/zinc20-ML/
- 主ファイル: `ZINC20-ML_smiles.tar.gz`  
  - 内容: 100 個の SMILES ファイル（各 1,000 万分子）  
  - 合計で ZINC20 分子空間のスナップショットを構成する 

## 利用目的

- OCP パイプラインにおける大規模バーチャルスクリーニング用の化合物ライブラリ
- AutoDock Vina → GROMACS の「入力候補集合」として利用
- DB (`ocp_results.sqlite`) と連携し、  
  - 既に使った ZINC ID を避ける
  - 新しいランごとに未使用 ZINC ID を割り当てる

## 運用ポリシー

1. **ローカル利用のみ**  
   - ZINC ライセンスに従い「検索結果やスクリーニング結果の共有」は可だが、  
     データベースそのものの再配布は行わない 。

2. **段階的展開**  
   - `ZINC20-ML_smiles.tar.gz` 全体（8.9GB）から .smi ファイルを一度に全展開するとストレージ負荷が大きい。  
   - 最初は数ファイルのみ展開し、必要に応じて追加する。

3. **ライブラリコード**  
   - OCP 内のライブラリコードは `zinc_2d_smi_v1` を継続利用する。  
   - 実体は ZINC20-ML 由来であることを本ドキュメントに明記。

4. **前処理**  
   - SMILES → 3D / pdbqt 変換は別途前処理ステージ（Kubernetes Job またはローカルバッチ）で実施。  
   - 変換済み pdbqt は  
     `catalog/libraries/zinc_2d_smi_v1/processed/pdbqt/`  
     に配置し、`vina-prep` が参照する。

EOPOL

echo "[INFO] library source policy written: $POLICY_FILE"

# cleanup
rm -f "$TMP_LIST" "$EXTRACT_LIST"

echo "[INFO] All done. ZINC20-ML subset ready under:"
echo "  $RAW_DIR"
echo "[INFO] Example files:"
ls -lh "$RAW_DIR" | head
