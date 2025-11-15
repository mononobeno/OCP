#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

REPORT_DIR="$ROOT_DIR/reports"
REPORT_FILE="$REPORT_DIR/$(date +%Y-%m-%d)_zinc20_2d_pipeline_status.md"
DB="$ROOT_DIR/catalog/db/ocp_results.sqlite"

echo "[INFO] ROOT_DIR   = $ROOT_DIR"
echo "[INFO] REPORT_DIR = $REPORT_DIR"
echo "[INFO] REPORT_FILE= $REPORT_FILE"

mkdir -p "$REPORT_DIR"

# --- 現状値を取得 ---
PDBQT_COUNT=$(ls "$ROOT_DIR"/catalog/libraries/zinc_2d_smi_v1/processed/pdbqt/*.pdbqt 2>/dev/null | wc -l || echo 0)
LIGANDS_COUNT=$(sqlite3 "$DB" "SELECT COUNT(*) FROM ligands WHERE library_id=1;" 2>/dev/null || echo 0)

# 直近 10 件くらいの ligands をダンプ（zinc_id + smiles）
LIGANDS_SAMPLE=$(sqlite3 "$DB" "SELECT zinc_id, substr(smiles,1,40) FROM ligands WHERE library_id=1 ORDER BY id DESC LIMIT 10;" 2>/dev/null || echo "")

cat > "$REPORT_FILE" << EOFMD
# OCP 化合物ライブラリ／2D 前処理パイプライン状況レポート

日付: $(date +%Y-%m-%d)
リポジトリ: OCP
対象ライブラリコード: \`zinc_2d_smi_v1\`

---

## 1. 概要

ZINC20-ML 由来の SMILES サブセットを取り込み、以下を実施した。

- SMILES → PDBQT 変換を **2D 構造（\`--gen3d\` なし）** で実施
- 100 化合物を約 10 秒で処理
- DB (\`ocp_results.sqlite\`) の \`ligands\` テーブルにメタ情報を登録
- Kubernetes 上の \`ligand-selector\` → \`vina-prep\` が、このライブラリを参照して正常に動作することを確認

本レポートでは、変更点と現状の状態、および今後の方針を整理する。

---

## 2. ライブラリ・DB の現状

### 2.1 PDBQT ファイル数

- ディレクトリ: \`catalog/libraries/zinc_2d_smi_v1/processed/pdbqt/\`
- PDBQT ファイル数: **${PDBQT_COUNT} 件**

（例）

\`\`\`bash
ls catalog/libraries/zinc_2d_smi_v1/processed/pdbqt | head
\`\`\`

### 2.2 DB 上の ligands レコード数

- DB: \`catalog/db/ocp_results.sqlite\`
- 対象テーブル: \`ligands\`
- \`library_id=1\` (\`zinc_2d_smi_v1\`) のレコード数: **${LIGANDS_COUNT} 件**

### 2.3 ligands テーブルのメタ情報

今回のマイグレーションにより、以下の列を追加・利用している。

- \`ligands.smiles\`      : SMILES 文字列（1 列目）
- \`ligands.source_file\`: 元になった .smi ファイル名（例: \`all_zinc20_ml_subset.smi\`）

サンプル（直近 10 件）:

\`\`\`text
zinc_id | smiles(先頭40文字)
--------------------------------
${LIGANDS_SAMPLE}
\`\`\`

---

## 3. スクリプトおよび構成変更点

### 3.1 DB マイグレーション

追加スクリプト:

- \`scripts/db_migrate_add_ligand_metadata.sh\`

主な内容:

- \`ligands.smiles TEXT\`
- \`ligands.source_file TEXT\`

を **存在しない場合のみ** 追加する。

### 3.2 2D-only 前処理スクリプト

修正済みスクリプト:

- \`scripts/preprocess_smiles_to_pdbqt.sh\`

主な変更点:

- \`--gen3d\` を削除し、2D 構造のまま \`.pdbqt\` を生成
- SMILES および source_file を DB (\`ligands\`) に登録／補完するように変更

ロジック概要:

1. 入力: \`.smi\` / \`.smi.gz\` （1列目 SMILES, 2列目 ID）
2. 各行について:
   - \`.pdbqt\` が未存在なら \`obabel -ismi\` で 2D のまま出力
   - \`ligands (zinc_id, library_id)\` を \`INSERT OR IGNORE\`
   - \`smiles\` と \`source_file\` を \`UPDATE\`（既存が NULL の場合のみ）

### 3.3 ZINC20-ML サブセット取り込み

新規スクリプト:

- \`scripts/fetch_zinc20_ml_smiles.sh\`

主な役割:

1. \`ZINC20-ML_smiles.tar.gz\`（約 9GB）を \`raw/\` に配置
   - 現状はブラウザでのダウンロード → WSL へコピーでも可
2. tarball 内から .smi ファイルを一部（環境変数 \`ZINC20_ML_NUM_FILES\` 個）展開
3. 展開した .smi を \`raw/\` にフラットに置き、\`all_zinc20_ml_subset.smi\` に統合
4. ライブラリ方針レポートを \`reports/policies/YYYY-MM-DD_library_source_policy.md\` として生成

### 3.4 ZINC20 サブセット用ブートストラップ

新規スクリプト（提案）:

- \`scripts/bootstrap_zinc20_subset_pipeline.sh\`

主な役割:

1. \`all_zinc20_ml_subset.smi\` から指定件数（デフォルト 1000 / 現状は 100 件など）だけ前処理
2. DB の \`ligands\` 件数を表示
3. \`scripts/run_ligand_selector_stage.sh\`
4. \`scripts/run_vina_prep_stage.sh\`

これにより、**「ライブラリ → DB → K8s パイプライン前半」** までを一括で検証できる。

---

## 4. パフォーマンス比較と設計判断

### 4.1 パフォーマンス比較（実測）

ユーザー実測値:

- \`--gen3d\` あり:
  - 1 化合物あたり ~20 分（推定）
  - 100 化合物 ≒ 33 時間
- \`--gen3d\` なし（2D のみ）:
  - 100 化合物を約 10 秒で処理
  - おおよそ **1.2 万倍高速**

この結果に基づき、現フェーズでは **2D 構造のまま 100 化合物程度でパイプラインを回す** 方針とした。

### 4.2 設計上の判断

- **DB を「唯一の真実」にする方針**

  - 「どの ZINC ID を計算対象にするか」は常に DB 由来とする。
  - \`ligand-selector\` Job は、DB 上で  
    - target (\`aps_apoh\`)  
    - library (\`zinc_2d_smi_v1\`)  
    - 過去の Vina 結果  
    を見ながら「未使用の ligand」を選定し、結果を PVC (\`ligands_raw/\`) に書き出す。

- **ファイルシステムは「実体置き場」**

  - \`processed/pdbqt\` は「DB で管理している ligands の実体」が置かれる場所。
  - DB には少なくとも \`zinc_id\`, \`library_id\`, \`smiles\`, \`source_file\` が残るため、
    ファイルが消えても「どの分子だったか」の再現は可能。

- **3D 化（\`--gen3d\`）は別フェーズで検討**

  - 現状は OCP/Kubernetes パイプラインの構築と検証が主目的。
  - 本格的なドッキングスコアの信頼性向上のための 3D 化は、
    別ライブラリコード（例: \`zinc_3d_v1\`）として切り出す想定。

---

## 5. Kubernetes パイプラインへの影響

### 5.1 ligand-selector

- 入力:
  - DB (\`ocp_results.sqlite\`)
  - \`targets\`, \`libraries\`, \`ligands\`, \`runs\` などのテーブル
- 出力 (PVC: \`/workspace/ligands_raw\`):
  - \`ligands_selected.tsv\`
  - \`ligands_zinc_ids.txt\`
  - \`run_info.json\`

→ **計算対象の選定ロジックは完全に DB 起点**。

### 5.2 vina-prep

- 入力:
  - \`ligands_zinc_ids.txt\`
  - \`/zinc-library/targets/aps_apoh/receptor.pdb\`
  - \`/zinc-library/libraries/zinc_2d_smi_v1/processed/pdbqt/\${ZINC_ID}.pdbqt\`
- 出力 (PVC: \`/workspace/vina_input\`):
  - \`receptor.pdb\`, \`receptor.pdbqt\`
  - \`ligands/\${ZINC_ID}.pdbqt\`

→ 「どの ZINC ID を使うか」は ligand-selector(=DB)が決めており、  
  vina-prep はそれに合わせてファイルを配置する役割に集中している。

---

## 6. 今後のタスク

1. **vina-runner Job の実装**
   - 入力: \`/workspace/vina_input\`（receptor + ligands）
   - 出力: \`/workspace/vina_output\`（Vina の out.pdbqt, log など）
   - Vina 実行結果を DB (\`vina_results\` 等) に登録

2. **vina2gromacs-prep → gmx-makebox 以降の GROMACS パイプライン実装**
   - Vina 出力構造を GROMACS 入力用に変換
   - CPU モードでの \`gmx grompp\`, \`gmx mdrun\` などを Job 化

3. **3D 構造ライブラリ（任意）**
   - 別ライブラリコード（例: \`zinc_3d_v1\`）として 3D 構造を管理
   - 2D ライブラリとは別枠で運用

---

## 7. まとめ

- ZINC20-ML ベースの化合物ライブラリを OCP に取り込み、
  2D のまま 100 化合物を高速前処理するパイプラインが構築できた。
- DB (\`ligands\`) には ZINC ID・SMILES・ソースファイルが記録され、
  ライブラリの管理は DB が「唯一の真実」として機能する。
- Kubernetes 上の \`ligand-selector\` → \`vina-prep\` は、
  この DB とライブラリを前提に正常動作している。
- 今後は \`vina-runner\` 以降のステージを実装しつつ、
  必要に応じて 3D 構造ライブラリを追加する。
EOFMD

echo "[INFO] report generated: $REPORT_FILE"
