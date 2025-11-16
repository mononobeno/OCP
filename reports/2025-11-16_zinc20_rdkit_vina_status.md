# ZINC20 + RDKit + Vina 状況レポート

日付: 2025-11-16  
担当: めけめけ もすもす

本ドキュメントは、ZINC20 ライブラリを RDKit で 3D 化し、AutoDock Vina でドッキングする
パイプラインの「現在の状況」と「今後の方針」を整理したステータスレポートである。

---

## 1. 対象ライブラリと入力データ

- 対象ライブラリ: `zinc_2d_smi_v1`
- ルートディレクトリ: `catalog/libraries/zinc_2d_smi_v1/`

### 1.1 raw データ

- `raw/ZINC20-ML_smiles.tar.gz`  
  ZINC20 ML 系サブセット（SMILES ベース、約 9 GB 規模）
- `raw/all_zinc20_ml_subset.smi`  
  展開済み / フィルター済みの SMILES 一覧

### 1.2 processed ディレクトリ構成

- `processed/rdkit_3d/`  
  RDKit により SMILES → 3D PDB に変換したファイル群
- `processed/pdbqt/`  
  Open Babel により 3D PDB → PDBQT に変換したファイル群（Vina 入力）

---

## 2. RDKit ベース 3D 生成の現状

### 2.1 使用スクリプト

- `scripts/preprocess_smiles_with_rdkit.sh`
- 内部で `tools/rdkit_gen3d_to_pdbqt.py` を使用
- 主なフロー:
  1. SMILES を N 件ずつバッチ読み込み
  2. RDKit (ETKDG) で 3D conformer 生成
  3. PDB 出力 (`processed/rdkit_3d/`)
  4. Open Babel で PDB → PDBQT 変換 (`processed/pdbqt/`)
  5. DB `ligands` テーブルの `has_3d`, `conformer_method` を更新

### 2.2 DB 上の状態（概念）

- `ligands.smiles`         … ZINC の SMILES  
- `ligands.has_3d`         … 3D 生成済みフラグ (0/1)  
- `ligands.conformer_method` … 例: `rdkit_etkdg`

RDKit による 3D 生成完了件数（目安）は以下で確認する:

```sql
SELECT COUNT(*) AS ligands_rdkit_etkdg
FROM ligands
WHERE has_3d = 1
  AND conformer_method = 'rdkit_etkdg';
```

---

## 3. AutoDock Vina 実行の現状

### 3.1 実行方法（Docker Compose 版）

スクリプト: `scripts/run_vina_docker_compose.sh`

役割:
- receptor / ligands の PDBQT を workspace に集約
- `images/vina-runner` イメージを起動
- Vina を実行し、log / pose / score を生成

### 3.2 現状の挙動と制約

ログ例:
- total ligands to process: N
- best score = 0（スコアが 0 のケースが多い）

RUN_ID / RUN_UUID:
- 現時点では Docker Compose 版に RUN 情報が連携されておらず、
  RUN_ID = <none> の状態で実行されている。
- そのため vina_results テーブルへの INSERT は防御的にスキップしている。

---

## 4. Receptor PDBQT の品質について

現在使用している receptor.pdbqt は簡易変換スクリプト由来であり、
AutoDock Vina 標準の MGLTools ワークフローを通していない。

想定される問題:
- 原子タイプ / 電荷 / プロトン化状態が不完全
- その結果、スコアが 0 kcal/mol 付近に張り付いている可能性が高い

### 4.1 対応方針

MGLTools (AutoDockTools) の導入
- prepare_receptor4.py による正式な receptor.pdbqt 生成
- 再度 Vina を実行し、スコアレンジが
  おおむね -5〜-10 kcal/mol に入るか確認（ターゲット依存）

---

## 5. 今後の対応計画（要約）

### RUN / RUN_LIGANDS を DB から生成するスクリプトの整備

`scripts/db_create_run_from_ids.sh` を作成し、
ZINC ID リストから runs / run_ligands を自動生成する。

RUN_ID / RUN_UUID を vina-runner に渡し、vina_results に結果を蓄積する。

### RDKit 3D 生成の進捗管理

`tools/sql/debug_vina_entities.sql` で
ligands.has_3d の状況を可視化する。

定期的に `preprocess_smiles_with_rdkit.sh` を回し、
3D 生成済みの割合を増やしていく。

### Receptor PDBQT の品質改善

MGLTools 導入と receptor 再生成により、
Vina スコアが有効レンジに入るか検証する。

---

## 6. メモ

本レポートは「現状のスナップショット」であり、
具体的な件数やターゲット名は SQL とログから埋めていく。

「ZINC20 → RDKit → Vina → DB」の流れを土台に、
今後 GROMACS や GitHub Pages による可視化を載せていく予定。
