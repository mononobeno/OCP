# Vina × DB 完全統合設計書

日付: 2025-11-16  
担当: めけめけ もすもす

本ドキュメントは、OCP リポジトリにおける
「AutoDock Vina を SQLite ベースのカタログ DB と完全に統合する」ための
設計および変更方針をまとめたものである。

---

## 1. 背景と目的

現在の Vina パイプラインには以下の問題がある:

### RUN 情報が DB と連携されていない

Docker Compose 版 Vina では RUN_ID = <none> の状態で実行されており、
Vina 結果が vina_results テーブルに記録されていない。

### ZINC ID のソースがファイルベース

ligands_zinc_ids.txt を直接読んでいるため、
run_ligands テーブルが活用されていない。

### DB を source of truth にしたい方針とズレている

OCP 全体としては「DB が唯一の真実で、ファイルは再生成可能なキャッシュ」という
ポリシーを採用したいが、Vina 周辺の実装はまだ途中である。

### 目的:

- Vina の実行〜結果保存までを **完全に DB 経由** にする。
- RUN / RUN_LIGANDS / VINA_RESULTS を一貫して使える状態にする。
- 将来の GROMACS や可視化 Job から見て扱いやすいデータモデルに整える。

---

## 2. スコープと前提

### 2.1 スコープ

本設計書がカバーする範囲:

- Docker Compose 版 Vina
- SQLite DB (`catalog/db/ocp_results.sqlite`)
- テーブル:
  - runs
  - run_ligands
  - ligands
  - targets
  - libraries
  - vina_results

### 2.2 前提（仮のスキーマ）

以下のようなカラム構成を前提とする（実際の schema に合わせて調整すること）:

- `libraries(id, code, name, ...)`
- `targets(id, code, name, pdb_id, ...)`
- `ligands(id, zinc_id, smiles, source_file, has_3d, conformer_method, ...)`
- `runs(id, uuid, name, library_id, target_id, status, created_at, ...)`
- `run_ligands(id, run_id, ligand_id, order_index, ...)`
- `vina_results(id, run_id, ligand_id, score, mode, log_path, pose_path, created_at, ...)`

---

## 3. Vina 実行フロー（目標像）

### ZINC ID リストの決定

例: `catalog/libraries/zinc_2d_smi_v1/processed/ligands_zinc_ids.txt`

ligand-selector や手動抽出で作成。

### RUN の登録

`scripts/db_create_run_from_ids.sh`（新規）で

- runs に 1 レコード INSERT
- run_ligands に (run_id, ligand_id, order_index) を INSERT
- RUN_ID / RUN_UUID を取得し、Vina 実行側へ渡す。

### Vina 実行（Docker Compose 版）

`scripts/run_vina_docker_compose.sh` を拡張し、
RUN_ID / RUN_UUID を環境変数として vina-runner コンテナに渡す。

`images/vina-runner/vina_runner.sh` は以下の挙動に変更:

- RUN_ID があれば DB から run_ligands → ligands → ZINC ID を取得
- RUN_ID がなければ従来どおり ligands_zinc_ids.txt を使う（フォールバック）

### 結果保存

Vina 実行後、vina_results に INSERT:

- run_id
- ligand_id
- score
- log_path / pose_path
- created_at

---

## 4. 新規スクリプト db_create_run_from_ids.sh の役割

### 4.1 インターフェース

```bash
scripts/db_create_run_from_ids.sh \
  <LIBRARY_CODE> \
  <TARGET_CODE> \
  <ZINC_ID_LIST_FILE>
```

例:

```bash
scripts/db_create_run_from_ids.sh \
  zinc20_ml_v1 \
  apoh \
  catalog/libraries/zinc_2d_smi_v1/processed/ligands_zinc_ids.txt
```

### 4.2 処理内容

1. libraries.code, targets.code から library_id, target_id を取得。
2. runs に 1 レコード INSERT（uuid, name, status=pending 等）。
3. ZINC ID リストを読み込み:
   - ligands.zinc_id から ligand_id を取得。
   - run_ligands に (run_id, ligand_id, order_index) を INSERT。
4. 実行結果として RUN_ID=..., RUN_UUID=... を echo する。

これにより:

```bash
eval "$(scripts/db_create_run_from_ids.sh ...)"
export RUN_ID RUN_UUID
```

のように、RUN_ID / RUN_UUID を環境変数として後続のシェルに引き継げる。

---

## 5. 既存スクリプトへの変更方針（概要）

### 5.1 scripts/run_vina_docker_compose.sh

やりたいこと:

- Docker Compose 実行前に db_create_run_from_ids.sh を呼び出す。
- そこで得た RUN_ID / RUN_UUID を環境変数として vina-runner に渡す。

イメージ:

```bash
eval "$(scripts/db_create_run_from_ids.sh "${LIB_CODE}" "${TARGET_CODE}" "${ZINC_LIST_FILE}")"
export RUN_ID RUN_UUID

docker compose -f docker-compose.vina.yml up --build
```

### 5.2 docker-compose.vina.yml（例）

vina-runner サービスに RUN_ID / RUN_UUID を渡す:

```yaml
services:
  vina-runner:
    environment:
      - RUN_ID=${RUN_ID}
      - RUN_UUID=${RUN_UUID}
      # DB パスや他の設定もここで渡す
```

### 5.3 images/vina-runner/vina_runner.sh

- RUN_ID が設定されている場合:
  - DB の run_ligands / ligands から ZINC ID を取得する。
- RUN_ID が無い場合:
  - 現状通り ligands_zinc_ids.txt を読むフォールバックパスを維持する。
- Vina 実行後は、RUN_ID がある場合のみ vina_results に INSERT する。

---

## 6. RDKit 3D 生成との関係

将来的には「PDBQT ファイルが存在する」ではなく、
DB 上の has_3d = 1 や conformer_method = 'rdkit_etkdg' を
Vina 対象選択の条件に加えることを想定。

まずは現在のファイルベースの判定を維持しつつ、
段階的に DB ベースのフラグへ移行する。

---

## 7. デバッグ用 SQL と運用

`tools/sql/debug_vina_entities.sql` に、以下のようなクエリをまとめる:

- ライブラリ / ターゲット一覧
- RDKit 3D 生成の進捗
- runs / run_ligands の確認
- vina_results のランキング

VSCode 等から:

```bash
sqlite3 catalog/db/ocp_results.sqlite < tools/sql/debug_vina_entities.sql
```

でまとめて実行できる。

---

## 8. 実装チェックリスト

- [ ] `scripts/db_create_run_from_ids.sh` を作成・テストする。
- [ ] `run_vina_docker_compose.sh` を拡張し、RUN_ID / RUN_UUID を vina-runner に渡す。
- [ ] `vina_runner.sh` を修正し、RUN_ID がある場合は DB から ZINC ID を取得するようにする。
- [ ] `vina_results` への INSERT を有効化し、実際に結果が貯まることを確認する。
- [ ] RDKit 3D と receptor 品質改善を並行しつつ、スコアが妥当なレンジになるか確認する。

以上。
