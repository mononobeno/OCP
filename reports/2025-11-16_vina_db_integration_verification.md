# Vina × DB 統合 動作確認手順書

日付: 2025-11-16  
担当: めけめけ もすもす

本ドキュメントは、Vina パイプラインを DB 駆動型に変更した後の動作確認手順をまとめたものです。

---

## 前提条件

### 1. DB スキーマの準備

以下のテーブルが存在していることを確認してください:

```bash
sqlite3 catalog/db/ocp_results.sqlite << 'EOF'
.schema libraries
.schema targets
.schema ligands
.schema runs
.schema run_ligands
.schema vina_results
EOF
```

### 2. 必要なマスターデータ

- `libraries` テーブルに `code='zinc20_ml_v1'` (または適切なコード) が登録されていること
- `targets` テーブルに `code='apoh'` (または適切なコード) が登録されていること
- `ligands` テーブルに ZINC ID が登録されていること (has_3d=1 のものが望ましい)

マスターデータ確認:

```bash
sqlite3 catalog/db/ocp_results.sqlite << 'EOF'
SELECT id, code, name FROM libraries;
SELECT id, code, name FROM targets;
SELECT COUNT(*) AS total, SUM(has_3d) AS with_3d FROM ligands;
EOF
```

---

## 動作確認フロー

### ステップ 1: ZINC ID リストの準備

Vina で実行したい ZINC ID のリストを用意します。

```bash
# 例: has_3d=1 の ligands から 5 件抽出
sqlite3 catalog/db/ocp_results.sqlite << 'EOF' > catalog/libraries/zinc_2d_smi_v1/processed/ligands_zinc_ids.txt
SELECT zinc_id 
FROM ligands 
WHERE has_3d = 1 
LIMIT 5;
EOF

# 内容確認
cat catalog/libraries/zinc_2d_smi_v1/processed/ligands_zinc_ids.txt
```

### ステップ 2: RUN の作成

`scripts/db_create_run_from_ids.sh` を使って RUN を DB に登録します。

```bash
# LIB_CODE, TARGET_CODE, ZINC_ID_LIST_FILE を指定
bash scripts/db_create_run_from_ids.sh \
  zinc20_ml_v1 \
  apoh \
  catalog/libraries/zinc_2d_smi_v1/processed/ligands_zinc_ids.txt
```

**期待される出力:**

```
[INFO] DB : catalog/db/ocp_results.sqlite
[INFO] LIB_CODE : zinc20_ml_v1 (id=1)
[INFO] TARGET : apoh (id=1)
[INFO] ZINC LIST: catalog/libraries/zinc_2d_smi_v1/processed/ligands_zinc_ids.txt
[INFO] RUN_NAME : run_zinc20_ml_v1_apoh_20251116_123456
[INFO] Created run: id=1, uuid=xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx
[INFO] Inserted 5 rows into run_ligands (skipped 0).
RUN_ID=1
RUN_UUID=xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx
```

### ステップ 3: DB 内容の確認

RUN と run_ligands が正しく登録されたか確認します。

```bash
sqlite3 catalog/db/ocp_results.sqlite << 'EOF'
-- 最新の RUN
SELECT id, uuid, name, library_id, target_id, status, created_at 
FROM runs 
ORDER BY id DESC 
LIMIT 1;

-- run_ligands の内容
SELECT rl.run_id, rl.order_index, l.zinc_id
FROM run_ligands rl
JOIN ligands l ON l.id = rl.ligand_id
WHERE rl.run_id = (SELECT MAX(id) FROM runs)
ORDER BY rl.order_index;
EOF
```

### ステップ 4: Vina 実行

DB 駆動型の Vina 実行スクリプトを実行します。

```bash
# リポジトリルートで実行
cd /home/dev/OCP
bash scripts/run_vina_docker_compose.sh
```

**期待される動作:**

1. `db_create_run_from_ids.sh` が自動的に呼ばれる
2. RUN_ID と RUN_UUID が取得される
3. Docker Compose が起動し、vina-runner コンテナに RUN_ID/RUN_UUID が渡される
4. vina-runner は RUN_ID から DB 経由で ZINC ID を取得する
5. Vina を実行し、結果を `vina_results` テーブルに保存する

**ログ確認:**

```bash
# コンテナログ確認
docker logs vina-runner-standalone

# 以下のようなメッセージが出ていれば成功
# [INFO] RUN_ID detected, fetching ZINC IDs from database...
# [INFO] Retrieved 5 ZINC IDs from DB for RUN_ID=1
# [INFO] [1/5] running vina for ZINC001241740817_1
# [INFO] [1/5] best score = -0.5
# [INFO] [1/5] saved to DB: run_id=1, ligand_id=123, score=-0.5
```

### ステップ 5: 結果の確認

#### 5.1 ファイルベースの結果

```bash
# pose ファイルの確認
ls -lh results/vina_output/poses/

# log ファイルの確認
ls -lh results/vina_output/logs/

# 特定の log を見る
head -20 results/vina_output/logs/ZINC001241740817_1.log
```

#### 5.2 DB ベースの結果

```bash
# vina_results テーブルの確認
sqlite3 catalog/db/ocp_results.sqlite << 'EOF'
SELECT COUNT(*) AS total_results FROM vina_results;

-- 最新 RUN の結果
SELECT 
  vr.id,
  vr.run_id,
  l.zinc_id,
  vr.score
FROM vina_results vr
JOIN ligands l ON l.id = vr.ligand_id
WHERE vr.run_id = (SELECT MAX(id) FROM runs)
ORDER BY vr.score ASC;

-- スコアランキング (上位 10 件)
SELECT 
  vr.run_id,
  l.zinc_id,
  vr.score,
  vr.created_at
FROM vina_results vr
JOIN ligands l ON l.id = vr.ligand_id
ORDER BY vr.score ASC
LIMIT 10;
EOF
```

#### 5.3 debug_vina_entities.sql による包括確認

```bash
sqlite3 catalog/db/ocp_results.sqlite < tools/sql/debug_vina_entities.sql
```

このスクリプトは以下の情報をまとめて出力します:

- libraries / targets の一覧
- RDKit 3D 生成状況
- 最近の runs 一覧
- 最新 RUN の run_ligands
- vina_results の概要とランキング

---

## トラブルシューティング

### 問題 1: `library not found for code=...`

**原因:** `libraries` テーブルに該当する `code` が存在しない

**対処:**

```bash
sqlite3 catalog/db/ocp_results.sqlite << 'EOF'
INSERT INTO libraries (code, name) VALUES ('zinc20_ml_v1', 'ZINC20 ML subset');
EOF
```

### 問題 2: `target not found for code=...`

**原因:** `targets` テーブルに該当する `code` が存在しない

**対処:**

```bash
sqlite3 catalog/db/ocp_results.sqlite << 'EOF'
INSERT INTO targets (code, name, pdb_id) VALUES ('apoh', 'APS ApoH', '1C1Z');
EOF
```

### 問題 3: `ligand not found for zinc_id=..., skipped.`

**原因:** ZINC ID リストに含まれる ID が `ligands` テーブルに存在しない

**対処:**

- ZINC ID が DB に登録されているか確認: `SELECT * FROM ligands WHERE zinc_id='ZINC...';`
- 登録されていない場合は ligand-selector や登録スクリプトを実行

### 問題 4: `No ligands found in run_ligands for RUN_ID=...`

**原因:** `run_ligands` テーブルに RUN_ID のレコードが無い

**対処:**

- `db_create_run_from_ids.sh` が正常に実行されたか確認
- `run_ligands` テーブルを直接確認: `SELECT * FROM run_ligands WHERE run_id=1;`

### 問題 5: vina_results に結果が保存されない

**原因 A:** RUN_ID が無い (フォールバックモードで実行されている)

**対処:**

- `run_vina_docker_compose.sh` が `db_create_run_from_ids.sh` を呼んでいるか確認
- Docker Compose のログで `RUN_ID = <none>` になっていないか確認

**原因 B:** ligand_id が取得できていない

**対処:**

- vina-runner のログで `[WARN] skipping DB insert: LIG_ID=, ...` が出ていないか確認
- ZINC ID と ligands テーブルの `zinc_id` カラムの形式が一致しているか確認

---

## 動作確認完了後

### DB の状態を確認

```bash
sqlite3 catalog/db/ocp_results.sqlite << 'EOF'
-- RUN の総数
SELECT COUNT(*) AS total_runs FROM runs;

-- run_ligands の総数
SELECT COUNT(*) AS total_run_ligands FROM run_ligands;

-- vina_results の総数
SELECT COUNT(*) AS total_vina_results FROM vina_results;

-- RUN ごとの結果数
SELECT 
  r.id,
  r.name,
  COUNT(vr.id) AS result_count,
  MIN(vr.score) AS best_score
FROM runs r
LEFT JOIN vina_results vr ON vr.run_id = r.id
GROUP BY r.id
ORDER BY r.id DESC;
EOF
```

### クリーンアップ (必要に応じて)

```bash
# テスト用の RUN を削除する場合
sqlite3 catalog/db/ocp_results.sqlite << 'EOF'
DELETE FROM vina_results WHERE run_id = 1;
DELETE FROM run_ligands WHERE run_id = 1;
DELETE FROM runs WHERE id = 1;
EOF

# 結果ファイルを削除する場合
rm -rf results/vina_output
```

---

## 次のステップ

1. **Receptor 品質の改善**
   - MGLTools の `prepare_receptor4.py` を使って正しい PDBQT を作成
   - スコアが妥当な範囲 (-5 ~ -10 kcal/mol 程度) に入るか確認

2. **RDKit 3D 生成の拡充**
   - `scripts/preprocess_smiles_with_rdkit.sh` を使って 3D conformer を追加生成
   - `has_3d=1` の ligands を増やす

3. **本番実行**
   - 50 ~ 100 ligands でテスト実行
   - スコアランキングを確認し、上位候補を GROMACS へ

4. **可視化 (GitHub Pages)**
   - `vina_results` から JSON/CSV を生成
   - Pages で結果を表示

---

以上で Vina × DB 統合の動作確認手順は完了です。

