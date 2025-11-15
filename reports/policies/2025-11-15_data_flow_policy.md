# OCP データフロー・方針書（DB / Kubernetes / GitHub Pages）

- 作成日: 2025-11-15
- 対象: ZINC ライブラリ / Vina / GROMACS / Pages 連携

---

## 1. 全体コンセプト

1. ZINC などの巨大ライブラリと受容体は **catalog ディレクトリ**で一元管理する。
2. ドッキング・MD の結果は **単一の SQLite DB (catalog/db/ocp_results.sqlite)** に集約する。
3. Kubernetes Job と GitHub Actions は **同じ DB ファイルを共有**し、
   それを「真実のソース（source of truth）」とする。
4. GitHub Pages に公開するレポートは、DB の内容から自動生成された静的ファイル（docs/ 以下）とする。

---

## 2. ディレクトリ構成ポリシー

- ZINC ライブラリ（入力）

  - `catalog/libraries/zinc_2d_smi_v1/raw/`
    - ZINC-downloader 由来の SMILES 等
  - `catalog/libraries/zinc_2d_smi_v1/processed/`
    - Vina 用に 3D 変換・pdbqt 化したもの（将来）

- ターゲット受容体（入力）

  - `catalog/targets/aps_apoh/receptor.pdb`
  - `catalog/targets/aps_prothrombin/receptor.pdb`
  - `catalog/targets/aps_tlr4_md2/receptor.pdb`

- 実行時ワークスペース（PVC）

  - `/workspace/` 以下に receptor / ligands_raw / vina_input / vina_output / gmx_* / reports などを配置。
  - 一回の run で使った一時ファイルは基本的に /workspace に置き、必要なメタ情報だけ DB に書き戻す。

- 結果 DB

  - `catalog/db/ocp_results.sqlite`
  - スキーマは `catalog/db/schema.sql` に定義。

---

## 3. Kubernetes から見たデータフロー

1. ligand-selector Job
   - /db/ocp_results.sqlite を sqlite3 で参照。
   - 「ある target + library に対してまだ実行していない ZINC ID」を SQL で選択。
   - 選んだ化合物のリストを `/workspace/ligands_raw/` に展開し、ligands テーブルに挿入。

2. vina-prep / vina-runner / vina2gromacs-prep
   - `/workspace` だけを主に参照し、必要に応じて catalog/targets を読む。
   - Vina のスコアは run_id / ligand_id 単位で `vina_results` に INSERT。

3. gmx-* Job
   - EM / EQ / PROD の計算を行い、RMSD やエネルギーなどの要約値を `md_results` に INSERT。
   - HTML / PNG などの解析レポートは `/workspace/reports/` に出力し、
     その相対パスを `md_results.report_relpath` に保存。

4. Pages-trigger Job（任意）
   - 新しい run が完了したタイミングで GitHub Actions の workflow_dispatch API を叩き、
     「DB から docs/ を再生成する」ワークフローを起動する。

---

## 4. GitHub Actions から見たデータフロー

1. `actions/checkout` によりリポジトリを取得すると、
   `catalog/db/ocp_results.sqlite` も同時に手元に来る。
2. `sqlite3 catalog/db/ocp_results.sqlite "SELECT ..."` で集計クエリを実行。
3. Python やスクリプト（`scripts/generate_pages_from_db.py`）で結果を HTML / Markdown に変換し、
   `docs/` 以下にファイルを生成する。
4. `docs/` を commit & push することで、GitHub Pages が自動更新される。
5. 更新した run_id と URL を `pages_publish` テーブルに INSERT しておくことで、
   「どの run がいつ公開されたか」を履歴として追跡できる。

---

## 5. 「かぶりのない化合物選択」のポリシー

- 同じ target × library 組み合わせで、同じ ZINC ID を二度使わないことを原則とする。
- ligand-selector Job では以下の条件で SQL を発行して化合物を選ぶ。

  - libraries.code = 対象ライブラリ（例: zinc_2d_smi_v1）
  - targets.code   = 対象ターゲット（例: aps_apoh）
  - その ligand_id について、過去の runs + vina_results にレコードが存在しないものだけを候補とする。

- この仕組みにより、GitHub Actions や別の Pod から参照しても、
  DB を見るだけで「どの化合物をいつ試したか」が再現可能になる。

---

## 6. 将来の拡張

- 必要になった時点で、SQLite から外部の RDB（PostgreSQL 等）へ移行することを視野に入れる。
- スキーマは極力シンプルに保ち、移行時はテーブル定義をほぼそのまま持ち込めるようにする。
- 解析メトリクス（エネルギー、相互作用解析結果など）が増えた場合、
  md_results にカラム追加、あるいは別テーブルを追加する形で拡張する。
