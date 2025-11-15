#!/usr/bin/env bash
set -euo pipefail

# === 0. 前提チェック ===
if ! command -v sqlite3 >/dev/null 2>&1; then
  echo "ERROR: sqlite3 コマンドが見つかりません。" >&2
  echo "sudo apt-get install -y sqlite3 などでインストールしてください。" >&2
  exit 1
fi

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

echo "[INFO] ROOT_DIR = $ROOT_DIR"

# === 1. ディレクトリ作成 ===
mkdir -p catalog/db
mkdir -p reports/policies

DB="$ROOT_DIR/catalog/db/ocp_results.sqlite"
SCHEMA="$ROOT_DIR/catalog/db/schema.sql"

echo "[INFO] DB file   = $DB"
echo "[INFO] SCHEMA    = $SCHEMA"

# === 2. スキーマ定義を書き出し ===
cat > "$SCHEMA" << 'EOSQL'
-- OCP results database schema

PRAGMA foreign_keys = ON;

CREATE TABLE IF NOT EXISTS targets (
  id            INTEGER PRIMARY KEY AUTOINCREMENT,
  code          TEXT NOT NULL UNIQUE,    -- "aps_apoh" など
  name          TEXT NOT NULL,           -- 表示名
  pdb_id        TEXT,                    -- 1C1Z, 6C2W ...
  receptor_path TEXT                     -- catalog/targets/.../receptor.pdb
);

CREATE TABLE IF NOT EXISTS libraries (
  id            INTEGER PRIMARY KEY AUTOINCREMENT,
  code          TEXT NOT NULL UNIQUE,    -- "zinc_2d_smi_v1" など
  description   TEXT
);

CREATE TABLE IF NOT EXISTS ligands (
  id            INTEGER PRIMARY KEY AUTOINCREMENT,
  zinc_id       TEXT NOT NULL,
  library_id    INTEGER NOT NULL,
  smiles_hash   TEXT,
  UNIQUE (zinc_id, library_id),
  FOREIGN KEY (library_id) REFERENCES libraries(id)
);

CREATE TABLE IF NOT EXISTS runs (
  id            INTEGER PRIMARY KEY AUTOINCREMENT,
  run_uuid      TEXT NOT NULL UNIQUE,
  started_at    TEXT NOT NULL,
  target_id     INTEGER NOT NULL,
  library_id    INTEGER NOT NULL,
  notes         TEXT,
  FOREIGN KEY (target_id)  REFERENCES targets(id),
  FOREIGN KEY (library_id) REFERENCES libraries(id)
);

CREATE TABLE IF NOT EXISTS vina_results (
  run_id        INTEGER NOT NULL,
  ligand_id     INTEGER NOT NULL,
  mode_rank     INTEGER NOT NULL,
  affinity_kcal REAL,
  rmsd_lb       REAL,
  rmsd_ub       REAL,
  out_relpath   TEXT,
  PRIMARY KEY (run_id, ligand_id, mode_rank),
  FOREIGN KEY (run_id)    REFERENCES runs(id),
  FOREIGN KEY (ligand_id) REFERENCES ligands(id)
);

CREATE TABLE IF NOT EXISTS md_results (
  run_id         INTEGER NOT NULL,
  ligand_id      INTEGER NOT NULL,
  stage          TEXT NOT NULL,   -- "em", "eq", "prod" など
  rmsd_bb        REAL,
  energy_min     REAL,
  stable_flag    INTEGER,
  report_relpath TEXT,
  PRIMARY KEY (run_id, ligand_id, stage),
  FOREIGN KEY (run_id)  REFERENCES runs(id),
  FOREIGN KEY (ligand_id) REFERENCES ligands(id)
);

CREATE TABLE IF NOT EXISTS pages_publish (
  id            INTEGER PRIMARY KEY AUTOINCREMENT,
  run_id        INTEGER NOT NULL,
  published_at  TEXT NOT NULL,
  commit_hash   TEXT,
  url           TEXT,
  FOREIGN KEY (run_id) REFERENCES runs(id)
);
EOSQL

echo "[INFO] applying schema.sql ..."
sqlite3 "$DB" < "$SCHEMA"

# === 3. 初期データ投入（targets / libraries） ===
sqlite3 "$DB" << 'EOSQL'
INSERT OR IGNORE INTO targets (code, name, pdb_id, receptor_path) VALUES
  ('aps_apoh',        'β2-glycoprotein I (APOH)', '1C1Z', 'catalog/targets/aps_apoh/receptor.pdb'),
  ('aps_prothrombin', 'Prothrombin (F2)',         '6C2W', 'catalog/targets/aps_prothrombin/receptor.pdb'),
  ('aps_tlr4_md2',    'TLR4/MD-2 complex',        '3FXI', 'catalog/targets/aps_tlr4_md2/receptor.pdb');

INSERT OR IGNORE INTO libraries (code, description) VALUES
  ('zinc_2d_smi_v1', 'ZINC downloader 2D SMILES (v1)');
EOSQL

echo "[INFO] DB initialized."

# === 4. APS 受容体選択レポートを書き出し ===
APS_REPORT="$ROOT_DIR/reports/2025-11-15_aps_receptor_selection.md"

cat > "$APS_REPORT" << 'EOF_MD'
# APS（抗リン脂質抗体症候群）ターゲット受容体選定レポート

- 作成日: 2025-11-15
- プロジェクト: OCP (Open Container Pipeline) – APS 解析系
- 対象疾患: 抗リン脂質抗体症候群（Antiphospholipid Syndrome; APS）

---

## 1. 背景

抗リン脂質抗体症候群（APS）は、抗リン脂質抗体（aPL）が持続的に陽性となり、
動静脈血栓や妊娠合併症を引き起こす自己免疫疾患である。
aPL は、陰性荷電リン脂質そのものだけでなく、リン脂質に結合した蛋白質
（β2-glycoprotein I, prothrombin, annexin など）を標的とすることが知られている。

本パイプラインでは「小分子ドッキング + MD 解析によって、APS 病態に関わる分子相互作用を
定量的に評価する」ことを目的として、以下の受容体を第 1 世代の標的セットとして採用する。

- β2-glycoprotein I（APOH）
- Prothrombin（F2）
- TLR4/MD-2 複合体

選定基準は以下の通り。

- APS 病態への関与が文献で繰り返し報告されている
- PDB に高品質な立体構造が存在し、ドッキンググリッドを定義しやすい
- 将来的に薬剤設計（阻害剤・修飾剤）が合理的に行えるドメイン構造を持つ

---

## 2. 個別ターゲットの選定理由

### 2.1 β2-glycoprotein I（APOH） – `aps_apoh`

- APS における最重要抗原の一つであり、抗 β2GPI 抗体だけで血栓形成を惹起しうる動物実験が多数報告されている。
- 第 V ドメインがリン脂質結合に重要であり、自己抗体はこの近傍を認識するとされる。
- 全長構造が解かれており、リン脂質結合部位や抗体エピトープ候補の立体配置を評価しやすい。
- 小分子によって β2GPI–リン脂質相互作用や β2GPI–自己抗体相互作用を遮断するという薬剤コンセプトを構築しやすい。

採用 PDB 構造:

- PDB ID: **1C1Z** — Human β2-glycoprotein I, full-length structure

OCP 内格納パス:

- `catalog/targets/aps_apoh/receptor.pdb`

---

### 2.2 Prothrombin（F2） – `aps_prothrombin`

- β2GPI と並ぶ APS の主要抗原であり、抗プロトロンビン抗体は血栓リスクと関連する。
- プロトロンビンは凝固カスケードの中心因子であり、異常な活性化は血栓形成と直結する。
- 全長プロトロンビンの closed コンフォメーション構造が報告されており、allosteric site を含む立体構造に基づいた設計が可能。

採用 PDB 構造:

- PDB ID: **6C2W** — Human prothrombin, closed conformation

OCP 内格納パス:

- `catalog/targets/aps_prothrombin/receptor.pdb`

---

### 2.3 TLR4/MD-2 複合体 – `aps_tlr4_md2`

- APS 患者由来 aPL が TLR4 経路を介して炎症・血栓促進シグナルを誘導することが報告されている。
- TLR4/MD-2 は LPS などの脂質リガンドを認識するシグナル受容体であり、既に小分子アンタゴニストの構造情報が豊富に存在する。
- aPL 依存的に活性化される TLR4 シグナルを小分子で抑制することで、炎症・血栓の二次的増悪を抑える戦略が立てやすい。

採用 PDB 構造:

- PDB ID: **3FXI** — Human TLR4–MD-2–LPS complex

OCP 内格納パス:

- `catalog/targets/aps_tlr4_md2/receptor.pdb`

---

## 3. 今後の拡張候補

- Annexin A5（ANXA5）: aPL による annexin A5 シールド破壊が APS 血栓形成に関与するとされる。
- 補体 C5 / C5aR: APS における補体活性化が病態形成に関与することが示されており、既に C5 阻害薬が臨床応用されている。

これらは現時点では OCP パイプラインの「主要 3 ターゲット」には含めていないが、
将来のモジュール追加候補として検討する。
EOF_MD

echo "[INFO] APS receptor selection report written: $APS_REPORT"

# === 5. データフロー方針書を書き出し ===
POLICY_REPORT="$ROOT_DIR/reports/policies/2025-11-15_data_flow_policy.md"

cat > "$POLICY_REPORT" << 'EOF_MD'
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
EOF_MD

echo "[INFO] Policy report written: $POLICY_REPORT"

echo "[INFO] All done. DB とレポート類の初期化が完了しました。"
