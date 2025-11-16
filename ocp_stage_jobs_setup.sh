#!/usr/bin/env bash
#
# ocp_stage_jobs_setup.sh
#
# OCP の「ステージごとの Pod / Job」用エントリポイントスクリプト群をまとめて生成する。
#
# ポリシー:
#   - すべてのステージは SQLite DB (catalog/db/ocp_results.sqlite) を参照して
#     「何を処理するか」を決定する（DB = source of truth）。
#   - GROMACS ステージは CUDA サポート付き GPU 計算を前提としたコマンド雛形を持つ。
#   - 生の化合物ライブラリは RDKit による前処理済みで DB に登録されている前提。
#
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS_DIR="${ROOT}/scripts"
REPORTS_DESIGN_DIR="${ROOT}/reports/design"

mkdir -p "${SCRIPTS_DIR}" "${REPORTS_DESIGN_DIR}"

DB_DEFAULT_PATH="${ROOT}/catalog/db/ocp_results.sqlite"

############################################################
# 共通ヘルパーコメント
############################################################
COMMON_DB_HEADER='#
# 共通仕様:
#   - RUN_ID は第1引数、または環境変数 RUN_ID から取得する。
#   - DB_PATH は環境変数 DB_PATH_OVERRIDE で上書き可能。
#   - DB は catalog/db/ocp_results.sqlite を前提とし、
#     ここから対象 ligand / target / run を決める（DB = source of truth）。
#'

############################################################
# ligand-selector Job エントリポイント
############################################################
cat > "${SCRIPTS_DIR}/k8s_job_ligand_selector.sh" <<EOF
#!/usr/bin/env bash
#
# k8s_job_ligand_selector.sh
#
# 役割:
#   - ZINC 化合物ライブラリ（RDKit 整形済み & DB 登録済み）から、
#     対象ライブラリ / RUN 用の ligand を選択し、
#     /workspace/ligands_raw にコピーするステージ。
#   - どの ligand を選ぶかは DB を見て決定する。
${COMMON_DB_HEADER}
set -euo pipefail

RUN_ID="\${1:-\${RUN_ID:-}}"
if [[ -z "\${RUN_ID}" ]]; then
  echo "Usage: \$0 <RUN_ID>  # or set RUN_ID env" >&2
  exit 1
fi

SCRIPT_DIR="\$(cd "\$(dirname "\${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="\${SCRIPT_DIR%/scripts}"
DB_PATH="\${DB_PATH_OVERRIDE:-${DB_DEFAULT_PATH}}"

WORKSPACE="\${WORKSPACE:-/workspace}"
LIGANDS_RAW_DIR="\${WORKSPACE}/ligands_raw"

mkdir -p "\${LIGANDS_RAW_DIR}"

if [[ ! -f "\${DB_PATH}" ]]; then
  echo "ERROR: DB not found: \${DB_PATH}" >&2
  exit 1
fi

if ! command -v sqlite3 >/dev/null 2>&1; then
  echo "ERROR: sqlite3 command not found." >&2
  exit 1
fi

echo "[INFO] ligand-selector: RUN_ID=\${RUN_ID}"
echo "[INFO] DB_PATH       : \${DB_PATH}"
echo "[INFO] LIGANDS_RAW   : \${LIGANDS_RAW_DIR}"

echo "[INFO] Selecting ligands for RUN_ID=\${RUN_ID} from DB..."

# 想定スキーマ:
#   runs(id, library_id, target_id, ...)
#   run_ligands(run_id, ligand_id, order_index, ...)
#   ligands(id, zinc_id, source_file, has_3d, conformer_method, ...)
#
# RDKit による整形済み ligand のみを対象:
#   - has_3d = 1
#   - conformer_method = 'rdkit_etkdg'
#
sqlite3 "\${DB_PATH}" <<SQL
.headers on
.mode column
SELECT
  rl.order_index,
  l.id AS ligand_id,
  l.zinc_id,
  l.source_file,
  l.has_3d,
  l.conformer_method
FROM run_ligands rl
JOIN ligands l ON l.id = rl.ligand_id
WHERE rl.run_id = \${RUN_ID}
  AND l.has_3d = 1
  AND l.conformer_method = 'rdkit_etkdg'
ORDER BY rl.order_index;
SQL

echo "[INFO] Copying selected ligands into \${LIGANDS_RAW_DIR} (if source_file exists)..."

sqlite3 -csv "\${DB_PATH}" <<SQL | while IFS=',' read -r order_index ligand_id zinc_id source_file has_3d conformer_method; do
SELECT
  rl.order_index,
  l.id AS ligand_id,
  l.zinc_id,
  l.source_file,
  l.has_3d,
  l.conformer_method
FROM run_ligands rl
JOIN ligands l ON l.id = rl.ligand_id
WHERE rl.run_id = \${RUN_ID}
  AND l.has_3d = 1
  AND l.conformer_method = 'rdkit_etkdg'
ORDER BY rl.order_index;
SQL
do
  if [[ -z "\${source_file}" || "\${source_file}" = "NULL" ]]; then
    echo "[WARN] ligand_id=\${ligand_id} zinc_id=\${zinc_id}: source_file is NULL, skipping." >&2
    continue
  fi
  SRC_PATH="\${REPO_ROOT}/\${source_file}"
  if [[ ! -f "\${SRC_PATH}" ]]; then
    echo "[WARN] ligand_id=\${ligand_id} zinc_id=\${zinc_id}: source file not found: \${SRC_PATH}, skipping." >&2
    continue
  fi
  BASENAME="\$(basename "\${SRC_PATH}")"
  DEST_PATH="\${LIGANDS_RAW_DIR}/\${zinc_id}_\${BASENAME}"
  echo "[INFO] copy: \${SRC_PATH} -> \${DEST_PATH}"
  cp "\${SRC_PATH}" "\${DEST_PATH}"
done

echo "[INFO] ligand-selector finished."
EOF

chmod +x "${SCRIPTS_DIR}/k8s_job_ligand_selector.sh"

############################################################
# vina-prep Job エントリポイント
############################################################
cat > "${SCRIPTS_DIR}/k8s_job_vina_prep.sh" <<EOF
#!/usr/bin/env bash
#
# k8s_job_vina_prep.sh
#
# 役割:
#   - Receptor PDB → PDBQT 変換（MGLTools を想定）
#   - RDKit 整形済み ligand PDB/SDF から ligand PDBQT 生成
#   - /workspace/vina_input に receptor.pdbqt + ligands PDBQT を並べる。
${COMMON_DB_HEADER}
set -euo pipefail

RUN_ID="\${1:-\${RUN_ID:-}}"
if [[ -z "\${RUN_ID}" ]]; then
  echo "Usage: \$0 <RUN_ID>  # or set RUN_ID env" >&2
  exit 1
fi

SCRIPT_DIR="\$(cd "\$(dirname "\${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="\${SCRIPT_DIR%/scripts}"
DB_PATH="\${DB_PATH_OVERRIDE:-${DB_DEFAULT_PATH}}"

WORKSPACE="\${WORKSPACE:-/workspace}"
VINA_INPUT_DIR="\${WORKSPACE}/vina_input"
LIGANDS_RAW_DIR="\${WORKSPACE}/ligands_raw"

mkdir -p "\${VINA_INPUT_DIR}"

if [[ ! -f "\${DB_PATH}" ]]; then
  echo "ERROR: DB not found: \${DB_PATH}" >&2
  exit 1
fi

echo "[INFO] vina-prep: RUN_ID=\${RUN_ID}"
echo "[INFO] DB_PATH   : \${DB_PATH}"
echo "[INFO] VINA_INPUT: \${VINA_INPUT_DIR}"

# Receptor 情報取得:
#   runs.target_id -> targets.id -> targets.code, targets.receptor_pdb_path などを想定。
sqlite3 "\${DB_PATH}" <<SQL
.headers on
.mode column
SELECT
  r.id         AS run_id,
  t.id         AS target_id,
  t.code       AS target_code,
  t.receptor_pdb_path,
  t.receptor_pdbqt_path,
  t.receptor_preparation_method
FROM runs r
JOIN targets t ON t.id = r.target_id
WHERE r.id = \${RUN_ID};
SQL

# 実際の receptor.pdb -> receptor.pdbqt の変換は MGLTools に委ねる。
# ここではパスの存在チェックと雛形コマンドのみ。
RECEPTOR_PDB_PATH="\$(sqlite3 "\${DB_PATH}" "SELECT t.receptor_pdb_path FROM runs r JOIN targets t ON t.id = r.target_id WHERE r.id = \${RUN_ID} LIMIT 1;")"

if [[ -z "\${RECEPTOR_PDB_PATH}" || "\${RECEPTOR_PDB_PATH}" = "NULL" ]]; then
  echo "[WARN] receptor_pdb_path is not set in DB; please update targets.receptor_pdb_path." >&2
else
  ABS_RECEPTOR_PDB="\${REPO_ROOT}/\${RECEPTOR_PDB_PATH}"
  echo "[INFO] receptor PDB: \${ABS_RECEPTOR_PDB}"
  if [[ -f "\${ABS_RECEPTOR_PDB}" ]]; then
    echo "[INFO] (NOTE) ここで MGLTools の prepare_receptor4.py を呼び出して receptor.pdbqt を生成する想定です。"
    echo "[INFO] 例:"
    echo "  prepare_receptor4.py -r \${ABS_RECEPTOR_PDB} -o \${VINA_INPUT_DIR}/receptor.pdbqt"
  else
    echo "[WARN] receptor PDB not found at: \${ABS_RECEPTOR_PDB}" >&2
  fi
fi

echo "[INFO] Preparing ligand PDBQT from \${LIGANDS_RAW_DIR} into \${VINA_INPUT_DIR} ..."

if ! command -v obabel >/dev/null 2>&1; then
  echo "[WARN] obabel not found; skipping actual PDBQT conversion. This is a skeleton." >&2
else
  for f in "\${LIGANDS_RAW_DIR}"/*; do
    [[ -e "\$f" ]] || continue
    base="\$(basename "\$f")"
    zinc_id="\${base%%_*}"
    out_pdbqt="\${VINA_INPUT_DIR}/\${zinc_id}.pdbqt"
    echo "[INFO] obabel: \$f -> \${out_pdbqt}"
    obabel "\$f" -O "\${out_pdbqt}" >/dev/null 2>&1 || echo "[WARN] obabel failed for \$f" >&2
  done
fi

echo "[INFO] vina-prep finished."
EOF

chmod +x "${SCRIPTS_DIR}/k8s_job_vina_prep.sh"

############################################################
# vina2gromacs-prep Job エントリポイント
############################################################
cat > "${SCRIPTS_DIR}/k8s_job_vina2gromacs_prep.sh" <<EOF
#!/usr/bin/env bash
#
# k8s_job_vina2gromacs_prep.sh
#
# 役割:
#   - Vina 結果 (vina_results) を DB から読み、
#     top ヒットの pose を GROMACS 用入力に変換する準備ステージ。
#   - 具体的には:
#       - top N (例: 1〜10) の pose PDBQT を取り出す
#       - receptor + ligand 複合体 PDB を生成
#       - ligand の topology (acpype 等) を作る
#   - ここでは DB 駆動 & ファイル配置の骨格のみ実装する。
${COMMON_DB_HEADER}
set -euo pipefail

RUN_ID="\${1:-\${RUN_ID:-}}"
TOP_N="\${TOP_N:-1}"  # デフォルトは top 1 ヒット
if [[ -z "\${RUN_ID}" ]]; then
  echo "Usage: \$0 <RUN_ID>  # or set RUN_ID env" >&2
  exit 1
fi

SCRIPT_DIR="\$(cd "\$(dirname "\${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="\${SCRIPT_DIR%/scripts}"
DB_PATH="\${DB_PATH_OVERRIDE:-${DB_DEFAULT_PATH}}"

WORKSPACE="\${WORKSPACE:-/workspace}"
VINA_OUTPUT_DIR="\${WORKSPACE}/vina_output"
GMX_INPUT_DIR="\${WORKSPACE}/gmx_input"

mkdir -p "\${GMX_INPUT_DIR}"

if [[ ! -f "\${DB_PATH}" ]]; then
  echo "ERROR: DB not found: \${DB_PATH}" >&2
  exit 1
fi

echo "[INFO] vina2gromacs-prep: RUN_ID=\${RUN_ID}"
echo "[INFO] DB_PATH          : \${DB_PATH}"
echo "[INFO] VINA_OUTPUT_DIR  : \${VINA_OUTPUT_DIR}"
echo "[INFO] GMX_INPUT_DIR    : \${GMX_INPUT_DIR}"
echo "[INFO] TOP_N            : \${TOP_N}"

echo "[INFO] Top \${TOP_N} hits from vina_results (by affinity_kcal ASC):"
sqlite3 "\${DB_PATH}" <<SQL
.headers on
.mode column
SELECT
  vr.id,
  vr.run_id,
  vr.ligand_id,
  l.zinc_id,
  vr.affinity_kcal,
  vr.pose_path
FROM vina_results vr
JOIN ligands l ON l.id = vr.ligand_id
WHERE vr.run_id = \${RUN_ID}
ORDER BY vr.affinity_kcal ASC
LIMIT \${TOP_N};
SQL

# ここで本来は:
#   - pose_path (PDBQT) を PDB に変換 (obabel)
#   - receptor PDB と merge して複合体 PDB を作成
#   - acpype で ligand トポロジー生成
#   を行う。ここではファイルパスを workspace にコピーする骨格のみ。

sqlite3 -csv "\${DB_PATH}" <<SQL | while IFS=',' read -r vr_id run_id ligand_id zinc_id affinity_kcal pose_path; do
SELECT
  vr.id,
  vr.run_id,
  vr.ligand_id,
  l.zinc_id,
  vr.affinity_kcal,
  vr.pose_path
FROM vina_results vr
JOIN ligands l ON l.id = vr.ligand_id
WHERE vr.run_id = \${RUN_ID}
ORDER BY vr.affinity_kcal ASC
LIMIT \${TOP_N};
SQL
do
  if [[ -z "\${pose_path}" || "\${pose_path}" = "NULL" ]]; then
    echo "[WARN] pose_path is NULL for zinc_id=\${zinc_id}, skipping." >&2
    continue
  fi
  SRC_PDBQT="\${VINA_OUTPUT_DIR}/\${pose_path##*/}"
  if [[ ! -f "\${SRC_PDBQT}" ]]; then
    # DB に相対パスが入っている想定 (vina_out/poses/...), REPO_ROOT 由来のケースも考慮
    ALT_PDBQT="\${REPO_ROOT}/\${pose_path}"
    if [[ -f "\${ALT_PDBQT}" ]]; then
      SRC_PDBQT="\${ALT_PDBQT}"
    else
      echo "[WARN] pose file not found for zinc_id=\${zinc_id}: \${SRC_PDBQT} or \${ALT_PDBQT}" >&2
      continue
    fi
  fi
  DEST_PDBQT="\${GMX_INPUT_DIR}/\${zinc_id}_pose.pdbqt"
  echo "[INFO] copy pose: \${SRC_PDBQT} -> \${DEST_PDBQT}"
  cp "\${SRC_PDBQT}" "\${DEST_PDBQT}"
done

echo "[INFO] (NOTE) ここで obabel により PDBQT -> PDB 変換、および acpype によるトポロジー生成を行う想定です。"
echo "[INFO] vina2gromacs-prep finished."
EOF

chmod +x "${SCRIPTS_DIR}/k8s_job_vina2gromacs_prep.sh"

############################################################
# gmx-makebox Job エントリポイント (GPU/CUDA 前提)
############################################################
cat > "${SCRIPTS_DIR}/k8s_job_gmx_makebox.sh" <<EOF
#!/usr/bin/env bash
#
# k8s_job_gmx_makebox.sh
#
# 役割:
#   - GROMACS によるシミュレーションボックス作成ステージ。
#   - 入力: /workspace/gmx_input 内の receptor + ligand 複合体 (例: complex.gro or .pdb)。
#   - 出力: /workspace/gmx_em (最小化前の tpr など)。
#   - GPU (CUDA) サポート付き GROMACS を前提とする。
${COMMON_DB_HEADER}
set -euo pipefail

RUN_ID="\${1:-\${RUN_ID:-}}"
if [[ -z "\${RUN_ID}" ]]; then
  echo "Usage: \$0 <RUN_ID>  # or set RUN_ID env" >&2
  exit 1
fi

WORKSPACE="\${WORKSPACE:-/workspace}"
GMX_INPUT_DIR="\${WORKSPACE}/gmx_input"
GMX_EM_DIR="\${WORKSPACE}/gmx_em"

mkdir -p "\${GMX_EM_DIR}"

echo "[INFO] gmx-makebox: RUN_ID=\${RUN_ID}"
echo "[INFO] GMX_INPUT : \${GMX_INPUT_DIR}"
echo "[INFO] GMX_EM    : \${GMX_EM_DIR}"

if ! command -v gmx_mpi >/dev/null 2>&1 && ! command -v gmx >/dev/null 2>&1; then
  echo "[WARN] gmx_mpi / gmx not found; this is a skeleton script." >&2
  exit 0
fi

COMPLEX_FILE="\${GMX_INPUT_DIR}/complex.gro"
if [[ ! -f "\${COMPLEX_FILE}" ]]; then
  echo "[WARN] complex.gro not found in \${GMX_INPUT_DIR}; please generate it in vina2gromacs-prep." >&2
fi

echo "[INFO] (NOTE) 以下は CUDA GPU を利用する GROMACS コマンドの例です。"
echo "[INFO]       実際の mdp ファイル名やオプションは適宜調整してください。"

cat <<'CMD'
# 例: ボックス作成
# gmx editconf -f complex.gro -o boxed.gro -c -d 1.0 -bt dodecahedron

# 例: 溶媒追加
# gmx solvate -cp boxed.gro -cs spc216.gro -o solv.gro -p topol.top

# 例: 最小化用 tpr 作成
# gmx grompp -f minim.mdp -c solv.gro -p topol.top -o em.tpr

# ここでは GPU を使うのは次ステージ(gmx-minimize)だが、
# コンテナイメージ自体は CUDA サポート付き GROMACS を前提とする。
CMD

echo "[INFO] gmx-makebox finished (skeleton)."
EOF

chmod +x "${SCRIPTS_DIR}/k8s_job_gmx_makebox.sh"

############################################################
# gmx-minimize Job エントリポイント (GPU/CUDA 前提)
############################################################
cat > "${SCRIPTS_DIR}/k8s_job_gmx_minimize.sh" <<EOF
#!/usr/bin/env bash
#
# k8s_job_gmx_minimize.sh
#
# 役割:
#   - GROMACS によるエネルギー最小化ステージ。
#   - 入力: /workspace/gmx_em/em.tpr
#   - 出力: em.gro, em.edr, em.log など。
#   - GPU (CUDA) サポート付き GROMACS を前提とし、GPU で mdrun を実行する。
${COMMON_DB_HEADER}
set -euo pipefail

RUN_ID="\${1:-\${RUN_ID:-}}"
if [[ -z "\${RUN_ID}" ]]; then
  echo "Usage: \$0 <RUN_ID>  # or set RUN_ID env" >&2
  exit 1
fi

WORKSPACE="\${WORKSPACE:-/workspace}"
GMX_EM_DIR="\${WORKSPACE}/gmx_em"

GMX_NTASKS="\${GMX_NTASKS:-4}"

echo "[INFO] gmx-minimize: RUN_ID=\${RUN_ID}"
echo "[INFO] GMX_EM      : \${GMX_EM_DIR}"
echo "[INFO] GMX_NTASKS  : \${GMX_NTASKS}"

if ! command -v gmx_mpi >/dev/null 2>&1 && ! command -v gmx >/dev/null 2>&1; then
  echo "[WARN] gmx_mpi / gmx not found; this is a skeleton script." >&2
  exit 0
fi

cd "\${GMX_EM_DIR}" || exit 1

echo "[INFO] (NOTE) CUDA GPU を利用した mdrun の例:"

cat <<CMD
# mpirun -np \${GMX_NTASKS} gmx_mpi mdrun \\
#   -deffnm em \\
#   -ntmpi \${GMX_NTASKS} \\
#   -nb gpu -pme gpu -pmefft gpu -bonded gpu -update gpu
CMD

echo "[INFO] gmx-minimize finished (skeleton)."
EOF

chmod +x "${SCRIPTS_DIR}/k8s_job_gmx_minimize.sh"

############################################################
# min-analysis-pages Job エントリポイント
############################################################
cat > "${SCRIPTS_DIR}/k8s_job_min_analysis_pages.sh" <<EOF
#!/usr/bin/env bash
#
# k8s_job_min_analysis_pages.sh
#
# 役割:
#   - 最小化結果を評価し、GitHub Pages 用レポートを生成するステージ。
#   - 入力: /workspace/gmx_em/em.edr, em.log など。
#   - 出力: docs/minimization/ 配下の Markdown/HTML/PNG、Git コミット & push（将来的）。
${COMMON_DB_HEADER}
set -euo pipefail

RUN_ID="\${1:-\${RUN_ID:-}}"
if [[ -z "\${RUN_ID}" ]]; then
  echo "Usage: \$0 <RUN_ID>  # or set RUN_ID env" >&2
  exit 1
fi

SCRIPT_DIR="\$(cd "\$(dirname "\${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="\${SCRIPT_DIR%/scripts}"
WORKSPACE="\${WORKSPACE:-/workspace}"
GMX_EM_DIR="\${WORKSPACE}/gmx_em"
DOCS_DIR="\${REPO_ROOT}/docs/minimization"

mkdir -p "\${DOCS_DIR}"

echo "[INFO] min-analysis-pages: RUN_ID=\${RUN_ID}"
echo "[INFO] GMX_EM  : \${GMX_EM_DIR}"
echo "[INFO] DOCS_DIR: \${DOCS_DIR}"

# ここで gmx energy / gmx rms などを用いて最小化結果を解析し、
# 画像 (PNG) やテキストを生成して docs/ に配置する想定。
REPORT_MD="\${DOCS_DIR}/run_\${RUN_ID}_min_report.md"

cat > "\${REPORT_MD}" <<RPT
# Minimization Report (RUN_ID=\${RUN_ID})

- Workspace: \${GMX_EM_DIR}
- This is a skeleton report. Use GROMACS tools to extract energy vs step,
  and embed plots here in future.

RPT

echo "[INFO] Generated report: \${REPORT_MD}"
echo "[INFO] min-analysis-pages finished (skeleton)."
EOF

chmod +x "${SCRIPTS_DIR}/k8s_job_min_analysis_pages.sh"

############################################################
# gmx-equilibration Job エントリポイント (GPU/CUDA 前提)
############################################################
cat > "${SCRIPTS_DIR}/k8s_job_gmx_equilibration.sh" <<EOF
#!/usr/bin/env bash
#
# k8s_job_gmx_equilibration.sh
#
# 役割:
#   - GROMACS による平衡化 (NVT / NPT) ステージ。
#   - 入力: 最小化後の構造 (em.gro 等)。
#   - 出力: /workspace/gmx_eq (eq.tpr, eq.trr, eq.edr, eq.gro 等)。
#   - CUDA GPU サポート付き GROMACS 前提。
${COMMON_DB_HEADER}
set -euo pipefail

RUN_ID="\${1:-\${RUN_ID:-}}"
if [[ -z "\${RUN_ID}" ]]; then
  echo "Usage: \$0 <RUN_ID>  # or set RUN_ID env" >&2
  exit 1
fi

WORKSPACE="\${WORKSPACE:-/workspace}"
GMX_EQ_DIR="\${WORKSPACE}/gmx_eq"

GMX_NTASKS="\${GMX_NTASKS:-4}"

mkdir -p "\${GMX_EQ_DIR}"

echo "[INFO] gmx-equilibration: RUN_ID=\${RUN_ID}"
echo "[INFO] GMX_EQ      : \${GMX_EQ_DIR}"
echo "[INFO] GMX_NTASKS  : \${GMX_NTASKS}"

if ! command -v gmx_mpi >/dev/null 2>&1 && ! command -v gmx >/dev/null 2>&1; then
  echo "[WARN] gmx_mpi / gmx not found; this is a skeleton script." >&2
  exit 0
fi

cd "\${GMX_EQ_DIR}" || exit 1

cat <<CMD
# 例: NVT / NPT 平衡化 (CUDA GPU 利用)
# mpirun -np \${GMX_NTASKS} gmx_mpi mdrun \\
#   -deffnm nvt \\
#   -ntmpi \${GMX_NTASKS} \\
#   -nb gpu -pme gpu -pmefft gpu -bonded gpu -update gpu
#
# mpirun -np \${GMX_NTASKS} gmx_mpi mdrun \\
#   -deffnm npt \\
#   -ntmpi \${GMX_NTASKS} \\
#   -nb gpu -pme gpu -pmefft gpu -bonded gpu -update gpu
CMD

echo "[INFO] gmx-equilibration finished (skeleton)."
EOF

chmod +x "${SCRIPTS_DIR}/k8s_job_gmx_equilibration.sh"

############################################################
# eq-analysis-pages Job エントリポイント
############################################################
cat > "${SCRIPTS_DIR}/k8s_job_eq_analysis_pages.sh" <<EOF
#!/usr/bin/env bash
#
# k8s_job_eq_analysis_pages.sh
#
# 役割:
#   - 平衡化結果 (NVT / NPT) を解析し、GitHub Pages 用レポートを生成するステージ。
${COMMON_DB_HEADER}
set -euo pipefail

RUN_ID="\${1:-\${RUN_ID:-}}"
if [[ -z "\${RUN_ID}" ]]; then
  echo "Usage: \$0 <RUN_ID>  # or set RUN_ID env" >&2
  exit 1
fi

SCRIPT_DIR="\$(cd "\$(dirname "\${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="\${SCRIPT_DIR%/scripts}"
WORKSPACE="\${WORKSPACE:-/workspace}"
GMX_EQ_DIR="\${WORKSPACE}/gmx_eq"
DOCS_DIR="\${REPO_ROOT}/docs/equilibration"

mkdir -p "\${DOCS_DIR}"

echo "[INFO] eq-analysis-pages: RUN_ID=\${RUN_ID}"
echo "[INFO] GMX_EQ  : \${GMX_EQ_DIR}"
echo "[INFO] DOCS_DIR: \${DOCS_DIR}"

REPORT_MD="\${DOCS_DIR}/run_\${RUN_ID}_eq_report.md"

cat > "\${REPORT_MD}" <<RPT
# Equilibration Report (RUN_ID=\${RUN_ID})

- Workspace: \${GMX_EQ_DIR}
- This is a skeleton report. Use GROMACS tools to extract temperature / pressure /
  energy vs time and embed plots here in future.

RPT

echo "[INFO] Generated report: \${REPORT_MD}"
echo "[INFO] eq-analysis-pages finished (skeleton)."
EOF

chmod +x "${SCRIPTS_DIR}/k8s_job_eq_analysis_pages.sh"

############################################################
# gmx-production Job エントリポイント (GPU/CUDA 前提)
############################################################
cat > "${SCRIPTS_DIR}/k8s_job_gmx_production.sh" <<EOF
#!/usr/bin/env bash
#
# k8s_job_gmx_production.sh
#
# 役割:
#   - GROMACS による本番 Production run ステージ。
#   - 入力: 平衡化後の構造 (eq.gro 等)。
#   - 出力: /workspace/gmx_md (md.xtc, md.edr, md.log 等)。
#   - CUDA GPU サポート付き GROMACS 前提。
${COMMON_DB_HEADER}
set -euo pipefail

RUN_ID="\${1:-\${RUN_ID:-}}"
if [[ -z "\${RUN_ID}" ]]; then
  echo "Usage: \$0 <RUN_ID>  # or set RUN_ID env" >&2
  exit 1
fi

WORKSPACE="\${WORKSPACE:-/workspace}"
GMX_MD_DIR="\${WORKSPACE}/gmx_md"

GMX_NTASKS="\${GMX_NTASKS:-4}"

mkdir -p "\${GMX_MD_DIR}"

echo "[INFO] gmx-production: RUN_ID=\${RUN_ID}"
echo "[INFO] GMX_MD      : \${GMX_MD_DIR}"
echo "[INFO] GMX_NTASKS  : \${GMX_NTASKS}"

if ! command -v gmx_mpi >/dev/null 2>&1 && ! command -v gmx >/dev/null 2>&1; then
  echo "[WARN] gmx_mpi / gmx not found; this is a skeleton script." >&2
  exit 0
fi

cd "\${GMX_MD_DIR}" || exit 1

cat <<CMD
# 例: Production run (CUDA GPU 利用)
# mpirun -np \${GMX_NTASKS} gmx_mpi mdrun \\
#   -deffnm md \\
#   -ntmpi \${GMX_NTASKS} \\
#   -nb gpu -pme gpu -pmefft gpu -bonded gpu -update gpu
CMD

echo "[INFO] gmx-production finished (skeleton)."
EOF

chmod +x "${SCRIPTS_DIR}/k8s_job_gmx_production.sh"

############################################################
# prod-analysis-pages Job エントリポイント
############################################################
cat > "${SCRIPTS_DIR}/k8s_job_prod_analysis_pages.sh" <<EOF
#!/usr/bin/env bash
#
# k8s_job_prod_analysis_pages.sh
#
# 役割:
#   - Production run の結果を解析し、GitHub Pages 用に
#     結合安定性評価やランキングなどを出力するステージ。
${COMMON_DB_HEADER}
set -euo pipefail

RUN_ID="\${1:-\${RUN_ID:-}}"
if [[ -z "\${RUN_ID}" ]]; then
  echo "Usage: \$0 <RUN_ID>  # or set RUN_ID env" >&2
  exit 1
fi

SCRIPT_DIR="\$(cd "\$(dirname "\${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="\${SCRIPT_DIR%/scripts}"
WORKSPACE="\${WORKSPACE:-/workspace}"
GMX_MD_DIR="\${WORKSPACE}/gmx_md"
DOCS_DIR="\${REPO_ROOT}/docs/production"

mkdir -p "\${DOCS_DIR}"

echo "[INFO] prod-analysis-pages: RUN_ID=\${RUN_ID}"
echo "[INFO] GMX_MD  : \${GMX_MD_DIR}"
echo "[INFO] DOCS_DIR: \${DOCS_DIR}"

REPORT_MD="\${DOCS_DIR}/run_\${RUN_ID}_prod_report.md"

cat > "\${REPORT_MD}" <<RPT
# Production MD Report (RUN_ID=\${RUN_ID})

- Workspace: \${GMX_MD_DIR}
- This is a skeleton report. Use GROMACS tools and analysis scripts
  to compute binding stability metrics and summarize them here.
- 最終的には、RUN_ID ごとのランキングや可視化を GitHub Pages
  トップページにリンクさせる想定。

RPT

echo "[INFO] Generated report: \${REPORT_MD}"
echo "[INFO] prod-analysis-pages finished (skeleton)."
EOF

chmod +x "${SCRIPTS_DIR}/k8s_job_prod_analysis_pages.sh"

############################################################
# デザインメモ（このフェーズの方針を固定）
############################################################
cat > "${REPORTS_DESIGN_DIR}/2025-11-16_stage_jobs_gpu_db_design.md" <<'EOF'
# ステージ Job / Pod 設計メモ（DB 駆動 + GPU GROMACS 方針）

日付: 2025-11-16  
担当: めけめて もすもす / ChatGPT 補助

## 基本方針

- すべてのステージ (ligand-selector / vina-prep / vina2gromacs-prep / gmx-* / *-analysis-pages) は、
  **SQLite DB (catalog/db/ocp_results.sqlite)** を source of truth として扱う。
  - 計算対象 (RUN / ligand / target) は DB から取得する。
  - ファイルは再生成可能なキャッシュとして扱う。

- GROMACS を使用するステージ (gmx-*) は、
  **必ず CUDA サポート付きの GPU 計算を前提としたイメージ** を使う。
  - スクリプトでは `gmx_mpi` / `mdrun` の GPU オプションを雛形として記述。

- 生の化合物ライブラリは、
  **RDKit による 3D 化・整形を通した上で DB に格納しておく**。
  - `ligands.has_3d = 1` と `ligands.conformer_method = 'rdkit_etkdg'` のようなフラグで管理。
  - ligand-selector はこのフラグを見て、3D 化済みのものだけを候補にする。

## このセットアップで作成されたスクリプト

- `scripts/k8s_job_ligand_selector.sh`
- `scripts/k8s_job_vina_prep.sh`
- `scripts/k8s_job_vina2gromacs_prep.sh`
- `scripts/k8s_job_gmx_makebox.sh`
- `scripts/k8s_job_gmx_minimize.sh`
- `scripts/k8s_job_min_analysis_pages.sh`
- `scripts/k8s_job_gmx_equilibration.sh`
- `scripts/k8s_job_eq_analysis_pages.sh`
- `scripts/k8s_job_gmx_production.sh`
- `scripts/k8s_job_prod_analysis_pages.sh`

これらはすべて **RUN_ID を引数または環境変数で受け取り、DB を参照する骨格** を持つ。
GROMACS 部分は CUDA GPU 前提のコマンドをコメントとして記述しており、
今後の実装で具体的な mdp / tpr ファイル名などを埋めていくことを想定している。

EOF

echo "[OK] Created stage job entrypoint skeletons under scripts/"
echo "[OK] Design memo: reports/design/2025-11-16_stage_jobs_gpu_db_design.md"
