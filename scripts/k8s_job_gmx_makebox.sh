#!/usr/bin/env bash
#
# k8s_job_gmx_makebox.sh
#
# 役割:
#   - GROMACS によるシミュレーションボックス作成ステージ。
#   - 入力: /workspace/gmx_input 内の receptor + ligand 複合体 (例: complex.gro or .pdb)。
#   - 出力: /workspace/gmx_em (最小化前の tpr など)。
#   - GPU (CUDA) サポート付き GROMACS を前提とする。
#
# 共通仕様:
#   - RUN_ID は第1引数、または環境変数 RUN_ID から取得する。
#   - DB_PATH は環境変数 DB_PATH_OVERRIDE で上書き可能。
#   - DB は catalog/db/ocp_results.sqlite を前提とし、
#     ここから対象 ligand / target / run を決める（DB = source of truth）。
#
set -euo pipefail

RUN_ID="${1:-${RUN_ID:-}}"
if [[ -z "${RUN_ID}" ]]; then
  echo "Usage: $0 <RUN_ID>  # or set RUN_ID env" >&2
  exit 1
fi

WORKSPACE="${WORKSPACE:-/workspace}"
GMX_INPUT_DIR="${WORKSPACE}/gmx_input"
GMX_EM_DIR="${WORKSPACE}/gmx_em"

mkdir -p "${GMX_EM_DIR}"

echo "[INFO] gmx-makebox: RUN_ID=${RUN_ID}"
echo "[INFO] GMX_INPUT : ${GMX_INPUT_DIR}"
echo "[INFO] GMX_EM    : ${GMX_EM_DIR}"

if ! command -v gmx_mpi >/dev/null 2>&1 && ! command -v gmx >/dev/null 2>&1; then
  echo "[WARN] gmx_mpi / gmx not found; this is a skeleton script." >&2
  exit 0
fi

COMPLEX_FILE="${GMX_INPUT_DIR}/complex.gro"
if [[ ! -f "${COMPLEX_FILE}" ]]; then
  echo "[WARN] complex.gro not found in ${GMX_INPUT_DIR}; please generate it in vina2gromacs-prep." >&2
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
