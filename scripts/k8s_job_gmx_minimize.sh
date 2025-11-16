#!/usr/bin/env bash
#
# k8s_job_gmx_minimize.sh
#
# 役割:
#   - GROMACS によるエネルギー最小化ステージ。
#   - 入力: /workspace/gmx_em/em.tpr
#   - 出力: em.gro, em.edr, em.log など。
#   - GPU (CUDA) サポート付き GROMACS を前提とし、GPU で mdrun を実行する。
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
GMX_EM_DIR="${WORKSPACE}/gmx_em"

GMX_NTASKS="${GMX_NTASKS:-4}"

echo "[INFO] gmx-minimize: RUN_ID=${RUN_ID}"
echo "[INFO] GMX_EM      : ${GMX_EM_DIR}"
echo "[INFO] GMX_NTASKS  : ${GMX_NTASKS}"

if ! command -v gmx_mpi >/dev/null 2>&1 && ! command -v gmx >/dev/null 2>&1; then
  echo "[WARN] gmx_mpi / gmx not found; this is a skeleton script." >&2
  exit 0
fi

cd "${GMX_EM_DIR}" || exit 1

echo "[INFO] (NOTE) CUDA GPU を利用した mdrun の例:"

cat <<CMD
# mpirun -np ${GMX_NTASKS} gmx_mpi mdrun \
#   -deffnm em \
#   -ntmpi ${GMX_NTASKS} \
#   -nb gpu -pme gpu -pmefft gpu -bonded gpu -update gpu
CMD

echo "[INFO] gmx-minimize finished (skeleton)."
