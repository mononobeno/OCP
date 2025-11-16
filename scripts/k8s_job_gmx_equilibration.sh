#!/usr/bin/env bash
#
# k8s_job_gmx_equilibration.sh
#
# 役割:
#   - GROMACS による平衡化 (NVT / NPT) ステージ。
#   - 入力: 最小化後の構造 (em.gro 等)。
#   - 出力: /workspace/gmx_eq (eq.tpr, eq.trr, eq.edr, eq.gro 等)。
#   - CUDA GPU サポート付き GROMACS 前提。
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
GMX_EQ_DIR="${WORKSPACE}/gmx_eq"

GMX_NTASKS="${GMX_NTASKS:-4}"

mkdir -p "${GMX_EQ_DIR}"

echo "[INFO] gmx-equilibration: RUN_ID=${RUN_ID}"
echo "[INFO] GMX_EQ      : ${GMX_EQ_DIR}"
echo "[INFO] GMX_NTASKS  : ${GMX_NTASKS}"

if ! command -v gmx_mpi >/dev/null 2>&1 && ! command -v gmx >/dev/null 2>&1; then
  echo "[WARN] gmx_mpi / gmx not found; this is a skeleton script." >&2
  exit 0
fi

cd "${GMX_EQ_DIR}" || exit 1

cat <<CMD
# 例: NVT / NPT 平衡化 (CUDA GPU 利用)
# mpirun -np ${GMX_NTASKS} gmx_mpi mdrun \
#   -deffnm nvt \
#   -ntmpi ${GMX_NTASKS} \
#   -nb gpu -pme gpu -pmefft gpu -bonded gpu -update gpu
#
# mpirun -np ${GMX_NTASKS} gmx_mpi mdrun \
#   -deffnm npt \
#   -ntmpi ${GMX_NTASKS} \
#   -nb gpu -pme gpu -pmefft gpu -bonded gpu -update gpu
CMD

echo "[INFO] gmx-equilibration finished (skeleton)."
