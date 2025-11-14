#!/usr/bin/env bash
set -euo pipefail

CLUSTER_NAME="ocp-kind"

echo "==============================="
echo "[STEP 1] kind クラスタ確認/作成"
echo "==============================="

if kind get clusters 2>/dev/null | grep -qx "${CLUSTER_NAME}"; then
  echo "[INFO] kind cluster '${CLUSTER_NAME}' は既に存在します。スキップ。"
else
  echo "[INFO] kind cluster '${CLUSTER_NAME}' が無いので作成します。"
  if [ -x scripts/create_kind_cluster.sh ]; then
    scripts/create_kind_cluster.sh
  else
    echo "[ERROR] scripts/create_kind_cluster.sh が見つからないか実行権限がありません。" >&2
    exit 1
  fi
fi

echo
echo "==============================="
echo "[STEP 2] GROMACS CUDA+MPI イメージ build & kind へロード"
echo "==============================="

if [ -x scripts/build_and_load_gromacs_image.sh ]; then
  scripts/build_and_load_gromacs_image.sh
else
  echo "[ERROR] scripts/build_and_load_gromacs_image.sh が見つからないか実行権限がありません。" >&2
  exit 1
fi

echo
echo "==============================="
echo "[STEP 3] GROMACS GPU テスト Job を apply"
echo "==============================="

if [ -x scripts/apply_gromacs_gpu_test_job.sh ]; then
  scripts/apply_gromacs_gpu_test_job.sh
else
  echo "[ERROR] scripts/apply_gromacs_gpu_test_job.sh が見つからないか実行権限がありません。" >&2
  exit 1
fi

echo
echo "==============================="
echo "[STEP 4] Job の完了待ち & ログ表示"
echo "==============================="

JOB_NAME="gromacs-gpu-test"
NS="drug-pipeline"

echo "[INFO] Job/${JOB_NAME} が完了するまで待機します..."
kubectl -n "${NS}" wait --for=condition=complete --timeout=600s job/${JOB_NAME} || {
  echo "[WARN] Job が完了しませんでした。Pod の状態を確認します。"
  kubectl -n "${NS}" get pods
  exit 1
}

POD_NAME=$(kubectl -n "${NS}" get pods -l job-name=${JOB_NAME} -o jsonpath='{.items[0].metadata.name}')

echo
echo "[INFO] Pod ログ (nvidia-smi / gmx_mpi --version など) を表示します:"
echo "-------------------------------------------------------"
kubectl -n "${NS}" logs "${POD_NAME}"
echo "-------------------------------------------------------"

echo
echo "[DONE] GROMACS GPU テスト Job 実行完了。"
