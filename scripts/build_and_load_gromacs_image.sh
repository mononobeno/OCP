#!/usr/bin/env bash
set -euo pipefail

IMAGE_NAME="gromacs-mpi-cuda:local"
CLUSTER_NAME="ocp-kind"

echo "[INFO] Build GROMACS (CUDA+MPI) image: ${IMAGE_NAME}"
docker build -t "${IMAGE_NAME}" images/gromacs-mpi-cuda

echo "[INFO] Load image into kind cluster: ${CLUSTER_NAME}"
kind load docker-image "${IMAGE_NAME}" --name "${CLUSTER_NAME}"

echo "[INFO] Done. You can now run the gromacs-gpu-test Job."
