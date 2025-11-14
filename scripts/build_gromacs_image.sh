#!/usr/bin/env bash
set -euo pipefail

IMAGE_NAME="ocp/gromacs-mpi-cuda:local"

echo "[INFO] Building GROMACS CUDA+MPI image: ${IMAGE_NAME}"
docker build -t "${IMAGE_NAME}" images/gromacs-mpi-cuda

echo "[INFO] Build finished."
echo "[INFO] Test run: gmx_mpi --version"
docker run --rm --gpus all "${IMAGE_NAME}"
