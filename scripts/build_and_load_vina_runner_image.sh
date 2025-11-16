#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IMG_DIR="$ROOT_DIR/images/vina-runner"

IMAGE_NAME="pipeline/vina-runner:local"
KIND_CLUSTER_NAME="${KIND_CLUSTER_NAME:-ocp-kind}"

echo "[INFO] building image: $IMAGE_NAME"
cd "$IMG_DIR"
docker build -t "$IMAGE_NAME" .

echo "[INFO] loading image into kind cluster: $KIND_CLUSTER_NAME"
kind load docker-image "$IMAGE_NAME" --name "$KIND_CLUSTER_NAME"

echo "[INFO] build_and_load_vina_runner_image.sh done."
