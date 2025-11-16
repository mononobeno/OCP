#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
KIND_CLUSTER_NAME="${KIND_CLUSTER_NAME:-ocp-kind}"
NS="drug-pipeline"

echo "[INFO] building & loading vina-runner image..."
"$ROOT_DIR/scripts/build_and_load_vina_runner_image.sh"

echo "[INFO] applying namespace/pvc/config (idempotent)"
kubectl apply -f "$ROOT_DIR/k8s/base/namespace.yaml"
kubectl apply -f "$ROOT_DIR/k8s/base/storage-pvc.yaml"
kubectl apply -f "$ROOT_DIR/k8s/config/pipeline-config.yaml"

echo "[INFO] deleting old vina-runner job if exists..."
kubectl -n "$NS" delete job vina-runner --ignore-not-found

echo "[INFO] applying job-vina-runner.yaml..."
kubectl apply -f "$ROOT_DIR/k8s/jobs/job-vina-runner.yaml"

echo "[INFO] waiting for job completion..."
kubectl -n "$NS" wait --for=condition=complete job/vina-runner --timeout=600s

echo "[INFO] job completed. logs:"
kubectl -n "$NS" logs job/vina-runner
