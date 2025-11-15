#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

IMAGE_NAME="pipeline/ligand-selector:local"
KIND_CLUSTER_NAME="${KIND_CLUSTER_NAME:-ocp-kind}"
NAMESPACE="drug-pipeline"

echo "[INFO] building image: $IMAGE_NAME"
docker build -t "$IMAGE_NAME" images/ligand-selector

echo "[INFO] loading image into kind (cluster: $KIND_CLUSTER_NAME)"
kind load docker-image "$IMAGE_NAME" --name "$KIND_CLUSTER_NAME"

echo "[INFO] applying namespace, pvc, config, job manifests"
kubectl apply -f k8s/base/namespace.yaml
kubectl apply -f k8s/base/storage-pvc.yaml
kubectl apply -f k8s/config/pipeline-config.yaml

# 既存 Job があれば消す
echo "[INFO] deleting old job if exists"
kubectl -n "$NAMESPACE" delete job ligand-selector --ignore-not-found

echo "[INFO] applying ligand-selector job"
kubectl apply -f k8s/jobs/job-ligand-selector.yaml

echo "[INFO] waiting for job to complete..."
kubectl -n "$NAMESPACE" wait --for=condition=complete job/ligand-selector --timeout=600s

echo "[INFO] job completed. logs:"
kubectl -n "$NAMESPACE" logs job/ligand-selector

echo
echo "[INFO] workspace files (ligands_raw) on PVC will contain:"
echo "  - ligands_selected.tsv"
echo "  - ligands_zinc_ids.txt"
echo "  - run_info.json"
