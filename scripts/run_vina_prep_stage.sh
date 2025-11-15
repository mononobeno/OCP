#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

IMAGE_NAME="pipeline/vina-prep:local"
KIND_CLUSTER_NAME="${KIND_CLUSTER_NAME:-ocp-kind}"
NAMESPACE="drug-pipeline"

echo "[INFO] Preparing catalog data for kind cluster..."
# カタログデータをtarアーカイブにして転送
CONTROL_PLANE_CONTAINER="${KIND_CLUSTER_NAME}-control-plane"
echo "[INFO] Creating temporary tar archive..."
tar -C catalog -czf /tmp/ocp-catalog-data.tar.gz targets libraries db
echo "[INFO] Copying to kind node..."
docker cp /tmp/ocp-catalog-data.tar.gz "$CONTROL_PLANE_CONTAINER:/tmp/"
docker exec "$CONTROL_PLANE_CONTAINER" sh -c 'mkdir -p /tmp/ocp-catalog && cd /tmp/ocp-catalog && tar xzf /tmp/ocp-catalog-data.tar.gz && rm /tmp/ocp-catalog-data.tar.gz'
echo "[INFO] Verifying copied data..."
docker exec "$CONTROL_PLANE_CONTAINER" ls -la /tmp/ocp-catalog/
rm -f /tmp/ocp-catalog-data.tar.gz
echo "[INFO] Catalog data prepared"

echo "[INFO] building image: $IMAGE_NAME"
docker build -t "$IMAGE_NAME" images/vina-prep

echo "[INFO] loading image into kind (cluster: $KIND_CLUSTER_NAME)"
kind load docker-image "$IMAGE_NAME" --name "$KIND_CLUSTER_NAME"

echo "[INFO] applying namespace, pvc, config, job manifests"
kubectl apply -f k8s/base/namespace.yaml
kubectl apply -f k8s/base/storage-pvc.yaml
kubectl apply -f k8s/config/pipeline-config.yaml

echo "[INFO] loading catalog data into PVC..."
kubectl -n "$NAMESPACE" delete job catalog-loader --ignore-not-found
kubectl apply -f k8s/jobs/job-catalog-loader.yaml
kubectl -n "$NAMESPACE" wait --for=condition=complete job/catalog-loader --timeout=300s
echo "[INFO] catalog data loaded"

echo "[INFO] deleting old job if exists"
kubectl -n "$NAMESPACE" delete job vina-prep --ignore-not-found

echo "[INFO] applying vina-prep job"
kubectl apply -f k8s/jobs/job-vina-prep.yaml

echo "[INFO] waiting for job to complete..."
kubectl -n "$NAMESPACE" wait --for=condition=complete job/vina-prep --timeout=600s

echo "[INFO] job completed. logs:"
kubectl -n "$NAMESPACE" logs job/vina-prep

echo
echo "[INFO] workspace files (vina_input) on PVC should now contain:"
echo "  - vina_input/receptor.pdb"
echo "  - vina_input/receptor.pdbqt (変換に成功していれば)"
echo "  - vina_input/ligands/*.pdbqt (preprocessed が存在するもののみ)"
