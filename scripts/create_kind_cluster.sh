#!/usr/bin/env bash
set -euo pipefail

CLUSTER_NAME="ocp-kind"
KIND_CONFIG="kind-config.yaml"

echo "[INFO] Creating kind config: ${KIND_CONFIG}"

cat > "${KIND_CONFIG}" << 'CFG'
kind: Cluster
apiVersion: kind.x-k8s.io/v1alpha4
nodes:
  - role: control-plane
    extraPortMappings:
      - containerPort: 30080
        hostPort: 30080
      - containerPort: 30443
        hostPort: 30443
  - role: worker
CFG

echo "[INFO] Creating kind cluster: ${CLUSTER_NAME}"
kind create cluster --name "${CLUSTER_NAME}" --config "${KIND_CONFIG}"

echo "[INFO] Current nodes:"
kubectl get nodes -o wide
