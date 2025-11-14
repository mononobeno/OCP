#!/usr/bin/env bash
set -euo pipefail

echo "[INFO] Apply namespace and GPU test Job"

kubectl apply -f k8s/base/namespace.yaml
kubectl apply -f k8s/jobs/job-gpu-nvidia-smi.yaml

echo "[INFO] Pods in drug-pipeline namespace:"
kubectl -n drug-pipeline get pods
