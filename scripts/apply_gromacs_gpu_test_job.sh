#!/usr/bin/env bash
set -euo pipefail

echo "[INFO] Applying namespace (drug-pipeline)"
kubectl apply -f k8s/base/namespace.yaml

echo "[INFO] Deleting existing gromacs-gpu-test Job if present"
kubectl -n drug-pipeline delete job gromacs-gpu-test --ignore-not-found=true

echo "[INFO] Applying GROMACS GPU test Job"
kubectl apply -f k8s/jobs/job-gromacs-gpu-test.yaml

echo "[INFO] Pods in drug-pipeline:"
kubectl -n drug-pipeline get pods
