#!/usr/bin/env bash
set -euo pipefail

echo "[INFO] Applying namespace (drug-pipeline)"
kubectl apply -f k8s/base/namespace.yaml

echo "[INFO] Applying GROMACS version test Job"
kubectl apply -f k8s/jobs/job-gmx-version-test.yaml

echo "[INFO] Pods in drug-pipeline:"
kubectl -n drug-pipeline get pods
