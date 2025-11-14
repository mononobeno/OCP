#!/usr/bin/env bash
set -euo pipefail

NS="drug-pipeline"

echo "[INFO] Namespace を作成/更新"
kubectl apply -f k8s/base/namespace.yaml

echo "[INFO] ConfigMap / PVC を apply"
kubectl apply -f k8s/config/pipeline-config.yaml
kubectl apply -f k8s/base/storage-pvc.yaml

echo "[INFO] パイプライン各 Job 定義を apply (まだ中身はスケルトン)"
kubectl apply -f k8s/jobs/job-ligand-selector.yaml
kubectl apply -f k8s/jobs/job-vina-prep.yaml
kubectl apply -f k8s/jobs/job-vina-runner.yaml
kubectl apply -f k8s/jobs/job-vina2gromacs-prep.yaml
kubectl apply -f k8s/jobs/job-gmx-makebox.yaml
kubectl apply -f k8s/jobs/job-gmx-minimize.yaml
kubectl apply -f k8s/jobs/job-min-analysis-pages.yaml
kubectl apply -f k8s/jobs/job-gmx-equilibration.yaml
kubectl apply -f k8s/jobs/job-eq-analysis-pages.yaml
kubectl apply -f k8s/jobs/job-gmx-production.yaml
kubectl apply -f k8s/jobs/job-prod-analysis-pages.yaml

echo "[INFO] drug-pipeline namespace の Job 一覧:"
kubectl -n "${NS}" get jobs
