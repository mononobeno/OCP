#!/usr/bin/env bash
set -euo pipefail

# カタログデータをPVCに直接コピーするスクリプト

NAMESPACE="drug-pipeline"
PVC_NAME="zinc-library-pvc"

echo "[INFO] Copying catalog data directly to PVC..."

# 一時的なPodを作成してデータをコピー
kubectl -n "$NAMESPACE" run catalog-copier --image=busybox --restart=Never \
  --overrides='{
    "spec": {
      "containers": [{
        "name": "catalog-copier",
        "image": "busybox",
        "command": ["sleep", "3600"],
        "volumeMounts": [{
          "name": "zinc-library",
          "mountPath": "/zinc-library"
        }]
      }],
      "volumes": [{
        "name": "zinc-library",
        "persistentVolumeClaim": {
          "claimName": "'"$PVC_NAME"'"
        }
      }]
    }
  }'

# Podが起動するまで待つ
kubectl -n "$NAMESPACE" wait --for=condition=Ready pod/catalog-copier --timeout=60s

# データをコピー
echo "[INFO] Copying targets..."
kubectl -n "$NAMESPACE" exec catalog-copier -- mkdir -p /zinc-library/targets
for target_dir in catalog/targets/*; do
  target_name=$(basename "$target_dir")
  echo "[INFO]   - $target_name"
  kubectl -n "$NAMESPACE" exec catalog-copier -- mkdir -p "/zinc-library/targets/$target_name"
  for file in "$target_dir"/*; do
    kubectl cp "$file" "$NAMESPACE/catalog-copier:/zinc-library/targets/$target_name/$(basename "$file")"
  done
done

echo "[INFO] Copying libraries..."
kubectl -n "$NAMESPACE" exec catalog-copier -- mkdir -p /zinc-library/libraries
for lib_dir in catalog/libraries/*; do
  lib_name=$(basename "$lib_dir")
  echo "[INFO]   - $lib_name"
  kubectl -n "$NAMESPACE" exec catalog-copier -- mkdir -p "/zinc-library/libraries/$lib_name"
  for item in "$lib_dir"/*; do
    if [ -f "$item" ]; then
      kubectl cp "$item" "$NAMESPACE/catalog-copier:/zinc-library/libraries/$lib_name/$(basename "$item")"
    elif [ -d "$item" ]; then
      subdir_name=$(basename "$item")
      kubectl -n "$NAMESPACE" exec catalog-copier -- mkdir -p "/zinc-library/libraries/$lib_name/$subdir_name"
      # サブディレクトリにファイルがある場合のみコピー
      if [ "$(ls -A "$item" 2>/dev/null)" ]; then
        for file in "$item"/*; do
          if [ -f "$file" ]; then
            kubectl cp "$file" "$NAMESPACE/catalog-copier:/zinc-library/libraries/$lib_name/$subdir_name/$(basename "$file")"
          fi
        done
      fi
    fi
  done
done

# マーカーファイルを作成
kubectl -n "$NAMESPACE" exec catalog-copier -- sh -c 'date > /zinc-library/.catalog_loaded'

echo "[INFO] Verifying..."
kubectl -n "$NAMESPACE" exec catalog-copier -- ls -laR /zinc-library/ | head -50

# Podを削除
kubectl -n "$NAMESPACE" delete pod catalog-copier

echo "[INFO] Catalog data copied successfully"
