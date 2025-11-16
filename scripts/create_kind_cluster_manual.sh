#!/usr/bin/env bash
set -euo pipefail

# kindの起動検出問題を回避して手動でクラスタを作成

CLUSTER_NAME="ocp-kind"
API_SERVER_PORT="45427"

echo "[INFO] Creating kind cluster with manual initialization..."

# Step 1: コンテナを起動（バックグラウンドで）
echo "[STEP 1] Starting kind containers..."
kind create cluster --name $CLUSTER_NAME --config kind-config.yaml 2>&1 | tee /tmp/kind-create.log &
KIND_PID=$!

# コンテナが作成されるまで待つ（最大90秒）
echo "[INFO] Waiting for containers to be created..."
for i in {1..90}; do
  if docker ps --filter "name=${CLUSTER_NAME}-control-plane" --format "{{.Names}}" | grep -q control-plane; then
    echo "[INFO] Containers created after $i seconds"
    break
  fi
  sleep 1
done

# kindプロセスがまだ動いていたらkill
if ps -p $KIND_PID > /dev/null 2>&1; then
  echo '[INFO] Killing hung kind process...'
  kill $KIND_PID 2>/dev/null || true
  wait $KIND_PID 2>/dev/null || true
fi

# Step 2: コンテナが起動しているか確認
echo "[STEP 2] Checking if containers are running..."
sleep 5

CONTROL_PLANE=$(docker ps --filter "name=${CLUSTER_NAME}-control-plane" --format "{{.Names}}" | head -1)
WORKER=$(docker ps --filter "name=${CLUSTER_NAME}-worker" --format "{{.Names}}" | head -1)

if [ -z "$CONTROL_PLANE" ]; then
  echo "[ERROR] Control plane container not found!"
  exit 1
fi

echo "[INFO] Found containers:"
echo "  - $CONTROL_PLANE"
if [ -n "$WORKER" ]; then
  echo "  - $WORKER"
fi

# Step 3: systemdが起動するまで待つ
echo "[STEP 3] Waiting for systemd to be ready..."
for i in {1..30}; do
  STATUS=$(docker exec $CONTROL_PLANE systemctl is-system-running 2>/dev/null || echo "not-ready")
  if [ "$STATUS" = "running" ] || [ "$STATUS" = "degraded" ]; then
    echo "[INFO] systemd is $STATUS"
    break
  fi
  echo "  Waiting... ($i/30)"
  sleep 2
done

# Step 4: Kubernetesコンポーネントが起動するまで待つ
echo "[STEP 4] Waiting for Kubernetes components..."
for i in {1..60}; do
  if docker exec $CONTROL_PLANE test -f /etc/kubernetes/admin.conf 2>/dev/null; then
    echo "[INFO] Kubernetes initialized!"
    break
  fi
  
  # kubeletとcontainerdのプロセス確認
  if [ $((i % 10)) -eq 0 ]; then
    echo "  Checking processes ($i/60)..."
    docker exec $CONTROL_PLANE ps aux | grep -E "kubelet|containerd|kube-apiserver" | grep -v grep || true
  fi
  sleep 3
done

# Step 5: kubeconfigを取得
echo "[STEP 5] Extracting kubeconfig..."
if docker exec $CONTROL_PLANE test -f /etc/kubernetes/admin.conf; then
  # kubeconfigを取得して設定
  docker exec $CONTROL_PLANE cat /etc/kubernetes/admin.conf > /tmp/kind-kubeconfig-$CLUSTER_NAME
  
  # ポート番号を修正
  ACTUAL_PORT=$(docker port $CONTROL_PLANE 6443 | cut -d: -f2)
  sed -i "s/6443/$ACTUAL_PORT/g" /tmp/kind-kubeconfig-$CLUSTER_NAME
  
  # kubeconfigをマージ
  export KUBECONFIG=/tmp/kind-kubeconfig-$CLUSTER_NAME:$HOME/.kube/config
  kubectl config view --flatten > /tmp/merged-config
  mv /tmp/merged-config $HOME/.kube/config
  kubectl config use-context kind-$CLUSTER_NAME
  
  echo "[INFO] Kubeconfig configured"
else
  echo "[WARN] Kubernetes not fully initialized, trying to trigger it..."
  
  # kubeadmで強制初期化を試みる
  docker exec $CONTROL_PLANE bash -c '
    if [ ! -f /etc/kubernetes/admin.conf ]; then
      echo "Running kubeadm init..."
      # kind用のデフォルト設定でkubeadm init
      cat > /tmp/kubeadm-config.yaml <<EOF
apiVersion: kubeadm.k8s.io/v1beta3
kind: InitConfiguration
localAPIEndpoint:
  advertiseAddress: 0.0.0.0
  bindPort: 6443
nodeRegistration:
  criSocket: unix:///run/containerd/containerd.sock
  kubeletExtraArgs:
    node-ip: $(hostname -i)
---
apiVersion: kubeadm.k8s.io/v1beta3
kind: ClusterConfiguration
networking:
  podSubnet: 10.244.0.0/16
  serviceSubnet: 10.96.0.0/12
EOF
      kubeadm init --config=/tmp/kubeadm-config.yaml --ignore-preflight-errors=all
    fi
  '
  
  sleep 10
  
  # 再試行
  if docker exec $CONTROL_PLANE test -f /etc/kubernetes/admin.conf; then
    docker exec $CONTROL_PLANE cat /etc/kubernetes/admin.conf > /tmp/kind-kubeconfig-$CLUSTER_NAME
    ACTUAL_PORT=$(docker port $CONTROL_PLANE 6443 | cut -d: -f2)
    sed -i "s/6443/$ACTUAL_PORT/g" /tmp/kind-kubeconfig-$CLUSTER_NAME
    export KUBECONFIG=/tmp/kind-kubeconfig-$CLUSTER_NAME:$HOME/.kube/config
    kubectl config view --flatten > /tmp/merged-config
    mv /tmp/merged-config $HOME/.kube/config
    kubectl config use-context kind-$CLUSTER_NAME
  else
    echo "[ERROR] Failed to initialize Kubernetes"
    exit 1
  fi
fi

# Step 6: ノードの確認
echo "[STEP 6] Verifying cluster..."
kubectl get nodes
kubectl get pods -A

echo ""
echo "[SUCCESS] Kind cluster '$CLUSTER_NAME' is ready!"
echo "  Context: kind-$CLUSTER_NAME"
echo ""
