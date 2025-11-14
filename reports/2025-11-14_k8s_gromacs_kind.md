# Kubernetes + GROMACS(kind) 動作確認レポート (2025-11-14)

本レポートは、kind ベースのローカル Kubernetes クラスタ上で GROMACS (CUDA+MPI) コンテナを実行し、
動作確認を行った結果をまとめたものである。

---

## 1. kind クラスタ構成

作成スクリプト: scripts/create_kind_cluster.sh

構成:
- クラスタ名: ocp-kind
- ノード:
  - control-plane: 1 ノード
  - worker: 1 ノード
- Kubernetes バージョン: v1.30.0
- コンテナランタイム: containerd (kind 標準)

確認コマンド:
- kind get clusters
- kubectl get nodes -o wide

結果:
- ocp-kind-control-plane: Ready
- ocp-kind-worker: NotReady → その後 Ready に遷移

---

## 2. GROMACS (CUDA+MPI) イメージ

Dockerfile: images/gromacs-mpi-cuda/Dockerfile

主なポイント:
- ベースイメージ:
  - nvidia/cuda:12.2.0-devel-ubuntu22.04
- ビルドオプション:
  - GMX_MPI=ON (MPI 対応)
  - GMX_GPU=CUDA (CUDA 対応)
  - FFTW3 を使用 (FFTW3 ライブラリ)
- インストール先:
  - /usr/local/gromacs
- デフォルト CMD:
  - nvidia-smi (あれば)
  - gmx_mpi --version もしくは gmx --version

ビルド:
- docker build -t gromacs-mpi-cuda:local images/gromacs-mpi-cuda

kind へのロード:
- kind load docker-image gromacs-mpi-cuda:local --name ocp-kind

---

## 3. GROMACS テスト Job (kind 上)

Job 定義: k8s/jobs/job-gromacs-gpu-test.yaml  
実行スクリプト: scripts/run_gromacs_gpu_test.sh

実行フロー:
1. kind クラスタが無ければ作成
2. GROMACS イメージをビルドして kind にロード
3. drug-pipeline namespace を apply
4. gromacs-gpu-test Job を apply
5. Job 完了を待機し、Pod ログを取得

Pod 内で実行される内容（要約）:
- nvidia-smi (GPU ランタイムが無い場合は not found)
- gmx_mpi --version

---

## 4. kind クラスタにおける GPU 利用の注意点

今回の調査で分かったこと:

- kind はローカルの containerd ベースクラスタであり、
  デフォルトではホストの NVIDIA GPU を Kubernetes リソースとして認識しない。
- 実際に GPU を K8s リソースとして使うには:
  - NVIDIA device plugin の導入
  - ノード側で GPU デバイスをコンテナに渡せる設定
  が必要になる。
- ローカル kind 環境では GPU ハードウェアアクセスに制限があるため、
  今回は **GPU リソース要求 (resources.limits.nvidia.com/gpu)** を Job から削除し、
  「GPU 付きでビルドされた GROMACS バイナリが正しく動くか」を中心に検証した。

その結果:

- Pod 内で nvidia-smi は `not found`（GPU ランタイム非有効のため）となるが、
- gmx_mpi --version の結果から:
  - GROMACS 2023.3
  - MPI 対応
  - CUDA サポート
  - AVX2 SIMD
  などが正しく有効化されていることを確認。

---

## 5. ログ要約 (gromacs-gpu-test Job)

Job 実行結果 (要約):

- nvidia-smi:
  - not found (GPU ランタイム非有効なコンテナを想定したメッセージ)
- gmx_mpi --version:
  - GROMACS 2023.3
  - Precision: mixed
  - MPI: enabled
  - OpenMP: enabled
  - GPU support: CUDA
  - SIMD: AVX2_256
  - CUDA runtime: 12.2
  などが表示され、ビルド自体は意図通り成功している。

---

## 6. 結論と今後の方針

- kind 上の Kubernetes 環境で:
  - GROMACS (CUDA+MPI) コンテナの起動
  - gmx_mpi バージョン情報の取得
  までを確認できた。
- ローカル kind 環境では GPU リソースを直接使った大規模計算は難しいため、
  今後の方針としては:

1. **ワークフロー設計用の K8s 環境として kind を利用**
   - Vina → GROMACS → 解析 → Pages 更新の「ジョブ連鎖」の動作確認
   - CPU 実行や短いテストランに利用

2. **本番レベルの GPU 計算は WSL+Docker 直、または別途 GPU 対応 K8s クラスタで実施**
   - 既に確認済みの `docker run --gpus all ...` を活用
   - 将来的に kubeadm 等で GPU ノードを持つクラスタを構築する選択肢あり

以上より、現時点では:
- 「GROMACS CUDA+MPI コンテナのビルドと kind 上での動作検証」は完了
- 次のステップとして、「実際の計算 (gmx mdrun) を行う Job」と、
  Vina とのパイプライン統合に進める状態となった。
