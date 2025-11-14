# プロジェクト構成 (OCP)

## images/

コンテナイメージ用の Dockerfile を配置するディレクトリ。

- images/gromacs-mpi-cuda/  
  - GROMACS (CUDA + OpenMPI 対応) コンテナ用 Dockerfile
- images/vina/  
  - AutoDock Vina 用コンテナイメージ
- images/common/  
  - 構造変換やレポート生成など、共通ツールをまとめたイメージ

## k8s/

Kubernetes マニフェストを配置するディレクトリ。

- k8s/base/  
  - namespace, PV/PVC など環境共通のベース定義
- k8s/jobs/  
  - ligand-selector, vina-runner, gmx-minimize など各 Job 定義
- k8s/config/  
  - ConfigMap, Secret など設定系リソース

## scripts/

クラスタ作成や一括実行用スクリプトを配置。

- 例: scripts/create_kind_cluster.sh
- 例: scripts/run_gpu_test_job.sh

## docs/

ドキュメントや GitHub Pages 連携を想定したファイル群。

- PROJECT_STRUCTURE.md (このファイル)
- 今後、解析レポートや設計書を配置していく。
