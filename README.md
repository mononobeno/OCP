# OCP - k8s + MD パイプライン環境

## 目的

WSL2 上の Kubernetes クラスタで、AutoDock Vina と GROMACS（CUDA + OpenMPI 対応）を使った  
**「ドッキング → MD シミュレーション → GitHub Pages レポート更新」** までを一気通貫で実行できる環境を作る。

---

## 全体アーキテクチャ概要

### ホスト環境

- Windows 11
  - WSL2: Ubuntu 20.04
  - GPU: NVIDIA（WSL から `nvidia-smi` が見える状態）
- コンテナランタイム
  - Docker + NVIDIA Container Toolkit（予定）
- Kubernetes
  - Kind もしくは kubeadm で構築したクラスタ
  - GPU 対応 worker ノードに `gpu=true` ラベルを付与し、GROMACS 用 Pod をスケジューリング

### Kubernetes リソース（論理構成）

- Namespace: `drug-pipeline`（仮）
- ストレージ:
  - `zinc-library-pv` / `zinc-library-pvc`  
    - ZINC 化合物ライブラリ用（read-only）
  - `work-pv` / `work-pvc`  
    - 中間ファイル・ログ・最終結果・レポート用
- 認証情報:
  - `Secret: github-credentials`（PAT or SSH key）
- 設定:
  - `ConfigMap: pipeline-config`（抽出件数、パス、GROMACS パラメータなど）

---

## パイプライン構成（Pod / Job 単位）

ZINC ライブラリから化合物を選び、AutoDock Vina → GROMACS → GitHub Pages 更新までを  
「ステージごとの Pod / Job」として分割する。

1. **ligand-selector Job**
   - ZINC 化合物ライブラリから任意件数の化合物を抽出して `/workspace/ligands_raw` にコピー
   - 入力: `zinc-library-pvc`（read-only）
   - 出力: `work-pvc`

2. **vina-prep Job**
   - Vina 計算前処理
     - Receptor PDB → PDBQT 変換
     - Ligand PDB/SDF → PDBQT 変換
   - 出力: `/workspace/vina_input`（receptor + ligands）

3. **vina-runner Job**
   - AutoDock Vina によるドッキング計算
   - 入力: `/workspace/vina_input`
   - 出力: `/workspace/vina_output`
     - 各 ligand の out.pdbqt
     - `summary.csv` など

4. **vina2gromacs-prep Job**
   - Vina 結果を GROMACS 用に整形
     - 最上位 pose を抽出
     - receptor + ligand の複合体 PDB 生成
     - ligand の topology 生成（acpype 等）
   - 出力: `/workspace/gmx_input`（`.gro`, `.top`, `.itp` など）

5. **gmx-makebox Job**（GROMACS, CUDA+MPI イメージ）
   - シミュレーションボックス作成
     - `gmx editconf`（box 定義）
     - `gmx solvate`（溶媒追加）
     - `gmx grompp`（最小化用 tpr）
   - 出力: `/workspace/gmx_em`（em.tpr 等）
   - GPU 必須ではないが、同一イメージを使う想定

6. **gmx-minimize Job**（GROMACS, CUDA + OpenMPI）
   - エネルギー最小化
   - 例: `mpirun -np 4 gmx_mpi mdrun -deffnm em`
   - 入力: `/workspace/gmx_em`
   - 出力: `em.gro`, `em.edr`, `em.log` など

7. **min-analysis-pages Job**
   - 最小化結果の評価 + GitHub Pages 更新
   - 解析:
     - エネルギー収束
     - RMSD など
   - GitHub Pages:
     - `docs/minimization/` 配下にレポート（Markdown/HTML/PNG）を生成
     - Git commit & push

8. **gmx-equilibration Job**（GROMACS, CUDA + OpenMPI）
   - 平衡化 (NVT / NPT)
   - 入力: 最小化最終構造
   - 出力: `/workspace/gmx_eq`（eq.tpr, eq.trr, eq.edr, eq.gro 等）

9. **eq-analysis-pages Job**
   - 平衡化結果の評価 + Pages 更新
   - 例: 温度・圧力・エネルギーの安定性の可視化

10. **gmx-production Job**（GROMACS, CUDA + OpenMPI）
    - 本番 Production run（ns オーダー）
    - 出力: `/workspace/gmx_md`（md.xtc, md.edr, md.log など）

11. **prod-analysis-pages Job**
    - Production run の結果を集計・評価
    - 最終的な結合安定性評価やランキングを GitHub Pages に反映

---

## ワークフロー（論理フロー図）

```mermaid
graph TD
    A[GitHub Actions<br/>select ligands] --> B[ligand-selector Job<br/>ZINC → ligands_raw]
    B --> C[vina-prep Job<br/>構造変換(PDB/SDF → PDBQT)]
    C --> D[vina-runner Job<br/>AutoDock Vina]
    D --> E[vina2gromacs-prep Job<br/>Vina → GROMACS入力]
    E --> F[gmx-makebox Job<br/>box & solvate]
    F --> G[gmx-minimize Job<br/>GPU+MPI]
    G --> H[min-analysis-pages Job<br/>Pages更新]
    H --> I[gmx-equilibration Job<br/>GPU+MPI]
    I --> J[eq-analysis-pages Job<br/>Pages更新]
    J --> K[gmx-production Job<br/>GPU+MPI]
    K --> L[prod-analysis-pages Job<br/>Pages更新]

