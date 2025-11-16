# ステージ Job / Pod 設計メモ（DB 駆動 + GPU GROMACS 方針）

日付: 2025-11-16  
担当: めけめて もすもす / ChatGPT 補助

## 基本方針

- すべてのステージ (ligand-selector / vina-prep / vina2gromacs-prep / gmx-* / *-analysis-pages) は、
  **SQLite DB (catalog/db/ocp_results.sqlite)** を source of truth として扱う。
  - 計算対象 (RUN / ligand / target) は DB から取得する。
  - ファイルは再生成可能なキャッシュとして扱う。

- GROMACS を使用するステージ (gmx-*) は、
  **必ず CUDA サポート付きの GPU 計算を前提としたイメージ** を使う。
  - スクリプトでは `gmx_mpi` / `mdrun` の GPU オプションを雛形として記述。

- 生の化合物ライブラリは、
  **RDKit による 3D 化・整形を通した上で DB に格納しておく**。
  - `ligands.has_3d = 1` と `ligands.conformer_method = 'rdkit_etkdg'` のようなフラグで管理。
  - ligand-selector はこのフラグを見て、3D 化済みのものだけを候補にする。

## このセットアップで作成されたスクリプト

- `scripts/k8s_job_ligand_selector.sh`
- `scripts/k8s_job_vina_prep.sh`
- `scripts/k8s_job_vina2gromacs_prep.sh`
- `scripts/k8s_job_gmx_makebox.sh`
- `scripts/k8s_job_gmx_minimize.sh`
- `scripts/k8s_job_min_analysis_pages.sh`
- `scripts/k8s_job_gmx_equilibration.sh`
- `scripts/k8s_job_eq_analysis_pages.sh`
- `scripts/k8s_job_gmx_production.sh`
- `scripts/k8s_job_prod_analysis_pages.sh`

これらはすべて **RUN_ID を引数または環境変数で受け取り、DB を参照する骨格** を持つ。
GROMACS 部分は CUDA GPU 前提のコマンドをコメントとして記述しており、
今後の実装で具体的な mdp / tpr ファイル名などを埋めていくことを想定している。

