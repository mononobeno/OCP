# GitHub Pages階層構造実装レポート

日付: 2025-11-17  
実装者: めけめて もすもす

## 概要

ターゲットタンパク質ベースの階層的ページ構成と、Plotly.jsによる対話的グラフ可視化を実装した。

---

## 新しいページ構成

### 階層構造

```
トップページ (index.md)
├── Kubernetes Job Status
├── Target Proteins一覧
│
└─→ ターゲット詳細ページ (targets/target_X.md)
    ├── タンパク質情報 (PDB ID, 3D構造リンク)
    ├── 化合物結果テーブル
    │
    └─→ 化合物詳細ページ (compounds/compound_X.md)
        ├── 化合物情報 (SMILES, Vina score)
        ├── MD結果サマリー
        ├── 平衡化評価グラフ (EM/NVT/NPT)
        └── Production MD解析グラフ (RMSD, H-bond)
```

### URL構造

- **トップ**: `/pages/index.html`
- **ターゲット**: `/pages/targets/target_1.html`
- **化合物詳細**: `/pages/compounds/compound_1.html`

---

## 実装した可視化グラフ

### 1. 平衡化品質評価

#### Energy Minimization (EM)
- **グラフ**: Potential Energy vs Step
- **評価基準**: エネルギーが収束しているか
- **データソース**: `energy_em_potential.xvg`
- **プロット色**: 赤

#### NVT Equilibration
- **グラフ**: Temperature vs Time
- **評価基準**: 300K付近で安定しているか
- **データソース**: `energy_nvt_temp.xvg`
- **プロット色**: 緑
- **特徴**: 300K基準線を赤破線で表示

#### NPT Equilibration
- **グラフ1**: Pressure vs Time
- **評価基準**: 1 bar付近で変動しているか
- **データソース**: `energy_npt_pressure.xvg`
- **プロット色**: 青

- **グラフ2**: Density vs Time
- **評価基準**: 密度が安定しているか
- **データソース**: `energy_npt_density.xvg`
- **プロット色**: 紫

### 2. Production MD解析

#### RMSD (Backbone)
- **グラフ**: RMSD vs Time
- **評価**: タンパク質構造の安定性
- **データソース**: `rmsd_backbone.xvg`
- **プロット色**: デフォルト青

#### Hydrogen Bonds
- **グラフ**: H-bond count vs Time
- **評価**: タンパク質内水素結合の変動
- **データソース**: `hbond.xvg`
- **プロット色**: オレンジ

---

## データベース拡張

### 新規テーブル

#### md_systems
システム構成と計算条件を保存

```sql
CREATE TABLE md_systems (
    id INTEGER PRIMARY KEY,
    vina_result_id INTEGER,
    total_atoms INTEGER,           -- 総原子数
    protein_atoms INTEGER,         -- タンパク質原子数
    ligand_atoms INTEGER,          -- リガンド原子数
    water_molecules INTEGER,       -- 水分子数
    ions_count INTEGER,            -- イオン数
    box_x/y/z REAL,               -- Box寸法
    box_volume REAL,              -- Box体積
    force_field TEXT,             -- 力場 (OPLSAA等)
    water_model TEXT,             -- 水モデル (SPC/E等)
    temperature REAL,             -- 温度 (K)
    pressure REAL,                -- 圧力 (bar)
    timestep_ps REAL,             -- タイムステップ (ps)
    performance_nsday REAL,       -- 性能 (ns/day)
    gpu_model TEXT                -- GPU型番
);
```

#### md_analysis
解析統計値を保存

```sql
CREATE TABLE md_analysis (
    id INTEGER PRIMARY KEY,
    vina_result_id INTEGER,
    rmsd_backbone_avg REAL,       -- RMSD平均
    rmsd_backbone_std REAL,       -- RMSD標準偏差
    rmsd_ligand_avg REAL,         -- Ligand RMSD平均
    hbond_avg REAL,               -- H-bond平均
    hbond_std REAL,               -- H-bond標準偏差
    hbond_max INTEGER,            -- H-bond最大値
    rg_avg REAL,                  -- 回転半径平均
    rmsf_avg REAL,                -- RMSF平均
    potential_energy_avg REAL,    -- ポテンシャルエネルギー平均
    analysis_json TEXT            -- 解析JSONパス
);
```

### 登録データ例

**md_systems (vina_result_id=1)**:
- total_atoms: 38,392
- protein_atoms: 1,960
- water_molecules: 12,144
- force_field: OPLSAA
- water_model: SPC/E
- temperature: 300.0 K
- pressure: 1.0 bar
- timestep_ps: 0.002
- performance_nsday: 408.7
- gpu_model: NVIDIA RTX 4070

**md_analysis (vina_result_id=1)**:
- rmsd_backbone_avg: 0.0702 nm
- hbond_avg: 89.5
- analysis_json: `/home/dev/OCP/results/md_gpu_test_1ns/analysis.json`

---

## 実装スクリプト

### 1. generate_hierarchical_pages_v2.sh
階層的ページ生成メインスクリプト

**機能**:
- トップページ生成 (K8s jobs + ターゲット一覧)
- ターゲット詳細ページ生成 (PDB情報 + 化合物テーブル)
- 化合物詳細ページ生成 (Plotlyグラフ埋め込み)

**使用法**:
```bash
bash scripts/generate_hierarchical_pages_v2.sh
```

### 2. generate_pages_with_equilibration.sh
平衡化グラフ追加版

**機能**:
- EM/NVT/NPT平衡化グラフ生成
- Production MD グラフ生成
- XVGデータをJavaScript配列に変換
- Plotly.js CDN読み込み

**使用法**:
```bash
bash scripts/generate_pages_with_equilibration.sh
```

### 3. analyze_md_trajectory.sh
MD軌跡解析実行

**機能**:
- RMSD計算 (backbone)
- 水素結合解析
- Radius of gyration計算
- RMSF計算
- エネルギー抽出 (EM/NVT/NPT)

**出力**: `analysis.json`

### 4. db_register_md_system.sh
システム情報DB登録

**機能**:
- TPRファイルから情報抽出
- md_systemsテーブルに登録
- md_analysisテーブルに統計値登録

---

## 技術詳細

### Plotly.js統合

#### CDN読み込み
```html
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
```

#### データ埋め込み形式
```javascript
const rmsdData = {
  x: [0.0, 0.01, 0.02, ...],  // XVGの1列目
  y: [0.0, 0.056, 0.062, ...]  // XVGの2列目
};

Plotly.newPlot('rmsd-plot', 
  [{x: rmsdData.x, y: rmsdData.y, type: 'scatter', mode: 'lines'}],
  {title: 'RMSD', xaxis: {title: 'Time (ns)'}, yaxis: {title: 'RMSD (nm)'}}
);
```

#### XVG → JavaScript変換
```bash
echo "const data = {x: [" >> page.md
awk '/^[^@#]/ {printf "%s,", $1}' file.xvg >> page.md
echo "], y: [" >> page.md
awk '/^[^@#]/ {printf "%s,", $2}' file.xvg >> page.md
echo "]};" >> page.md
```

### SMILES取得修正

**問題**: ligandsテーブルからSMILESが取得できなかった

**原因**: テスト用ligandのsmiles列が空

**修正**:
```sql
UPDATE ligands 
SET smiles = 'CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O' 
WHERE id = 1;
```

### targets.description追加

**追加理由**: ターゲット詳細ページで説明文表示

**実装**:
```sql
ALTER TABLE targets ADD COLUMN description TEXT;
UPDATE targets SET description = 'Glycoprotein receptor protein' WHERE id = 1;
```

---

## 生成結果

### ファイル構成

```
docs/pages/
├── index.md                        # トップページ
├── targets/
│   ├── target_1.md                # β2-glycoprotein I
│   ├── target_2.md                # Prothrombin
│   └── target_3.md                # TLR4/MD-2
└── compounds/
    └── compound_1.md              # ZINC00000001詳細
```

### トップページ内容

- Pipeline Architecture図
- Kubernetes Jobs状態
- Target Proteins一覧テーブル
  - PDB ID, Name, Vina Results, MD Done, Link

### ターゲットページ内容

- Protein Info (PDB ID, Description, RCSB link)
- Compounds Table
  - ZINC ID, SMILES (30文字省略), Vina Score, MD Status, Link

### 化合物詳細ページ内容

1. **Compound Information**
   - ZINC ID, SMILES, Target, Vina Affinity

2. **MD Simulation Summary**
   - Status, Time, RMSD, Performance

3. **Equilibration Quality Assessment**
   - EM graph
   - NVT Temperature graph
   - NPT Pressure/Density graphs

4. **Production MD Analysis**
   - RMSD graph
   - H-bond graph

---

## パフォーマンス考察

### 課題
"パフォーマンスはそれだけだと意味なくて。別テーブルに系の情報、原子総数とか入れないとむずいかな"

### 解決策
**md_systemsテーブル**で以下を管理:

- **原子数情報**: total_atoms, protein_atoms, water_molecules
- **システムサイズ**: box dimensions, volume
- **計算条件**: force_field, water_model, timestep
- **性能指標**: performance_nsday + gpu_model

### 性能評価指標

| 系サイズ | 原子数 | GPU | 性能 (ns/day) | 評価 |
|---------|--------|-----|---------------|------|
| Small   | ~40k   | RTX 4070 | 408.7 | ✅ Good |
| Medium  | ~100k  | RTX 4070 | ~150-200 (推定) | - |
| Large   | ~500k  | RTX 4070 | ~30-50 (推定) | - |

---

## GitHub Pages公開

### URL
- **Base**: https://mononobeno.github.io/OCP/
- **Pages**: https://mononobeno.github.io/OCP/pages/index.html
- **Target**: https://mononobeno.github.io/OCP/pages/targets/target_1.html
- **Compound**: https://mononobeno.github.io/OCP/pages/compounds/compound_1.html

### 自動更新フロー

```
MD完了 → analyze_md_trajectory.sh → XVG生成
       → db_register_md_result.sh → DB更新
       → db_register_md_system.sh → システム情報登録
       → generate_pages_with_equilibration.sh → HTML生成
       → git commit & push → GitHub Pages更新
```

---

## 次のステップ

### 短期 (immediate)
- [ ] Ligand RMSD追加 (タンパク質-リガンド複合体)
- [ ] 距離解析グラフ (Protein-Ligand距離)
- [ ] RMSF per-residue heatmap

### 中期 (1-2週間)
- [ ] 複数化合物での一括解析
- [ ] 比較グラフ (複数化合物のRMSD重ね合わせ)
- [ ] エネルギー分解解析

### 長期 (1ヶ月)
- [ ] MM-PBSA結合自由エネルギー
- [ ] Principal Component Analysis (PCA)
- [ ] 3D構造ビューア統合 (Mol*, NGL Viewer)

---

## コミット履歴

- **c54b174**: 階層的ページ構成とPlotlyグラフ実装
- **d2349fb**: 平衡化評価グラフとシステム情報テーブル追加

---

**レポート作成日**: 2025-11-17  
**ステータス**: ✅ 完了
