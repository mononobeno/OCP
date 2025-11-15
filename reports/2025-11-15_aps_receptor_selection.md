# APS（抗リン脂質抗体症候群）ターゲット受容体選定レポート

- 作成日: 2025-11-15
- プロジェクト: OCP (Open Container Pipeline) – APS 解析系
- 対象疾患: 抗リン脂質抗体症候群（Antiphospholipid Syndrome; APS）

---

## 1. 背景

抗リン脂質抗体症候群（APS）は、抗リン脂質抗体（aPL）が持続的に陽性となり、
動静脈血栓や妊娠合併症を引き起こす自己免疫疾患である。
aPL は、陰性荷電リン脂質そのものだけでなく、リン脂質に結合した蛋白質
（β2-glycoprotein I, prothrombin, annexin など）を標的とすることが知られている。

本パイプラインでは「小分子ドッキング + MD 解析によって、APS 病態に関わる分子相互作用を
定量的に評価する」ことを目的として、以下の受容体を第 1 世代の標的セットとして採用する。

- β2-glycoprotein I（APOH）
- Prothrombin（F2）
- TLR4/MD-2 複合体

選定基準は以下の通り。

- APS 病態への関与が文献で繰り返し報告されている
- PDB に高品質な立体構造が存在し、ドッキンググリッドを定義しやすい
- 将来的に薬剤設計（阻害剤・修飾剤）が合理的に行えるドメイン構造を持つ

---

## 2. 個別ターゲットの選定理由

### 2.1 β2-glycoprotein I（APOH） – `aps_apoh`

- APS における最重要抗原の一つであり、抗 β2GPI 抗体だけで血栓形成を惹起しうる動物実験が多数報告されている。
- 第 V ドメインがリン脂質結合に重要であり、自己抗体はこの近傍を認識するとされる。
- 全長構造が解かれており、リン脂質結合部位や抗体エピトープ候補の立体配置を評価しやすい。
- 小分子によって β2GPI–リン脂質相互作用や β2GPI–自己抗体相互作用を遮断するという薬剤コンセプトを構築しやすい。

採用 PDB 構造:

- PDB ID: **1C1Z** — Human β2-glycoprotein I, full-length structure

OCP 内格納パス:

- `catalog/targets/aps_apoh/receptor.pdb`

---

### 2.2 Prothrombin（F2） – `aps_prothrombin`

- β2GPI と並ぶ APS の主要抗原であり、抗プロトロンビン抗体は血栓リスクと関連する。
- プロトロンビンは凝固カスケードの中心因子であり、異常な活性化は血栓形成と直結する。
- 全長プロトロンビンの closed コンフォメーション構造が報告されており、allosteric site を含む立体構造に基づいた設計が可能。

採用 PDB 構造:

- PDB ID: **6C2W** — Human prothrombin, closed conformation

OCP 内格納パス:

- `catalog/targets/aps_prothrombin/receptor.pdb`

---

### 2.3 TLR4/MD-2 複合体 – `aps_tlr4_md2`

- APS 患者由来 aPL が TLR4 経路を介して炎症・血栓促進シグナルを誘導することが報告されている。
- TLR4/MD-2 は LPS などの脂質リガンドを認識するシグナル受容体であり、既に小分子アンタゴニストの構造情報が豊富に存在する。
- aPL 依存的に活性化される TLR4 シグナルを小分子で抑制することで、炎症・血栓の二次的増悪を抑える戦略が立てやすい。

採用 PDB 構造:

- PDB ID: **3FXI** — Human TLR4–MD-2–LPS complex

OCP 内格納パス:

- `catalog/targets/aps_tlr4_md2/receptor.pdb`

---

## 3. 今後の拡張候補

- Annexin A5（ANXA5）: aPL による annexin A5 シールド破壊が APS 血栓形成に関与するとされる。
- 補体 C5 / C5aR: APS における補体活性化が病態形成に関与することが示されており、既に C5 阻害薬が臨床応用されている。

これらは現時点では OCP パイプラインの「主要 3 ターゲット」には含めていないが、
将来のモジュール追加候補として検討する。
