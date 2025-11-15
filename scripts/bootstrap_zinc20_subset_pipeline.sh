#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RAW_DIR="$ROOT_DIR/catalog/libraries/zinc_2d_smi_v1/raw"
DB="$ROOT_DIR/catalog/db/ocp_results.sqlite"

INPUT_SMI_DEFAULT="$RAW_DIR/all_zinc20_ml_subset.smi"
INPUT_SMI="${1:-$INPUT_SMI_DEFAULT}"
MAX_COUNT="${2:-1000}"   # とりあえず 1000 件だけ前処理（0 なら全件）

echo "[INFO] ROOT_DIR = $ROOT_DIR"
echo "[INFO] RAW_DIR  = $RAW_DIR"
echo "[INFO] DB       = $DB"
echo "[INFO] INPUT_SMI = $INPUT_SMI"
echo "[INFO] MAX_COUNT = $MAX_COUNT"

if [ ! -f "$INPUT_SMI" ]; then
  echo "[ERROR] INPUT_SMI not found: $INPUT_SMI" >&2
  echo "  先に scripts/fetch_zinc20_ml_smiles.sh を完了させてください。" >&2
  exit 1
fi

if ! command -v obabel >/dev/null 2>&1; then
  echo "[ERROR] obabel が見つかりません。" >&2
  echo "  conda か apt で Open Babel をインストールしてください。" >&2
  exit 1
fi

if [ ! -x "$ROOT_DIR/scripts/preprocess_smiles_to_pdbqt.sh" ]; then
  echo "[ERROR] scripts/preprocess_smiles_to_pdbqt.sh が見つからないか実行権限がありません。" >&2
  exit 1
fi

# 1. SMILES -> pdbqt 前処理
echo "[STEP 1] SMILES -> pdbqt 前処理を開始します"
"$ROOT_DIR/scripts/preprocess_smiles_to_pdbqt.sh" \
  "$INPUT_SMI" \
  zinc_2d_smi_v1 \
  "$MAX_COUNT"

echo "[STEP 1] 前処理完了。生成された pdbqt の一部を確認:"
ls -lh "$ROOT_DIR"/catalog/libraries/zinc_2d_smi_v1/processed/pdbqt | head

# 2. DB 上の ligands 件数を確認
echo "[STEP 2] DB の ligands 件数確認 (library=zinc_2d_smi_v1)"
sqlite3 "$DB" "SELECT COUNT(*) FROM ligands WHERE library_id = (SELECT id FROM libraries WHERE code='zinc_2d_smi_v1');"

# 3. Kubernetes パイプライン: ligand-selector → vina-prep
echo "[STEP 3] Kubernetes 上で ligand-selector を実行します (cluster=ocp-kind)"
"$ROOT_DIR/scripts/run_ligand_selector_stage.sh"

echo "[STEP 4] Kubernetes 上で vina-prep を実行します (cluster=ocp-kind)"
"$ROOT_DIR/scripts/run_vina_prep_stage.sh"

echo "[DONE] ZINC20 subset を使ったパイプライン前半 (ligand-selector → vina-prep) が完了しました。"
echo "[INFO] /workspace/vina_input 内の中身を確認すると、実際に選ばれた ligand の pdbqt が見えるはずです。"
