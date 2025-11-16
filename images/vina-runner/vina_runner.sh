#!/usr/bin/env bash
set -euo pipefail

echo "[INFO] vina-runner starting..."

# ---- 環境変数とデフォルト ----
WORKSPACE="${WORKSPACE:-/workspace}"
DB_PATH="${DB_PATH:-/workspace/db/ocp_results.sqlite}"
TARGET_CODE="${TARGET_CODE:-aps_apoh}"
LIB_CODE="${LIB_CODE:-zinc_2d_smi_v1}"

VINA_EXHAUSTIVENESS="${VINA_EXHAUSTIVENESS:-8}"
VINA_NUM_MODES="${VINA_NUM_MODES:-9}"
VINA_CPU="${VINA_CPU:-0}"   # 0 = all

# Box parameters (環境変数でオーバーライド可能)
CENTER_X="${CENTER_X:-0.0}"
CENTER_Y="${CENTER_Y:-0.0}"
CENTER_Z="${CENTER_Z:-0.0}"
SIZE_X="${SIZE_X:-20}"
SIZE_Y="${SIZE_Y:-20}"
SIZE_Z="${SIZE_Z:-20}"

# ディレクトリ構造
RAW_DIR="$WORKSPACE/ligands"
LIGAND_DIR="$RAW_DIR/pdbqt"
RECEPTOR_PDBQT="$WORKSPACE/receptor/receptor.pdbqt"
VINA_OUT_DIR="$WORKSPACE/vina_out"

mkdir -p "$VINA_OUT_DIR/poses" "$VINA_OUT_DIR/logs"

echo "[INFO] DB_PATH    = $DB_PATH"
echo "[INFO] WORKSPACE  = $WORKSPACE"
echo "[INFO] TARGET_CODE= $TARGET_CODE"
echo "[INFO] LIB_CODE   = $LIB_CODE"
echo "[INFO] RECEPTOR   = $RECEPTOR_PDBQT"
echo "[INFO] LIGAND_DIR = $LIGAND_DIR"
echo "[INFO] VINA_OUT_DIR= $VINA_OUT_DIR"

if [ ! -f "$DB_PATH" ]; then
  echo "[ERROR] DB not found: $DB_PATH" >&2
  exit 1
fi

if [ ! -f "$RECEPTOR_PDBQT" ]; then
  echo "[ERROR] receptor.pdbqt not found: $RECEPTOR_PDBQT" >&2
  exit 1
fi

if [ ! -d "$LIGAND_DIR" ]; then
  echo "[ERROR] ligand dir not found: $LIGAND_DIR" >&2
  exit 1
fi

if ! command -v vina >/dev/null 2>&1; then
  echo "[ERROR] vina command not found in container." >&2
  exit 1
fi

if ! command -v sqlite3 >/dev/null 2>&1; then
  echo "[ERROR] sqlite3 command not found in container." >&2
  exit 1
fi

if ! command -v jq >/dev/null 2>&1; then
  echo "[ERROR] jq command not found in container." >&2
  exit 1
fi

# ---- vina_results テーブルは既に存在するはず（スキーマ準拠）----
# 既存スキーマ: vina_results(run_id, ligand_id, mode_rank, affinity_kcal, rmsd_lb, rmsd_ub, out_relpath)
# 注: このスクリプトでは mode_rank=1 (best mode) のみ記録する

# ---- RUN 情報を環境変数または run_info.json から取得 ----
RUN_ID="${RUN_ID:-}"
RUN_UUID="${RUN_UUID:-}"

# 環境変数になければ run_info.json から読む (後方互換性)
if [ -z "$RUN_ID" ] || [ -z "$RUN_UUID" ]; then
  RUN_INFO_JSON="$RAW_DIR/run_info.json"
  if [ -f "$RUN_INFO_JSON" ]; then
    echo "[INFO] reading run_info.json: $RUN_INFO_JSON"
    RUN_ID="$(jq -r '.run_id // empty' "$RUN_INFO_JSON" || true)"
    RUN_UUID="$(jq -r '.run_uuid // empty' "$RUN_INFO_JSON" || true)"
  fi
fi

echo "[INFO] RUN_ID   = ${RUN_ID:-<none>}"
echo "[INFO] RUN_UUID = ${RUN_UUID:-<none>}"

# library_id の取得
LIB_ID="$(sqlite3 "$DB_PATH" "SELECT id FROM libraries WHERE code='$LIB_CODE';" || true)"
if [ -z "$LIB_ID" ]; then
  echo "[ERROR] library code '$LIB_CODE' not found in DB." >&2
  exit 1
fi
echo "[INFO] LIB_ID = $LIB_ID"

# ---- 処理対象 ZINC ID の決定 ----
ZINC_ID_LIST_FILE="$RAW_DIR/ligands_zinc_ids.txt"
USE_FILE=""

# NOTE: 現在のスキーマには run_ligands テーブルが存在しないため、
#       常にファイルから ZINC ID を読み込む
if [ -f "$ZINC_ID_LIST_FILE" ]; then
  echo "[INFO] using ZINC IDs from file: $ZINC_ID_LIST_FILE"
  USE_FILE="$ZINC_ID_LIST_FILE"
else
  echo "[ERROR] ZINC ID list file not found: $ZINC_ID_LIST_FILE" >&2
  exit 1
fi

TOTAL=$(wc -l < "$USE_FILE" || echo 0)
echo "[INFO] total ligands to process: $TOTAL"

if [ "$TOTAL" -eq 0 ]; then
  echo "[WARN] nothing to do, exiting."
  exit 0
fi

# ---- Vina 実行ループ ----
i=0
while IFS=$'\t ' read -r ZINC_ID _REST; do
  ZINC_ID_TRIM="$(echo "$ZINC_ID" | tr -d '[:space:]')"
  [ -z "$ZINC_ID_TRIM" ] && continue

  i=$((i+1))
  LIG_PDBQT="$LIGAND_DIR/${ZINC_ID_TRIM}.pdbqt"
  OUT_PDBQT="$VINA_OUT_DIR/poses/${ZINC_ID_TRIM}.pdbqt"
  LOG_FILE="$VINA_OUT_DIR/logs/${ZINC_ID_TRIM}.log"
  CONFIG_FILE="$VINA_OUT_DIR/config_${ZINC_ID_TRIM}.txt"

  if [ ! -f "$LIG_PDBQT" ]; then
    echo "[WARN] [$i/$TOTAL] ligand pdbqt missing, skip: $LIG_PDBQT"
    continue
  fi

  cat > "$CONFIG_FILE" << EOF_CFG
receptor = $RECEPTOR_PDBQT
ligand   = $LIG_PDBQT

center_x = $CENTER_X
center_y = $CENTER_Y
center_z = $CENTER_Z

size_x = $SIZE_X
size_y = $SIZE_Y
size_z = $SIZE_Z

exhaustiveness = ${VINA_EXHAUSTIVENESS:-8}
num_modes      = ${VINA_NUM_MODES:-9}
EOF_CFG

  echo "[INFO] [$i/$TOTAL] running vina for $ZINC_ID_TRIM"
  # Vina 1.2.3は--logオプションがないので、標準出力をリダイレクト
  if vina --config "$CONFIG_FILE" --out "$OUT_PDBQT" --cpu "$VINA_CPU" > "$LOG_FILE" 2>&1; then
    # Vinaのログからベストスコアを抽出（mode 1の行）
    BEST_SCORE="$(grep '   1 ' "$LOG_FILE" | head -1 | awk '{print $2}' || echo "")"
    
    if [ -z "$BEST_SCORE" ]; then
      # 別のフォーマットを試す
      BEST_SCORE="$(grep -E '^\s+1\s+' "$OUT_PDBQT" | head -1 | awk '{print $4}' || echo "")"
    fi
    
    echo "[INFO] [$i/$TOTAL] best score = ${BEST_SCORE:-N/A}"

    # ligand_id を取得
    LIG_ID="$(sqlite3 "$DB_PATH" "SELECT id FROM ligands WHERE zinc_id='$ZINC_ID_TRIM' AND library_id=$LIB_ID;" 2>/dev/null || true)"

    # DB に保存 (実際のスキーマに合わせる)
    # vina_results(run_id, ligand_id, mode_rank, affinity_kcal, rmsd_lb, rmsd_ub, out_relpath)
    if [ -n "$LIG_ID" ] && [ -n "$RUN_ID" ] && [ -n "$BEST_SCORE" ]; then
      POSE_RELPATH="vina_out/poses/${ZINC_ID_TRIM}.pdbqt"
      sqlite3 "$DB_PATH" << EOSQL
INSERT OR REPLACE INTO vina_results (run_id, ligand_id, mode_rank, affinity_kcal, out_relpath)
VALUES ($RUN_ID, $LIG_ID, 1, $BEST_SCORE, '$POSE_RELPATH');
EOSQL
      echo "[INFO] [$i/$TOTAL] saved to DB: run_id=$RUN_ID, ligand_id=$LIG_ID, score=$BEST_SCORE"
    else
      echo "[WARN] [$i/$TOTAL] skipping DB insert: LIG_ID=$LIG_ID, RUN_ID=$RUN_ID, SCORE=$BEST_SCORE"
    fi
  else
    echo "[WARN] [$i/$TOTAL] vina failed for $ZINC_ID_TRIM"
  fi

done < "$USE_FILE"

echo "[INFO] vina-runner finished."
