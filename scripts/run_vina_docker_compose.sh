#!/usr/bin/env bash
set -euo pipefail

# Docker ComposeでVina-runnerを実行（Kubernetes不要）
# DB駆動モード: runs / run_ligands から ZINC ID を取得

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJ_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

cd "$PROJ_ROOT"

# ========================================
# 設定パラメータ (必要に応じて変更)
# ========================================
LIB_CODE="zinc_2d_smi_v1"
TARGET_CODE="aps_apoh"
ZINC_IDS_FILE="catalog/libraries/zinc_2d_smi_v1/processed/ligands_zinc_ids.txt"
DB_PATH="catalog/db/ocp_results.sqlite"

echo "[INFO] ========================================"
echo "[INFO] Vina DB-Driven Execution"
echo "[INFO] ========================================"
echo "[INFO] LIB_CODE     : $LIB_CODE"
echo "[INFO] TARGET_CODE  : $TARGET_CODE"
echo "[INFO] ZINC_IDS_FILE: $ZINC_IDS_FILE"
echo "[INFO] DB_PATH      : $DB_PATH"
echo ""

# 前回の結果ディレクトリがあれば確認
if [ -d "results/vina_output" ]; then
  echo "[WARN] Previous vina_output exists. Contents will be preserved."
  echo "[INFO] To start fresh: rm -rf results/vina_output"
fi

# 必要なディレクトリ作成
mkdir -p results/vina_output/{poses,logs}

# receptor.pdbqtの存在確認
RECEPTOR_PATH="catalog/targets/aps_apoh/receptor.pdbqt"
if [ ! -f "$RECEPTOR_PATH" ]; then
  echo "[ERROR] Receptor file not found: $RECEPTOR_PATH"
  echo "[INFO] Please prepare receptor.pdbqt for your target"
  exit 1
fi

# ligands_zinc_ids.txtの存在確認
if [ ! -f "$ZINC_IDS_FILE" ]; then
  echo "[WARN] ligands_zinc_ids.txt not found, creating sample list..."
  ls catalog/libraries/zinc_2d_smi_v1/processed/pdbqt/*.pdbqt | \
    head -10 | \
    xargs -n1 basename | \
    sed 's/\.pdbqt$//' > "$ZINC_IDS_FILE"
  echo "[INFO] Created sample list with $(wc -l < "$ZINC_IDS_FILE") ligands"
fi

# ========================================
# DB に RUN を登録
# ========================================
echo "[INFO] Creating RUN in database..."
if [ ! -f "scripts/db_create_run_from_ids.sh" ]; then
  echo "[ERROR] scripts/db_create_run_from_ids.sh not found"
  exit 1
fi

# db_create_run_from_ids.sh を実行して RUN_ID / RUN_UUID を取得
eval "$(bash scripts/db_create_run_from_ids.sh "$LIB_CODE" "$TARGET_CODE" "$ZINC_IDS_FILE")"

if [ -z "${RUN_ID:-}" ] || [ -z "${RUN_UUID:-}" ]; then
  echo "[ERROR] Failed to create RUN in database"
  exit 1
fi

echo "[INFO] ========================================"
echo "[INFO] RUN created successfully"
echo "[INFO] RUN_ID   : $RUN_ID"
echo "[INFO] RUN_UUID : $RUN_UUID"
echo "[INFO] ========================================"
echo ""

# run_info.json の作成 (後方互換性のため残す)
cat > "results/vina_output/run_info.json" << EOF
{
  "run_id": $RUN_ID,
  "run_uuid": "$RUN_UUID",
  "library_code": "$LIB_CODE",
  "target_code": "$TARGET_CODE",
  "timestamp": "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
}
EOF

echo "[INFO] Starting vina-runner with Docker Compose..."
echo "[INFO] RUN_ID      : $RUN_ID"
echo "[INFO] RUN_UUID    : $RUN_UUID"
echo "[INFO] Target      : $TARGET_CODE"
echo "[INFO] Library     : $LIB_CODE"
echo ""

# Docker Compose実行 (環境変数を渡す)
export RUN_ID
export RUN_UUID
export OCP_DB_PATH="$DB_PATH"
docker compose -f docker-compose.vina.yml up --build --abort-on-container-exit

echo ""
echo "[INFO] Vina run completed!"
echo "[INFO] Check results:"
echo "  ls -lh results/vina_output/poses/"
echo "  sqlite3 catalog/db/ocp_results.sqlite 'SELECT COUNT(*) FROM vina_results WHERE run_id=$RUN_ID;'"
echo "  sqlite3 catalog/db/ocp_results.sqlite 'SELECT * FROM vina_results WHERE run_id=$RUN_ID ORDER BY affinity_kcal ASC;'"
echo "  sqlite3 catalog/db/ocp_results.sqlite < tools/sql/debug_vina_entities.sql"
