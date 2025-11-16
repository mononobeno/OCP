#!/usr/bin/env bash
set -euo pipefail

# ローカルでvina-runnerコンテナをテストするスクリプト
# (Kubernetes不要、Dockerだけで動作確認)

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJ_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

IMAGE_NAME="vina-runner:local"
CONTAINER_NAME="vina-test-$(date +%s)"

echo "[TEST] Building vina-runner image..."
docker build -t "$IMAGE_NAME" "$PROJ_ROOT/images/vina-runner/"

echo "[TEST] Creating test workspace..."
TEST_WS="/tmp/vina-test-workspace"
rm -rf "$TEST_WS"
mkdir -p "$TEST_WS"/{receptor,ligands/pdbqt,vina_out/{logs,poses},db}

# サンプルreceptor（ダミー）
cat > "$TEST_WS/receptor/receptor.pdbqt" << 'EOF'
REMARK  dummy receptor for testing
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00     0.000 C
ROOT
ENDROOT
TORSDOF 0
EOF

# サンプルligand（ダミー）
cat > "$TEST_WS/ligands/pdbqt/ZINC000000000001_1.pdbqt" << 'EOF'
REMARK  dummy ligand for testing
ATOM      1  C   LIG A   1       1.000   1.000   1.000  1.00  0.00     0.000 C
ROOT
ENDROOT
TORSDOF 0
EOF

# ZINC IDリスト
echo "ZINC000000000001_1" > "$TEST_WS/ligands/ligands_zinc_ids.txt"

# DB作成
DB_PATH="$TEST_WS/db/ocp_results.sqlite"
sqlite3 "$DB_PATH" << EOSQL
CREATE TABLE IF NOT EXISTS libraries (
  id INTEGER PRIMARY KEY,
  code TEXT UNIQUE NOT NULL,
  name TEXT
);
INSERT OR IGNORE INTO libraries (id, code, name) VALUES (1, 'test_lib', 'Test Library');

CREATE TABLE IF NOT EXISTS ligands (
  id INTEGER PRIMARY KEY,
  library_id INTEGER NOT NULL,
  zinc_id TEXT,
  FOREIGN KEY (library_id) REFERENCES libraries(id)
);
INSERT OR IGNORE INTO ligands (id, library_id, zinc_id) VALUES (1, 1, 'ZINC000000000001');

CREATE TABLE IF NOT EXISTS targets (
  id INTEGER PRIMARY KEY,
  code TEXT UNIQUE NOT NULL,
  name TEXT
);
INSERT OR IGNORE INTO targets (id, code, name) VALUES (1, 'test_target', 'Test Target');

CREATE TABLE IF NOT EXISTS runs (
  id INTEGER PRIMARY KEY,
  run_uuid TEXT UNIQUE NOT NULL,
  target_id INTEGER,
  library_id INTEGER,
  FOREIGN KEY (target_id) REFERENCES targets(id),
  FOREIGN KEY (library_id) REFERENCES libraries(id)
);
INSERT OR IGNORE INTO runs (id, run_uuid, target_id, library_id) 
VALUES (1, 'test-run-uuid-12345', 1, 1);

CREATE TABLE IF NOT EXISTS vina_results (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  run_id INTEGER NOT NULL,
  ligand_id INTEGER NOT NULL,
  zinc_id TEXT,
  library_id INTEGER,
  score REAL,
  created_at TEXT DEFAULT (datetime('now')),
  FOREIGN KEY (run_id) REFERENCES runs(id),
  FOREIGN KEY (ligand_id) REFERENCES ligands(id),
  FOREIGN KEY (library_id) REFERENCES libraries(id)
);
EOSQL

# run_info.json
cat > "$TEST_WS/run_info.json" << 'EOF'
{
  "RUN_ID": 1,
  "RUN_UUID": "test-run-uuid-12345",
  "TARGET_ID": 1,
  "LIB_ID": 1
}
EOF

echo "[TEST] Running vina-runner container..."
docker run --rm --name "$CONTAINER_NAME" \
  -v "$TEST_WS:/workspace:rw" \
  -e WORKSPACE=/workspace \
  -e DB_PATH=/workspace/db/ocp_results.sqlite \
  -e TARGET_CODE=test_target \
  -e LIB_CODE=test_lib \
  -e CENTER_X=0.0 \
  -e CENTER_Y=0.0 \
  -e CENTER_Z=0.0 \
  -e SIZE_X=10 \
  -e SIZE_Y=10 \
  -e SIZE_Z=10 \
  -e VINA_CPU=2 \
  "$IMAGE_NAME"

echo ""
echo "[TEST] Checking results in DB..."
sqlite3 "$DB_PATH" "SELECT * FROM vina_results;"

echo ""
echo "[TEST] Done. Test workspace: $TEST_WS"
