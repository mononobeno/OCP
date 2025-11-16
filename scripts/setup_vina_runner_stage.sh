#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IMG_DIR="$ROOT_DIR/images/vina-runner"
K8S_JOB_DIR="$ROOT_DIR/k8s/jobs"

mkdir -p "$IMG_DIR" "$K8S_JOB_DIR"

echo "[INFO] ROOT_DIR = $ROOT_DIR"
echo "[INFO] IMG_DIR  = $IMG_DIR"

########################################
# 1. Dockerfile（vina-runner イメージ）
########################################
cat > "$IMG_DIR/Dockerfile" << 'EOF_DF'
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      ca-certificates \
      wget \
      sqlite3 \
      jq \
      autodock-vina && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY vina_runner.sh /app/vina_runner.sh

RUN chmod +x /app/vina_runner.sh

# デフォルトのエントリポイント
ENTRYPOINT ["/app/vina_runner.sh"]
EOF_DF

########################################
# 2. コンテナ内のメインスクリプト
########################################
cat > "$IMG_DIR/vina_runner.sh" << 'EOF_RUN'
#!/usr/bin/env bash
set -euo pipefail

echo "[INFO] vina-runner starting..."

# ---- 環境変数とデフォルト ----
WORKSPACE="${WORKSPACE:-/workspace}"
DB_PATH="${DB_PATH:-/zinc-library/db/ocp_results.sqlite}"
TARGET_CODE="${TARGET_CODE:-aps_apoh}"
LIB_CODE="${LIB_CODE:-zinc_2d_smi_v1}"

VINA_EXHAUSTIVENESS="${VINA_EXHAUSTIVENESS:-8}"
VINA_NUM_MODES="${VINA_NUM_MODES:-9}"
VINA_CPU="${VINA_CPU:-0}"   # 0 = all

VINA_IN_DIR="$WORKSPACE/vina_input"
VINA_OUT_DIR="$WORKSPACE/vina_output"
RAW_DIR="$WORKSPACE/ligands_raw"

RECEPTOR_PDBQT="$VINA_IN_DIR/receptor.pdbqt"
LIGAND_DIR="$VINA_IN_DIR/ligands"

mkdir -p "$VINA_OUT_DIR/poses" "$VINA_OUT_DIR/logs"

echo "[INFO] DB_PATH    = $DB_PATH"
echo "[INFO] WORKSPACE  = $WORKSPACE"
echo "[INFO] TARGET_CODE= $TARGET_CODE"
echo "[INFO] LIB_CODE   = $LIB_CODE"
echo "[INFO] VINA_IN_DIR= $VINA_IN_DIR"
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

# ---- vina_results テーブルを準備（なければ作る）----
sqlite3 "$DB_PATH" << 'EOSQL'
CREATE TABLE IF NOT EXISTS vina_results (
  id INTEGER PRIMARY KEY,
  run_id INTEGER,
  ligand_id INTEGER,
  zinc_id TEXT,
  library_id INTEGER,
  score REAL,
  created_at TEXT DEFAULT (datetime('now'))
);
EOSQL

# ---- RUN 情報（run_info.json）を読む ----
RUN_INFO_JSON="$RAW_DIR/run_info.json"
RUN_ID=""
RUN_UUID=""

if [ -f "$RUN_INFO_JSON" ]; then
  echo "[INFO] reading run_info.json: $RUN_INFO_JSON"
  RUN_ID="$(jq -r '.run_id // empty' "$RUN_INFO_JSON" || true)"
  RUN_UUID="$(jq -r '.run_uuid // empty' "$RUN_INFO_JSON" || true)"
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
TMP_ZINC_LIST="/tmp/zinc_ids_from_db.txt"

if [ -n "$RUN_ID" ]; then
  echo "[INFO] fetching zinc_id list from DB for run_id=$RUN_ID"
  sqlite3 "$DB_PATH" << EOSQL > "$TMP_ZINC_LIST"
SELECT l.zinc_id
FROM run_ligands rl
JOIN ligands l   ON l.id = rl.ligand_id
JOIN runs   r    ON r.id = rl.run_id
WHERE r.id = $RUN_ID
  AND l.library_id = $LIB_ID
  AND l.zinc_id IS NOT NULL
ORDER BY rl.id;
EOSQL

  if [ -s "$TMP_ZINC_LIST" ]; then
    echo "[INFO] got \$(wc -l < "$TMP_ZINC_LIST") ligands from DB."
  else
    echo "[WARN] no ligands found in DB for run_id=$RUN_ID, fallback to ligands_zinc_ids.txt" >&2
  fi
else
  echo "[WARN] RUN_ID is empty, fallback to ligands_zinc_ids.txt" >&2
fi

USE_FILE=""
if [ -s "$TMP_ZINC_LIST" ]; then
  USE_FILE="$TMP_ZINC_LIST"
elif [ -f "$ZINC_ID_LIST_FILE" ]; then
  echo "[INFO] using ZINC IDs from $ZINC_ID_LIST_FILE"
  USE_FILE="$ZINC_ID_LIST_FILE"
else
  echo "[ERROR] no ZINC ID list found (DB and file both missing/empty)." >&2
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

# box parameters (暫定値・必要なら後で ConfigMap から差し替え)
center_x = 0
center_y = 0
center_z = 0

size_x = 20
size_y = 20
size_z = 20

exhaustiveness = $VINA_EXHAUSTIVENESS
num_modes      = $VINA_NUM_MODES
EOF_CFG

  echo "[INFO] [$i/$TOTAL] running vina for $ZINC_ID_TRIM"
  if vina --config "$CONFIG_FILE" --out "$OUT_PDBQT" --log "$LOG_FILE" --cpu "$VINA_CPU" >/dev/null 2>&1; then
    # Vina ログからスコア抽出
    BEST_SCORE="$(grep -m1 '^REMARK VINA RESULT:' "$LOG_FILE" | awk '{print $4}' || echo "")"
    echo "[INFO] [$i/$TOTAL] best score = ${BEST_SCORE:-N/A}"

    # ligand_id を取得
    LIG_ID="$(sqlite3 "$DB_PATH" "SELECT id FROM ligands WHERE zinc_id='$ZINC_ID_TRIM' AND library_id=$LIB_ID;" || true)"

    if [ -n "$LIG_ID" ] && [ -n "$RUN_ID" ] && [ -n "$BEST_SCORE" ]; then
      sqlite3 "$DB_PATH" << EOSQL
INSERT INTO vina_results (run_id, ligand_id, zinc_id, library_id, score)
VALUES ($RUN_ID, $LIG_ID, '$ZINC_ID_TRIM', $LIB_ID, $BEST_SCORE);
EOSQL
    fi
  else
    echo "[WARN] [$i/$TOTAL] vina failed for $ZINC_ID_TRIM"
  fi

done < "$USE_FILE"

echo "[INFO] vina-runner finished."
EOF_RUN

chmod +x "$IMG_DIR/vina_runner.sh"

########################################
# 3. K8s Job マニフェスト
########################################
cat > "$K8S_JOB_DIR/job-vina-runner.yaml" << 'EOF_JOB'
apiVersion: batch/v1
kind: Job
metadata:
  name: vina-runner
  namespace: drug-pipeline
spec:
  backoffLimit: 0
  template:
    spec:
      restartPolicy: Never
      containers:
        - name: vina-runner
          image: pipeline/vina-runner:local
          imagePullPolicy: IfNotPresent
          envFrom:
            - configMapRef:
                name: pipeline-config
          env:
            - name: WORKSPACE
              value: /workspace
            - name: DB_PATH
              value: /zinc-library/db/ocp_results.sqlite
          volumeMounts:
            - name: zinc-library-pvc
              mountPath: /zinc-library
            - name: work-pvc
              mountPath: /workspace
      volumes:
        - name: zinc-library-pvc
          persistentVolumeClaim:
            claimName: zinc-library-pvc
        - name: work-pvc
          persistentVolumeClaim:
            claimName: work-pvc
EOF_JOB

########################################
# 4. イメージ build & kind へ load
########################################
cat > "$ROOT_DIR/scripts/build_and_load_vina_runner_image.sh" << 'EOF_BUILD'
#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IMG_DIR="$ROOT_DIR/images/vina-runner"

IMAGE_NAME="pipeline/vina-runner:local"
KIND_CLUSTER_NAME="${KIND_CLUSTER_NAME:-ocp-kind}"

echo "[INFO] building image: $IMAGE_NAME"
cd "$IMG_DIR"
docker build -t "$IMAGE_NAME" .

echo "[INFO] loading image into kind cluster: $KIND_CLUSTER_NAME"
kind load docker-image "$IMAGE_NAME" --name "$KIND_CLUSTER_NAME"

echo "[INFO] build_and_load_vina_runner_image.sh done."
EOF_BUILD

chmod +x "$ROOT_DIR/scripts/build_and_load_vina_runner_image.sh"

########################################
# 5. Job 実行用スクリプト
########################################
cat > "$ROOT_DIR/scripts/run_vina_runner_stage.sh" << 'EOF_RUNSTAGE'
#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
KIND_CLUSTER_NAME="${KIND_CLUSTER_NAME:-ocp-kind}"
NS="drug-pipeline"

echo "[INFO] building & loading vina-runner image..."
"$ROOT_DIR/scripts/build_and_load_vina_runner_image.sh"

echo "[INFO] applying namespace/pvc/config (idempotent)"
kubectl apply -f "$ROOT_DIR/k8s/base/namespace.yaml"
kubectl apply -f "$ROOT_DIR/k8s/base/storage-pvc.yaml"
kubectl apply -f "$ROOT_DIR/k8s/config/pipeline-config.yaml"

echo "[INFO] deleting old vina-runner job if exists..."
kubectl -n "$NS" delete job vina-runner --ignore-not-found

echo "[INFO] applying job-vina-runner.yaml..."
kubectl apply -f "$ROOT_DIR/k8s/jobs/job-vina-runner.yaml"

echo "[INFO] waiting for job completion..."
kubectl -n "$NS" wait --for=condition=complete job/vina-runner --timeout=600s

echo "[INFO] job completed. logs:"
kubectl -n "$NS" logs job/vina-runner
EOF_RUNSTAGE

chmod +x "$ROOT_DIR/scripts/run_vina_runner_stage.sh"

echo "[INFO] setup_vina_runner_stage.sh done."
