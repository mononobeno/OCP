#!/usr/bin/env bash
#
# k8s_job_min_analysis_pages.sh
#
# 役割:
#   - 最小化結果を評価し、GitHub Pages 用レポートを生成するステージ。
#   - 入力: /workspace/gmx_em/em.edr, em.log など。
#   - 出力: docs/minimization/ 配下の Markdown/HTML/PNG、Git コミット & push（将来的）。
#
# 共通仕様:
#   - RUN_ID は第1引数、または環境変数 RUN_ID から取得する。
#   - DB_PATH は環境変数 DB_PATH_OVERRIDE で上書き可能。
#   - DB は catalog/db/ocp_results.sqlite を前提とし、
#     ここから対象 ligand / target / run を決める（DB = source of truth）。
#
set -euo pipefail

RUN_ID="${1:-${RUN_ID:-}}"
if [[ -z "${RUN_ID}" ]]; then
  echo "Usage: $0 <RUN_ID>  # or set RUN_ID env" >&2
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="${SCRIPT_DIR%/scripts}"
WORKSPACE="${WORKSPACE:-/workspace}"
GMX_EM_DIR="${WORKSPACE}/gmx_em"
DOCS_DIR="${REPO_ROOT}/docs/minimization"

mkdir -p "${DOCS_DIR}"

echo "[INFO] min-analysis-pages: RUN_ID=${RUN_ID}"
echo "[INFO] GMX_EM  : ${GMX_EM_DIR}"
echo "[INFO] DOCS_DIR: ${DOCS_DIR}"

# ここで gmx energy / gmx rms などを用いて最小化結果を解析し、
# 画像 (PNG) やテキストを生成して docs/ に配置する想定。
REPORT_MD="${DOCS_DIR}/run_${RUN_ID}_min_report.md"

cat > "${REPORT_MD}" <<RPT
# Minimization Report (RUN_ID=${RUN_ID})

- Workspace: ${GMX_EM_DIR}
- This is a skeleton report. Use GROMACS tools to extract energy vs step,
  and embed plots here in future.

RPT

echo "[INFO] Generated report: ${REPORT_MD}"
echo "[INFO] min-analysis-pages finished (skeleton)."
