#!/usr/bin/env bash
#
# k8s_job_eq_analysis_pages.sh
#
# 役割:
#   - 平衡化結果 (NVT / NPT) を解析し、GitHub Pages 用レポートを生成するステージ。
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
GMX_EQ_DIR="${WORKSPACE}/gmx_eq"
DOCS_DIR="${REPO_ROOT}/docs/equilibration"

mkdir -p "${DOCS_DIR}"

echo "[INFO] eq-analysis-pages: RUN_ID=${RUN_ID}"
echo "[INFO] GMX_EQ  : ${GMX_EQ_DIR}"
echo "[INFO] DOCS_DIR: ${DOCS_DIR}"

REPORT_MD="${DOCS_DIR}/run_${RUN_ID}_eq_report.md"

cat > "${REPORT_MD}" <<RPT
# Equilibration Report (RUN_ID=${RUN_ID})

- Workspace: ${GMX_EQ_DIR}
- This is a skeleton report. Use GROMACS tools to extract temperature / pressure /
  energy vs time and embed plots here in future.

RPT

echo "[INFO] Generated report: ${REPORT_MD}"
echo "[INFO] eq-analysis-pages finished (skeleton)."
