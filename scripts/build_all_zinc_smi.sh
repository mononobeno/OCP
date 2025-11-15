#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RAW_DIR="$ROOT_DIR/catalog/libraries/zinc_2d_smi_v1/raw"

echo "[INFO] ROOT_DIR = $ROOT_DIR"
echo "[INFO] RAW_DIR  = $RAW_DIR"

if ! command -v aria2c >/dev/null 2>&1; then
  echo "ERROR: aria2c が見つかりません。" >&2
  echo "  sudo apt-get install -y aria2 などでインストールしてください。" >&2
  exit 1
fi

cd "$RAW_DIR"

# ① .uri を自動検出
URI_FILE="$(ls *smi*.uri 2>/dev/null | head -n1 || true)"

if [ -z "$URI_FILE" ]; then
  echo "[ERROR] *.uri ファイルが RAW_DIR に見つかりません: $RAW_DIR" >&2
  ls -l
  exit 1
fi

echo "[INFO] 使用する URI ファイル: $URI_FILE"

# 古い all_zinc.smi があれば削除（空ファイル対策）
if [ -f all_zinc.smi ]; then
  echo "[INFO] 既存 all_zinc.smi を削除します"
  rm -f all_zinc.smi
fi

# ② aria2c で .smi.gz をダウンロード
echo "[INFO] aria2c で .smi.gz をダウンロード中..."
aria2c -i "$RAW_DIR/$URI_FILE" -x 16 -s 16 -k 1M

# ③ .smi.gz を解凍
if ls "$RAW_DIR"/*.smi.gz >/dev/null 2>&1; then
  echo "[INFO] .smi.gz を解凍中..."
  gunzip -f "$RAW_DIR"/*.smi.gz
else
  echo "[WARN] .smi.gz が1つもありません。ダウンロードが失敗している可能性があります。" >&2
fi

# ④ .smi を統合
if ls "$RAW_DIR"/*.smi >/dev/null 2>&1; then
  echo "[INFO] .smi を結合して all_zinc.smi を作成します..."
  cat "$RAW_DIR"/*.smi > "$RAW_DIR/all_zinc.smi"
  echo "[OK] all_zinc.smi 作成完了:"
  ls -lh "$RAW_DIR/all_zinc.smi"
else
  echo "[WARN] .smi ファイルが1つもありません。all_zinc.smi は作成されませんでした。" >&2
fi
