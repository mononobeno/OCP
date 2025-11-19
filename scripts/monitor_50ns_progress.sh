#!/bin/bash

echo "=========================================="
echo "50ns MD Simulation Progress Monitor"
echo "=========================================="
echo "現在時刻: $(date '+%Y-%m-%d %H:%M:%S JST')"
echo ""

# Check processes
if pgrep -f "run_3_compounds_50ns" > /dev/null; then
    echo "✅ メインスクリプト: 実行中"
else
    echo "⚠️  メインスクリプト: 停止"
fi

if pgrep -f "gmx.*mdrun" > /dev/null; then
    echo "✅ GROMACS mdrun: 実行中"
else
    echo "⚠️  GROMACS mdrun: 停止"
fi

echo ""
echo "=========================================="
echo "最新の進捗:"
echo "=========================================="

# Get latest step
LATEST=$(grep "step.*will finish" /tmp/50ns_md_run.log 2>/dev/null | tail -1)
if [ -n "$LATEST" ]; then
    STEP=$(echo "$LATEST" | awk '{print $2}' | tr -d ',')
    PROGRESS=$(awk "BEGIN {printf \"%.1f\", ($STEP / 25000000) * 100}")
    echo "ステップ: $STEP / 25,000,000 ($PROGRESS%)"
    echo ""
    
    # Extract time: "Wed Nov 19 08:49:14 2025"
    UTC_STR=$(echo "$LATEST" | grep -oP '[A-Z][a-z]{2}\s+[A-Z][a-z]{2}\s+\d{2}\s+\d{2}:\d{2}:\d{2}\s+\d{4}')
    if [ -n "$UTC_STR" ]; then
        # Parse and add 9 hours for JST
        JST_STR=$(date -d "$UTC_STR UTC + 9 hours" '+%Y年%m月%d日 %H:%M:%S' 2>/dev/null)
        if [ -n "$JST_STR" ]; then
            echo "完了予定 (JST): $JST_STR"
        else
            echo "完了予定 (UTC): $UTC_STR"
        fi
    fi
else
    echo "ステップ情報なし"
fi

echo ""
echo "現在の化合物: $(grep -E 'Compound [0-9]: ZINC' /tmp/50ns_md_run.log 2>/dev/null | tail -1 | grep -oP 'ZINC\S+')"
echo "完了: $(grep '✅ MD simulation completed' /tmp/50ns_md_run.log 2>/dev/null | wc -l) / 3 化合物"

echo ""
echo "=========================================="
