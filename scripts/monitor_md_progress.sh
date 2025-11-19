#!/bin/bash
# Monitor real 10ns MD simulation progress

LOG_FILE="/tmp/real_md_run.log"

echo "=========================================="
echo "Real 10ns GPU MD Simulation Monitor"
echo "=========================================="
echo "Start time: $(date '+%Y-%m-%d %H:%M:%S')"
echo ""

# Check if MD is running
if pgrep -f "gmx.*mdrun" > /dev/null; then
    echo "✅ GROMACS mdrun is RUNNING"
    echo "PID: $(pgrep -f 'gmx.*mdrun')"
else
    echo "⚠️  GROMACS mdrun process not found"
fi

echo ""
echo "----------------------------------------"
echo "Latest Progress:"
echo "----------------------------------------"
grep "step.*will finish" "$LOG_FILE" | tail -1

echo ""
echo "----------------------------------------"
echo "PME Grid Optimization:"
echo "----------------------------------------"
grep "optimal pme grid" "$LOG_FILE" | tail -1

echo ""
echo "----------------------------------------"
echo "Recent Output (last 15 lines):"
echo "----------------------------------------"
tail -15 "$LOG_FILE"

echo ""
echo "=========================================="
echo "Monitor command: watch -n 10 bash scripts/monitor_md_progress.sh"
echo "Log file: tail -f $LOG_FILE"
echo "=========================================="
