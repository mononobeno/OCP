#!/usr/bin/env python3
"""
50ns MD Simulation Progress Monitor
Displays current progress with JST time conversion
"""

import re
import subprocess
from datetime import datetime, timedelta

print("=" * 50)
print("50ns MD Simulation Progress Monitor")
print("=" * 50)
print(f"現在時刻: {datetime.now().strftime('%Y-%m-%d %H:%M:%S JST')}")
print()

# Check if main script is running
main_running = subprocess.run(['pgrep', '-f', 'run_3_compounds_50ns'], 
                              capture_output=True).returncode == 0
print(f"{'✅' if main_running else '⚠️ '} メインスクリプト: {'実行中' if main_running else '停止'}")

# Check if GROMACS mdrun is running
mdrun_running = subprocess.run(['pgrep', '-f', 'gmx.*mdrun'], 
                               capture_output=True).returncode == 0
print(f"{'✅' if mdrun_running else '⚠️ '} GROMACS mdrun: {'実行中' if mdrun_running else '停止'}")

print()
print("=" * 50)
print("最新の進捗:")
print("=" * 50)

# Read log file
try:
    with open('/tmp/50ns_md_run.log', 'r') as f:
        lines = f.readlines()
    
    # Find latest "step X will finish ..." line
    step_lines = [l for l in lines if 'step' in l and 'will finish' in l]
    
    if step_lines:
        latest = step_lines[-1]
        
        # Extract step number
        step_match = re.search(r'step\s+(\d+)', latest)
        if step_match:
            step = int(step_match.group(1))
            progress = (step / 25_000_000) * 100
            print(f"ステップ: {step:,} / 25,000,000 ({progress:.1f}%)")
            print()
        
        # Extract UTC time: "Wed Nov 19 08:49:14 2025"
        time_match = re.search(r'([A-Z][a-z]{2}\s+[A-Z][a-z]{2}\s+\d{2}\s+\d{2}:\d{2}:\d{2}\s+\d{4})', latest)
        if time_match:
            utc_str = time_match.group(1)
            
            # Parse UTC time
            try:
                utc_time = datetime.strptime(utc_str, '%a %b %d %H:%M:%S %Y')
                # Convert to JST (UTC + 9 hours)
                jst_time = utc_time + timedelta(hours=9)
                print(f"完了予定 (JST): {jst_time.strftime('%Y年%m月%d日 %H:%M:%S')}")
                print(f"完了予定 (UTC): {utc_str}")
            except ValueError:
                print(f"完了予定 (UTC): {utc_str}")
    else:
        print("ステップ情報なし")
    
    # Count completed compounds
    completed = len([l for l in lines if 'Finished mdrun' in l])
    print()
    print(f"現在の化合物: {completed + 1} / 3")
    print(f"完了: {completed} / 3 化合物")

except FileNotFoundError:
    print("ログファイルが見つかりません")

print()
print("=" * 50)
