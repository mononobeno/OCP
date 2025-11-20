#!/usr/bin/env python3
"""
Extract GROMACS XVG data and convert to Plotly JSON format
"""
import sys
import json
import argparse

def parse_xvg(filepath):
    """Parse XVG file and return time, value arrays"""
    time_data = []
    value_data = []
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            # Skip comments and empty lines
            if line.startswith('#') or line.startswith('@') or not line:
                continue
            
            parts = line.split()
            if len(parts) >= 2:
                try:
                    # XVG time is already in ns
                    time_data.append(round(float(parts[0]), 3))
                    value_data.append(float(parts[1]))
                except ValueError:
                    continue
    
    return time_data, value_data

def main():
    parser = argparse.ArgumentParser(description='Convert XVG to JSON for Plotly')
    parser.add_argument('--input', required=True, help='Input XVG file')
    parser.add_argument('--output', required=True, help='Output JSON file')
    parser.add_argument('--name', default='Data', help='Trace name')
    
    args = parser.parse_args()
    
    time_data, value_data = parse_xvg(args.input)
    
    plotly_data = {
        "x": time_data,
        "y": value_data,
        "type": "scatter",
        "mode": "lines",
        "name": args.name
    }
    
    with open(args.output, 'w') as f:
        json.dump(plotly_data, f, indent=2)
    
    print(f"✅ Converted {len(time_data)} data points")
    print(f"   Time range: {time_data[0]:.2f} - {time_data[-1]:.2f} ns")
    print(f"   Value range: {min(value_data):.4f} - {max(value_data):.4f}")

if __name__ == '__main__':
    main()
