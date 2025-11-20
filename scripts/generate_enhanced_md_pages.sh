#!/bin/bash
# Generate enhanced pages with real plot data for 5 compounds

set -e

COMPOUNDS=(
  "ZINC001241740817_1"
  "ZINC001241740915_1"
  "ZINC001241740923_1"
  "ZINC001241740924_1"
  "ZINC001241740926_1"
)

OUTPUT_DIR="docs/pages/compounds"
IMAGE_DIR="docs/images/compounds"

mkdir -p "$OUTPUT_DIR" "$IMAGE_DIR"

echo "📊 Generating enhanced MD result pages..."
echo ""

for COMPOUND in "${COMPOUNDS[@]}"; do
  MD_DIR="results/md_output/${COMPOUND}_10ns"
  
  if [ ! -d "$MD_DIR" ]; then
    echo "❌ [$COMPOUND] MD directory not found"
    continue
  fi
  
  echo "📄 [$COMPOUND] Processing..."
  
  # Extract data from XVG files to JSON
  python3 tools/extract_xvg_to_json.py \
    --input "$MD_DIR/rmsd_backbone.xvg" \
    --output "/tmp/${COMPOUND}_rmsd_bb.json" \
    --name "Backbone RMSD" 2>&1 | grep -E "✅|⚠️" || true
  
  python3 tools/extract_xvg_to_json.py \
    --input "$MD_DIR/rmsd_calpha.xvg" \
    --output "/tmp/${COMPOUND}_rmsd_ca.json" \
    --name "C-alpha RMSD" 2>&1 | grep -E "✅|⚠️" || true
  
  python3 tools/extract_xvg_to_json.py \
    --input "$MD_DIR/gyrate.xvg" \
    --output "/tmp/${COMPOUND}_gyrate.json" \
    --name "Radius of Gyration" 2>&1 | grep -E "✅|⚠️" || true
  
  python3 tools/extract_xvg_to_json.py \
    --input "$MD_DIR/hbond.xvg" \
    --output "/tmp/${COMPOUND}_hbond.json" \
    --name "Hydrogen Bonds" 2>&1 | grep -E "✅|⚠️" || true
  
  # Read JSON data
  RMSD_BB_DATA=$(cat "/tmp/${COMPOUND}_rmsd_bb.json")
  RMSD_CA_DATA=$(cat "/tmp/${COMPOUND}_rmsd_ca.json")
  GYRATE_DATA=$(cat "/tmp/${COMPOUND}_gyrate.json")
  HBOND_DATA=$(cat "/tmp/${COMPOUND}_hbond.json")
  
  # Get final values
  FINAL_RMSD_BB=$(awk 'END {print $2}' "$MD_DIR/rmsd_backbone.xvg")
  FINAL_RMSD_CA=$(awk 'END {print $2}' "$MD_DIR/rmsd_calpha.xvg")
  FINAL_RG=$(awk 'END {print $2}' "$MD_DIR/gyrate.xvg")
  FINAL_HBOND=$(awk 'END {print $2}' "$MD_DIR/hbond.xvg")
  
  # Calculate averages
  AVG_RMSD_BB=$(awk '!/^[@#]/ && NF>=2 {sum+=$2; n++} END {print sum/n}' "$MD_DIR/rmsd_backbone.xvg")
  AVG_RG=$(awk '!/^[@#]/ && NF>=2 {sum+=$2; n++} END {print sum/n}' "$MD_DIR/gyrate.xvg")
  AVG_HBOND=$(awk '!/^[@#]/ && NF>=2 {sum+=$2; n++} END {print sum/n}' "$MD_DIR/hbond.xvg")
  
  # Generate Markdown page
  cat > "$OUTPUT_DIR/${COMPOUND}_10ns.md" << 'EOFPAGE'
---
layout: default
title: COMPOUND_ID - 10ns MD
parent: Compounds
nav_order: 10
---

# COMPOUND_ID - 10ns MD Simulation

## Basic Information
- **Compound ID**: COMPOUND_ID
- **Receptor**: Prothrombin (6C2W)
- **Simulation Time**: 10 ns
- **Status**: ✅ Completed

## Analysis Summary

| Metric | Average | Final Value |
|--------|---------|-------------|
| Backbone RMSD (nm) | AVG_RMSD_BB | FINAL_RMSD_BB |
| C-alpha RMSD (nm) | AVG_RMSD_CA | FINAL_RMSD_CA |
| Radius of Gyration (nm) | AVG_RG | FINAL_RG |
| Hydrogen Bonds | AVG_HBOND | FINAL_HBOND |

## Visualizations

### RMSD Analysis

<div id="rmsd-plot"></div>

<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
<script>
var data = [
  RMSD_BB_DATA_JSON,
  RMSD_CA_DATA_JSON
];

var layout = {
  title: 'RMSD vs Time',
  xaxis: {title: 'Time (ns)'},
  yaxis: {title: 'RMSD (nm)'},
  hovermode: 'closest'
};

Plotly.newPlot('rmsd-plot', data, layout);
</script>

### Radius of Gyration

<div id="gyrate-plot"></div>

<script>
var data = [GYRATE_DATA_JSON];

var layout = {
  title: 'Radius of Gyration vs Time',
  xaxis: {title: 'Time (ns)'},
  yaxis: {title: 'Rg (nm)'},
  hovermode: 'closest'
};

Plotly.newPlot('gyrate-plot', data, layout);
</script>

### Hydrogen Bonds

<div id="hbond-plot"></div>

<script>
var data = [HBOND_DATA_JSON];

var layout = {
  title: 'Hydrogen Bonds vs Time',
  xaxis: {title: 'Time (ns)'},
  yaxis: {title: 'Number of H-bonds'},
  hovermode: 'closest'
};

Plotly.newPlot('hbond-plot', data, layout);
</script>

## Files
- Trajectory: `results/md_output/COMPOUND_ID_10ns/md.xtc`
- Final Structure: `results/md_output/COMPOUND_ID_10ns/complex_final.pdb`
- RMSD Data: `results/md_output/COMPOUND_ID_10ns/rmsd_backbone.xvg`
- Gyration Data: `results/md_output/COMPOUND_ID_10ns/gyrate.xvg`
- H-bond Data: `results/md_output/COMPOUND_ID_10ns/hbond.xvg`

EOFPAGE

  # Replace placeholders
  sed -i "s/COMPOUND_ID/$COMPOUND/g" "$OUTPUT_DIR/${COMPOUND}_10ns.md"
  sed -i "s/AVG_RMSD_BB/$AVG_RMSD_BB/g" "$OUTPUT_DIR/${COMPOUND}_10ns.md"
  sed -i "s/FINAL_RMSD_BB/$FINAL_RMSD_BB/g" "$OUTPUT_DIR/${COMPOUND}_10ns.md"
  sed -i "s/AVG_RMSD_CA/$AVG_RMSD_CA/g" "$OUTPUT_DIR/${COMPOUND}_10ns.md"
  sed -i "s/FINAL_RMSD_CA/$FINAL_RMSD_CA/g" "$OUTPUT_DIR/${COMPOUND}_10ns.md"
  sed -i "s/AVG_RG/$AVG_RG/g" "$OUTPUT_DIR/${COMPOUND}_10ns.md"
  sed -i "s/FINAL_RG/$FINAL_RG/g" "$OUTPUT_DIR/${COMPOUND}_10ns.md"
  sed -i "s/AVG_HBOND/$AVG_HBOND/g" "$OUTPUT_DIR/${COMPOUND}_10ns.md"
  sed -i "s/FINAL_HBOND/$FINAL_HBOND/g" "$OUTPUT_DIR/${COMPOUND}_10ns.md"
  
  # Embed JSON data (escape for sed)
  python3 -c "
import json, sys
with open('/tmp/${COMPOUND}_rmsd_bb.json') as f:
    data = json.load(f)
    data['name'] = 'Backbone RMSD'
    data['line'] = {'color': 'blue'}
    print(json.dumps(data))
" > /tmp/${COMPOUND}_rmsd_bb_styled.json
  
  python3 -c "
import json, sys
with open('/tmp/${COMPOUND}_rmsd_ca.json') as f:
    data = json.load(f)
    data['name'] = 'C-alpha RMSD'
    data['line'] = {'color': 'orange', 'dash': 'dash'}
    print(json.dumps(data))
" > /tmp/${COMPOUND}_rmsd_ca_styled.json
  
  # Use Python to properly embed JSON (avoiding sed escaping issues)
  python3 << EOFPYTHON
import re

with open('$OUTPUT_DIR/${COMPOUND}_10ns.md', 'r') as f:
    content = f.read()

# Read JSON data
with open('/tmp/${COMPOUND}_rmsd_bb_styled.json') as f:
    rmsd_bb = f.read().strip()
with open('/tmp/${COMPOUND}_rmsd_ca_styled.json') as f:
    rmsd_ca = f.read().strip()
with open('/tmp/${COMPOUND}_gyrate.json') as f:
    gyrate = f.read().strip()
with open('/tmp/${COMPOUND}_hbond.json') as f:
    hbond = f.read().strip()

# Replace placeholders
content = content.replace('RMSD_BB_DATA_JSON', rmsd_bb)
content = content.replace('RMSD_CA_DATA_JSON', rmsd_ca)
content = content.replace('GYRATE_DATA_JSON', gyrate)
content = content.replace('HBOND_DATA_JSON', hbond)

with open('$OUTPUT_DIR/${COMPOUND}_10ns.md', 'w') as f:
    f.write(content)
EOFPYTHON

  echo "   ✅ Page generated: ${COMPOUND}_10ns.md"
done

echo ""
echo "✅ All pages generated with interactive plots"
echo "📂 Output directory: $OUTPUT_DIR"
