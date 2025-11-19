#!/bin/bash
set -e

echo "=========================================="
echo "Generating GitHub Pages for 50ns MD Results"
echo "=========================================="

DOCS_DIR="/home/dev/OCP/docs/pages/compounds"
DB_FILE="/home/dev/OCP/catalog/db/ocp_results.sqlite"

# Get compounds with 50ns MD completed
COMPOUNDS=$(sqlite3 "$DB_FILE" << 'SQL'
SELECT DISTINCT l.zinc_id
FROM vina_results v
JOIN ligands l ON v.ligand_id = l.id
WHERE v.md_status LIKE '%50ns%'
ORDER BY v.affinity ASC;
SQL
)

if [ -z "$COMPOUNDS" ]; then
  echo "No 50ns MD results found in database"
  exit 1
fi

echo "Found $(echo "$COMPOUNDS" | wc -l) compounds with 50ns MD results"
echo ""

# Function to convert XVG to JavaScript arrays
xvg_to_arrays() {
  local file=$1
  local time_var=$2
  local val_var=$3
  
  tail -n +18 "$file" | awk -v tv="$time_var" -v vv="$val_var" '
    BEGIN {t=""; v=""}
    {
      if (t != "") {t = t ","}
      if (v != "") {v = v ","}
      t = t $1
      v = v $2
    }
    END {
      print "var " tv " = [" t "];"
      print "var " vv " = [" v "];"
    }
  '
}

# Generate page for each compound
COMPOUND_NUM=0
for ZINC_ID in $COMPOUNDS; do
  COMPOUND_NUM=$((COMPOUND_NUM + 1))
  
  echo "Generating page $COMPOUND_NUM: $ZINC_ID"
  
  # Get compound data
  COMPOUND_DATA=$(sqlite3 "$DB_FILE" << SQL
SELECT 
  l.smiles,
  v.affinity,
  v.md_rmsd_avg,
  v.md_performance_nsday,
  v.md_simulation_time_ns
FROM vina_results v
JOIN ligands l ON v.ligand_id = l.id
WHERE l.zinc_id = '$ZINC_ID'
  AND v.md_status LIKE '%50ns%'
ORDER BY v.run_id DESC
LIMIT 1;
SQL
)
  
  IFS='|' read -r SMILES AFFINITY RMSD_AVG PERF SIM_TIME <<< "$COMPOUND_DATA"
  
  OUTPUT_DIR="/home/dev/OCP/results/md_output/${ZINC_ID}_50ns_md"
  PAGE_FILE="$DOCS_DIR/compound_50ns_${COMPOUND_NUM}.md"
  
  if [ ! -d "$OUTPUT_DIR" ]; then
    echo "  ⚠️  Output directory not found: $OUTPUT_DIR"
    continue
  fi
  
  cd "$OUTPUT_DIR"
  
  # Generate markdown page
  cat > "$PAGE_FILE" << MDPAGE
---
layout: compound
title: "50ns MD - $ZINC_ID"
---

# 50ns GPU MD Simulation Results

**Compound**: $ZINC_ID  
**Target**: Prothrombin (F2)  
**Vina Affinity**: $AFFINITY kcal/mol  
**Simulation Time**: ${SIM_TIME} ns  
**Performance**: $PERF ns/day (RTX 4070)  

---

## RMSD Evolution (Multi-Scale)

<div id="rmsd_plot" style="width:100%;height:500px;"></div>

<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
<script>
MDPAGE

  # Add RMSD data
  for rmsd_type in backbone calpha mainchain protein; do
    if [ -f "rmsd_${rmsd_type}.xvg" ]; then
      xvg_to_arrays "rmsd_${rmsd_type}.xvg" "time_${rmsd_type}" "rmsd_${rmsd_type}" >> "$PAGE_FILE"
    fi
  done

  cat >> "$PAGE_FILE" << 'MDPAGE'

var rmsd_traces = [
  {name: 'Backbone', x: time_backbone, y: rmsd_backbone, mode: 'lines', line: {color: 'rgb(31,119,180)', width: 2}},
  {name: 'C-alpha', x: time_calpha, y: rmsd_calpha, mode: 'lines', line: {color: 'rgb(255,127,14)', width: 2}},
  {name: 'MainChain', x: time_mainchain, y: rmsd_mainchain, mode: 'lines', line: {color: 'rgb(44,160,44)', width: 2}},
  {name: 'Protein', x: time_protein, y: rmsd_protein, mode: 'lines', line: {color: 'rgb(214,39,40)', width: 2}}
];

Plotly.newPlot('rmsd_plot', rmsd_traces, {
  title: 'RMSD Evolution (50ns Real MD)',
  xaxis: {title: 'Time (ns)'},
  yaxis: {title: 'RMSD (nm)'}
});
</script>

---

## RMSF (Per-Residue Fluctuation)

<div id="rmsf_plot" style="width:100%;height:400px;"></div>

<script>
MDPAGE

  if [ -f "rmsf.xvg" ]; then
    xvg_to_arrays "rmsf.xvg" "res_num" "rmsf_val" >> "$PAGE_FILE"
  fi

  cat >> "$PAGE_FILE" << 'MDPAGE'

Plotly.newPlot('rmsf_plot', [{
  x: res_num,
  y: rmsf_val,
  mode: 'lines',
  line: {color: 'rgb(148,103,189)', width: 2}
}], {
  title: 'Per-Residue Fluctuation',
  xaxis: {title: 'Residue Number'},
  yaxis: {title: 'RMSF (nm)'}
});
</script>

---

## Radius of Gyration

<div id="gyrate_plot" style="width:100%;height:400px;"></div>

<script>
MDPAGE

  if [ -f "gyrate.xvg" ]; then
    tail -n +18 gyrate.xvg | awk '
      BEGIN {t=""; rg=""}
      {
        if (t != "") {t = t ","}
        if (rg != "") {rg = rg ","}
        t = t $1
        rg = rg $2
      }
      END {
        print "var time_rg = [" t "];"
        print "var rg_val = [" rg "];"
      }
    ' >> "$PAGE_FILE"
  fi

  cat >> "$PAGE_FILE" << 'MDPAGE'

Plotly.newPlot('gyrate_plot', [{
  x: time_rg,
  y: rg_val,
  mode: 'lines',
  line: {color: 'rgb(140,86,75)', width: 2}
}], {
  title: 'Radius of Gyration',
  xaxis: {title: 'Time (ns)'},
  yaxis: {title: 'Rg (nm)'}
});
</script>

---

## Hydrogen Bonds

<div id="hbond_plot" style="width:100%;height:400px;"></div>

<script>
MDPAGE

  if [ -f "hbond.xvg" ]; then
    xvg_to_arrays "hbond.xvg" "time_hb" "hbond_num" >> "$PAGE_FILE"
  fi

  cat >> "$PAGE_FILE" << MDPAGE

Plotly.newPlot('hbond_plot', [{
  x: time_hb,
  y: hbond_num,
  mode: 'lines',
  line: {color: 'rgb(227,119,194)', width: 2}
}], {
  title: 'Hydrogen Bonds (Intra-protein)',
  xaxis: {title: 'Time (ns)'},
  yaxis: {title: 'Number of H-bonds'}
});
</script>

---

## Summary

| Metric | Value |
|--------|-------|
| **ZINC ID** | $ZINC_ID |
| **SMILES** | $SMILES |
| **Vina Affinity** | $AFFINITY kcal/mol |
| **RMSD (Backbone)** | $RMSD_AVG nm |
| **Performance** | $PERF ns/day |
| **Simulation Time** | ${SIM_TIME} ns |
| **Force Field** | OPLSAA |
| **Water Model** | SPC/E |
| **GPU** | NVIDIA RTX 4070 |

---

**This is real 50ns MD calculation data.**

[Back to Compound Index](index.html)
MDPAGE

  echo "  ✅ Generated: $PAGE_FILE"
done

# Update index page
cat > "$DOCS_DIR/index.md" << 'EOF'
---
layout: default
title: "MD Simulation Results"
---

# Molecular Dynamics Simulation Results

## 50ns MD Simulations (Real GPU Calculations)

EOF

COMPOUND_NUM=0
for ZINC_ID in $COMPOUNDS; do
  COMPOUND_NUM=$((COMPOUND_NUM + 1))
  
  COMPOUND_DATA=$(sqlite3 "$DB_FILE" << SQL
SELECT v.affinity, v.md_rmsd_avg
FROM vina_results v
JOIN ligands l ON v.ligand_id = l.id
WHERE l.zinc_id = '$ZINC_ID'
  AND v.md_status LIKE '%50ns%'
ORDER BY v.run_id DESC
LIMIT 1;
SQL
)
  
  IFS='|' read -r AFFINITY RMSD_AVG <<< "$COMPOUND_DATA"
  
  cat >> "$DOCS_DIR/index.md" << EOF
### [Compound ${COMPOUND_NUM} - $ZINC_ID](compound_50ns_${COMPOUND_NUM}.html)
- **Affinity**: $AFFINITY kcal/mol
- **RMSD**: $RMSD_AVG nm
- **Simulation**: 50 ns (Real GPU MD)

EOF
done

cat >> "$DOCS_DIR/index.md" << 'EOF'

---

## 10ns MD Simulations

### [Real 10ns MD - ZINC001241750201_1](real_md_compound_1.html)
- **RMSD**: 0.1018 nm
- **Performance**: 346.6 ns/day
- **Status**: ✅ First validation run

---

[Back to Main Index](../../index.html)
EOF

echo ""
echo "=========================================="
echo "✅ GitHub Pages Generation Complete!"
echo "=========================================="
echo "Generated pages: $COMPOUND_NUM"
echo "Index updated: $DOCS_DIR/index.md"
echo ""
echo "Pages will be available at:"
for i in $(seq 1 $COMPOUND_NUM); do
  echo "  - compound_50ns_${i}.html"
done
echo ""
