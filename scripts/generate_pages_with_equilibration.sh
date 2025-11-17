#!/bin/bash
# 平衡化グラフを含む完全版ページ生成
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DB="${ROOT}/catalog/db/ocp_results.sqlite"
PAGES_DIR="${ROOT}/docs/pages"

echo "Generating pages with equilibration graphs..."

mkdir -p "${PAGES_DIR}/compounds"

# 化合物詳細ページ更新（平衡化グラフ追加）
sqlite3 -separator '|' "${DB}" "SELECT v.rowid, l.zinc_id, l.smiles, v.affinity_kcal, v.md_output_dir, v.md_rmsd_avg, v.md_simulation_time_ns, v.md_performance_nsday, t.pdb_id, t.name FROM vina_results v JOIN ligands l ON v.ligand_id=l.id JOIN runs r ON v.run_id=r.id JOIN targets t ON r.target_id=t.id WHERE v.md_status='completed';" | while IFS='|' read -r vrid zid smi aff mddir rmsd simtime perf pdb tname; do
    
    COMPOUND_PAGE="${PAGES_DIR}/compounds/compound_${vrid}.md"
    
    # ヘッダー部分（既存と同じ）
    cat > "${COMPOUND_PAGE}" <<HEADER
---
layout: default
title: ${zid}
---

# ${zid} - MD Analysis

[← Back to ${tname}](../targets/target_$(sqlite3 "${DB}" "SELECT r.target_id FROM vina_results v JOIN runs r ON v.run_id=r.id WHERE v.rowid=${vrid};").html)

## Compound Information

- **ZINC ID**: ${zid}
- **SMILES**: \`${smi}\`
- **Target**: ${tname} (${pdb})
- **Vina Affinity**: **${aff} kcal/mol**

## MD Simulation Summary

- **Status**: ✅ Completed
- **Simulation Time**: ${simtime} ns
- **RMSD (Backbone)**: **${rmsd} nm**
- **Performance**: ${perf} ns/day on RTX 4070

<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>

---

## Equilibration Quality Assessment

### 1. Energy Minimization

Potential energy should converge to a stable minimum.

<div id="em-plot"></div>

<script>
HEADER

    # EM data
    if [[ -f "${mddir}/energy_em_potential.xvg" ]]; then
        echo "const emData = {x: [" >> "${COMPOUND_PAGE}"
        awk '/^[^@#]/ {printf "%s,", $1}' "${mddir}/energy_em_potential.xvg" >> "${COMPOUND_PAGE}"
        echo "], y: [" >> "${COMPOUND_PAGE}"
        awk '/^[^@#]/ {printf "%s,", $2}' "${mddir}/energy_em_potential.xvg" >> "${COMPOUND_PAGE}"
        echo "]};" >> "${COMPOUND_PAGE}"
        echo "Plotly.newPlot('em-plot', [{x: emData.x, y: emData.y, type: 'scatter', mode: 'lines', name: 'Potential Energy', line: {color: 'red'}}], {title: 'Energy Minimization', xaxis: {title: 'Step'}, yaxis: {title: 'Potential (kJ/mol)'}});" >> "${COMPOUND_PAGE}"
    fi

    cat >> "${COMPOUND_PAGE}" <<'NVT_SECTION'
</script>

### 2. NVT Equilibration (Temperature)

Temperature should stabilize around 300K.

<div id="nvt-temp-plot"></div>

<script>
NVT_SECTION

    # NVT temperature
    if [[ -f "${mddir}/energy_nvt_temp.xvg" ]]; then
        echo "const nvtData = {x: [" >> "${COMPOUND_PAGE}"
        awk '/^[^@#]/ {printf "%s,", $1}' "${mddir}/energy_nvt_temp.xvg" >> "${COMPOUND_PAGE}"
        echo "], y: [" >> "${COMPOUND_PAGE}"
        awk '/^[^@#]/ {printf "%s,", $2}' "${mddir}/energy_nvt_temp.xvg" >> "${COMPOUND_PAGE}"
        echo "]};" >> "${COMPOUND_PAGE}"
        echo "Plotly.newPlot('nvt-temp-plot', [{x: nvtData.x, y: nvtData.y, type: 'scatter', mode: 'lines', name: 'Temperature', line: {color: 'green'}}], {title: 'NVT Temperature', xaxis: {title: 'Time (ps)'}, yaxis: {title: 'Temperature (K)'}, shapes: [{type: 'line', x0: 0, x1: 100, y0: 300, y1: 300, line: {color: 'red', dash: 'dash'}}]});" >> "${COMPOUND_PAGE}"
    fi

    cat >> "${COMPOUND_PAGE}" <<'NPT_SECTION'
</script>

### 3. NPT Equilibration (Pressure & Density)

Pressure should fluctuate around 1 bar, density should stabilize.

<div id="npt-pressure-plot"></div>
<div id="npt-density-plot"></div>

<script>
NPT_SECTION

    # NPT pressure
    if [[ -f "${mddir}/energy_npt_pressure.xvg" ]]; then
        echo "const nptPressData = {x: [" >> "${COMPOUND_PAGE}"
        awk '/^[^@#]/ {printf "%s,", $1}' "${mddir}/energy_npt_pressure.xvg" >> "${COMPOUND_PAGE}"
        echo "], y: [" >> "${COMPOUND_PAGE}"
        awk '/^[^@#]/ {printf "%s,", $2}' "${mddir}/energy_npt_pressure.xvg" >> "${COMPOUND_PAGE}"
        echo "]};" >> "${COMPOUND_PAGE}"
        echo "Plotly.newPlot('npt-pressure-plot', [{x: nptPressData.x, y: nptPressData.y, type: 'scatter', mode: 'lines', name: 'Pressure', line: {color: 'blue'}}], {title: 'NPT Pressure', xaxis: {title: 'Time (ps)'}, yaxis: {title: 'Pressure (bar)'}});" >> "${COMPOUND_PAGE}"
    fi

    # NPT density
    if [[ -f "${mddir}/energy_npt_density.xvg" ]]; then
        echo "const nptDensData = {x: [" >> "${COMPOUND_PAGE}"
        awk '/^[^@#]/ {printf "%s,", $1}' "${mddir}/energy_npt_density.xvg" >> "${COMPOUND_PAGE}"
        echo "], y: [" >> "${COMPOUND_PAGE}"
        awk '/^[^@#]/ {printf "%s,", $2}' "${mddir}/energy_npt_density.xvg" >> "${COMPOUND_PAGE}"
        echo "]};" >> "${COMPOUND_PAGE}"
        echo "Plotly.newPlot('npt-density-plot', [{x: nptDensData.x, y: nptDensData.y, type: 'scatter', mode: 'lines', name: 'Density', line: {color: 'purple'}}], {title: 'NPT Density', xaxis: {title: 'Time (ps)'}, yaxis: {title: 'Density (kg/m³)'}});" >> "${COMPOUND_PAGE}"
    fi

    echo "</script>" >> "${COMPOUND_PAGE}"

    # Production MDセクション
    cat >> "${COMPOUND_PAGE}" <<'PROD_SECTION'

---

## Production MD Analysis

### RMSD (Backbone Stability)

<div id="rmsd-plot"></div>

<script>
PROD_SECTION

    # RMSD
    if [[ -f "${mddir}/rmsd_backbone.xvg" ]]; then
        echo "const rmsdData = {x: [" >> "${COMPOUND_PAGE}"
        awk '/^[^@#]/ {printf "%s,", $1}' "${mddir}/rmsd_backbone.xvg" >> "${COMPOUND_PAGE}"
        echo "], y: [" >> "${COMPOUND_PAGE}"
        awk '/^[^@#]/ {printf "%s,", $2}' "${mddir}/rmsd_backbone.xvg" >> "${COMPOUND_PAGE}"
        echo "]};" >> "${COMPOUND_PAGE}"
        echo "Plotly.newPlot('rmsd-plot', [{x: rmsdData.x, y: rmsdData.y, type: 'scatter', mode: 'lines', name: 'RMSD'}], {title: 'Backbone RMSD', xaxis: {title: 'Time (ns)'}, yaxis: {title: 'RMSD (nm)'}});" >> "${COMPOUND_PAGE}"
    fi

    cat >> "${COMPOUND_PAGE}" <<'HBOND_SECTION'
</script>

### Hydrogen Bonds (Protein-Protein)

<div id="hbond-plot"></div>

<script>
HBOND_SECTION

    # H-bonds
    if [[ -f "${mddir}/hbond.xvg" ]]; then
        echo "const hbondData = {x: [" >> "${COMPOUND_PAGE}"
        awk '/^[^@#]/ {printf "%s,", $1}' "${mddir}/hbond.xvg" >> "${COMPOUND_PAGE}"
        echo "], y: [" >> "${COMPOUND_PAGE}"
        awk '/^[^@#]/ {printf "%s,", $2}' "${mddir}/hbond.xvg" >> "${COMPOUND_PAGE}"
        echo "]};" >> "${COMPOUND_PAGE}"
        echo "Plotly.newPlot('hbond-plot', [{x: hbondData.x, y: hbondData.y, type: 'scatter', mode: 'lines', name: 'H-bonds', line: {color: 'orange'}}], {title: 'Hydrogen Bonds', xaxis: {title: 'Time (ns)'}, yaxis: {title: 'Count'}});" >> "${COMPOUND_PAGE}"
    fi

    echo "</script>" >> "${COMPOUND_PAGE}"
    echo "" >> "${COMPOUND_PAGE}"
    echo "_Generated: $(date '+%Y-%m-%d %H:%M:%S')_" >> "${COMPOUND_PAGE}"

    echo "  ✓ ${zid} (with equilibration graphs)"
done

echo "✅ Pages with equilibration graphs generated"
