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

### Multi-Scale RMSD Analysis

タンパク質の動きを複数のスケールで評価し、結合ポケット領域の安定性を確認します。

<div id="multi-rmsd-plot" style="height:500px;"></div>

<script>
PROD_SECTION

    # 複数のRMSDデータを収集
    echo "const rmsdTraces = [];" >> "${COMPOUND_PAGE}"
    
    # Backbone RMSD
    if [[ -f "${mddir}/rmsd_backbone.xvg" ]]; then
        echo "rmsdTraces.push({x: [" >> "${COMPOUND_PAGE}"
        awk '/^[^@#]/ {printf "%s,", $1}' "${mddir}/rmsd_backbone.xvg" >> "${COMPOUND_PAGE}"
        echo "], y: [" >> "${COMPOUND_PAGE}"
        awk '/^[^@#]/ {printf "%s,", $2}' "${mddir}/rmsd_backbone.xvg" >> "${COMPOUND_PAGE}"
        echo "], type: 'scatter', mode: 'lines', name: 'Backbone', line: {color: '#2E86AB', width: 2}});" >> "${COMPOUND_PAGE}"
    fi
    
    # C-alpha RMSD
    if [[ -f "${mddir}/rmsd_calpha.xvg" ]]; then
        echo "rmsdTraces.push({x: [" >> "${COMPOUND_PAGE}"
        awk '/^[^@#]/ {printf "%s,", $1}' "${mddir}/rmsd_calpha.xvg" >> "${COMPOUND_PAGE}"
        echo "], y: [" >> "${COMPOUND_PAGE}"
        awk '/^[^@#]/ {printf "%s,", $2}' "${mddir}/rmsd_calpha.xvg" >> "${COMPOUND_PAGE}"
        echo "], type: 'scatter', mode: 'lines', name: 'C-alpha', line: {color: '#00A878', width: 2}});" >> "${COMPOUND_PAGE}"
    fi
    
    # MainChain RMSD
    if [[ -f "${mddir}/rmsd_mainchain.xvg" ]]; then
        echo "rmsdTraces.push({x: [" >> "${COMPOUND_PAGE}"
        awk '/^[^@#]/ {printf "%s,", $1}' "${mddir}/rmsd_mainchain.xvg" >> "${COMPOUND_PAGE}"
        echo "], y: [" >> "${COMPOUND_PAGE}"
        awk '/^[^@#]/ {printf "%s,", $2}' "${mddir}/rmsd_mainchain.xvg" >> "${COMPOUND_PAGE}"
        echo "], type: 'scatter', mode: 'lines', name: 'MainChain', line: {color: '#F77F00', width: 2}});" >> "${COMPOUND_PAGE}"
    fi
    
    # Pocket region RMSD
    if [[ -f "${mddir}/rmsd_pocket.xvg" ]]; then
        echo "rmsdTraces.push({x: [" >> "${COMPOUND_PAGE}"
        awk '/^[^@#]/ {printf "%s,", $1}' "${mddir}/rmsd_pocket.xvg" >> "${COMPOUND_PAGE}"
        echo "], y: [" >> "${COMPOUND_PAGE}"
        awk '/^[^@#]/ {printf "%s,", $2}' "${mddir}/rmsd_pocket.xvg" >> "${COMPOUND_PAGE}"
        echo "], type: 'scatter', mode: 'lines', name: 'Pocket Region', line: {color: '#C1121F', width: 3, dash: 'dot'}});" >> "${COMPOUND_PAGE}"
    fi
    
    echo "Plotly.newPlot('multi-rmsd-plot', rmsdTraces, {title: 'RMSD Multi-Scale Analysis', xaxis: {title: 'Time (ns)'}, yaxis: {title: 'RMSD (nm)'}, hovermode: 'closest'});" >> "${COMPOUND_PAGE}"

    cat >> "${COMPOUND_PAGE}" <<'POCKET_SECTION'
</script>

### Residue Flexibility (RMSF)

各残基の揺らぎを示すRMSF。結合ポケット領域の柔軟性を評価できます。

<div id="rmsf-plot" style="height:400px;"></div>

<script>
POCKET_SECTION

    if [[ -f "${mddir}/rmsf.xvg" ]]; then
        echo "const rmsfData = {x: [" >> "${COMPOUND_PAGE}"
        awk '/^[^@#]/ {printf "%s,", $1}' "${mddir}/rmsf.xvg" >> "${COMPOUND_PAGE}"
        echo "], y: [" >> "${COMPOUND_PAGE}"
        awk '/^[^@#]/ {printf "%s,", $2}' "${mddir}/rmsf.xvg" >> "${COMPOUND_PAGE}"
        echo "]};" >> "${COMPOUND_PAGE}"
        echo "Plotly.newPlot('rmsf-plot', [{x: rmsfData.x, y: rmsfData.y, type: 'scatter', mode: 'lines', name: 'RMSF', line: {color: '#9D4EDD', width: 2}, fill: 'tozeroy', fillcolor: 'rgba(157,78,221,0.2)'}], {title: 'Residue Flexibility (RMSF)', xaxis: {title: 'Residue Number'}, yaxis: {title: 'RMSF (nm)'}});" >> "${COMPOUND_PAGE}"
    fi

    cat >> "${COMPOUND_PAGE}" <<'HBOND_SECTION'
</script>

### Radius of Gyration

タンパク質のコンパクトさを示す回転半径。

<div id="gyrate-plot" style="height:400px;"></div>

<script>
HBOND_SECTION

    if [[ -f "${mddir}/gyrate.xvg" ]]; then
        echo "const gyrateData = {x: [" >> "${COMPOUND_PAGE}"
        awk '/^[^@#]/ {printf "%s,", $1}' "${mddir}/gyrate.xvg" >> "${COMPOUND_PAGE}"
        echo "], y: [" >> "${COMPOUND_PAGE}"
        awk '/^[^@#]/ {printf "%s,", $2}' "${mddir}/gyrate.xvg" >> "${COMPOUND_PAGE}"
        echo "]};" >> "${COMPOUND_PAGE}"
        echo "Plotly.newPlot('gyrate-plot', [{x: gyrateData.x, y: gyrateData.y, type: 'scatter', mode: 'lines', name: 'Rg', line: {color: '#06A77D', width: 2}}], {title: 'Protein Compactness (Radius of Gyration)', xaxis: {title: 'Time (ns)'}, yaxis: {title: 'Rg (nm)'}});" >> "${COMPOUND_PAGE}"
    fi

    cat >> "${COMPOUND_PAGE}" <<'HBOND_FINAL_SECTION'
</script>

### Hydrogen Bonds (Protein Stability)

<div id="hbond-plot" style="height:400px;"></div>

<script>
HBOND_FINAL_SECTION

    # H-bonds
    if [[ -f "${mddir}/hbond.xvg" ]]; then
        echo "const hbondData = {x: [" >> "${COMPOUND_PAGE}"
        awk '/^[^@#]/ {printf "%s,", $1}' "${mddir}/hbond.xvg" >> "${COMPOUND_PAGE}"
        echo "], y: [" >> "${COMPOUND_PAGE}"
        awk '/^[^@#]/ {printf "%s,", $2}' "${mddir}/hbond.xvg" >> "${COMPOUND_PAGE}"
        echo "]};" >> "${COMPOUND_PAGE}"
        echo "Plotly.newPlot('hbond-plot', [{x: hbondData.x, y: hbondData.y, type: 'scatter', mode: 'lines', name: 'H-bonds', line: {color: '#FF6B35', width: 2}}], {title: 'Hydrogen Bonds', xaxis: {title: 'Time (ns)'}, yaxis: {title: 'Count'}});" >> "${COMPOUND_PAGE}"
    fi

    echo "</script>" >> "${COMPOUND_PAGE}"
    echo "" >> "${COMPOUND_PAGE}"
    echo "_Generated: $(date '+%Y-%m-%d %H:%M:%S')_" >> "${COMPOUND_PAGE}"

    echo "  ✓ ${zid} (with equilibration graphs)"
done

echo "✅ Pages with equilibration graphs generated"
