#!/usr/bin/env bash
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
energy_src="$GROUP_HOME/sherry/sarcomere/sarcomere/analysis/visualizations/gaphunt-rig/status2_gate_20261028"
energy_out="$SCRATCH/sarcomere/interactive_energy_20261028"
energy_durable="$energy_src/results/interactive_energy"
"$GROUP_HOME/sherry/mamba_root/envs/sarcomere/bin/python" "$energy_src/build_interactive_energy.py" --output "$energy_out" --plotly "$energy_out/plotly-3.1.0.min.js"
node "$energy_src/check_interactive_energy.cjs" "$energy_out/connector_energy_3d.html" "$energy_out/js_validation.json"
mkdir -p "$energy_durable"
cp "$energy_out/connector_energy_3d.html" "$energy_out/connector_energy_3d_preview.png" "$energy_out/connector_energy_examples.json" "$energy_out/validation.json" "$energy_out/js_validation.json" "$energy_durable/"
sha256sum "$energy_out/plotly-3.1.0.min.js" "$energy_out/connector_energy_3d.html" > "$energy_durable/SHA256SUMS"
ls -lh "$energy_out/connector_energy_3d.html"
