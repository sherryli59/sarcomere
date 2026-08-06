#!/bin/bash
cd /home/groups/rotskoff/sherry/sarcomere/sarcomere
EXE=cpp/build_release/sarcomere
PY=/home/groups/rotskoff/sherry/mamba_root/envs/torch_env/bin/python
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
OUT=/home/groups/rotskoff/sherry/sarcomere/sarcomere/bondsweep_out; mkdir -p $OUT
run () {
  local tag=$1 kaa=$2 tether=$3 steric=$4
  local f=$OUT/${tag}.h5
  $EXE --initial_structure=shared_thick_thin_lattice --six_ring_prebind_aa     --dt=1e-5 --nsteps=2000 --save_every=100 --temperature_kelvin=310 --thermal_noise_scale=1     --n_actins=16 --n_myosins=4 --Lx=2.06 --Ly=0.190525589 --Lz=0.165 --periodic_mask=111 --actin_length=1 --myosin_radius=0.025 --actin_radius=0.015 --am_surface_gap=0.015 --am_heads_per_half=141 --am_cycle_period=0.05     --k_aa=$kaa --am_tether_k=$tether --actin_steric_k=$steric     --filename=$f > $OUT/${tag}.log 2>&1
  if [ $? -ne 0 ]; then echo "$tag RUN_FAIL"; return; fi
  $PY analysis/analyze_shared_lattice.py $f --output $OUT/${tag}.csv > /dev/null 2>&1
  if [ $? -ne 0 ]; then echo "$tag ANALYZE_FAIL"; return; fi
  $PY -c "import pandas as pd;d=pd.read_csv('$OUT/${tag}.csv');print('%-14s kaa=%-6s tether=%-6s steric=%-6s | caps=%-8d psi6=%.3f occ=%.2f'%('$tag','$kaa','$tether','$steric',int(d['cumulative_cap_events'].iloc[-1]),d['mean_psi6'].iloc[-1],d['expected_aa_edge_occupancy'].iloc[-1]))"
}
echo '=== dt=1e-5 : soften BONDS, hold STERIC stiff (10000) ==='
run bond_hard   3000  3000  10000
run bond_mid    300   300   10000
run bond_soft   30    30    10000
run bond_vsoft  3     3     10000
echo '=== controls ==='
run all_soft    3     3     100
run steric_only 3000  3000  100
echo DONE
