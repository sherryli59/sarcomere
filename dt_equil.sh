#!/bin/bash
cd /home/groups/rotskoff/sherry/sarcomere/sarcomere
EXE=cpp/build_release/sarcomere
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
OUT=dt_equil_out; mkdir -p $OUT
BOX='--n_actins=16 --n_myosins=4 --Lx=2.06 --Ly=0.190525589 --Lz=0.165 --periodic_mask=111 --actin_length=1 --myosin_radius=0.025 --actin_radius=0.015 --am_surface_gap=0.015'
CAL='--k_aa=3000 --am_tether_k=3000 --actin_steric_k=10000 --am_heads_per_half=141 --am_cycle_period=0.05 --temperature_kelvin=310 --thermal_noise_scale=1'
run () { local tag=$1 dt=$2 nsteps=$3 se=$4
  $EXE --initial_structure=shared_thick_thin_lattice --six_ring_prebind_aa $BOX $CAL     --dt=$dt --nsteps=$nsteps --save_every=$se --filename=$OUT/${tag}.h5 > $OUT/${tag}.log 2>&1
  echo "$tag exit=$? (dt=$dt nsteps=$nsteps -> physical=$(python3 -c "print($dt*$nsteps)")s)"
}
run ref_1e6  1e-6 40000 400
run test_1e5 1e-5 4000  40
echo DONE
