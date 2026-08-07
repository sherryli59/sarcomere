#!/bin/bash
#SBATCH --partition=rotskoff
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH --time=01:00:00
#SBATCH --job-name=am_cycle_sweep
#SBATCH --output=/scratch/users/xsli/sarcomere/cycle_sweep/sweep_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/cycle_sweep/sweep_%j.err

set -euo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH
export OMP_NUM_THREADS=1 OMP_DYNAMIC=FALSE HDF5_USE_FILE_LOCKING=FALSE

ROOT=/scratch/users/xsli/sarcomere/cycle_sweep
mkdir -p "${ROOT}"
EXE=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp/build_crossover_20260806/sarcomere

# Baseline AM cycle rates. Scaling all of them plus the stroke-execution speed
# by the same factor f is a pure time rescaling of the motor: every branching
# ratio, and therefore duty ratio and force per head, is unchanged.
K_ON=100; K_WEAK_OFF=200; K_PS=200; K_STRONG_OFF=160; V0=5

run_one () {
  f=$1
  kon=$(python3 -c "print($K_ON*$f)")
  kwo=$(python3 -c "print($K_WEAK_OFF*$f)")
  kps=$(python3 -c "print($K_PS*$f)")
  kso=$(python3 -c "print($K_STRONG_OFF*$f)")
  v=$(python3 -c "print($V0*$f)")
  for PT in 0 1; do
    extra=""
    [ "$PT" = "1" ] && extra="--am_crossover_pull_through=true"
    ${EXE} \
      --n_actins=1 --n_myosins=1 --n_fixed_myosins=1 \
      --Lx=10 --Ly=2 --Lz=2 --dimension=3 --pbc_mask=111 \
      --initial_structure=single_am \
      --actin_length=1 --myosin_length=1.5 \
      --actin_radius=0.015 --myosin_radius=0.025 \
      --pre_phase_steps=0 --interaction_ramp_steps=0 \
      --physical_mobility=true --temperature_kelvin=310 \
      --thermal_noise_scale=0 --solvent_viscosity=0.001 \
      --sterics=true --actin_steric_k=2000 --myosin_steric_k=1000 --am_steric_k=1000 \
      --am_cutoff=0.10 --am_optimal=0.055 --am_surface_gap=0.015 \
      --k_am=0 --kappa_am=0 \
      --dynamic_am=true --am_spatial_heads=true --am_heads_per_half=9 \
      --am_sites_per_actin=182 --am_max_heads_per_actin=24 \
      --am_k_on=${kon} --am_k_weak_off=${kwo} \
      --am_k_powerstroke=${kps} --am_k_strong_off=${kso} \
      --am_unloaded_velocity=${v} \
      --am_motor_k=1000 --am_power_stroke=0.005 --am_stall_force=4.5 \
      --am_tether_k=3000 --am_alignment_k=120 \
      --am_thin_activation=true --am_initialize_q_calcium_equilibrium=true \
      --directional=true --k_aa=1000 --kappa_aa=100 \
      --dt=5e-7 --nsteps=2000000 --save_every=10000 --seed=3 \
      ${extra} \
      --filename=${ROOT}/f${f}_pt${PT}.h5 > ${ROOT}/f${f}_pt${PT}.log 2>&1 &
  done
}

for f in 1 2 4 8 16; do run_one $f; done
wait
echo "sweep complete"
ls -la ${ROOT}/*.h5
