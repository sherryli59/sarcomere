#!/bin/bash
# Locate the unstable torque mode.
#  arms 0-3: soft springs, dt bisection -> where do torque caps first appear
#  arms 4-8: at the failing dt, soften ONE torque term at a time -> attribution
# All arms: 20 ms of physical time, A160/M80, seed 20260806.
#SBATCH --partition=owners
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8GB
#SBATCH --time=01:00:00
#SBATCH --job-name=kappa_ladder
#SBATCH --no-requeue
#SBATCH --array=0-4
#SBATCH --output=/scratch/users/xsli/sarcomere/kappa_ladder/logs/%A_%a.out
#SBATCH --error=/scratch/users/xsli/sarcomere/kappa_ladder/logs/%A_%a.err

set -uo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export OMP_DYNAMIC=FALSE HDF5_USE_FILE_LOCKING=FALSE

ROOT=/scratch/users/xsli/sarcomere/kappa_ladder
mkdir -p "${ROOT}"/{data,logs}
EXE=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp/build_crossover_20260806/sarcomere

# torque terms, defaults
KAPPA_AA=100; ALIGN=120; PERP=50; HEADT=100

case ${SLURM_ARRAY_TASK_ID} in
  0) TAG=k1_dt2e-6;  DT=2e-6;  N=10000; SAVE=500; KAPPA_AA=1 ;;
  1) TAG=k1_dt5e-6;  DT=5e-6;  N=4000;  SAVE=200; KAPPA_AA=1 ;;
  2) TAG=k1_dt1e-5;  DT=1e-5;  N=2000;  SAVE=100; KAPPA_AA=1 ;;
  3) TAG=k1_dt2e-5;  DT=2e-5;  N=1000;  SAVE=50;  KAPPA_AA=1 ;;
  4) TAG=k10_dt1e-5; DT=1e-5;  N=2000;  SAVE=100; KAPPA_AA=10 ;;
esac

srun --kill-on-bad-exit=1 "${EXE}" \
  --n_actins=160 --n_myosins=80 --Lx=3.86 --Ly=1 --Lz=1 \
  --dimension=3 --pbc_mask=111 \
  --initial_structure=myosin_x_noise \
  --actin_direction_noise=-1 --myosin_direction_noise=1 \
  --actin_length=1 --myosin_length=1.5 \
  --actin_radius=0.015 --myosin_radius=0.025 \
  --pre_phase_steps=0 --interaction_ramp_steps=0 \
  --physical_mobility=true --temperature_kelvin=310 \
  --thermal_noise_scale=1 --solvent_viscosity=0.001 \
  --max_actin_displacement=0.005 --max_myosin_displacement=0.005 \
  --max_actin_rotation=0.005 --max_myosin_rotation=0.005 \
  --sterics=true \
  --am_tether_k=285 --actin_steric_k=285 --am_steric_k=285 \
  --myosin_steric_k=285 --am_motor_k=285 --k_aa=285 \
  --am_power_stroke=0.0175 \
  --am_cutoff=0.10 --am_optimal=0.055 --am_surface_gap=0.015 \
  --k_am=0 --kappa_am=0 \
  --dynamic_am=true --am_offcenter_torque=true --am_spatial_heads=true \
  --am_heads_per_half=9 --am_sites_per_actin=182 --am_max_heads_per_actin=24 \
  --am_k_on=400 --am_k_weak_off=800 --am_k_powerstroke=800 --am_k_strong_off=640 \
  --am_unloaded_velocity=20 --am_stall_force=4.5 \
  --am_alignment_k=${ALIGN} --am_head_transverse_k=${HEADT} \
  --am_strain_dependent_detachment=true --am_strain_detachment_theta=0.1 \
  --am_thin_activation=true --am_initialize_q_calcium_equilibrium=true \
  --am_crossover_pull_through=true \
  --directional=true --aa_local_multilink=true --aa_barbed_zone=0.20 \
  --kappa_aa=${KAPPA_AA} --aa_perpendicular_k=${PERP} --aa_status1_radial_k=75 \
  --aa_preferred_centerline=0.055 --aa_capture_margin=0.045 \
  --aa_antiparallel_overlap=0.10 \
  --aa_status2_registry_on=0.040 --aa_status2_registry_off=0.040 \
  --k_on=2500 --k_off=1 --max_aa_bonds_per_actin=4 \
  --aa_link_rupture_length=0 \
  --titin_mode=explicit_dynamic --titin_k=100 --titin_rest_length=0.665 \
  --titin_capture_distance=0.72 --titin_k_on=5000 --titin_k_off=0.01 \
  --titin_transition_theta=1.0 \
  --dt=${DT} --nsteps=${N} --save_every=${SAVE} --seed=20260806 \
  --filename="${ROOT}/data/${TAG}.h5"
echo "TAG=${TAG} dt=${DT} kappa_aa=${KAPPA_AA} align=${ALIGN} perp=${PERP} headt=${HEADT} rc=$?"
