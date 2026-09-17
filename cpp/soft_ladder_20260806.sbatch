#!/bin/bash
# Uniform spring softening + dt ladder.
# All stiff translational springs -> 285 pN/um so tau = 1/(k*mu) = 5.0e-5 s with
# mu_actin = 70.2 um/(pN*s); dt=1e-5 then sits at 0.2*tau. am_power_stroke is
# scaled up so motor force k*stroke stays at ~5 pN (stall is 4.5).
# Arm 0 is the unsoftened reference at the current dt.
#SBATCH --partition=owners
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8GB
#SBATCH --time=01:00:00
#SBATCH --job-name=soft_ladder
#SBATCH --no-requeue
#SBATCH --array=0-3
#SBATCH --output=/scratch/users/xsli/sarcomere/soft_ladder/logs/%A_%a.out
#SBATCH --error=/scratch/users/xsli/sarcomere/soft_ladder/logs/%A_%a.err

set -uo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export OMP_DYNAMIC=FALSE HDF5_USE_FILE_LOCKING=FALSE

ROOT=/scratch/users/xsli/sarcomere/soft_ladder
mkdir -p "${ROOT}"/{data,logs}
EXE=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp/build_crossover_20260806/sarcomere

MS=20   # same physical time in every arm
case ${SLURM_ARRAY_TASK_ID} in
  0) TAG=ref_stiff_dt5e-7; DT=5e-7;  N=40000; SAVE=2000; K=stiff ;;
  1) TAG=soft_dt5e-6;      DT=5e-6;  N=4000;  SAVE=200;  K=soft ;;
  2) TAG=soft_dt1e-5;      DT=1e-5;  N=2000;  SAVE=100;  K=soft ;;
  3) TAG=soft_dt2e-5;      DT=2e-5;  N=1000;  SAVE=50;   K=soft ;;
esac

if [ "$K" = "soft" ]; then
  SPR="--am_tether_k=285 --actin_steric_k=285 --am_steric_k=285 \
       --myosin_steric_k=285 --am_motor_k=285 --k_aa=285 \
       --am_power_stroke=0.0175"
else
  SPR="--am_tether_k=3000 --actin_steric_k=2000 --am_steric_k=1000 \
       --myosin_steric_k=1000 --am_motor_k=1000 --k_aa=1000 \
       --am_power_stroke=0.005"
fi

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
  --sterics=true ${SPR} \
  --am_cutoff=0.10 --am_optimal=0.055 --am_surface_gap=0.015 \
  --k_am=0 --kappa_am=0 \
  --dynamic_am=true --am_offcenter_torque=true --am_spatial_heads=true \
  --am_heads_per_half=9 --am_sites_per_actin=182 --am_max_heads_per_actin=24 \
  --am_k_on=400 --am_k_weak_off=800 --am_k_powerstroke=800 --am_k_strong_off=640 \
  --am_unloaded_velocity=20 --am_stall_force=4.5 \
  --am_alignment_k=120 --am_head_transverse_k=100 \
  --am_strain_dependent_detachment=true --am_strain_detachment_theta=0.1 \
  --am_thin_activation=true --am_initialize_q_calcium_equilibrium=true \
  --am_crossover_pull_through=true \
  --directional=true --aa_local_multilink=true --aa_barbed_zone=0.20 \
  --kappa_aa=100 --aa_perpendicular_k=50 --aa_status1_radial_k=75 \
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
echo "TAG=${TAG} dt=${DT} springs=${K} rc=$?"
