#!/bin/bash
# dt ladder judged on OBSERVABLES, with the displacement cap scaled to dt.
#
# F_max = max_disp / (mu * dt), so holding max_disp fixed while raising dt
# tightens the force cap as 1/dt and manufactures cap events that are not
# instability. Scaling max_disp with dt keeps the cap's physical meaning
# ("do not jump more than X per step") constant across arms, so the counters
# are comparable and the observables are the real test.
#   dt=5e-7 -> 0.005 um   (the historical value)
#   dt=2e-6 -> 0.020
#   dt=5e-6 -> 0.050
#   dt=1e-5 -> 0.100
# Rotation cap scaled the same way.
#SBATCH --partition=owners
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8GB
#SBATCH --time=04:00:00
#SBATCH --job-name=dt_orig
#SBATCH --no-requeue
#SBATCH --array=0-3
#SBATCH --output=/scratch/users/xsli/sarcomere/dt_orig/logs/%A_%a.out
#SBATCH --error=/scratch/users/xsli/sarcomere/dt_orig/logs/%A_%a.err

set -uo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export OMP_DYNAMIC=FALSE HDF5_USE_FILE_LOCKING=FALSE

ROOT=/scratch/users/xsli/sarcomere/dt_orig
mkdir -p "${ROOT}"/{data,logs}
EXE=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp/build_crossover_20260806/sarcomere

# 100 ms of physical time in every arm; same seed.
SPR="--am_tether_k=285 --actin_steric_k=285 --am_steric_k=285 --myosin_steric_k=285 --am_motor_k=285 --k_aa=285 --am_power_stroke=0.0175"
KAPPA=1; ALIGN=30
case ${SLURM_ARRAY_TASK_ID} in
  # do the ORIGINAL stiffnesses survive dt=1e-5 once the cap is scaled?
  0) TAG=orig_all_dt1e-5; DT=1e-5; N=10000; SAVE=200; DISP=0.100; ROT=0.100
     SPR="--am_tether_k=3000 --actin_steric_k=2000 --am_steric_k=1000 --myosin_steric_k=1000 --am_motor_k=1000 --k_aa=1000 --am_power_stroke=0.005"
     KAPPA=100; ALIGN=120 ;;
  1) TAG=orig_kappa_dt1e-5; DT=1e-5; N=10000; SAVE=200; DISP=0.100; ROT=0.100; KAPPA=100 ;;
  2) TAG=orig_align_dt1e-5; DT=1e-5; N=10000; SAVE=200; DISP=0.100; ROT=0.100; ALIGN=120 ;;
  3) TAG=orig_springs_dt1e-5; DT=1e-5; N=10000; SAVE=200; DISP=0.100; ROT=0.100
     SPR="--am_tether_k=3000 --actin_steric_k=2000 --am_steric_k=1000 --myosin_steric_k=1000 --am_motor_k=1000 --k_aa=1000 --am_power_stroke=0.005" ;;
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
  --max_actin_displacement=${DISP} --max_myosin_displacement=${DISP} \
  --max_actin_rotation=${ROT} --max_myosin_rotation=${ROT} \
  --sterics=true \
  ${SPR} \
  --am_cutoff=0.10 --am_optimal=0.055 --am_surface_gap=0.015 \
  --k_am=0 --kappa_am=0 \
  --dynamic_am=true --am_offcenter_torque=true --am_spatial_heads=true \
  --am_heads_per_half=9 --am_sites_per_actin=182 --am_max_heads_per_actin=24 \
  --am_k_on=400 --am_k_weak_off=800 --am_k_powerstroke=800 --am_k_strong_off=640 \
  --am_unloaded_velocity=20 --am_stall_force=4.5 \
  --am_alignment_k=${ALIGN} --am_head_transverse_k=100 \
  --am_strain_dependent_detachment=true --am_strain_detachment_theta=0.1 \
  --am_thin_activation=true --am_initialize_q_calcium_equilibrium=true \
  --am_crossover_pull_through=true \
  --directional=true --aa_local_multilink=true --aa_barbed_zone=0.20 \
  --kappa_aa=${KAPPA} --aa_perpendicular_k=50 --aa_status1_radial_k=75 \
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
echo "TAG=${TAG} dt=${DT} disp_cap=${DISP} rc=$?"
