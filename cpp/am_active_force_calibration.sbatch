#!/bin/bash
#SBATCH --partition=owners,hns,normal
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=24GB
#SBATCH --time=04:00:00
#SBATCH --job-name=am_force_cal
#SBATCH --array=0-14
#SBATCH --open-mode=append
#SBATCH --no-requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/logs_am_force_calibration/%A_%a.out
#SBATCH --error=/scratch/users/xsli/sarcomere/logs_am_force_calibration/%A_%a.err

set -euo pipefail

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export HDF5_USE_FILE_LOCKING=FALSE

CPP_DIR=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp
EXE="${CPP_DIR}/build_release/sarcomere"
ROOT=/scratch/users/xsli/sarcomere
OUT="${AM_FORCE_CAL_OUTPUT_DIR:-${ROOT}/data_am_force_calibration_20260721}"
mkdir -p "${OUT}"

conditions=(k10_f1 k200_f1 k1000_f1 k1000_f4p5 k1000_f4p5_no_thermal)
index=$((SLURM_ARRAY_TASK_ID / 3))
rep=$((SLURM_ARRAY_TASK_ID % 3))
condition="${conditions[index]}"
# Matched seeds across force conditions expose the same initial geometry and
# Brownian number stream before the mechanics diverge.
seed=$((213001 + rep))

motor_k=10
stall_force=1
noise_scale=1
case "${condition}" in
    k10_f1) ;;
    k200_f1) motor_k=200 ;;
    k1000_f1) motor_k=1000 ;;
    k1000_f4p5) motor_k=1000; stall_force=4.5 ;;
    k1000_f4p5_no_thermal) motor_k=1000; stall_force=4.5; noise_scale=0 ;;
    *) echo "Unknown condition ${condition}" >&2; exit 2 ;;
esac

tag="${condition}_seed${seed}"
output="${OUT}/${tag}.h5"
scontrol update jobid="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID}" name="afc_${condition}_${rep}" 2>/dev/null || true
exec 9>"${output}.lock"
flock -n 9 || exit 0
if [ -e "${output}" ]; then
    echo "Output already exists: ${output}"
    exit 0
fi

srun --kill-on-bad-exit=1 "${EXE}" \
    --n_actins=500 --n_myosins=250 \
    --Lx=6 --Ly=1.41421356237 --Lz=1.41421356237 \
    --dimension=3 --pbc_mask=111 \
    --initial_structure=myosin_x_noise \
    --actin_direction_noise=0 --myosin_direction_noise=0 \
    --myosin_direction_tilt_deg=20 \
    --pre_phase_steps=20000 --interaction_ramp_steps=10000 \
    --physical_mobility=true --temperature_kelvin=310 \
    --thermal_noise_scale="${noise_scale}" \
    --actin_radius=0.015 --myosin_radius=0.025 --am_surface_gap=0.015 \
    --actin_steric_k=100 --myosin_steric_k=1000 --am_steric_k=1000 \
    --dynamic_am=true --am_offcenter_torque=true \
    --am_spatial_heads=true --am_heads_per_half=6 \
    --am_tether_k=1000 --am_tether_slack=0.005 \
    --am_binding_azimuth_selection=true \
    --am_binding_azimuth_half_width_deg=30 \
    --k_am=3000 --kappa_am=0 \
    --am_k_on=1000 --am_k_weak_off=200 --am_k_powerstroke=200 \
    --am_k_strong_off=80 --am_force_scale=1 \
    --am_motor_k="${motor_k}" --am_power_stroke=0.005 \
    --am_stall_force="${stall_force}" \
    --k_on=0 --k_aa=0 --kappa_aa=0 \
    --k_mm=0 --k_bundle_max=0 --bundle_ramp_steps=0 \
    --directional=true --sterics=true \
    --dt=1e-6 --nsteps=20000 --save_every=200 \
    --seed="${seed}" --filename="${output}"
