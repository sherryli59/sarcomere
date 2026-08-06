#!/bin/bash
#SBATCH --partition=owners,hns,normal
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4GB
#SBATCH --time=02:00:00
#SBATCH --job-name=am_fv_xlink
#SBATCH --array=0-19
#SBATCH --open-mode=append
#SBATCH --no-requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/logs_am_force_velocity_crosslink/%A_%a.out
#SBATCH --error=/scratch/users/xsli/sarcomere/logs_am_force_velocity_crosslink/%A_%a.err

set -euo pipefail

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH
export OMP_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

CPP_DIR=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp
EXE="${CPP_DIR}/build_release/sarcomere"
ROOT=/scratch/users/xsli/sarcomere
OUT="${AM_FV_XLINK_OUTPUT_DIR:-${ROOT}/data_am_force_velocity_crosslink_overlap_20260721}"
mkdir -p "${OUT}"

conditions=(free_noaa free_persistentaa isometric_noaa isometric_persistentaa)
condition_index=$((SLURM_ARRAY_TASK_ID / 5))
rep=$((SLURM_ARRAY_TASK_ID % 5))
condition="${conditions[condition_index]}"
seed=$((215001 + rep))

aa_k_on=0
base_lifetime=0
isometric=false
case "${condition}" in
    free_noaa) ;;
    free_persistentaa) aa_k_on=5000; base_lifetime=1 ;;
    isometric_noaa) isometric=true ;;
    isometric_persistentaa) aa_k_on=5000; base_lifetime=1; isometric=true ;;
    *) echo "Unknown condition ${condition}" >&2; exit 2 ;;
esac

printf -v rep_label 'r%02d' "${rep}"
output="${OUT}/am_xlink_${condition}_${rep_label}.h5"
scontrol update jobid="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID}" name="fvx_${condition}_${rep}" 2>/dev/null || true
exec 9>"${output}.lock"
flock -n 9 || exit 0
if [ -e "${output}" ]; then
    echo "Output already exists: ${output}"
    exit 0
fi

srun --kill-on-bad-exit=1 "${EXE}" \
    --dynamic_am=true --initial_structure=cb --pre_phase_steps=0 \
    --n_actins=2 --n_myosins=4 --n_fixed_myosins=4 \
    --Lx=4 --Ly=1 --Lz=1 --dimension=2 --pbc_mask=100 \
    --physical_mobility=false \
    --actin_diff_coeff_trans=0.2 --actin_diff_coeff_rot=0 \
    --myosin_diff_coeff_trans=0.05 --myosin_diff_coeff_rot=0.05 \
    --thermal_noise_scale=1 --isometric_actin_axial="${isometric}" \
    --am_cutoff=0.08 --am_optimal=0.06 \
    --am_spatial_heads=true --am_heads_per_half=6 \
    --am_tether_k=1000 --am_tether_slack=0.005 \
    --am_binding_azimuth_selection=true \
    --am_k_on=1000 --am_k_weak_off=200 --am_k_powerstroke=200 \
    --am_k_strong_off=80 --am_force_scale=1 \
    --am_motor_k=1000 --am_power_stroke=0.005 \
    --am_stall_force=4.5 --am_unloaded_velocity=5 \
    --k_on="${aa_k_on}" --k_off=1 --base_lifetime="${base_lifetime}" \
    --directional_base_lifetime=0 --lifetime_coeff=0 \
    --k_aa=1000 --kappa_aa=100 --aa_transition_state_theta=0.1 \
    --aa_antiparallel_overlap=0.1 \
    --max_myosin_bonds=2 --interaction_ramp_steps=5000 \
    --k_mm=0 \
    --dt=1e-6 --nsteps=50000 --save_every=5000 \
    --seed="${seed}" --deterministic \
    --filename="${output}"
