#!/bin/bash
#SBATCH --partition=owners,hns,normal
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4GB
#SBATCH --time=01:00:00
#SBATCH --job-name=six_barbed_exp
#SBATCH --array=0-5
#SBATCH --open-mode=append
#SBATCH --no-requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/logs_six_actin_barbed_explicit/%A_%a.out
#SBATCH --error=/scratch/users/xsli/sarcomere/logs_six_actin_barbed_explicit/%A_%a.err

set -euo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH
export OMP_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

CPP_DIR=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp
EXE="${CPP_DIR}/build_release/sarcomere"
ROOT=/scratch/users/xsli/sarcomere
OUT="${SIX_BARBED_OUTPUT_DIR:-${ROOT}/data_six_actin_overlap_explicit_20260721}"
mkdir -p "${OUT}" /scratch/users/xsli/sarcomere/logs_six_actin_barbed_explicit

rep=$((SLURM_ARRAY_TASK_ID % 3))
if [ "${SLURM_ARRAY_TASK_ID}" -lt 3 ]; then
    condition=perfect
    angle_noise=0
    radial_noise=0
    prebind=true
    seed=$((221000 + rep))
else
    condition=anneal
    angle_noise=20
    radial_noise=0.005
    prebind=false
    seed=$((222000 + rep))
fi
output="${OUT}/${condition}_rep-${rep}.h5"
exec 9>"${output}.lock"
flock -n 9 || exit 0
[ ! -e "${output}" ] || exit 0

srun --kill-on-bad-exit=1 "${EXE}" \
    --nsteps=200000 --dt=1e-6 --save_every=2000 \
    --n_actins=6 --n_myosins=1 --n_fixed_myosins=1 \
    --Lx=2 --Ly=2 --Lz=2 --pbc_mask=000 \
    --initial_structure=six_actin_ring --pre_phase_steps=0 \
    --six_ring_angle_noise_deg="${angle_noise}" \
    --six_ring_radial_noise="${radial_noise}" \
    --six_ring_prebind_aa="${prebind}" \
    --physical_mobility=true --temperature_kelvin=310 \
    --isometric_actin_axial=true --directional=false \
    --myosin_radius=0.025 --actin_radius=0.015 --am_surface_gap=0.015 \
    --aa_spacing_ratio=0.92 --aa_capture_margin=0.02 \
    --aa_antiparallel_overlap=0.1 \
    --k_am=3000 --kappa_am=50 --k_aa=3000 --kappa_aa=50 \
    --dynamic_am=true --am_spatial_heads=true --am_heads_per_half=6 \
    --am_tether_k=1000 --am_tether_slack=0.005 \
    --am_binding_azimuth_selection=true \
    --am_k_on=1000 --am_k_weak_off=200 --am_k_powerstroke=200 \
    --am_k_strong_off=80 --am_force_scale=1 \
    --am_motor_k=1000 --am_power_stroke=0.005 \
    --am_stall_force=4.5 --am_unloaded_velocity=5 \
    --k_on=5000 --base_lifetime=10 --directional_base_lifetime=0 \
    --lifetime_coeff=0 --tau_rec=0 --max_aa_bonds_per_actin=2 \
    --aa_energy_gated_association=true --aa_transition_state_theta=0.1 \
    --aa_passive_detailed_balance=false \
    --max_myosin_bonds=6 --interaction_ramp_steps=5000 \
    --seed="${seed}" --deterministic --filename="${output}"
