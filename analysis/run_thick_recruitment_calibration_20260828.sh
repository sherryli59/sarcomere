#!/usr/bin/env bash
set -euo pipefail

export HDF5_USE_FILE_LOCKING=FALSE
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MPLCONFIGDIR=/tmp/sarcomere-mpl-calibration
mkdir -p "${MPLCONFIGDIR}"

repo_dir=/home/groups/rotskoff/sherry/sarcomere/sarcomere
exe="${repo_dir}/cpp/build_release/sarcomere"
output_dir="${1:-${repo_dir}/analysis/thick_recruitment_calibration_20260828_v2}"
mkdir -p "${output_dir}"

dt=1e-4
nsteps=10000
save_every=100
for scale in 1 3 6; do
    inputs=()
    for tension in 0 1 3 6; do
        tag="scale-${scale}_tension-${tension}"
        output="${output_dir}/${tag}.h5"
        inputs+=("${output}")
        if [[ -s "${output}" ]]; then
            continue
        fi
        "${exe}" \
            --n_actins=1 --n_myosins=64 --Lx=10 --Ly=10 --Lz=10 \
            --dimension=3 --pbc_mask=111 --initial_structure=myosin_x_noise \
            --pre_phase_steps=0 --interaction_ramp_steps=0 \
            --dynamic_am=true --am_heads_per_half=9 --am_k_on=0 \
            --am_thin_activation=false --am_thin_feedback=false \
            --am_activation_floor=1 --am_cycle_period=0 \
            --am_thick_feedback=true --am_thick_recruit_rate=2 \
            --am_thick_fold_rate=18 --am_thick_initial_on_fraction=0.1 \
            --am_thick_tension_scale="${scale}" \
            --am_thick_force_clamp="${tension}" \
            --sterics=false --thermal_noise_scale=0 --n_fixed_myosins=64 \
            --aa_local_multilink=true --k_on=0 --base_lifetime=0 \
            --directional_base_lifetime=0 --lifetime_coeff=0 \
            --dt="${dt}" --nsteps="${nsteps}" --save_every="${save_every}" \
            --seed="$((628000 + 100 * scale + tension))" --deterministic \
            --filename="${output}" >"${output%.h5}.log" 2>&1
    done
    python "${repo_dir}/analysis/check_thick_recruitment.py" \
        "${inputs[@]}" --dt="${dt}" --tail-fraction=0.5 \
        --output-prefix="${output_dir}/scale-${scale}_audit"
done

