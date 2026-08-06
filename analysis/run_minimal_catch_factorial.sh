#!/usr/bin/env bash
set -euo pipefail

export HDF5_USE_FILE_LOCKING=FALSE
export OMP_NUM_THREADS=1

repo_dir=/home/groups/rotskoff/sherry/sarcomere/sarcomere
exe="${repo_dir}/cpp/build_release/sarcomere"
output_dir="${1:-${repo_dir}/analysis/minimal_catch_factorial_20260722}"
mkdir -p "${output_dir}"

read -r -a variants <<< "${MINIMAL_VARIANTS:-perfect anneal one_motor no_motor parallel wrong_polarity}"
replicas="${MINIMAL_REPLICAS:-3}"
aa_perpendicular_k="${MINIMAL_AA_PERPENDICULAR_K:-50}"
aa_stiffness_reference_links="${MINIMAL_AA_STIFFNESS_REFERENCE_LINKS:-3}"
status2_memory="${MINIMAL_STATUS2_MEMORY:-0.0005}"
lifetime_coeff="${MINIMAL_LIFETIME_COEFF:-0.4}"
aa_links_per_pair="${MINIMAL_AA_LINKS_PER_PAIR:-2}"
aa_k_on="${MINIMAL_AA_K_ON:-1000}"
aa_base_lifetime="${MINIMAL_AA_BASE_LIFETIME:-0.1}"
aa_directional_lifetime="${MINIMAL_AA_DIRECTIONAL_LIFETIME:-0.05}"
aa_status2_force_on="${MINIMAL_AA_STATUS2_FORCE_ON:-1.0}"
aa_status2_force_off="${MINIMAL_AA_STATUS2_FORCE_OFF:-0.5}"
aa_status_mode="${MINIMAL_AA_STATUS_MODE:-productive_register}"
am_force_scale="${MINIMAL_AM_FORCE_SCALE:-1}"
am_cycle_period="${MINIMAL_AM_CYCLE_PERIOD:-0.05}"
am_activation_rise="${MINIMAL_AM_ACTIVATION_RISE:-0.002}"
am_activation_decay="${MINIMAL_AM_ACTIVATION_DECAY:-0.008}"
am_activation_active_fraction="${MINIMAL_AM_ACTIVATION_ACTIVE_FRACTION:-0.5}"
am_heads_per_half="${MINIMAL_AM_HEADS_PER_HALF:-141}"
am_max_heads_per_actin_half="${MINIMAL_AM_MAX_HEADS_PER_ACTIN_HALF:-24}"
am_strong_off_floor_fraction="${MINIMAL_AM_STRONG_OFF_FLOOR_FRACTION:-1}"
am_k_strong_off="${MINIMAL_AM_K_STRONG_OFF:-160}"
am_k_on="${MINIMAL_AM_K_ON:-1000}"
am_weak_binding_floor="${MINIMAL_AM_WEAK_BINDING_FLOOR:-0.1}"
am_motor_k="${MINIMAL_AM_MOTOR_K:-1000}"
am_tether_k="${MINIMAL_AM_TETHER_K:-1000}"
am_alignment_k="${MINIMAL_AM_ALIGNMENT_K:-30}"
am_alignment_slack_deg="${MINIMAL_AM_ALIGNMENT_SLACK_DEG:-2}"
thermal_noise_scale="${MINIMAL_THERMAL_NOISE_SCALE:-1}"
aa_transition_state_theta="${MINIMAL_AA_TRANSITION_STATE_THETA:-0.1}"
minimal_dt="${MINIMAL_DT:-1e-6}"
minimal_save_every="${MINIMAL_SAVE_EVERY:-200}"
minimal_nsteps="${MINIMAL_NSTEPS:-200000}"
for variant in "${variants[@]}"; do
    for ((replica=0; replica<replicas; ++replica)); do
        seed=$((227100 + 10 * replica + ${#variant}))
        structure=minimal_catch
        if [[ "${variant}" != perfect ]]; then
            structure="minimal_catch_${variant}"
        fi
        n_actins=2
        lx=5
        pbc_mask=000
        if [[ "${variant}" == balanced ]]; then
            n_actins=4
            lx=5.32
            pbc_mask=100
        fi
        output="${output_dir}/${variant}_rep${replica}.h5"
        if [[ -s "${output}" ]]; then
            echo "skip existing ${output}"
            continue
        fi
        echo "run ${variant} replica ${replica}"
        "${exe}" \
            --n_actins="${n_actins}" --n_myosins=2 \
            --Lx="${lx}" --Ly=3 --Lz=3 --dimension=3 --pbc_mask="${pbc_mask}" \
            --initial_structure="${structure}" --pre_phase_steps=0 \
            --interaction_ramp_steps=5000 \
            --physical_mobility=true --temperature_kelvin=310 \
            --thermal_noise_scale="${thermal_noise_scale}" \
            --myosin_radius=0.025 --actin_radius=0.015 \
            --am_surface_gap=0.015 --sterics=true \
            --actin_steric_k=100 --myosin_steric_k=1000 --am_steric_k=1000 \
            --dynamic_am=true --am_spatial_heads=true --am_heads_per_half="${am_heads_per_half}" \
            --am_max_heads_per_actin_half="${am_max_heads_per_actin_half}" \
            --am_binding_azimuth_selection=true --am_tether_k="${am_tether_k}" \
            --am_tether_slack=0.005 --am_tether_saturation_heads=6 \
            --am_alignment_k="${am_alignment_k}" \
            --am_alignment_slack_deg="${am_alignment_slack_deg}" \
            --am_rupture_distance=0.1 \
            --am_strain_dependent_detachment=true \
            --am_strain_detachment_theta=0.1 \
            --am_strain_detachment_max_exponent=20 \
            --am_strong_off_floor_fraction="${am_strong_off_floor_fraction}" \
            --am_k_on="${am_k_on}" --am_weak_binding_floor="${am_weak_binding_floor}" \
            --am_k_weak_off=200 --am_k_powerstroke=200 \
            --am_k_strong_off="${am_k_strong_off}" --am_force_scale="${am_force_scale}" \
            --am_motor_k="${am_motor_k}" --am_power_stroke=0.005 \
            --am_stall_force=4.5 --am_unloaded_velocity=5 \
            --am_cycle_period="${am_cycle_period}" --am_activation_rise="${am_activation_rise}" \
            --am_activation_decay="${am_activation_decay}" --am_activation_floor=0 \
            --am_activation_active_fraction="${am_activation_active_fraction}" \
            --am_thin_activation=true --am_thin_feedback=true \
            --am_q_calcium_open_rate=100 --am_q_close_rate=100 \
            --am_thick_feedback=true --am_thick_recruit_rate=20 \
            --am_thick_fold_rate=20 --am_thick_tension_scale=0.1 \
            --am_thick_initial_on_fraction=0.5 \
            --aa_local_multilink=true --aa_prebundle_links_per_pair="${aa_links_per_pair}" \
            --aa_spacing_ratio=0.92 --aa_capture_margin=0.02 \
            --aa_antiparallel_overlap=0.1 --aa_perpendicular_k="${aa_perpendicular_k}" \
            --aa_stiffness_reference_links="${aa_stiffness_reference_links}" \
            --aa_status_mode="${aa_status_mode}" \
            --aa_status2_force_on="${aa_status2_force_on}" \
            --aa_status2_force_off="${aa_status2_force_off}" \
            --am_status2_memory="${status2_memory}" \
            --k_on="${aa_k_on}" --k_off=1 --base_lifetime="${aa_base_lifetime}" \
            --directional_base_lifetime="${aa_directional_lifetime}" --lifetime_coeff="${lifetime_coeff}" \
            --stabilize_status1=false \
            --k_aa=3000 --kappa_aa=50 \
            --aa_energy_gated_association=true \
            --aa_transition_state_theta="${aa_transition_state_theta}" --aa_encounter_k=10 \
            --aa_encounter_kappa=0 --aa_passive_detailed_balance=false \
            --tau_rec=0.002 --max_aa_bonds_per_actin=2 \
            --max_myosin_bonds=2 --directional=true \
            --dt="${minimal_dt}" --nsteps="${minimal_nsteps}" --save_every="${minimal_save_every}" \
            --seed="${seed}" --deterministic --filename="${output}" \
            >"${output%.h5}.log" 2>&1
    done
done

python "${repo_dir}/analysis/analyze_minimal_catch_factorial.py" \
    "${output_dir}" --output "${output_dir}/summary.csv" \
    --plot "${output_dir}/heartbeat_and_persistence.png"
