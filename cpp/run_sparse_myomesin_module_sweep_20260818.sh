#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd "${script_dir}/.." && pwd)"
output_dir="${1:-${script_dir}/diagnostics/sparse_myomesin_module_20260818}"
binary="${repo_dir}/cpp/build/sarcomere"
validator="${repo_dir}/analysis/check_myomesin_annealing.py"

mkdir -p "${output_dir}"

run_case() {
    local axial_k="$1"
    local force_max="$2"
    local stem="module_k${axial_k}_cap${force_max}_rad300_50ms"
    OMP_NUM_THREADS=1 /usr/bin/time -f \
        'WALL_SECONDS=%e MAX_RSS_KB=%M' \
        "${binary}" \
        --filename="${output_dir}/${stem}.h5" \
        --initial_structure=minimal_catch_balanced_double_myosin_anneal \
        --n_actins=4 --n_myosins=4 --Lx=5.32 --Ly=1 --Lz=1 --pbc_mask=100 \
        --nsteps=50000 --dt=1e-6 --save_every=250 --seed=831704 --deterministic \
        --physical_mobility=true --temperature_kelvin=310 --thermal_noise_scale=1 \
        --actin_crowding_parallel=10 --actin_crowding_perpendicular=10 \
        --actin_crowding_tumble=160 --actin_crowding_roll=160 \
        --myosin_crowding_parallel=20 --myosin_crowding_perpendicular=20 \
        --myosin_crowding_tumble=64 --myosin_crowding_roll=64 \
        --actin_steric_k=100 --myosin_steric_k=1000 --myosin_exclusion_distance=0 \
        --dynamic_am=true --am_k_on=1000 --am_k_weak_off=200 \
        --am_k_powerstroke=200 --am_k_strong_off=160 --am_motor_k=1000 \
        --am_stall_force=4.5 --am_tether_k=1000 --am_tether_slack=0.005 \
        --am_activation_floor=1 --am_cycle_period=0.05 \
        --am_activation_active_fraction=1 --am_thin_activation=false \
        --am_thick_feedback=false --am_thick_initial_on_fraction=1 \
        --aa_local_multilink=true --aa_prebundle_links_per_pair=1 \
        --aa_stiffness_reference_links=1 --k_aa=3000 --kappa_aa=50 \
        --base_lifetime=0.1 --directional_base_lifetime=0.05 \
        --lifetime_coeff=0.4 --tau_rec=0.002001 \
        --aa_status_mode=motor_transmitted_register \
        --aa_status2_force_on=1 --aa_status2_force_off=0.5 \
        --aa_motor_pair_support_window=0.005 \
        --aa_motor_pair_commit_time=0.010 \
        --aa_motor_pair_exchange_time=0.001 \
        --aa_motor_pair_dropout_time=0.005 \
        --aa_motor_pair_switch_ratio=1.5 --max_aa_bonds_per_actin=1 \
        --aa_energy_gated_association=true --k_mm="${axial_k}" \
        --myomesin_sites_per_filament=2 --myomesin_k_on=5000 \
        --myomesin_k_off=1 --myomesin_axial_capture=0.25 \
        --myomesin_radial_rest=0.11 --myomesin_radial_capture=0.04 \
        --myomesin_radial_k=300 --myomesin_alignment_cos=0.94 \
        --myomesin_bell_distance=0.0005 \
        --myomesin_force_max="${force_max}" \
        --titin_k=0 --k_bundle_max=0 --pre_phase_steps=0
    python "${validator}" "${output_dir}/${stem}.h5" \
        --pairs 0-2,1-3 --minimum-reduction-fraction 0.10 \
        --json "${output_dir}/${stem}_annealing.json"
}

run_case 120 20
run_case 240 40
