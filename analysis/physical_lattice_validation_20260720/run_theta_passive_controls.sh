#!/usr/bin/env bash
set -euo pipefail

run_dir="${1:-analysis/physical_lattice_validation_20260720/theta_passive_db}"
binary="${SARCOMERE_BINARY:-cpp/build_release/sarcomere}"
mkdir -p "$run_dir"

common=(
  --nsteps 200000 --dt 1e-6 --save_every 2000
  --n_actins 6 --n_myosins 1
  --Lx 2 --Ly 2 --Lz 2 --pbc_mask 000
  --initial_structure six_actin_ring --pre_phase_steps 0
  --physical_mobility=true --temperature_kelvin 310
  --dynamic_am=false --v_am 0 --n_fixed_myosins 1
  --directional=false --isometric_actin_axial=true
  --k_am 3000 --kappa_am 50 --k_aa 3000 --kappa_aa 50
  --aa_spacing_ratio 0.92 --aa_capture_margin 0.02
  --k_on 5000 --base_lifetime 10
  --directional_base_lifetime 0 --lifetime_coeff 0 --tau_rec 0
  --max_aa_bonds_per_actin 2
  --aa_energy_gated_association=true
  --aa_passive_detailed_balance=true
  --deterministic=true
)

for theta in 0 0.1 0.25 0.5 1; do
  label=${theta/./p}
  for rep in 0 1 2; do
    "$binary" "${common[@]}" \
      --seed "$((10000 + rep))" \
      --aa_transition_state_theta "$theta" \
      --six_ring_angle_noise_deg 20 --six_ring_radial_noise 0.005 \
      --six_ring_prebind_aa=false \
      --filename "$run_dir/anneal_theta${label}_rep-${rep}.h5" \
      >"$run_dir/anneal_theta${label}_rep-${rep}.log" 2>&1

    "$binary" "${common[@]}" \
      --seed "$((20000 + rep))" \
      --aa_transition_state_theta "$theta" \
      --six_ring_angle_noise_deg 0 --six_ring_radial_noise 0 \
      --six_ring_prebind_aa=true \
      --filename "$run_dir/perfect_theta${label}_rep-${rep}.h5" \
      >"$run_dir/perfect_theta${label}_rep-${rep}.log" 2>&1
  done
done
