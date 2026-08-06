#!/usr/bin/env bash
set -euo pipefail

repo_root=/home/groups/rotskoff/sherry/sarcomere/sarcomere
steps=${1:-1000000}
output=${2:-$repo_root/data/production_titin_like_20260731/production_1s_h141_seed1001.h5}
resume_value=${3:-false}
heads_per_half=${4:-141}

export OMP_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

exec "$repo_root/cpp/build_titin_production_20260731/sarcomere" \
  --nsteps "$steps" --seed 1001 --dt 1e-6 --save_every 5000 \
  --n_actins 4 --n_myosins 2 --Lx 4.7 --Ly 3 --Lz 3 --pbc_mask 100 \
  --actin_length 1 --myosin_length 1.5 --myosin_radius 0.025 \
  --am_cutoff 0.05 --am_optimal 0.055 --actin_radius 0.015 \
  --am_surface_gap 0.015 --am_steric_k 1000 \
  --actin_steric_k 2000 --myosin_steric_k 1000 \
  --aa_preferred_centerline 0.0506 --aa_capture_margin 0.02 \
  --aa_antiparallel_overlap 0.10 --init_aa_overlap 0.16 \
  --aa_barbed_zone 0.10 --aa_axial_register_ratio 1 \
  --aa_axial_register_force_max 12 --aa_register_all_statuses=true \
  --aa_status_mode productive_register --aa_status2_force_on 1 \
  --aa_status2_force_off 0.5 --aa_status2_registry_on 0.1 \
  --aa_status2_registry_off 0.15 --aa_register_saturation_error 0.05 \
  --aa_local_multilink=true --aa_prebundle_links_per_pair 3 \
  --aa_site_exclusion 0.0192 --aa_overlap_capture 0.1 \
  --aa_perpendicular_k 50 --aa_stiffness_reference_links 3 \
  --aa_neighbor_slot_memory 0 --aa_status2_slot_memory 0 --tau_rec 0 \
  --k_on 1000 --k_off 1 --base_lifetime 0.1 \
  --directional_base_lifetime 0.05 --lifetime_coeff 0.4 \
  --k_aa 1000 --kappa_aa 100 --dynamic_am=true \
  --am_offcenter_torque=true --am_spatial_heads=true \
  --am_heads_per_half "$heads_per_half" --am_max_heads_per_actin_half 24 \
  --max_myosin_bonds 24 --am_k_on 100 --am_weak_binding_floor 0.1 \
  --am_k_weak_off 200 --am_k_powerstroke 200 --am_k_strong_off 160 \
  --am_force_scale 1 --am_motor_k 1000 --am_power_stroke 0.005 \
  --am_stall_force 4.5 --am_unloaded_velocity 5 \
  --am_weak_stiffness_ratio 0.1 --am_tether_k 3000 \
  --am_tether_slack 0.002 --am_alignment_k 30 \
  --am_alignment_slack_deg 2 --am_tether_saturation_heads 6 \
  --am_binding_azimuth_selection=true \
  --am_binding_azimuth_half_width_deg 30 --am_rupture_distance 0.1 \
  --am_strain_dependent_detachment=true --am_strain_detachment_theta 0.1 \
  --am_strain_detachment_max_exponent 20 \
  --am_strong_off_floor_fraction 1 --physical_mobility=true \
  --temperature_kelvin 310 --thermal_noise_scale 1 \
  --solvent_viscosity 0.001 --actin_hydrodynamic_diameter 0.010 \
  --myosin_hydrodynamic_diameter 0.035 --actin_crowding_parallel 10 \
  --actin_crowding_perpendicular 10 --actin_crowding_tumble 160 \
  --actin_crowding_roll 160 --myosin_crowding_parallel 20 \
  --myosin_crowding_perpendicular 20 --myosin_crowding_tumble 64 \
  --myosin_crowding_roll 64 --max_actin_displacement 0.005 \
  --max_myosin_displacement 0.005 --max_actin_rotation 0.005 \
  --max_myosin_rotation 0.005 --titin_k 220 \
  --titin_rest_length 1.125 --titin_mode balanced_owning \
  --fixed_lattice_preequil_steps 50000 --n_fixed_myosins 0 \
  --initial_structure minimal_catch_balanced --dimension 3 \
  --resume="$resume_value" --filename "$output"
