#!/bin/bash
#SBATCH --partition=hns,normal
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --time=02:00:00
#SBATCH --job-name=titin_half_contact
#SBATCH --open-mode=append
#SBATCH --no-requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/titin_pilot_iterations_20260801/logs/half_contact_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/titin_pilot_iterations_20260801/logs/half_contact_%j.err

set -euo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export OMP_DYNAMIC=FALSE OMP_SCHEDULE=static
export OPENBLAS_NUM_THREADS=1 GOTO_NUM_THREADS=1 MKL_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE
JEMALLOC=/home/groups/rotskoff/sherry/mamba_root/envs/sarcomere/lib/libjemalloc.so
test -r "${JEMALLOC}"
export LD_PRELOAD="${JEMALLOC}"

ROOT=/scratch/users/xsli/sarcomere/titin_pilot_iterations_20260801
EXE=${ROOT}/bin/sarcomere_3a5061f9a6e2
EXPECTED=3a5061f9a6e20993424f87409315645ed22acd74f125ad7c3de6c181828fb28a
OUTPUT=${ROOT}/half_fullstack_rotcontact_v5_dt0p5us_pre5ms_ramp20ms_seed2080101.h5
test "$(sha256sum "${EXE}" | awk '{print $1}')" = "${EXPECTED}"

resume_flag=()
if [ -s "${OUTPUT}" ]; then resume_flag=(--resume=true); fi

cmd=(srun --kill-on-bad-exit=1 "${EXE}"
 --n_actins=640 --n_myosins=320 --Lx=6.127 --Ly=1.5874 --Lz=1.5874 --dimension=3 --pbc_mask=111
 --actin_length=1 --myosin_length=1.5 --actin_radius=0.015 --myosin_radius=0.025
 --initial_structure=myosin_x_noise --actin_direction_noise=-1 --myosin_direction_noise=1 --myosin_direction_tilt_deg=0
 --pre_phase_steps=10000 --interaction_ramp_steps=40000
 --physical_mobility=true --temperature_kelvin=310 --thermal_noise_scale=1 --solvent_viscosity=0.001
 --actin_hydrodynamic_diameter=0.010 --myosin_hydrodynamic_diameter=0.035
 --actin_crowding_parallel=10 --actin_crowding_perpendicular=10 --actin_crowding_tumble=160 --actin_crowding_roll=160
 --myosin_crowding_parallel=20 --myosin_crowding_perpendicular=20 --myosin_crowding_tumble=64 --myosin_crowding_roll=64
 --max_actin_displacement=0.005 --max_myosin_displacement=0.005 --max_actin_rotation=0.005 --max_myosin_rotation=0.005
 --sterics=true --actin_steric_k=2000 --myosin_steric_k=1000 --am_steric_k=1000 --am_surface_gap=0.015
 --k_am=10 --kappa_am=100 --dynamic_am=true --am_offcenter_torque=true --am_spatial_heads=true
 --am_heads_per_half=141 --am_max_heads_per_actin_half=24 --am_k_on=100 --am_weak_binding_floor=0.1
 --am_k_weak_off=200 --am_k_powerstroke=200 --am_k_strong_off=160 --am_force_scale=1
 --am_motor_k=1000 --am_power_stroke=0.005 --am_stall_force=4.5 --am_unloaded_velocity=5
 --am_weak_stiffness_ratio=0.1 --am_tether_k=3000 --am_tether_slack=0.002
 --am_alignment_k=30 --am_alignment_slack_deg=2 --am_tether_saturation_heads=6
 --am_binding_azimuth_selection=true --am_binding_azimuth_half_width_deg=30
 --am_rupture_distance=0.1 --am_strain_dependent_detachment=true --am_strain_detachment_theta=0.1
 --am_strain_detachment_max_exponent=20 --am_strong_off_floor_fraction=1
 --am_cycle_period=0 --am_activation_floor=1 --am_activation_active_fraction=1
 --am_thin_activation=true --am_thin_feedback=true --am_q_calcium_open_rate=20 --am_q_calcium_half=0.5
 --am_q_calcium_hill=3 --am_q_strong_feedback_rate=100 --am_q_strong_half=0.02 --am_q_strong_hill=2
 --am_q_close_rate=20 --am_thick_feedback=false --am_status2_memory=0
 --aa_local_multilink=true --aa_spacing_ratio=1 --aa_capture_margin=0.02 --aa_antiparallel_overlap=0.1
 --aa_barbed_zone=0.1 --aa_axial_register_ratio=3 --aa_axial_register_force_max=36
 --aa_register_all_statuses=false --aa_prebundle_links_per_pair=2 --aa_site_exclusion=0.0192
 --aa_overlap_capture=0.1 --k_aa=1000 --kappa_aa=100 --aa_perpendicular_k=50
 --aa_stiffness_reference_links=3 --aa_total_link_pool=0 --aa_status_mode=productive_register
 --aa_status2_force_on=1 --aa_status2_force_off=0.5 --k_on=5000 --k_off=1 --base_lifetime=1
 --directional_base_lifetime=0.5 --lifetime_coeff=0.4 --stabilize_status1=false
 --aa_energy_gated_association=false --aa_encounter_k=10 --aa_encounter_kappa=0
 --aa_passive_detailed_balance=false --tau_rec=0 --max_aa_bonds_per_actin=2
 --aa_neighbor_identity=true --aa_neighbor_slot_memory=0.02 --aa_status2_slot_memory=0.20
 --aa_association_max_direction_dot=-0.5 --aa_link_rupture_length=0.12
 --directional=true --max_myosin_bonds=24 --titin_mode=explicit_dynamic
 --titin_k=40 --titin_rest_length=0.5 --titin_links_per_half=3 --titin_capture_distance=0.70
 --titin_k_on=5000 --titin_k_off=0.5 --titin_strain_gated_kinetics=true --titin_transition_theta=0.5
 --k_mm=0 --k_bundle_max=0 --bundle_ramp_steps=0
 --dt=5e-7 --nsteps=40000 --save_every=2000 --seed=2080101 --filename="${OUTPUT}")
if [ "${#resume_flag[@]}" -gt 0 ]; then
    cmd+=("${resume_flag[@]}")
fi

printf 'command='; printf '%q ' "${cmd[@]}"; printf '\n'
exec "${cmd[@]}"
