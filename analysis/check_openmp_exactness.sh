#!/bin/bash
set -euo pipefail

if [ "$#" -lt 2 ] || [ "$#" -gt 3 ]; then
    echo "usage: $0 REFERENCE_EXE CANDIDATE_EXE [OUTPUT_DIR]" >&2
    exit 2
fi

reference_exe=$1
candidate_exe=$2
output_dir=${3:-/tmp/sarcomere_openmp_exactness}
mkdir -p "${output_dir}"
exact_dt=${SARC_EXACT_DT:-1e-6}
exact_nsteps=${SARC_EXACT_NSTEPS:-60}
exact_save_every=${SARC_EXACT_SAVE_EVERY:-1}

export OMP_DYNAMIC=FALSE
export OPENBLAS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export MKL_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

run_case() {
    local executable=$1
    local output=$2
    local threads=$3
    local schedule=$4
    env OMP_NUM_THREADS="${threads}" OMP_SCHEDULE="${schedule}" "${executable}" \
        --n_actins=64 --n_myosins=32 \
        --Lx=3.86 --Ly=1 --Lz=1 --dimension=3 --pbc_mask=111 \
        --actin_length=1 --myosin_length=1.5 \
        --actin_radius=0.015 --myosin_radius=0.025 \
        --initial_structure=myosin_x_noise \
        --actin_direction_noise=-1 \
        --myosin_direction_noise=1 --myosin_direction_tilt_deg=0 \
        --pre_phase_steps=10 --interaction_ramp_steps=20 \
        --physical_mobility=true --temperature_kelvin=310 \
        --thermal_noise_scale=1 --solvent_viscosity=0.001 \
        --actin_crowding_parallel=10 --actin_crowding_perpendicular=10 \
        --actin_crowding_tumble=160 --actin_crowding_roll=160 \
        --myosin_crowding_parallel=20 --myosin_crowding_perpendicular=20 \
        --myosin_crowding_tumble=64 --myosin_crowding_roll=64 \
        --max_actin_displacement=0.005 --max_myosin_displacement=0.005 \
        --max_actin_rotation=0.005 --max_myosin_rotation=0.005 \
        --sterics=true \
        --actin_steric_k=10000 --myosin_steric_k=1000 --am_steric_k=1000 \
        --am_surface_gap=0.015 --k_am=10 --kappa_am=100 \
        --dynamic_am=true --am_offcenter_torque=true --am_spatial_heads=true \
        --am_heads_per_half=141 --am_max_heads_per_actin_half=24 \
        --am_k_on=100 --am_weak_binding_floor=0.1 \
        --am_k_weak_off=200 --am_k_powerstroke=200 \
        --am_k_strong_off=160 --am_force_scale=1 \
        --am_motor_k=1000 --am_power_stroke=0.005 \
        --am_stall_force=4.5 --am_unloaded_velocity=5 \
        --am_weak_stiffness_ratio=0.1 \
        --am_tether_k=3000 --am_tether_slack=0.002 \
        --am_alignment_k=30 --am_alignment_slack_deg=2 \
        --am_tether_saturation_heads=6 \
        --am_binding_azimuth_selection=true \
        --am_binding_azimuth_half_width_deg=30 \
        --am_rupture_distance=0.1 \
        --am_strain_dependent_detachment=true \
        --am_strain_detachment_theta=0.1 \
        --am_strain_detachment_max_exponent=20 \
        --am_strong_off_floor_fraction=1 \
        --am_cycle_period=0.05 \
        --am_activation_rise=0.02 --am_activation_decay=0.2 \
        --am_activation_floor=0 --am_activation_active_fraction=0.5 \
        --am_thin_activation=true --am_thin_feedback=true \
        --am_q_calcium_open_rate=20 --am_q_calcium_half=0.5 \
        --am_q_calcium_hill=3 \
        --am_q_strong_feedback_rate=100 --am_q_strong_half=0.02 \
        --am_q_strong_hill=2 --am_q_close_rate=20 \
        --am_thick_feedback=false --am_status2_memory=0 \
        --aa_local_multilink=true \
        --aa_spacing_ratio=0.92 --aa_capture_margin=0.02 \
        --aa_antiparallel_overlap=0.1 \
        --aa_prebundle_links_per_pair=2 \
        --aa_site_exclusion=0.0192 --aa_overlap_capture=0.1 \
        --k_aa=3000 --kappa_aa=50 --aa_perpendicular_k=50 \
        --aa_stiffness_reference_links=3 --aa_total_link_pool=0 \
        --aa_status_mode=arm_load_productive \
        --aa_status2_force_on=1 --aa_status2_force_off=0.5 \
        --k_on=1000 --base_lifetime=0.1 \
        --directional_base_lifetime=0.05 --lifetime_coeff=0.4 \
        --stabilize_status1=false \
        --aa_energy_gated_association=true \
        --aa_transition_state_theta=0.1 \
        --aa_encounter_k=10 --aa_encounter_kappa=0 \
        --aa_passive_detailed_balance=false \
        --tau_rec=0.05 --max_aa_bonds_per_actin=3 \
        --aa_neighbor_identity=true \
        --aa_neighbor_slot_memory=0.02 --aa_status2_slot_memory=0.20 \
        --directional=true --max_myosin_bonds=10 \
        --k_mm=0 --titin_k=0 --k_bundle_max=0 --bundle_ramp_steps=0 \
        --dt="${exact_dt}" --nsteps="${exact_nsteps}" \
        --save_every="${exact_save_every}" \
        --seed=2072301 --filename="${output}" >/dev/null
}

reference_output="${output_dir}/reference.h5"
candidate_output="${output_dir}/candidate.h5"
rm -f "${reference_output}" "${candidate_output}"
reference_threads=${SARC_REFERENCE_THREADS:-${SARC_EXACT_THREADS:-8}}
candidate_threads=${SARC_CANDIDATE_THREADS:-${SARC_EXACT_THREADS:-8}}
reference_schedule=${SARC_REFERENCE_SCHEDULE:-static}
candidate_schedule=${SARC_CANDIDATE_SCHEDULE:-static}
run_case "${reference_exe}" "${reference_output}" \
    "${reference_threads}" "${reference_schedule}"
run_case "${candidate_exe}" "${candidate_output}" \
    "${candidate_threads}" "${candidate_schedule}"

python - "${reference_output}" "${candidate_output}" <<'PY'
import sys
import h5py
import numpy as np

left_name, right_name = sys.argv[1:]
mismatches = []

def values_equal(left, right):
    left = np.asarray(left)
    right = np.asarray(right)
    return (
        left.shape == right.shape
        and left.dtype == right.dtype
        and np.array_equal(left, right, equal_nan=True)
    )

with h5py.File(left_name, "r", locking=False) as left, \
     h5py.File(right_name, "r", locking=False) as right:
    left_objects = {}
    right_objects = {}
    left.visititems(lambda name, obj: left_objects.setdefault(name, obj))
    right.visititems(lambda name, obj: right_objects.setdefault(name, obj))
    if set(left_objects) != set(right_objects):
        mismatches.append("object names differ")
    for name in sorted(set(left_objects) & set(right_objects)):
        lobj = left_objects[name]
        robj = right_objects[name]
        if type(lobj) is not type(robj):
            mismatches.append(f"{name}: object types differ")
            continue
        if isinstance(lobj, h5py.Dataset) and not values_equal(lobj[...], robj[...]):
            mismatches.append(f"{name}: dataset differs")
        if set(lobj.attrs) != set(robj.attrs):
            mismatches.append(f"{name}: attribute names differ")
        for attr in set(lobj.attrs) & set(robj.attrs):
            if not values_equal(lobj.attrs[attr], robj.attrs[attr]):
                mismatches.append(f"{name}@{attr}: attribute differs")

if mismatches:
    print("EXACTNESS FAIL")
    print("\n".join(mismatches[:50]))
    raise SystemExit(1)
print("EXACTNESS PASS: all HDF5 datasets and attributes are bitwise equal")
PY
