#!/bin/bash
#SBATCH --job-name=alignv2
#SBATCH --partition=rotskoff,owners
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=4:00:00
#SBATCH --requeue
#SBATCH --array=0-13
#SBATCH --output=/scratch/users/xsli/sarcomere/align_probe_v2_20260903/logs/%x_%A_%a.out
#SBATCH --error=/scratch/users/xsli/sarcomere/align_probe_v2_20260903/logs/%x_%A_%a.err
set -euo pipefail
export OMP_NUM_THREADS=1        # 1-4 filaments: parallelism is pure overhead

root=/scratch/users/xsli/sarcomere/align_probe_v2_20260903
exe=${root}/bin/sarcomere_alignprobe_20260903
mkdir -p "${root}"/{data,params,logs}
test -x "${exe}"

configs=(
    "tethered_cap1|tethered|1|0|true|100"
    "tethered_cap2|tethered|2|0|true|100"
    "tethered_cap3|tethered|3|0|true|100"
    "tethered_cap6|tethered|6|0|true|100"
    "tethered_cap6_noofftorque|tethered|6|0|false|100"
    "tethered_cap6_notransverse|tethered|6|0|true|0"
    "tethered_cap6_alignON|tethered|6|120|true|100"
    "free_cap1|free|1|0|true|100"
    "free_cap2|free|2|0|true|100"
    "free_cap3|free|3|0|true|100"
    "free_cap6|free|6|0|true|100"
    "free_cap6_noofftorque|free|6|0|false|100"
    "free_cap6_notransverse|free|6|0|true|0"
    "free_cap6_alignON|free|6|120|true|100"
)

IFS='|' read -r label rig cap align_k offtorque transverse_k \
    <<< "${configs[$SLURM_ARRAY_TASK_ID]}"

case "${rig}" in
    tethered) init=minimal_catch;    na=2; nm=2; kaa=1000 ;;
    free)     init=single_am_tilted; na=1; nm=1; kaa=0 ;;
    *) exit 2 ;;
esac

nsteps=${ALIGN_PROBE_NSTEPS:-50000}
save_every=${ALIGN_PROBE_SAVE_EVERY:-50}

for seed in 1 2 3; do
    tag=${label}_seed${seed}
    out=${root}/data/${tag}.h5
    rm -f "${out}" "${out}.lock"
    printf '%s\n' "label=${label}" "rig=${rig}" "init=${init}" "cap=${cap}" \
        "am_alignment_k=${align_k}" "am_offcenter_torque=${offtorque}" \
        "am_head_transverse_k=${transverse_k}" "k_aa=${kaa}" "seed=${seed}" \
        >"${root}/params/${tag}.params"

    "${exe}" \
      --n_actins=${na} --n_myosins=${nm} --Lx=4 --Ly=1 --Lz=1 \
      --dimension=3 --pbc_mask=111 --actin_length=1 --myosin_length=1.5 \
      --actin_radius=0.015 --myosin_radius=0.05 \
      --initial_structure=${init} \
      --actin_direction_noise=-1 --myosin_direction_noise=0 \
      --pre_phase_steps=0 --interaction_ramp_steps=0 \
      --physical_mobility=true --temperature_kelvin=310 \
      --thermal_noise_scale=1 --solvent_viscosity=0.001 \
      --actin_hydrodynamic_diameter=0.010 --myosin_hydrodynamic_diameter=0.035 \
      --actin_crowding_parallel=10 --actin_crowding_perpendicular=10 \
      --actin_crowding_tumble=160 --actin_crowding_roll=160 \
      --myosin_crowding_parallel=20 --myosin_crowding_perpendicular=20 \
      --myosin_crowding_tumble=64 --myosin_crowding_roll=64 \
      --max_actin_displacement=0.005 --max_myosin_displacement=0.005 \
      --max_actin_rotation=0.005 --max_myosin_rotation=0.005 \
      --sterics=true --actin_steric_k=2000 --myosin_steric_k=1000 \
      --am_steric_k=1000 \
      --am_cutoff=0.05 --am_optimal=0.055 --am_surface_gap=0.015 \
      --k_am=0 --kappa_am=0 --dynamic_am=true \
      --am_offcenter_torque=${offtorque} --am_spatial_heads=true \
      --am_heads_per_half=9 \
      --am_max_heads_per_actin_half=${cap} --am_max_heads_per_actin=48 \
      --am_sites_per_actin=182 --am_k_on=100 --am_weak_binding_floor=1 \
      --am_k_weak_off=50 --am_k_powerstroke=200 --am_k_strong_off=160 \
      --am_force_scale=1 --am_motor_k=1000 --am_power_stroke=0.005 \
      --am_stall_force=4.5 --am_unloaded_velocity=5 \
      --am_weak_stiffness_ratio=0.1 \
      --am_tether_k=3000 --am_tether_slack=0.002 \
      --am_head_transverse_k=${transverse_k} \
      --am_alignment_k=${align_k} --am_alignment_slack_deg=2 \
      --am_tether_saturation_heads=6 \
      --am_binding_azimuth_selection=true \
      --am_binding_azimuth_half_width_deg=30 \
      --am_rupture_distance=0.1 --am_strain_dependent_detachment=true \
      --am_strain_detachment_theta=0.1 \
      --am_strain_detachment_max_exponent=20 \
      --am_strong_off_floor_fraction=1 \
      --am_cycle_period=0 --am_activation_floor=1 \
      --am_activation_active_fraction=1 --am_thin_activation=true \
      --am_thin_feedback=false --am_thick_feedback=false \
      --actin_turnover=false \
      --aa_barbed_lateral_range=0 --aa_barbed_lateral_k=0 \
      --aa_local_multilink=true --aa_preferred_centerline=0.0736 \
      --aa_capture_margin=0.045 --aa_antiparallel_overlap=0.10 \
      --aa_barbed_zone=0.20 --k_aa=${kaa} --kappa_aa=100 \
      --aa_perpendicular_k=50 --aa_status1_radial_k=75 \
      --aa_axial_register_ratio=3 --aa_axial_register_force_max=36 \
      --aa_status_mode=productive_register \
      --aa_status2_force_on=0.25 --aa_status2_force_off=0.10 \
      --aa_status2_registry_on=0.040 --aa_status2_registry_off=0.040 \
      --k_on=2500 --k_off=1 --base_lifetime=0.1 \
      --directional_base_lifetime=0.05 --lifetime_coeff=0.4 \
      --max_aa_bonds_per_actin=4 --aa_neighbor_identity=true \
      --aa_association_max_direction_dot=-0.5 \
      --directional=true --max_myosin_bonds=24 \
      --titin_mode=legacy_crosslinked --titin_k=0 --k_mm=0 --k_bundle_max=0 \
      --dt=1e-5 --nsteps=${nsteps} --save_every=${save_every} \
      --seed=${seed} --filename="${out}" \
      >"${root}/logs/${tag}.log" 2>&1
    echo "done ${tag}"
done
