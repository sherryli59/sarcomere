#!/bin/bash
#SBATCH --job-name=confprobe
#SBATCH --partition=rotskoff,owners
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=4:00:00
#SBATCH --requeue
#SBATCH --array=0-11
#SBATCH --output=/scratch/users/xsli/sarcomere/confined_probe_20260910/logs/%x_%A_%a.out
#SBATCH --error=/scratch/users/xsli/sarcomere/confined_probe_20260910/logs/%x_%A_%a.err
set -euo pipefail
export OMP_NUM_THREADS=1        # 2 filaments: threading is pure overhead, and
                                # 1 thread makes the model bitwise reproducible

root=/scratch/users/xsli/sarcomere/confined_probe_20260910
exe=${root}/bin/sarcomere_confined_20260910
mkdir -p "${root}"/{data,params,logs}
test -x "${exe}"
# Refuse to run a stale binary under a fresh campaign name.
case "$(basename "${exe}")" in *_20260910) ;; *) echo "exe/campaign date mismatch"; exit 3 ;; esac

# CONFINED ALIGNMENT PROBE
#
# Question: with NO alignment potential (am_alignment_k=0), does multi-head
# binding align the actin to the myosin, and is the operative variable the
# AXIAL SPAN of the bound heads rather than their number?
#
# The free probe could not answer this: with k_aa=0 and titin_k=0 nothing held
# the actin, so it escaped to 0.7-0.8 um against a 0.13 um retention reach and
# sat unbound 67-87% of frames.  Here:
#   * the myosin is frozen by setting its crowding factors to 1e9.  Crowding
#     divides the mobility (rod_mobility.h:77) and the integrator uses
#     sqrt(2 D dt) noise, so drift AND noise vanish together -- measured
#     residual motion 0.03 nm and 0.003 deg over 20k steps.
#   * the actin's centre is confined by flooring the occupancy multiplier on the
#     pair radial tether (am_tether_occupancy_floor).  That tether is harmonic
#     on the pair centerline distance with rest am_preferred_centerline and is
#     applied as a PURE CENTRE FORCE with no torque, so it fixes separation
#     without touching orientation -- leaving theta the only free rotational
#     degree of freedom.  k=100 pN/um gives perp median 76 nm (target 80),
#     max 247 nm, 86% of frames bound, 70% with >=2 heads, and force-cap
#     events on only 0.4% of steps.
#
# ARMS (predictions):
#   span3         hph=9  -> 3 crowns, spans 0/225/450 nm.  theta should relax.
#   span1         hph=3  -> crown_count=(3+2)/3=1, so ALL heads share one axial
#                           station and the span is identically 0.  If span is
#                           the operative variable, theta should NOT relax.
#                           Note this also cuts the head budget 9->3, so the
#                           analysis must condition on bound-head count.
#   notransverse  hph=9, am_head_transverse_k=0 -> if the couple comes from the
#                           radial head spring, the effect should vanish.
#   alignON       hph=9, am_alignment_k=120 -> positive control.
#
# Three initial tilts test convergence: theta should reach the same steady
# distribution from 10 deg (below) and 60 deg (above) if it is a real fixed point.

configs=(
  "span3_tilt10|9|100|0|10"    "span3_tilt30|9|100|0|30"    "span3_tilt60|9|100|0|60"
  "span1_tilt10|3|100|0|10"    "span1_tilt30|3|100|0|30"    "span1_tilt60|3|100|0|60"
  "notrans_tilt10|9|0|0|10"    "notrans_tilt30|9|0|0|30"    "notrans_tilt60|9|0|0|60"
  "alignON_tilt10|9|100|120|10" "alignON_tilt30|9|100|120|30" "alignON_tilt60|9|100|120|60"
)
IFS='|' read -r label hph transverse_k align_k tilt \
    <<< "${configs[$SLURM_ARRAY_TASK_ID]}"

nsteps=${CONF_PROBE_NSTEPS:-200000}          # 2 s at dt=1e-5
save_every=${CONF_PROBE_SAVE_EVERY:-50}      # 0.5 ms per frame -> 4001 frames

for seed in 1 2 3 4 5 6; do
    tag=${label}_seed${seed}
    out=${root}/data/${tag}.h5
    rm -f "${out}" "${out}.lock"
    printf '%s\n' "label=${label}" "am_heads_per_half=${hph}" \
        "am_head_transverse_k=${transverse_k}" "am_alignment_k=${align_k}" \
        "probe_tilt_deg=${tilt}" "confinement=am_tether_occupancy_floor=1.0,am_tether_k=100" \
        "myosin=frozen_via_crowding_1e9" "seed=${seed}" \
        >"${root}/params/${tag}.params"

    "${exe}" \
      --n_actins=1 --n_myosins=1 --Lx=4 --Ly=1 --Lz=1 \
      --dimension=3 --pbc_mask=111 --actin_length=1 --myosin_length=1.5 \
      --actin_radius=0.015 --myosin_radius=0.05 \
      --initial_structure=single_am_tilted --probe_tilt_deg=${tilt} \
      --actin_direction_noise=-1 --myosin_direction_noise=0 \
      --pre_phase_steps=0 --interaction_ramp_steps=0 \
      --physical_mobility=true --temperature_kelvin=310 \
      --thermal_noise_scale=1 --solvent_viscosity=0.001 \
      --actin_hydrodynamic_diameter=0.010 --myosin_hydrodynamic_diameter=0.035 \
      --actin_crowding_parallel=10 --actin_crowding_perpendicular=10 \
      --actin_crowding_tumble=160 --actin_crowding_roll=160 \
      --myosin_crowding_parallel=1e9 --myosin_crowding_perpendicular=1e9 \
      --myosin_crowding_tumble=1e9 --myosin_crowding_roll=1e9 \
      --max_actin_displacement=0.005 --max_myosin_displacement=0.005 \
      --max_actin_rotation=0.005 --max_myosin_rotation=0.005 \
      --sterics=true --actin_steric_k=2000 --myosin_steric_k=1000 --am_steric_k=1000 \
      --am_cutoff=0.05 --am_optimal=0.055 --am_surface_gap=0.015 \
      --k_am=0 --kappa_am=0 --dynamic_am=true \
      --am_offcenter_torque=true --am_spatial_heads=true \
      --am_heads_per_half=${hph} \
      --am_max_heads_per_actin_half=6 --am_max_heads_per_actin=48 \
      --am_sites_per_actin=182 --am_k_on=100 --am_weak_binding_floor=1 \
      --am_k_weak_off=50 --am_k_powerstroke=200 --am_k_strong_off=160 \
      --am_force_scale=1 --am_motor_k=1000 --am_power_stroke=0.005 \
      --am_stall_force=4.5 --am_unloaded_velocity=5 --am_weak_stiffness_ratio=0.1 \
      --am_tether_k=100 --am_tether_slack=0.002 \
      --am_tether_occupancy_floor=1.0 --am_tether_occupancy_floor_reach=1.0 \
      --am_head_transverse_k=${transverse_k} \
      --am_alignment_k=${align_k} --am_alignment_slack_deg=2 \
      --am_tether_saturation_heads=6 \
      --am_binding_azimuth_selection=true --am_binding_azimuth_half_width_deg=30 \
      --am_rupture_distance=0.1 --am_strain_dependent_detachment=true \
      --am_strain_detachment_theta=0.1 --am_strain_detachment_max_exponent=20 \
      --am_strong_off_floor_fraction=1 \
      --am_cycle_period=0 --am_activation_floor=1 \
      --am_activation_active_fraction=1 --am_thin_activation=true \
      --am_thin_feedback=false --am_thick_feedback=false --actin_turnover=false \
      --aa_barbed_lateral_range=0 --aa_barbed_lateral_k=0 \
      --aa_local_multilink=true --aa_preferred_centerline=0.0736 \
      --aa_capture_margin=0.045 --aa_antiparallel_overlap=0.10 \
      --aa_barbed_zone=0.20 --k_aa=0 --kappa_aa=100 \
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
