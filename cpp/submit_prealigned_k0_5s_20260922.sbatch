#!/bin/bash
# normal,hns rather than rotskoff: rotskoff had 152/160 CPUs allocated (8 idle)
# against a 12 x 32 = 384 CPU ask, while normal+hns had ~800 idle.  Both allow
# 7 days, same as rotskoff, and NEITHER preempts -- so this is also safer than
# adding owners back, where a requeue restarts from step 0 and overwrites the
# incrementally written trajectory.
#SBATCH --partition=normal,hns
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
# 16GB, not 64GB.  MEASURED peak RSS of the sarcomere step (.0, not .batch)
# in previous A1600/M320 production runs: 0.64 GB and 3.4 GB over a 21.6 h
# run.  The old 64 GB ask was ~19x the observed peak and was what stopped the
# job being placed -- owners had 4259 idle CPUs but few nodes with 32 free
# cores AND 64 GB free at once.  The torque breakdown adds only
# 12 sources x 32 threads x 1600 actins x 3 doubles = ~15 MB.
#SBATCH --mem=16GB
# owners caps wall time at 2 days; ~15 h expected for 500k steps at A1600.
#SBATCH --time=2-00:00:00
#SBATCH --open-mode=append
#SBATCH --no-requeue
# 120 s warning to the batch script before preemption/wall clock,
# so the self-resubmit trap below has time to run.
#SBATCH --signal=B:USR1@120
#SBATCH --array=0-11
#SBATCH --job-name=prealign
#SBATCH --output=/scratch/users/xsli/sarcomere/prealigned_k0_5s_20260922/logs/%x_%A_%a.out
#SBATCH --error=/scratch/users/xsli/sarcomere/prealigned_k0_5s_20260922/logs/%x_%A_%a.err

#
# CORRECTED 2026-09-22: the first attempt seeded every filament about +x,
# which is POLAR order.  MEASURED consequence: polar fraction 1.000, ZERO
# antiparallel pairs, and only 0.059 ACTN2 links/actin where an isotropic
# start reaches 4.05.  ACTN2 only links pairs with direction_dot <= -0.5
# (aa_association_max_direction_dot), so there was nothing for it to bind and
# the runs could not form a Z-disc at all.  Nematic order S is sign-blind, so
# polar_fraction = 0.5 gives the SAME S -- verified 0.9537 either way -- with
# the antiparallel population restored (0.0000 -> 0.5007).  Those 12 runs were
# cancelled and their trajectories deleted.
#
# PRE-ALIGNED START, am_alignment_k = 0.  Can the intrinsic couple HOLD
# order it did not have to create?
#
# WHY.  MEASURED: from an isotropic start with k=0, S reaches only 0.169 after
# 5 s (still climbing at +0.030/s, never crossing 0.9), while k=120 reaches
# 0.97 in 0.25 s.  That says the couple cannot BUILD order on this timescale.
# It does not say the couple cannot MAINTAIN it -- a different question, and
# the one that matters if alignment is established by some other route.
#
# THE KNOB.  actin_direction_noise / myosin_direction_noise are NOT angles:
# set_direction_axis_noise adds N(0, sigma) to each Cartesian component of
# (1,0,0) and renormalises.  Calibrated by direct sampling (200k draws):
#     sigma   rms angle    initial S
#     0.125     10.2 deg      0.954
#     0.300     24.7 deg      0.757
#     1.000     54.0 deg      0.173   <- production's myosin value
# Note the last line: production initialises myosin at S ~ 0.17, and actin
# ISOTROPICALLY (actin_direction_noise = -1, S ~ 0.011).
#
# THE ARMS.
#   sigma 0.125 -> S0 = 0.954, deliberately AT the k=120 steady state
#                  (0.957-0.984).  If the couple can hold order, S stays here;
#                  if not, it decays toward the 0.169 the isotropic k=0 run
#                  reached.  That is a clean binary outcome.
#   sigma 0.300 -> S0 = 0.757, partially ordered.  Together with the first arm
#                  this tests whether there is a threshold order above which
#                  the couple is self-sustaining and below which it is not.
#
# Both actin AND myosin are seeded about +x, so this is not a myosin-scaffold
# experiment: the whole system starts aligned.
#
# 2 sigmas x 6 seeds = 12 tasks, titin ON, am_alignment_k = 0 throughout.
#
# PRODUCTION on the current binary (commit e0f003f), with corrected titin
# release, titin ON vs OFF.
#
# TITIN.  titin_k_off_unsupported was 0, which per its own help text falls back
# to the single global titin_k_off = 0.01/s -- a 100 s mean lifetime, so titin
# never released and pinned every thin filament after its Z-disc died (measured
# bound in 99.9% of frames at 8-15 pN against 10-13 pN of motor force).
#   k_off_unsupported = 10 /s.  Not tuned: the parameter models the serial
#   actin-ACTN2-titin anchor, "whose composite lifetime is set by its weakest
#   element", and 10 /s = 1/base_lifetime, the ACTN2 link lifetime titin is in
#   series with.  The resulting 100 ms lifetime also matches the MEASURED
#   65-135 ms Z-disc death window.
#   titin_rebind_requires_actn2 false -> true.  REQUIRED, or the fix is inert:
#   the same-partner recapture path (sarcomere.cpp:4390) is deliberately exempt
#   from titin_require_live_actn2 and is gated by this flag instead, so a
#   released titin re-grabs the same dead-Z-disc actin inside the 20 ms
#   rebind_memory window.
#
# NOTE, not changed here: titin_strain_gated_kinetics=true is a NO-OP at
# titin_transition_theta=1, because dissociation_rate uses
# exponent = beta*(1-theta)*E, which is identically 0 at theta=1.  Load-
# dependent release would need theta < 1.  Left alone so this campaign changes
# one thing about titin chemistry, not two.
#
# MODEL CHANGES since the last production run:
#   am_head_azimuth=false   azimuth-off radial channel; the head anchors on the
#                           myosin AXIS against am_preferred_centerline, so the
#                           azimuth can no longer strain a head by geometry.
#                           The two old radial springs were in incompatible
#                           frames: net +0.15 pN outward in 6/6 runs.
#   am_head_radial_force_max  DERIVED from the displacement clamp; at dt=1e-5
#                           that is 1.92 pN/head (it scales as 1/dt, so it is
#                           HALF the 3.85 pN of the dt=5e-6 rig it was
#                           validated on -- watch the radial cap counters).
#   actin torque clamp      now budgets the component perpendicular to the
#                           filament axis.  The old clamp fired on 99.9% of
#                           actin-steps and discarded 50.2% of the real
#                           rotational drive.
#   kappa_aa 100 -> 3       CAVEAT: measured at am_alignment_k = 0; this run
#                           has am_alignment_k = 120, a different angular
#                           landscape, so that evidence does not transfer
#                           cleanly.  Flagged, not resolved.
#   rest 73.6 -> 79 nm      ACTN2 fought aa_barbed_lateral_range, a one-sided
#                           wall at 80 nm and 6x stiffer.  aa_capture_margin is
#                           trimmed 0.045 -> 0.0396 so aa_cutoff stays 0.1186
#                           and capture range is NOT a confound.
#
# Lx = 6 is MANDATORY and is honoured here (A1600/M320).
#
# Partition rotskoff only and --no-requeue: a preempted owners job restarts
# from step 0 and overwrites the incrementally written trajectory.
#
# COST, measured from the previous production runs: 25.8 h and 15.7 GB per
# seed at 32 cores.  2 arms x 6 seeds = 12 runs, ~9,900 core-hours.
#
set -euo pipefail
# Lx=6 is mandatory.  The model selects a 2.00 um sarcomere period (dominant
# harmonic k=3 in Lx=6).  Lx=3 would need k=1.5, which PBC cannot host, so the
# 3x2x2 box shows a flat band spectrum and no axial order.  Density is held at
# 66.7 actin/um^3 and 13.3 myosin/um^3, so Lx=6 implies A1600/M320.
export HDF5_USE_FILE_LOCKING=FALSE
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-32}

repo=/home/groups/rotskoff/sherry/sarcomere/sarcomere
root=${ASSEMBLY_PILOT_ROOT:-/scratch/users/xsli/sarcomere/prealigned_k0_5s_20260922}
exe=${root}/bin/sarcomere_prealign_20260922
monitor=${root}/source_snapshot/monitor_live_production_hdf5.py
mkdir -p "${root}"/{data,params,analysis,logs}

# ACTN2 rest length raised 0.055 -> 0.0736 um = 0.92 * am_preferred_centerline
# (lattice_spacing.h physiological_aa_am_ratio).  aa_capture_margin stays
# 0.045, so aa_cutoff grows 0.100 -> 0.1186 as an intended consequence.
# Rationale: the lateral repulsion stalls at ~46 nm because ACTN2 pulls back
# to its 55 nm rest length with ~40 pN, exactly matching the lateral cap.
# `arm` selects THICK-FILAMENT REGULATION and is consumed by the
# `case "${arm}"` dispatch below.  It is held at the winning factorial cell.
# Titin is an ORTHOGONAL factor and must not be folded into it: doing so made
# every task fall through that case to `exit 2`, and a permissive default
# would have been worse -- it would have silently run the fullon regulation.
arms=(intermediate)
arm=${arms[0]}
titin_arms=(on)
sigmas=(0.125 0.30)
seeds=(20260805 20260815 20260825 20260835 20260845 20260855)
task=${SLURM_ARRAY_TASK_ID}
titin_arm=on
sigma=${sigmas[$(( task / 6 ))]}
seed=${seeds[$(( task % 6 ))]}
# Force-only titin ablation: titin_k=0 and titin_links_per_half=0 are both
# rejected by the validator, so the off arm keeps the kinetics and the code
# path live and zeroes the force.  That makes this a one-parameter A/B.
if [ "${titin_arm}" = "off" ]; then titin_k_arm=1e-9; else titin_k_arm=100; fi

# Background held fixed at the winning factorial cell: ACTN2 rest length
# 73.6 nm (= 0.92 * am_preferred_centerline) plus barbed-end lateral repulsion.
# Only the thick-filament regulation differs between arms.
turnover=false
turnover_rate=0
clearance_ratio=0.80
lateral_range=0.080
lateral_k=2000
lateral_fmax=40

# am_thick_feedback=false initialises every head to `detached` and promotes any
# head found in off_srx, so it is full-ON and recruit/fold/initial_on are inert
# for that arm (sarcomere.cpp ~2331 and ~2922).  Verified, not assumed.
thick_feedback=false
recruit=2
fold=18
initial_on=0.1
# Nominal steady-state ON fraction is recruit/(recruit+fold): at zero tension
# the recruitment rate is `recruit` (am_motor_kinetics.h thick_recruitment_rate,
# exponentially accelerated by tension) and the folding rate is `fold`
# (load-slowed).  recruit+fold is held at 20 /s across every regulated arm, as
# in the original ladder, so the arms differ in availability and not in
# switching speed.
#
#   hcm          10/(10+10) = 0.50
#   intermediate 15/(15+ 5) = 0.75   <- midway between hcm and full-ON
#   fullon                    1.00
#
# The first version of this file reused the old ladder's intermediate
# (6/14 = 0.30), which sits BELOW hcm rather than between hcm and full-ON, and
# left the 0.50-1.00 interval unsampled -- the interval where full-ON and HCM
# were measured to differ (S_myo 0.947 vs 0.786, 5.78 vs 3.05 bound heads per
# filament).
case "${arm}" in
    fullon) ;;
    intermediate) thick_feedback=true; recruit=15; fold=5;  initial_on=0.75 ;;
    hcm)          thick_feedback=true; recruit=10; fold=10; initial_on=0.5 ;;
    *) echo "FATAL: unknown regulation arm '${arm}'; the case dispatch sets thick_feedback/recruit/fold/initial_on and falling through would run the wrong regulation" >&2; exit 2 ;;
esac

dt=1e-5
# Production duration per the plan: 5 s at dt=1e-5, 1000 saved frames,
# matching the A2400 Lx=6 reference cadence (~12 GB per trajectory).
nsteps=${ASSEMBLY_PILOT_NSTEPS:-500000}
save_every=${ASSEMBLY_PILOT_SAVE_EVERY:-500}
turnover_start=${ASSEMBLY_TURNOVER_START:-0.2}
duration_ms=$((nsteps / 100))
tag=A1600_M320_Lx6_${arm}_sig${sigma}_seed${seed}_${duration_ms}ms
output=${root}/data/${tag}.h5
params=${root}/params/${tag}.params
live_log=${root}/analysis/${tag}.live.jsonl

test -x "${exe}"
test -r "${monitor}"
# A task preempted before its first checkpoint leaves a file with zero
# committed frames.  Resuming from that aborts the binary with
#   "No frames available in actin center dataset for resume."
# so decide on committed frames, not on mere file existence.
committed=0
if [[ -e "${output}" ]]; then
    committed=$(python - "${output}" <<'PYFRAMES'
import sys
try:
    import h5py
    with h5py.File(sys.argv[1], "r", locking=False) as handle:
        n = int(handle["meta/committed_frames"][0]) if "meta/committed_frames" in handle else 0
        m = handle["actin/center"].shape[0] if "actin/center" in handle else 0
        print(min(n, m))
except Exception:
    print(0)
PYFRAMES
)
fi
if [[ -e "${output}" && "${committed}" -gt 0 && "${SLURM_RESTART_COUNT:-0}" -eq 0 && "${ASSEMBLY_PILOT_RESUME_EXISTING:-false}" != true ]]; then
    echo "refusing to overwrite ${output}" >&2
    exit 1
fi
resume=()
if [[ "${committed}" -gt 0 ]]; then
    resume=(--resume=true)
elif [[ -e "${output}" ]]; then
    echo "discarding zero-frame ${output} (preempted before first checkpoint)" >&2
    rm -f "${output}"
fi
exec 9>"${output}.lock"
flock -n 9

printf '%s\n' \
    "arm=${arm}" "seed=${seed}" "turnover_rate=${turnover_rate}" \
    "thick_feedback=${thick_feedback}" "recruit=${recruit}" "fold=${fold}" \
    "initial_on=${initial_on}" "tension_scale=3" \
    "clearance_ratio=${clearance_ratio}" \
    "lateral_range=${lateral_range}" "lateral_k=${lateral_k}" "lateral_fmax=${lateral_fmax}" \
    "turnover_start=${turnover_start}" >"${params}"

cmd=(srun --kill-on-bad-exit=1 "${exe}"
    --n_actins=1600 --n_myosins=320 --Lx=6 --Ly=2 --Lz=2
    --dimension=3 --pbc_mask=111 --actin_length=1 --myosin_length=1.5
    --actin_radius=0.015 --myosin_radius=0.05
    --initial_structure=myosin_x_noise --actin_direction_noise="${sigma}"
    --myosin_direction_noise="${sigma}" --myosin_direction_tilt_deg=0
    --actin_direction_polar_fraction=0.5
    --myosin_direction_polar_fraction=0.5
    --pre_phase_steps=500 --interaction_ramp_steps=2000
    --fixed_lattice_preequil_steps=0 --physical_mobility=true
    --temperature_kelvin=310 --thermal_noise_scale=1 --solvent_viscosity=0.001
    --actin_hydrodynamic_diameter=0.010 --myosin_hydrodynamic_diameter=0.035
    --actin_crowding_parallel=10 --actin_crowding_perpendicular=10
    --actin_crowding_tumble=160 --actin_crowding_roll=160
    --myosin_crowding_parallel=20 --myosin_crowding_perpendicular=20
    --myosin_crowding_tumble=64 --myosin_crowding_roll=64
    --max_actin_displacement=0.005 --max_myosin_displacement=0.005
    --max_actin_rotation=0.005 --max_myosin_rotation=0.005
    --sterics=true --actin_steric_k=2000 --myosin_steric_k=1000
    --am_steric_k=1000 --myosin_exclusion_distance=0
    --am_cutoff=0.05 --am_optimal=0.055 --am_surface_gap=0.015
    --k_am=0 --kappa_am=0 --dynamic_am=true --am_offcenter_torque=true
    --am_head_azimuth=false
    # NO torque breakdown and NO force breakdown in production.  The
    # breakdowns are diagnostics: they add ~37 extra per-actin datasets
    # (12 sources x 3 components + tb_total_precap) and are for probe
    # rigs, not production trajectories.  Both default to off.
    --am_spatial_heads=true --am_heads_per_half=9 --am_crossover_pull_through=true
    --am_max_heads_per_actin_half=24 --am_sites_per_actin=182
    --am_max_heads_per_actin=24 --am_k_on=100 --am_weak_binding_floor=0.1
    --am_k_weak_off=200 --am_k_powerstroke=200 --am_k_strong_off=160
    --am_force_scale=1 --am_motor_k=1000 --am_power_stroke=0.005
    --am_stall_force=4.5 --am_unloaded_velocity=5 --am_weak_stiffness_ratio=0.1
    --am_tether_k=3000 --am_tether_slack=0.002 --am_head_transverse_k=100
    --am_transverse_strain_detachment_theta=0 --am_alignment_k=0
    --am_alignment_slack_deg=2 --am_tether_saturation_heads=6
    --am_binding_azimuth_selection=true --am_binding_azimuth_half_width_deg=30
    --am_rupture_distance=0.1 --am_strain_dependent_detachment=true
    --am_strain_detachment_theta=0.1 --am_strain_detachment_max_exponent=20
    --am_strong_off_floor_fraction=1 --am_cycle_period=0 --am_activation_floor=0
    --am_activation_active_fraction=0.5 --am_thin_activation=true
    --am_thin_feedback=true --am_q_calcium_open_rate=20 --am_q_calcium_half=0.5
    --am_q_calcium_hill=3 --am_q_strong_feedback_rate=100
    --am_q_strong_half=0.02 --am_q_strong_hill=2 --am_q_close_rate=20
    --am_initialize_q_calcium_equilibrium=true --am_thick_feedback="${thick_feedback}"
    --am_thick_recruit_rate="${recruit}" --am_thick_fold_rate="${fold}"
    --am_thick_tension_scale=3 --am_thick_initial_on_fraction="${initial_on}"
    --am_status2_memory=0 --am_bidirectional_crossbridge=false
    --am_overlap_target_ratio=0 --actin_turnover="${turnover}"
    --actin_turnover_rate="${turnover_rate}" --actin_turnover_start_time="${turnover_start}"
    --actin_turnover_max_placement_attempts=5000
    --actin_turnover_clearance_ratio="${clearance_ratio}"
    --aa_barbed_lateral_range="${lateral_range}"
    --aa_barbed_lateral_k="${lateral_k}"
    --aa_barbed_lateral_force_max="${lateral_fmax}"
    --aa_barbed_lateral_axial_window=0.250
    --aa_local_multilink=true --aa_spacing_ratio=0.9875 --aa_preferred_centerline=0.079
    --aa_capture_margin=0.0396 --aa_antiparallel_overlap=0.10 --aa_barbed_zone=0.20
    --aa_axial_register_ratio=3 --aa_axial_register_force_max=36
    --aa_register_all_statuses=false --aa_status2_registry_on=0.040
    --aa_status2_registry_off=0.040 --aa_status2_support_window=0.010
    --aa_register_saturation_error=0.05 --aa_prebundle_links_per_pair=3
    --aa_site_exclusion=0.0192 --aa_overlap_capture=0.10 --k_aa=1000
    --kappa_aa=3 --aa_perpendicular_k=50 --aa_status1_radial_k=75
    --aa_status1_site_diffusion=0.010 --aa_arm_state_model=false
    --aa_actn2_contact_radius=0.005 --aa_myosin_contact_epsilon=0.00428
    --aa_myosin_contact_force_cap=40 --aa_myosin_contact_shell_width=0.015
    --aa_myosin_contact_shell_energy=0.0214 --aa_occlusion_rupture_penetration=0.010
    --aa_zdisc_tip_contact=true --aa_zdisc_tip_onset=0.050
    --aa_zdisc_tip_stiffness=1000 --aa_zdisc_tip_force_max=20
    --aa_zdisc_tip_transverse_flat=0.100 --aa_zdisc_tip_transverse_cutoff=0.130
    --aa_status2_locked_until_detach=true --aa_stiffness_reference_links=3
    --aa_total_link_pool=0 --aa_status_mode=productive_register
    --aa_status2_force_on=0.25 --aa_status2_force_off=0.10 --k_on=2500
    --k_off=1 --base_lifetime=0.1 --directional_base_lifetime=0.05
    --lifetime_coeff=0.4 --stabilize_status1=false
    --aa_energy_gated_association=false --aa_transition_state_theta=0.1
    --aa_encounter_k=10 --aa_encounter_kappa=0 --aa_passive_detailed_balance=false
    --tau_rec=0 --max_aa_bonds_per_actin=4 --aa_neighbor_identity=true
    --aa_neighbor_slot_memory=0.02 --aa_status2_slot_memory=0.20
    --aa_association_max_direction_dot=-0.5 --aa_link_rupture_length=0.12
    --aa_double_overlap_contact=true --aa_double_overlap_contact_diameter=0.012
    --aa_double_overlap_contact_epsilon=0.00428 --aa_double_overlap_contact_force_cap=0
    --directional=true --max_myosin_bonds=24 --titin_mode=explicit_dynamic
    --titin_k="${titin_k_arm}" --titin_rest_length=0.665 --titin_links_per_half=3
    --titin_capture_distance=0.72 --titin_k_on=5000 --titin_k_off=0.01
    --titin_strain_gated_kinetics=true --titin_transition_theta=1
    --titin_rebind_memory=0.020 --titin_rebind_rate_multiplier=100
    --titin_require_live_actn2=true --titin_k_off_unsupported=10
    --titin_rebind_requires_actn2=true --k_mm=0 --k_bundle_max=0
    --bundle_ramp_steps=0 --dt="${dt}" --nsteps="${nsteps}"
    --save_every="${save_every}" --seed="${seed}" --filename="${output}"
    ${resume[@]+"${resume[@]}"})

"${cmd[@]}" &
simulation_pid=$!
python "${monitor}" "${output}" "${live_log}" --pid "${simulation_pid}" \
    --initial-seconds 60 --interval-seconds 300 \
    --actin-force-cap-fraction-max 0.75 --myosin-force-cap-fraction-max 0.50 &
monitor_pid=$!
# --- self-resubmitting chain --------------------------------------------
# Sherlock rejects sleeper jobs, so a polling watcher is not permitted.  Each
# task supervises ITSELF instead: --signal=B:USR1@120 gives the batch script
# 120 s of warning before preemption or wall clock, and the trap below asks
# whether this task's OWN trajectory has all its frames.  Completeness is read
# from committed frames in the HDF5, never from the SLURM state, because a
# preempted task usually holds a good partial file and the launcher resumes
# from committed frames rather than restarting.  BABYSIT_RETRY rides along in
# the environment so a task that keeps dying stops after BABYSIT_MAX_RETRY
# instead of cycling on the allocation forever.
BABYSIT_RETRY=${BABYSIT_RETRY:-0}
BABYSIT_MAX_RETRY=${BABYSIT_MAX_RETRY:-5}
# The requeue pool is defined HERE, not passed through --export: a
# comma-separated value inside --export is split on its commas and each piece
# is treated as another variable name to export, which corrupts the submission.
BABYSIT_PARTITION=${BABYSIT_PARTITION:-normal,hns,rotskoff,owners}
CHAIN_LOG=${root}/logs/chain.log
CHAIN_SELF=${root}/source_snapshot/$(basename "$0")

trajectory_complete() {
    python3 - "$output" <<'PYDONE'
import os, sys
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
try:
    import h5py
    with h5py.File(sys.argv[1], "r", locking=False) as fh:
        cf = fh["meta/committed_frames"][()]
        cf = int(cf.ravel()[0]) if hasattr(cf, "ravel") else int(cf)
    sys.exit(0 if cf >= 1001 else 1)
except Exception:
    sys.exit(1)
PYDONE
}

resubmit_if_needed() {
    local tag="task ${SLURM_ARRAY_TASK_ID} (${arm}/${titin_arm} seed ${seed})"
    if trajectory_complete; then
        echo "[$(date '+%F %T')] ${tag}: COMPLETE" >> "$CHAIN_LOG"; return 0
    fi
    if [ "${BABYSIT_RETRY}" -ge "${BABYSIT_MAX_RETRY}" ]; then
        echo "[$(date '+%F %T')] ${tag}: INCOMPLETE, retry budget ${BABYSIT_MAX_RETRY} exhausted -- needs a human" >> "$CHAIN_LOG"
        return 0
    fi
    local nxt=$(( BABYSIT_RETRY + 1 ))
    local new
    # Requeue into the FULL pool, not just the partition that evicted us.
    local parts=${BABYSIT_PARTITION:-${SLURM_JOB_PARTITION}}
    new=$(sbatch --parsable --partition="${parts}" --no-requeue \
            --time=2-00:00:00 --mem=16G --cpus-per-task=32 \
            --array="${SLURM_ARRAY_TASK_ID}" \
            --export=ALL,ASSEMBLY_PILOT_ROOT="${root}",BABYSIT_RETRY=${nxt},BABYSIT_MAX_RETRY=${BABYSIT_MAX_RETRY},BABYSIT_PARTITION="${parts}" \
            --output="${root}/logs/%x_%A_%a.out" \
            --error="${root}/logs/%x_%A_%a.err" \
            "${CHAIN_SELF}" 2>>"$CHAIN_LOG")
    echo "[$(date '+%F %T')] ${tag}: INCOMPLETE on partition ${SLURM_JOB_PARTITION} -> resubmitted as ${new:-SBATCH_FAILED} to [${parts}] (retry ${nxt}/${BABYSIT_MAX_RETRY})" >> "$CHAIN_LOG"
}

on_usr1() {
    echo "[$(date '+%F %T')] task ${SLURM_ARRAY_TASK_ID}: USR1, preempt/timeout in ~120 s" >> "$CHAIN_LOG"
    kill -TERM "${simulation_pid}" 2>/dev/null
    wait "${simulation_pid}" 2>/dev/null
    resubmit_if_needed
    exit 0
}
trap on_usr1 USR1

set +e
wait "${simulation_pid}"
status=$?
wait "${monitor_pid}"
set -e
resubmit_if_needed
exit "${status}"

