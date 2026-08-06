#!/bin/bash
# Guarded runs for the explicit-titin / turnover-capable ACTN2 model with the
# ACTN2-support capture gate (--titin_require_live_actn2, titin schema 3).
#
# This is the v4 production script with three changes:
#   1. --titin_require_live_actn2=true, the formation-only capture specificity
#      rule that replaces nonspecific slack-region capture;
#   2. system size, box length, seed and step count are parameters, so the
#      density-matched occupancy diagnostics and the 1,280-actin pilot are the
#      same code path with the same frozen binary; and
#   3. the frozen executable hash points at the schema-3 build.
#
# Density is held at the production value of 41.451 actins/um^3 by scaling Lx
# with the actin count at the fixed 2 x 2 um cross-section. It still refuses to
# extend itself: PROD_NSTEPS is capped at 40000 steps = 20 ms at dt = 0.5 us.
#SBATCH --partition=hns,normal,owners
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32GB
#SBATCH --time=04:00:00
#SBATCH --job-name=titin_serial
#SBATCH --open-mode=append
#SBATCH --no-requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/prod_titin_serial_20260802/logs/%x_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/prod_titin_serial_20260802/logs/%x_%j.err

set -euo pipefail

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export OMP_DYNAMIC=FALSE
export OMP_SCHEDULE=static
export OPENBLAS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export MKL_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE
JEMALLOC=/home/groups/rotskoff/sherry/mamba_root/envs/sarcomere/lib/libjemalloc.so
test -r "${JEMALLOC}" || { echo "missing jemalloc: ${JEMALLOC}" >&2; exit 1; }
export LD_PRELOAD="${JEMALLOC}"

ROOT="${PROD_ROOT:-/scratch/users/xsli/sarcomere/prod_titin_serial_20260802}"
EXE="${PROD_EXE:-${ROOT}/bin/sarcomere_5374a633b1cb}"
EXPECTED_EXE_SHA256=5374a633b1cb6151c676f9e868409510f004f5da664800368164e1d5d04edea0

N_ACTINS="${PROD_N_ACTINS:-1280}"
N_MYOSINS="${PROD_N_MYOSINS:-640}"
LX="${PROD_LX:-7.72}"
SEED="${PROD_SEED:-2080201}"
NSTEPS="${PROD_NSTEPS:-40000}"
SAVE_EVERY="${PROD_SAVE_EVERY:-2000}"
TAG="${PROD_TAG:-A${N_ACTINS}_M${N_MYOSINS}_serial_seed${SEED}}"
# Dissociation rate applied while the serial actin-ACTN2-titin anchor has no
# live middle element. Order 100 /s puts the composite lifetime near the
# measured ~6 ms ACTN2 encounter lifetime instead of titin's own 2 s.
KOFFU="${PROD_KOFF_UNSUPPORTED:-100}"
# Sliding throttle. With it disabled the engaged A-M overlap grows without
# bound, myosin occupies the whole thin filament, and the crosslinkable
# segment left for ACTN2 collapses to zero -- measured median 0 by 15 ms.
OVLTGT="${PROD_OVERLAP_TARGET:-0}"
# Motor drive. The passing minimal cell ran 35 heads/half at activation
# floor 0 and active fraction 0.5; production drifted to 141 at floor 1 and
# fraction 1, roughly 8x the drive. Measured consequence: net axial actin
# drift 2.18 um/s in production vs 0.011 um/s in the minimal cell.
HEADS="${PROD_HEADS:-141}"
ACTFLOOR="${PROD_ACT_FLOOR:-1}"
ACTFRAC="${PROD_ACT_FRAC:-1}"
# AM nematic angular tether. Reciprocal (torques actin AND myosin) and
# polarity-blind, so it aligns axes without forcing filaments parallel --
# antiparallel Z-disc geometry survives. Its magnitude is scaled by
# min(1, bound_heads / am_tether_saturation_heads), so cutting heads 141->35
# silently cut this ensemble stiffness too; raising k restores it.
ALIGNK="${PROD_ALIGN_K:-30}"
# Myomesin M-band spring between transversely adjacent thick filaments.
# Draws their centres into a common axial plane. Was dead code until
# 2026-08-03 (forces computed, never applied), so every earlier run had
# k_mm effectively 0 regardless of what was passed.
KMM="${PROD_K_MM:-0}"
OUTPUT="${PROD_OUTPUT:-${ROOT}/data/${TAG}.h5}"
PARAMS="${PROD_PARAMS:-${ROOT}/params/${TAG}.params}"
mkdir -p "${ROOT}/data" "${ROOT}/logs" "${ROOT}/params" "${ROOT}/monitor"

test -x "${EXE}" || { echo "missing executable ${EXE}" >&2; exit 1; }
exe_sha256=$(sha256sum "${EXE}" | awk '{print $1}')
if [ "${exe_sha256}" != "${EXPECTED_EXE_SHA256}" ]; then
    echo "frozen executable hash mismatch: ${exe_sha256}" >&2
    exit 1
fi
if [ "${NSTEPS}" -gt 40000 ]; then
    echo "guarded run refuses PROD_NSTEPS=${NSTEPS}; maximum is 40000" >&2
    exit 1
fi
if [ -e "${OUTPUT}.failed_gate" ]; then
    echo "checkpoint is marked failed; refusing to resume ${OUTPUT}" >&2
    exit 1
fi
# The ratio and density are the scientific content of the diagnostic ladder, so
# a mistyped Lx must stop the job rather than silently produce an off-density
# trajectory that would be compared against the pilot.
python - "${N_ACTINS}" "${N_MYOSINS}" "${LX}" <<'PY'
import sys
n_actin, n_myosin, lx = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3])
if n_actin != 2 * n_myosin:
    raise SystemExit(f"actin:myosin must be 2:1, got {n_actin}:{n_myosin}")
density = n_actin / (lx * 2.0 * 2.0)
if abs(density - 41.4507772) > 1e-3:
    raise SystemExit(
        f"actin density {density:.4f}/um^3 does not match the production "
        f"41.4508/um^3 at Lx={lx}")
# Minimum image: the neighbour cutoff is max(actin,myosin length) + 0.10 plus a
# 15% skin = 1.84 um, and the titin capture shell is 0.85 um.
if lx <= 3.68:
    raise SystemExit(f"Lx={lx} is below the 3.68 um minimum-image floor")
PY

exec 9>"${OUTPUT}.lock"
flock -n 9 || { echo "another writer owns ${OUTPUT}" >&2; exit 0; }

resume_flag=()
if [ -s "${OUTPUT}" ] && python - "${OUTPUT}" <<'PY'
import h5py
import sys
with h5py.File(sys.argv[1], "r", locking=False) as h:
    ok = (
        ("/resume/current_step" in h or "/resume_prev/current_step" in h)
        and "/actin/center" in h
        and h["/actin/center"].shape[0] > 0
        and int(h["/meta"].attrs.get("titin_mode", -1)) == 2
        # Schema 2 trajectories were produced without the capture gate and
        # saturated; the binary rejects them too, but fail here with a clear
        # message rather than inside the resume guard.
        and int(h["/meta"].attrs.get("titin_schema_version", 0)) == 4
        and int(h["/meta"].attrs.get("titin_require_live_actn2", 0)) == 1
    )
raise SystemExit(0 if ok else 1)
PY
then
    echo "Found a valid schema-3 gated checkpoint; resuming toward ${NSTEPS} total steps."
    resume_flag=(--resume=true)
elif [ -e "${OUTPUT}" ]; then
    preserved="${OUTPUT}.incomplete_job${SLURM_JOB_ID}"
    echo "Preserving non-resumable output as ${preserved}" >&2
    mv "${OUTPUT}" "${preserved}"
fi

git_head=$(git -C /home/groups/rotskoff/sherry/sarcomere/sarcomere rev-parse HEAD 2>/dev/null || echo unknown)
cat >"${PARAMS}" <<EOF
model=explicit_crossbridge_relaxed_actn2_explicit_titin_serial_anchor_v7
purpose=capture_specificity_occupancy_gate
n_actins=${N_ACTINS}
n_myosins=${N_MYOSINS}
box=${LX}x2x2
actin_density_per_um3=$(python -c "print(${N_ACTINS}/(${LX}*4.0))")
seed=${SEED}
dt=5e-7
nsteps=${NSTEPS}
save_every=${SAVE_EVERY}
pre_phase_steps=10000
interaction_ramp_steps=40000
activation=tonic
actn2_links_per_pair=2
actn2_target_overlap_um=0.1
actn2_barbed_association_zone_um=0.20
actn2_status1_capture_distance_um=0.10
actn2_status1_radial_k_pN_per_um=75
actn2_status2_centerline_um=0.055
actn2_contour_length_um=0.12
titin_topology=explicit_thick_tip_domains_to_dynamic_actin_barbed_ends
titin_domains_per_half=3
titin_rest_length_um=0.75
titin_capture_distance_um=0.85
titin_k_pN_per_um=40
titin_k_on_per_s=5000
titin_k_off_per_s=0.5
titin_force_law=tension_only_flat_bottom
titin_transition_theta=1
titin_rebind_memory_s=0.020
titin_rebind_rate_multiplier=100
titin_require_live_actn2=true
titin_k_off_unsupported_per_s=${KOFFU}
titin_anchor_topology=serial_actin_actn2_titin_composite_lifetime
titin_rebind_requires_actn2=true
am_overlap_target_ratio=${OVLTGT}
am_heads_per_half=${HEADS}
am_activation_floor=${ACTFLOOR}
am_activation_active_fraction=${ACTFRAC}
am_alignment_k=${ALIGNK}
k_mm=${KMM}
titin_capture_gate=formation_requires_live_actn2_on_barbed_end_any_status
titin_schema_version=4
executable=${EXE}
executable_sha256=${exe_sha256}
git_head=${git_head}
slurm_job_id=${SLURM_JOB_ID}
omp_num_threads=${OMP_NUM_THREADS}
omp_schedule=${OMP_SCHEDULE}
openblas_num_threads=${OPENBLAS_NUM_THREADS}
ld_preload=${LD_PRELOAD}
EOF

cmd=(srun --kill-on-bad-exit=1 "${EXE}"
    --n_actins="${N_ACTINS}" --n_myosins="${N_MYOSINS}"
    --Lx="${LX}" --Ly=2 --Lz=2 --dimension=3 --pbc_mask=111
    --actin_length=1 --myosin_length=1.5
    --actin_radius=0.015 --myosin_radius=0.025
    --initial_structure=myosin_x_noise
    --actin_direction_noise=-1
    --myosin_direction_noise=1 --myosin_direction_tilt_deg=0
    --pre_phase_steps=10000 --interaction_ramp_steps=40000
    --physical_mobility=true --temperature_kelvin=310
    --thermal_noise_scale=1 --solvent_viscosity=0.001
    --actin_hydrodynamic_diameter=0.010
    --myosin_hydrodynamic_diameter=0.035
    --actin_crowding_parallel=10 --actin_crowding_perpendicular=10
    --actin_crowding_tumble=160 --actin_crowding_roll=160
    --myosin_crowding_parallel=20 --myosin_crowding_perpendicular=20
    --myosin_crowding_tumble=64 --myosin_crowding_roll=64
    --max_actin_displacement=0.005 --max_myosin_displacement=0.005
    --max_actin_rotation=0.005 --max_myosin_rotation=0.005
    --sterics=true
    --actin_steric_k=2000 --myosin_steric_k=1000 --am_steric_k=1000
    --am_surface_gap=0.015
    --k_am=10 --kappa_am=100
    --dynamic_am=true --am_offcenter_torque=true --am_spatial_heads=true
    --am_heads_per_half="${HEADS}" --am_max_heads_per_actin_half=24
    --am_k_on=100 --am_weak_binding_floor=0.1
    --am_k_weak_off=200 --am_k_powerstroke=200
    --am_k_strong_off=160 --am_force_scale=1
    --am_motor_k=1000 --am_power_stroke=0.005
    --am_stall_force=4.5 --am_unloaded_velocity=5
    --am_weak_stiffness_ratio=0.1
    --am_tether_k=3000 --am_tether_slack=0.002
    --am_alignment_k="${ALIGNK}" --am_alignment_slack_deg=2
    --am_tether_saturation_heads=6
    --am_binding_azimuth_selection=true
    --am_binding_azimuth_half_width_deg=30
    --am_rupture_distance=0.1
    --am_strain_dependent_detachment=true
    --am_strain_detachment_theta=0.1
    --am_strain_detachment_max_exponent=20
    --am_strong_off_floor_fraction=1
    --am_cycle_period=0
    --am_activation_floor="${ACTFLOOR}" --am_activation_active_fraction="${ACTFRAC}"
    --am_thin_activation=true --am_thin_feedback=true
    --am_q_calcium_open_rate=20 --am_q_calcium_half=0.5
    --am_q_calcium_hill=3
    --am_q_strong_feedback_rate=100 --am_q_strong_half=0.02
    --am_q_strong_hill=2 --am_q_close_rate=20
    --am_thick_feedback=false --am_status2_memory=0
    --aa_local_multilink=true
    --aa_spacing_ratio=1 --aa_preferred_centerline=0.055
    --aa_capture_margin=0.045
    --aa_antiparallel_overlap=0.1 --aa_barbed_zone=0.2
    --aa_axial_register_ratio=3 --aa_axial_register_force_max=36
    --aa_register_all_statuses=false
    --aa_prebundle_links_per_pair=2
    --aa_site_exclusion=0.0192 --aa_overlap_capture=0.1
    --k_aa=1000 --kappa_aa=100 --aa_perpendicular_k=50
    --aa_status1_radial_k=75
    --aa_stiffness_reference_links=3 --aa_total_link_pool=0
    --aa_status_mode=productive_register
    --aa_status2_force_on=1 --aa_status2_force_off=0.5
    --k_on=5000 --k_off=1 --base_lifetime=1
    --directional_base_lifetime=0.5 --lifetime_coeff=0.4
    --stabilize_status1=false
    --aa_energy_gated_association=false
    --aa_encounter_k=10 --aa_encounter_kappa=0
    --aa_passive_detailed_balance=false
    --tau_rec=0 --max_aa_bonds_per_actin=2
    --aa_neighbor_identity=true
    --aa_neighbor_slot_memory=0.02 --aa_status2_slot_memory=0.20
    --aa_association_max_direction_dot=-0.5
    --aa_link_rupture_length=0.12
    --directional=true --max_myosin_bonds=24
    --titin_mode=explicit_dynamic
    --titin_k=40 --titin_rest_length=0.75
    --titin_links_per_half=3 --titin_capture_distance=0.85
    --titin_k_on=5000 --titin_k_off=0.5
    --titin_strain_gated_kinetics=true --titin_transition_theta=1
    --titin_rebind_memory=0.020 --titin_rebind_rate_multiplier=100
    --titin_require_live_actn2=true
    --titin_k_off_unsupported="${KOFFU}"
    --titin_rebind_requires_actn2=true
    --am_overlap_target_ratio="${OVLTGT}"
    --k_mm="${KMM}" --k_bundle_max=0 --bundle_ramp_steps=0
    --dt=5e-7 --nsteps="${NSTEPS}" --save_every="${SAVE_EVERY}"
    --seed="${SEED}" --filename="${OUTPUT}")
if [ "${#resume_flag[@]}" -gt 0 ]; then
    cmd+=("${resume_flag[@]}")
fi

echo "started_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
printf 'command='
printf '%q ' "${cmd[@]}"
printf '\n'
/usr/bin/time -v "${cmd[@]}"
echo "finished_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
