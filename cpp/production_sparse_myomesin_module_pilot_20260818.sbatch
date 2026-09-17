#!/bin/bash
# A800/M320, 300,000-step (750-ms) production pilot at converged dt=2.5 us.
# Array tasks 0-2 are matched no-myomesin controls; 3-5 use k=240/cap=40.
# Concurrency is capped at two tasks (16 CPUs total).
#SBATCH --partition=owners
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH --time=12:00:00
#SBATCH --job-name=sarc_sparse_mm_pilot
#SBATCH --array=0-5%2
#SBATCH --open-mode=append
#SBATCH --requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/prod_sparse_myomesin_module_A800_M320_750ms_dt2p5us_300k_20260818/logs/%A_%a.out
#SBATCH --error=/scratch/users/xsli/sarcomere/prod_sparse_myomesin_module_A800_M320_750ms_dt2p5us_300k_20260818/logs/%A_%a.err

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

jemalloc=/home/groups/rotskoff/sherry/mamba_root/envs/sarcomere/lib/libjemalloc.so
test -r "${jemalloc}" || { echo "missing jemalloc: ${jemalloc}" >&2; exit 1; }
export LD_PRELOAD="${jemalloc}"

root="${SPARSE_MM_PILOT_ROOT:-/scratch/users/xsli/sarcomere/prod_sparse_myomesin_module_A800_M320_750ms_dt2p5us_300k_20260818}"
exe="${SPARSE_MM_PILOT_EXE:-${root}/bin/sarcomere_sparse_myomesin_module_20260818}"
expected_exe_sha256=2a133f81b13fcbc89f869eb1067f13984ab8be6aa03544244b1d776dbe85e9b6
analyzer="${root}/source_snapshot/check_production_myomesin_pilot.py"
live_checker="${root}/source_snapshot/check_live_production_hdf5.py"
live_monitor="${root}/source_snapshot/monitor_live_production_hdf5.py"

n_actins=800
n_myosins=320
lx=4.825
ly=2
lz=2
dt="${SPARSE_MM_PILOT_DT:-2.5e-6}"
nsteps="${SPARSE_MM_PILOT_NSTEPS:-300000}"
save_every="${SPARSE_MM_PILOT_SAVE_EVERY:-2000}"
pre_phase_steps="${SPARSE_MM_PILOT_PRE_PHASE_STEPS:-2000}"
interaction_ramp_steps="${SPARSE_MM_PILOT_RAMP_STEPS:-8000}"
am_tether_k="${SPARSE_MM_PILOT_AM_TETHER_K:-3000}"
am_head_transverse_k="${SPARSE_MM_PILOT_AM_HEAD_TRANSVERSE_K:-100}"
am_force_scale="${SPARSE_MM_PILOT_AM_FORCE_SCALE:-1}"
am_steric_k="${SPARSE_MM_PILOT_AM_STERIC_K:-1000}"
actin_steric_k="${SPARSE_MM_PILOT_ACTIN_STERIC_K:-2000}"
aa_k_on="${SPARSE_MM_PILOT_AA_K_ON:-2500}"
max_actin_displacement="${SPARSE_MM_PILOT_MAX_ACTIN_DISPLACEMENT:-0.005}"
max_myosin_displacement="${SPARSE_MM_PILOT_MAX_MYOSIN_DISPLACEMENT:-0.005}"
duration_ms=$(awk -v steps="${nsteps}" -v timestep="${dt}" \
    'BEGIN { printf "%.9g", steps * timestep * 1000.0 }')
seeds=(831704 831705 831706)
arms=(control module40)
k_axial=(0 240)

task=${SLURM_ARRAY_TASK_ID:-0}
if (( task < 0 || task >= 6 )); then
    echo "array task ${task} is outside 0-5" >&2
    exit 1
fi
arm_index=$((task / 3))
seed_index=$((task % 3))
arm=${arms[$arm_index]}
seed=${seeds[$seed_index]}
k_mm=${k_axial[$arm_index]}
tag="A800_M320_${arm}_seed${seed}_${duration_ms}ms"
output=${root}/data/${tag}.h5
params=${root}/params/${tag}.params
summary=${root}/analysis/${tag}.json
live_log=${root}/analysis/${tag}.live.jsonl
live_final=${root}/analysis/${tag}.live_final.json

mkdir -p "${root}"/{data,logs,params,analysis}
test -x "${exe}" || { echo "missing executable: ${exe}" >&2; exit 1; }
actual_sha256=$(sha256sum "${exe}" | awk '{print $1}')
if [ "${actual_sha256}" != "${expected_exe_sha256}" ]; then
    echo "frozen executable hash mismatch: ${actual_sha256}" >&2
    exit 1
fi
test -r "${analyzer}" || { echo "missing analyzer: ${analyzer}" >&2; exit 1; }
test -r "${live_checker}" || { echo "missing live checker: ${live_checker}" >&2; exit 1; }
test -r "${live_monitor}" || { echo "missing live monitor: ${live_monitor}" >&2; exit 1; }

python - "${n_actins}" "${n_myosins}" "${lx}" "${ly}" "${lz}" <<'PY'
import sys
n_actin, n_myosin = map(int, sys.argv[1:3])
lx, ly, lz = map(float, sys.argv[3:6])
if (n_actin, n_myosin) != (800, 320):
    raise SystemExit("pilot must retain the A800/M320 production scale")
density = n_actin / (lx * ly * lz)
if abs(density - 41.4507772) > 1e-3:
    raise SystemExit(f"unexpected actin density {density:.7f}/um^3")
if lx <= 3.68 or ly != 2.0 or lz != 2.0:
    raise SystemExit("production geometry or minimum-image guard failed")
PY

resume_requested=false
if [ -e "${output}" ]; then
    if [ "${SLURM_RESTART_COUNT:-0}" -gt 0 ]; then
        echo "resuming preempted trajectory from its latest committed frame"
        resume_requested=true
    elif [ "${SPARSE_MM_PILOT_RESUME_EXISTING:-false}" = true ]; then
        echo "resuming explicitly authorized existing trajectory from its latest committed frame"
        resume_requested=true
    else
        echo "refusing to overwrite existing trajectory: ${output}" >&2
        exit 1
    fi
fi
exec 9>"${output}.lock"
flock -n 9 || { echo "another writer owns ${output}" >&2; exit 1; }

git_head=$(git -C /home/groups/rotskoff/sherry/sarcomere/sarcomere rev-parse HEAD 2>/dev/null || echo unknown)
python - "${params}" <<PY
from pathlib import Path
Path("${params}").write_text("""model=sparse_local_dynamic_finite_valency_bare_zone_modules
stage=A800_M320_${duration_ms}ms_three_seed_matched_pilot
arm=${arm}
seed=${seed}
n_actins=${n_actins}
n_myosins=${n_myosins}
box_um=${lx}x${ly}x${lz}_periodic
dt_s=${dt}
nsteps=${nsteps}
duration_ms=${duration_ms}
save_every=${save_every}
k_mm_module_axial_pN_per_um=${k_mm}
myomesin_module_force_max_pN=40
myomesin_module_radial_k_pN_per_um=300
myomesin_module_radial_rest_um=0.11
myomesin_axial_capture_um=0.25
myomesin_radial_capture_um=0.04
myomesin_sites_per_filament=2
myosin_exclusion_distance_um=0.0953
k_aa_pN_per_um=3000
aa_stiffness_reference_links=3
aa_status_mode=motor_transmitted_register
am_tether_k_pN_per_um=${am_tether_k}
am_head_transverse_k_pN_per_um=${am_head_transverse_k}
am_force_scale=${am_force_scale}
am_steric_k_pN_per_um=${am_steric_k}
actin_steric_k_pN_per_um=${actin_steric_k}
aa_k_on_per_s=${aa_k_on}
max_actin_displacement_um=${max_actin_displacement}
max_myosin_displacement_um=${max_myosin_displacement}
executable=${exe}
executable_sha256=${actual_sha256}
git_head=${git_head}
slurm_job_id=${SLURM_ARRAY_JOB_ID:-local}
slurm_array_task_id=${task}
""")
PY

cmd=(srun --kill-on-bad-exit=1 "${exe}"
    --n_actins=${n_actins} --n_myosins=${n_myosins}
    --Lx=${lx} --Ly=${ly} --Lz=${lz} --dimension=3 --pbc_mask=111
    --actin_length=1 --myosin_length=1.5
    --actin_radius=0.015 --myosin_radius=0.025
    --initial_structure=myosin_x_noise
    --actin_direction_noise=-1
    --myosin_direction_noise=1 --myosin_direction_tilt_deg=0
    --pre_phase_steps=${pre_phase_steps}
    --interaction_ramp_steps=${interaction_ramp_steps}
    --fixed_lattice_preequil_steps=0
    --physical_mobility=true --temperature_kelvin=310
    --thermal_noise_scale=1 --solvent_viscosity=0.001
    --actin_hydrodynamic_diameter=0.010
    --myosin_hydrodynamic_diameter=0.035
    --actin_crowding_parallel=10 --actin_crowding_perpendicular=10
    --actin_crowding_tumble=160 --actin_crowding_roll=160
    --myosin_crowding_parallel=20 --myosin_crowding_perpendicular=20
    --myosin_crowding_tumble=64 --myosin_crowding_roll=64
    --max_actin_displacement=${max_actin_displacement} --max_myosin_displacement=${max_myosin_displacement}
    --max_actin_rotation=0.005 --max_myosin_rotation=0.005
    --sterics=true
    --actin_steric_k=${actin_steric_k} --myosin_steric_k=300 --am_steric_k=${am_steric_k}
    --myosin_exclusion_distance=0.0953
    --am_cutoff=0.05 --am_optimal=0.055 --am_surface_gap=0.015
    --k_am=0 --kappa_am=0
    --dynamic_am=true --am_offcenter_torque=true --am_spatial_heads=true
    --am_heads_per_half=9 --am_crossover_pull_through=true
    --am_max_heads_per_actin_half=24 --am_sites_per_actin=182
    --am_max_heads_per_actin=24 --am_k_on=100 --am_weak_binding_floor=0.1
    --am_k_weak_off=200 --am_k_powerstroke=200 --am_k_strong_off=160
    --am_force_scale=${am_force_scale} --am_motor_k=1000 --am_power_stroke=0.005
    --am_stall_force=4.5 --am_unloaded_velocity=5
    --am_weak_stiffness_ratio=0.1 --am_tether_k=${am_tether_k}
    --am_tether_slack=0.002 --am_head_transverse_k=${am_head_transverse_k}
    --am_transverse_strain_detachment_theta=0
    --am_alignment_k=120 --am_alignment_slack_deg=2
    --am_tether_saturation_heads=6 --am_binding_azimuth_selection=true
    --am_binding_azimuth_half_width_deg=30 --am_rupture_distance=0.1
    --am_strain_dependent_detachment=true --am_strain_detachment_theta=0.1
    --am_strain_detachment_max_exponent=20 --am_strong_off_floor_fraction=1
    --am_cycle_period=0 --am_activation_floor=0
    --am_activation_active_fraction=0.5 --am_thin_activation=true
    --am_thin_feedback=true --am_q_calcium_open_rate=20
    --am_q_calcium_half=0.5 --am_q_calcium_hill=3
    --am_q_strong_feedback_rate=100 --am_q_strong_half=0.02
    --am_q_strong_hill=2 --am_q_close_rate=20
    --am_initialize_q_calcium_equilibrium=true
    --am_thick_feedback=false --am_status2_memory=0
    --am_bidirectional_crossbridge=false --am_overlap_target_ratio=0
    --aa_local_multilink=true --aa_spacing_ratio=1
    --aa_preferred_centerline=0.055 --aa_capture_margin=0.045
    --aa_antiparallel_overlap=0.10 --aa_barbed_zone=0.20
    --aa_axial_register_ratio=3 --aa_axial_register_force_max=36
    --aa_register_all_statuses=false --aa_status2_registry_on=0.040
    --aa_status2_registry_off=0.040 --aa_status2_support_window=0.010
    --aa_register_saturation_error=0.05 --aa_prebundle_links_per_pair=3
    --aa_site_exclusion=0.0192 --aa_overlap_capture=0.10
    --k_aa=3000 --kappa_aa=100 --aa_perpendicular_k=50
    --aa_status1_radial_k=75 --aa_status1_site_diffusion=0.010
    --aa_status2_locked_until_detach=true --aa_stiffness_reference_links=3
    --aa_total_link_pool=0 --aa_status_mode=motor_transmitted_register
    --aa_status2_force_on=1 --aa_status2_force_off=0.5
    --aa_motor_pair_support_window=0.005
    --aa_motor_pair_commit_time=0.010
    --aa_motor_pair_exchange_time=0.001
    --aa_motor_pair_dropout_time=0.005
    --aa_motor_pair_switch_ratio=1.5
    --k_on=${aa_k_on} --k_off=1 --base_lifetime=0.1
    --directional_base_lifetime=0.05 --lifetime_coeff=0.4
    --stabilize_status1=false --aa_energy_gated_association=true
    --aa_transition_state_theta=0.1 --aa_encounter_k=10
    --aa_encounter_kappa=0 --aa_passive_detailed_balance=false
    --tau_rec=0 --max_aa_bonds_per_actin=2 --aa_neighbor_identity=true
    --aa_neighbor_slot_memory=0.02 --aa_status2_slot_memory=0.20
    --aa_association_max_direction_dot=-0.5 --aa_link_rupture_length=0.12
    --aa_double_overlap_contact=true
    --aa_double_overlap_contact_diameter=0.012
    --aa_double_overlap_contact_epsilon=0.00428
    --aa_double_overlap_contact_force_cap=0
    --directional=true --max_myosin_bonds=24
    --titin_mode=explicit_dynamic --titin_k=100 --titin_rest_length=0.665
    --titin_links_per_half=3 --titin_capture_distance=0.72
    --titin_k_on=5000 --titin_k_off=0.01
    --titin_strain_gated_kinetics=true --titin_transition_theta=1
    --titin_rebind_memory=0.020 --titin_rebind_rate_multiplier=100
    --titin_require_live_actn2=true --titin_k_off_unsupported=0
    --titin_rebind_requires_actn2=false
    --k_mm=${k_mm} --myomesin_sites_per_filament=2
    --myomesin_k_on=5000 --myomesin_k_off=1
    --myomesin_axial_capture=0.25 --myomesin_radial_rest=0.11
    --myomesin_radial_capture=0.04 --myomesin_radial_k=300
    --myomesin_alignment_cos=0.94 --myomesin_bell_distance=0.0005
    --myomesin_force_max=40 --k_bundle_max=0 --bundle_ramp_steps=0
    --dt=${dt} --nsteps=${nsteps} --save_every=${save_every}
    --seed=${seed} --filename=${output})
if [ "${resume_requested}" = true ]; then
    cmd+=(--resume=true)
fi

echo "arm=${arm} seed=${seed} k_mm=${k_mm}"
echo "started_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
printf 'command='
printf '%q ' "${cmd[@]}"
printf '\n'

set +e
/usr/bin/time -v "${cmd[@]}" &
simulation_pid=$!
python "${live_monitor}" "${output}" "${live_log}" \
    --pid "${simulation_pid}" --initial-seconds 60 --interval-seconds 300 &
monitor_pid=$!
wait "${simulation_pid}"
simulation_status=$?
kill "${monitor_pid}" 2>/dev/null
wait "${monitor_pid}" 2>/dev/null
set -e
if [ "${simulation_status}" -ne 0 ]; then
    echo "simulation exited with status ${simulation_status}" >&2
    exit "${simulation_status}"
fi

python "${live_checker}" "${output}" --require-healthy > "${live_final}"

python - "${output}" "${arm}" <<'PY'
import h5py
import math
import sys

path, arm = sys.argv[1:]
with h5py.File(path, "r", locking=False) as handle:
    meta = handle["meta"].attrs
    expected_k = 240.0 if arm == "module40" else 0.0
    model = meta.get("myomesin_model", b"")
    if isinstance(model, bytes):
        model = model.decode()
    checks = {
        "actin_count": handle["actin/center"].shape[1] == 800,
        "myosin_count": handle["myosin/center"].shape[1] == 320,
        "k_mm": math.isclose(float(meta.get("k_mm", -1)), expected_k),
        "k_aa": math.isclose(float(meta.get("k_aa", -1)), 3000.0),
        "module_force": math.isclose(
            float(meta.get("myomesin_force_max", -1)), 40.0
        ),
        "module_radial_k": math.isclose(
            float(meta.get("myomesin_radial_k", -1)), 300.0
        ),
        "module_model": (
            model == "sparse_local_dynamic_finite_valency_bare_zone_modules"
            if arm == "module40" else model == "disabled"
        ),
        "double_overlap_contact": int(
            meta.get("aa_double_overlap_contact_enabled", -1)
        ) == 1,
    }
    failed = [name for name, passed in checks.items() if not passed]
    if failed:
        raise SystemExit("metadata validation failed: " + ", ".join(failed))
print("completed_metadata_validation=PASS")
PY

python "${analyzer}" "${output}" --output "${summary}"
echo "finished_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
