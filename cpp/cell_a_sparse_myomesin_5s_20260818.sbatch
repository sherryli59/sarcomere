#!/bin/bash
# Cell A reproduced with the current binary, paired against Cell A plus only
# the sparse finite-valency myomesin module.
#SBATCH --partition=rotskoff
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64GB
#SBATCH --time=48:00:00
#SBATCH --array=0-1%2
#SBATCH --job-name=sarc_cellA_sparse_mm
#SBATCH --open-mode=append
#SBATCH --requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/cellA_sparse_myomesin_A800_M320_5s_dt10us_20260818/slurm/%x_%A_%a.out
#SBATCH --error=/scratch/users/xsli/sarcomere/cellA_sparse_myomesin_A800_M320_5s_dt10us_20260818/slurm/%x_%A_%a.err

set -euo pipefail

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-32}"
export OMP_DYNAMIC=FALSE
export OMP_SCHEDULE=static
export OPENBLAS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export MKL_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

jemalloc=/home/groups/rotskoff/sherry/mamba_root/envs/sarcomere/lib/libjemalloc.so
test -r "${jemalloc}" || { echo "missing jemalloc: ${jemalloc}" >&2; exit 1; }
export LD_PRELOAD="${jemalloc}"

run_root_base=/scratch/users/xsli/sarcomere/cellA_sparse_myomesin_A800_M320_5s_dt10us_20260818
run_kind=${CELL_A_RUN_KIND:-full5s}
root=${run_root_base}/${run_kind}
exe=${CELL_A_EXE:-/scratch/users/xsli/sarcomere/prod_sparse_myomesin_module_A800_M320_750ms_dt2p5us_300k_20260818/bin/sarcomere_sparse_myomesin_module_20260818}
expected_exe_sha256=${CELL_A_EXE_SHA256:-2a133f81b13fcbc89f869eb1067f13984ab8be6aa03544244b1d776dbe85e9b6}
repo=/home/groups/rotskoff/sherry/sarcomere/sarcomere
snapshot=${CELL_A_SNAPSHOT:-${run_root_base}/source_snapshot}
analyzer=${snapshot}/check_production_myomesin_pilot.py
live_checker=${snapshot}/check_live_production_hdf5.py
live_monitor=${snapshot}/monitor_live_production_hdf5.py

n_actins=800
n_myosins=320
lx=4.825
ly=2
lz=2
seed=20260805
dt=1e-5
nsteps=${CELL_A_NSTEPS:-500000}
save_every=${CELL_A_SAVE_EVERY:-500}
myosin_radius=${CELL_A_MYOSIN_RADIUS:-0.025}
am_cutoff=${CELL_A_AM_CUTOFF:-0.05}
am_optimal=${CELL_A_AM_OPTIMAL:-0.055}
myomesin_radial_rest=${CELL_A_MYOMESIN_RADIAL_REST:-0.11}
myomesin_sites_per_filament=${CELL_A_MYOMESIN_SITES_PER_FILAMENT:-2}
myosin_hard_core_diameter=$(awk -v r="${myosin_radius}" 'BEGIN {print 2*r}')
pre_phase_steps=500
interaction_ramp_steps=2000
actin_force_cap_fraction_max=0.75
myosin_force_cap_fraction_max=0.50

task=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}
case "${task}" in
    0) arm=control; k_mm=0 ;;
    1) arm=module40; k_mm=240 ;;
    *) echo "unexpected array task ${task}" >&2; exit 2 ;;
esac

duration_ms=$(awk -v n="${nsteps}" -v d="${dt}" 'BEGIN { printf "%.6g", n*d*1000 }')
tag=A800_M320_cellA_${arm}_seed${seed}_${duration_ms}ms
output=${root}/data/${tag}.h5
params=${root}/params/${tag}.params
summary=${root}/analysis/${tag}.json
live_log=${root}/analysis/${tag}.live.jsonl
live_final=${root}/analysis/${tag}.live_final.json

mkdir -p "${root}"/{data,params,analysis,logs}
test -x "${exe}" || { echo "missing executable: ${exe}" >&2; exit 1; }
actual_sha256=$(sha256sum "${exe}" | awk '{print $1}')
test "${actual_sha256}" = "${expected_exe_sha256}" || {
    echo "frozen executable hash mismatch: ${actual_sha256}" >&2; exit 1;
}
for file in "${analyzer}" "${live_checker}" "${live_monitor}"; do
    test -r "${file}" || { echo "missing analysis file: ${file}" >&2; exit 1; }
done

python - "${nsteps}" "${save_every}" "${dt}" \
    "${myosin_radius}" "${am_cutoff}" "${am_optimal}" \
    "${myomesin_radial_rest}" "${myomesin_sites_per_filament}" <<'PY'
import sys
nsteps, stride = map(int, sys.argv[1:3])
dt = float(sys.argv[3])
myosin_radius, am_cutoff, am_optimal, myomesin_rest = map(float, sys.argv[4:8])
myomesin_sites = int(sys.argv[8])
if nsteps <= 0 or stride <= 0 or nsteps % stride:
    raise SystemExit("nsteps must be positive and divisible by save_every")
if 500 % stride:
    raise SystemExit(
        "save_every must divide the 500-step prephase so saved physical steps "
        "remain aligned in this frozen executable"
    )
if dt != 1e-5:
    raise SystemExit("this comparison requires the original Cell A 10-us timestep")
if not (myosin_radius > 0 and am_cutoff > 0 and am_optimal > 0):
    raise SystemExit("geometric radii and AM distances must be positive")
if 2 * myosin_radius > myomesin_rest + 1e-12:
    raise SystemExit(
        "myomesin radial rest must not lie inside the myosin hard core"
    )
if myomesin_sites <= 0:
    raise SystemExit("myomesin sites per filament must be positive")
PY

resume_requested=false
if [ -e "${output}" ]; then
    if [ "${SLURM_RESTART_COUNT:-0}" -gt 0 ] || [ "${CELL_A_RESUME_EXISTING:-false}" = true ]; then
        echo "resuming from the latest committed frame"
        resume_requested=true
    else
        echo "refusing to overwrite existing trajectory: ${output}" >&2
        exit 1
    fi
fi
exec 9>"${output}.lock"
flock -n 9 || { echo "another writer owns ${output}" >&2; exit 1; }

git_head=$(git -C "${repo}" rev-parse HEAD 2>/dev/null || echo unknown)
python - "${params}" <<PY
from pathlib import Path
Path("${params}").write_text("""experiment=cell_A_plus_sparse_finite_valency_myomesin
run_kind=${run_kind}
arm=${arm}
seed=${seed}
n_actins=${n_actins}
n_myosins=${n_myosins}
box_um=${lx}x${ly}x${lz}_periodic
dt_s=${dt}
nsteps=${nsteps}
duration_ms=${duration_ms}
save_every=${save_every}
cell_A_reference_job=38011216
cell_A_reference_seed=20260805
cell_A_reference_duration_ms=5000
myosin_exclusion_distance_um=0
myosin_steric_k_pN_per_um=1000
myosin_radius_um=${myosin_radius}
myosin_hard_core_diameter_um=${myosin_hard_core_diameter}
am_cutoff_um=${am_cutoff}
am_optimal_um=${am_optimal}
k_aa_pN_per_um=1000
aa_status_mode=productive_register
aa_status2_force_on_pN=0.25
aa_status2_force_off_pN=0.10
aa_energy_gated_association=false
k_mm_module_axial_pN_per_um=${k_mm}
myomesin_module_force_max_pN=40
myomesin_module_radial_k_pN_per_um=300
myomesin_module_radial_rest_um=${myomesin_radial_rest}
myomesin_axial_capture_um=0.25
myomesin_radial_capture_um=0.04
myomesin_sites_per_filament=${myomesin_sites_per_filament}
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
    --actin_radius=0.015 --myosin_radius=${myosin_radius}
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
    --max_actin_displacement=0.005 --max_myosin_displacement=0.005
    --max_actin_rotation=0.005 --max_myosin_rotation=0.005
    --sterics=true
    --actin_steric_k=2000 --myosin_steric_k=1000 --am_steric_k=1000
    --myosin_exclusion_distance=0
    --am_cutoff=${am_cutoff} --am_optimal=${am_optimal} --am_surface_gap=0.015
    --k_am=0 --kappa_am=0
    --dynamic_am=true --am_offcenter_torque=true --am_spatial_heads=true
    --am_heads_per_half=9 --am_crossover_pull_through=true
    --am_max_heads_per_actin_half=24 --am_sites_per_actin=182
    --am_max_heads_per_actin=24 --am_k_on=100 --am_weak_binding_floor=0.1
    --am_k_weak_off=200 --am_k_powerstroke=200 --am_k_strong_off=160
    --am_force_scale=1 --am_motor_k=1000 --am_power_stroke=0.005
    --am_stall_force=4.5 --am_unloaded_velocity=5
    --am_weak_stiffness_ratio=0.1 --am_tether_k=3000
    --am_tether_slack=0.002 --am_head_transverse_k=100
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
    --k_aa=1000 --kappa_aa=100 --aa_perpendicular_k=50
    --aa_status1_radial_k=75 --aa_status1_site_diffusion=0.010
    --aa_status2_locked_until_detach=true --aa_stiffness_reference_links=3
    --aa_total_link_pool=0 --aa_status_mode=productive_register
    --aa_status2_force_on=0.25 --aa_status2_force_off=0.10
    --k_on=2500 --k_off=1 --base_lifetime=0.1
    --directional_base_lifetime=0.05 --lifetime_coeff=0.4
    --stabilize_status1=false --aa_energy_gated_association=false
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
    --k_mm=${k_mm} --myomesin_sites_per_filament=${myomesin_sites_per_filament}
    --myomesin_k_on=5000 --myomesin_k_off=1
    --myomesin_axial_capture=0.25 --myomesin_radial_rest=${myomesin_radial_rest}
    --myomesin_radial_capture=0.04 --myomesin_radial_k=300
    --myomesin_alignment_cos=0.94 --myomesin_bell_distance=0.0005
    --myomesin_force_max=40 --k_bundle_max=0 --bundle_ramp_steps=0
    --dt=${dt} --nsteps=${nsteps} --save_every=${save_every}
    --seed=${seed} --filename=${output})
if [ "${resume_requested}" = true ]; then
    cmd+=(--resume=true)
fi

echo "arm=${arm} seed=${seed} k_mm=${k_mm} run_kind=${run_kind}"
echo "started_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
printf 'command='; printf '%q ' "${cmd[@]}"; printf '\n'

set +e
/usr/bin/time -v "${cmd[@]}" &
simulation_pid=$!
python "${live_monitor}" "${output}" "${live_log}" \
    --pid "${simulation_pid}" --initial-seconds 60 --interval-seconds 300 \
    --actin-force-cap-fraction-max "${actin_force_cap_fraction_max}" \
    --myosin-force-cap-fraction-max "${myosin_force_cap_fraction_max}" &
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

python "${live_checker}" "${output}" --require-healthy \
    --actin-force-cap-fraction-max "${actin_force_cap_fraction_max}" \
    --myosin-force-cap-fraction-max "${myosin_force_cap_fraction_max}" \
    > "${live_final}"

python - "${output}" "${arm}" "${nsteps}" "${save_every}" <<'PY'
import h5py
import math
import sys
path, arm, expected_step, stride = sys.argv[1:]
expected_step, stride = int(expected_step), int(stride)
with h5py.File(path, "r", locking=False) as handle:
    meta = handle["meta"].attrs
    steps = handle["state/current_step"][:, 0]
    model = meta.get("myomesin_model", b"")
    if isinstance(model, bytes):
        model = model.decode()
    expected_k = 240.0 if arm == "module40" else 0.0
    expected_model = (
        "sparse_local_dynamic_finite_valency_bare_zone_modules"
        if arm == "module40" else "disabled"
    )
    checks = {
        "final_step": int(steps[-1]) == expected_step,
        "frame_count": len(steps) == expected_step // stride + 1,
        "actin_count": handle["actin/center"].shape[1] == 800,
        "myosin_count": handle["myosin/center"].shape[1] == 320,
        "dt": math.isclose(float(meta.get("dt", -1)), 1e-5),
        "k_mm": math.isclose(float(meta.get("k_mm", -1)), expected_k),
        "k_aa": math.isclose(float(meta.get("k_aa", -1)), 1000.0),
        "myosin_steric_k": math.isclose(float(meta.get("myosin_steric_k", -1)), 1000.0),
        "myosin_exclusion_distance": math.isclose(float(meta.get("myosin_exclusion_distance", -1)), 0.0),
        "status_mode": int(meta.get("aa_status_mode", -1)) == 7,
        "force_on": math.isclose(float(meta.get("aa_status2_force_on", -1)), 0.25),
        "force_off": math.isclose(float(meta.get("aa_status2_force_off", -1)), 0.10),
        "energy_gate_off": int(meta.get("aa_energy_gated_association", -1)) == 0,
        "module_model": model == expected_model,
        "module_force": math.isclose(float(meta.get("myomesin_force_max", -1)), 40.0),
        "module_radial_k": math.isclose(float(meta.get("myomesin_radial_k", -1)), 300.0),
    }
    failed = [name for name, passed in checks.items() if not passed]
    if failed:
        raise SystemExit("metadata validation failed: " + ", ".join(failed))
print("completed_metadata_validation=PASS")
PY

python "${analyzer}" "${output}" --output "${summary}"
echo "finished_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
