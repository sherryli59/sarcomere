#!/bin/bash
# Three matched A2100/M1000 interventions testing whether sarcomere assembly
# survives removal of the empirical AM overlap arrest and replacement of the
# coarse motor/crosslink rules by explicit molecular mechanics.
#
# task 0 / run 2: prescribed legacy motor, overlap throttle disabled, pair AA
# task 1 / run 3: explicit heads, overlap throttle disabled, same pair AA
# task 2 / run 4: explicit heads, local sliding ACTN2, motor-transmitted status 2
#
# The first two tasks retain the weak historical bundle field so run 3 isolates
# the motor replacement.  Task 2 removes that field together with pair-level AA
# registry.  All tasks share the same seed, physical mobility, full sterics,
# geometry, pre-relaxation, timestep, and zero directional/catch lifetime bonus.
#SBATCH --partition=owners,hns,normal
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48GB
#SBATCH --time=2-00:00:00
#SBATCH --job-name=sarc_overlap_mech
#SBATCH --array=0-2
#SBATCH --open-mode=append
#SBATCH --requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/overlap_mechanism_controls_A2100_20260819/logs/%A_%a.out
#SBATCH --error=/scratch/users/xsli/sarcomere/overlap_mechanism_controls_A2100_20260819/logs/%A_%a.err

set -euo pipefail

eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export OMP_DYNAMIC=FALSE
export OMP_SCHEDULE=static
export OPENBLAS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export MKL_NUM_THREADS=1
export HDF5_USE_FILE_LOCKING=FALSE

jemalloc=/home/groups/rotskoff/sherry/mamba_root/envs/sarcomere/lib/libjemalloc.so
test -r "${jemalloc}" || { echo "missing jemalloc: ${jemalloc}" >&2; exit 1; }
export LD_PRELOAD="${jemalloc}"

root="${OVERLAP_CONTROL_ROOT:-/scratch/users/xsli/sarcomere/overlap_mechanism_controls_A2100_20260819}"
exe="${OVERLAP_CONTROL_EXE:-${root}/bin/sarcomere_overlap_controls_20260819}"
expected_exe_sha256=e2e6db8d13f88cb935bbc9021911b53e9be0cbabc452eb8f5398b54a49d33431

n_actins=2100
n_myosins=1000
lx=12
ly=2
lz=2
seed=97500
dt="${OVERLAP_CONTROL_DT:-2.5e-6}"
nsteps="${OVERLAP_CONTROL_NSTEPS:-300000}"
save_every="${OVERLAP_CONTROL_SAVE_EVERY:-2000}"
pre_phase_steps="${OVERLAP_CONTROL_PRE_PHASE_STEPS:-2000}"
interaction_ramp_steps="${OVERLAP_CONTROL_RAMP_STEPS:-8000}"
resume_target=$((pre_phase_steps + nsteps))

arms=(legacy_no_throttle explicit_heads_pair_aa explicit_heads_local_actn2)
task=${SLURM_ARRAY_TASK_ID:-0}
if (( task < 0 || task >= 3 )); then
    echo "array task ${task} is outside 0-2" >&2
    exit 1
fi
arm=${arms[$task]}
tag="A${n_actins}_M${n_myosins}_${arm}_lambda0_seed${seed}"
output="${root}/data/${tag}.h5"
params="${root}/params/${tag}.params"
mkdir -p "${root}"/{data,logs,params,analysis}

test -x "${exe}" || { echo "missing executable: ${exe}" >&2; exit 1; }
actual_sha256=$(sha256sum "${exe}" | awk '{print $1}')
if [ "${actual_sha256}" != "${expected_exe_sha256}" ]; then
    echo "frozen executable hash mismatch: ${actual_sha256}" >&2
    exit 1
fi

dynamic_am=false
aa_local_multilink=false
aa_status_mode=legacy_motor_identity
k_am=10
kappa_am=100
k_bundle_max=0.1
aa_perpendicular_all_statuses=true
aa_status2_locked_until_detach=false
if [ "${arm}" = explicit_heads_pair_aa ]; then
    dynamic_am=true
    k_am=0
    kappa_am=0
elif [ "${arm}" = explicit_heads_local_actn2 ]; then
    dynamic_am=true
    aa_local_multilink=true
    aa_status_mode=motor_transmitted_register
    k_am=0
    kappa_am=0
    k_bundle_max=0
    aa_perpendicular_all_statuses=false
    aa_status2_locked_until_detach=true
fi

exec 9>"${output}.lock"
flock -n 9 || { echo "another writer owns ${output}" >&2; exit 1; }

resume_flag=()
run_nsteps=${nsteps}
checkpoint_step=-1
if [ -s "${output}" ]; then
    if [ "${SLURM_RESTART_COUNT:-0}" -gt 0 ] ||
       [ "${OVERLAP_CONTROL_RESUME_EXISTING:-false}" = true ]; then
        checkpoint_step=$(python - "${output}" "${n_actins}" "${n_myosins}" <<'PY'
import h5py
import sys

path, n_actins, n_myosins = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
with h5py.File(path, "r", locking=False) as handle:
    if handle["/actin/center"].shape[1] != n_actins:
        raise SystemExit("actin count mismatch")
    if handle["/myosin/center"].shape[1] != n_myosins:
        raise SystemExit("myosin count mismatch")
    candidates = []
    for group in ("/resume", "/resume_prev"):
        key = f"{group}/current_step"
        if key in handle:
            candidates.append(int(handle[key][...].reshape(-1)[-1]))
    if not candidates:
        raise SystemExit("no promoted resume checkpoint")
    print(max(candidates))
PY
        )
        if [ "${checkpoint_step}" -ge "${resume_target}" ]; then
            echo "${arm} already complete at step ${checkpoint_step}/${resume_target}"
            exit 0
        fi
        resume_flag=(--resume=true)
        run_nsteps=${resume_target}
        echo "resuming ${arm} from ${checkpoint_step} to ${resume_target}"
    else
        echo "refusing to overwrite existing trajectory: ${output}" >&2
        exit 1
    fi
elif [ -e "${output}" ]; then
    echo "refusing to overwrite empty trajectory path: ${output}" >&2
    exit 1
fi

python - "${params}" <<PY
from pathlib import Path

Path("${params}").write_text("""experiment=overlap_mechanism_controls_A2100_20260819
arm=${arm}
scientific_question=does_assembly_survive_without_empirical_overlap_arrest
n_actins=${n_actins}
n_myosins=${n_myosins}
box_um=${lx}x${ly}x${lz}_periodic
seed=${seed}
dt_s=${dt}
nsteps=${nsteps}
duration_s=$(awk -v n="${nsteps}" -v dt="${dt}" 'BEGIN { print n*dt }')
pre_phase_steps=${pre_phase_steps}
interaction_ramp_steps=${interaction_ramp_steps}
dynamic_am=${dynamic_am}
am_overlap_target_ratio=0
aa_local_multilink=${aa_local_multilink}
aa_status_mode=${aa_status_mode}
aa_perpendicular_all_statuses=${aa_perpendicular_all_statuses}
aa_status2_locked_until_detach=${aa_status2_locked_until_detach}
k_am=${k_am}
kappa_am=${kappa_am}
k_bundle_max=${k_bundle_max}
base_lifetime_s=1e-3
directional_base_lifetime_s=0
lifetime_coeff=0
physical_mobility=true
full_sterics=true
executable=${exe}
executable_sha256=${actual_sha256}
source_snapshot=${root}/source_snapshot
git_head=486fe5c76b6bb9867e0f97ca9444f46f5c0dec68
slurm_job_id=${SLURM_ARRAY_JOB_ID:-local}
slurm_array_task_id=${task}
resume_checkpoint_step=${checkpoint_step}
simulator_nsteps_argument=${run_nsteps}
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
    --max_actin_displacement=0.005 --max_myosin_displacement=0.005
    --max_actin_rotation=0.005 --max_myosin_rotation=0.005
    --sterics=true
    --actin_steric_k=2000 --myosin_steric_k=300 --am_steric_k=1000
    --myosin_exclusion_distance=0.0953
    --am_cutoff=0.05 --am_optimal=0.055 --am_surface_gap=0.015
    --k_am=${k_am} --kappa_am=${kappa_am}
    --dynamic_am=${dynamic_am}
    --am_offcenter_torque=true --am_spatial_heads=true
    --am_heads_per_half=9 --am_crossover_pull_through=true
    --am_max_heads_per_actin_half=24 --am_sites_per_actin=182
    --am_max_heads_per_actin=24
    --am_k_on=100 --am_weak_binding_floor=0.1
    --am_k_weak_off=200 --am_k_powerstroke=200 --am_k_strong_off=160
    --am_force_scale=1 --am_motor_k=1000 --am_power_stroke=0.005
    --am_stall_force=4.5 --am_unloaded_velocity=5
    --am_weak_stiffness_ratio=0.1
    --am_tether_k=3000 --am_tether_slack=0.002
    --am_head_transverse_k=100 --am_transverse_strain_detachment_theta=0
    --am_alignment_k=120 --am_alignment_slack_deg=2
    --am_tether_saturation_heads=6
    --am_binding_azimuth_selection=true
    --am_binding_azimuth_half_width_deg=30
    --am_rupture_distance=0.1
    --am_strain_dependent_detachment=true
    --am_strain_detachment_theta=0.1
    --am_strain_detachment_max_exponent=20
    --am_strong_off_floor_fraction=1
    --am_cycle_period=0 --am_activation_floor=0
    --am_thin_activation=false --am_thin_feedback=false
    --am_thick_feedback=false --am_status2_memory=0
    --am_bidirectional_crossbridge=false --am_overlap_target_ratio=0
    --aa_local_multilink=${aa_local_multilink}
    --aa_spacing_ratio=1 --aa_preferred_centerline=0.055
    --aa_capture_margin=0.045
    --aa_antiparallel_overlap=0.10 --aa_barbed_zone=0.20
    --aa_axial_register_ratio=1 --aa_axial_register_force_max=36
    --aa_register_all_statuses=false
    --aa_perpendicular_all_statuses=${aa_perpendicular_all_statuses}
    --aa_status2_registry_on=0.040 --aa_status2_registry_off=0.040
    --aa_status2_support_window=0.010
    --aa_register_saturation_error=0.05
    --aa_prebundle_links_per_pair=3
    --aa_site_exclusion=0.0192 --aa_overlap_capture=0.10
    --k_aa=10 --kappa_aa=100 --aa_perpendicular_k=50
    --aa_status1_radial_k=-1 --aa_status1_site_diffusion=0.010
    --aa_status2_locked_until_detach=${aa_status2_locked_until_detach}
    --aa_stiffness_reference_links=3 --aa_total_link_pool=0
    --aa_status_mode=${aa_status_mode}
    --aa_status2_force_on=1 --aa_status2_force_off=0.5
    --aa_motor_pair_support_window=0.005
    --aa_motor_pair_commit_time=0.010
    --aa_motor_pair_exchange_time=0.001
    --aa_motor_pair_dropout_time=0.005
    --aa_motor_pair_switch_ratio=1.5
    --k_on=1000 --k_off=1 --base_lifetime=1e-3
    --directional_base_lifetime=0 --lifetime_coeff=0
    --stabilize_status1=false
    --aa_energy_gated_association=false
    --aa_transition_state_theta=0.1
    --aa_encounter_k=10 --aa_encounter_kappa=0
    --aa_passive_detailed_balance=false
    --tau_rec=0 --max_aa_bonds_per_actin=2
    --aa_neighbor_identity=true
    --aa_neighbor_slot_memory=0.02 --aa_status2_slot_memory=0.20
    --aa_association_max_direction_dot=-0.5
    --aa_link_rupture_length=0.12
    --directional=true --max_myosin_bonds=6
    --titin_k=0 --k_mm=0
    --k_bundle_max=${k_bundle_max} --bundle_ramp_steps=5000
    --dt=${dt} --nsteps=${run_nsteps} --save_every=${save_every}
    --seed=${seed} --filename=${output})
if [ "${#resume_flag[@]}" -gt 0 ]; then
    cmd+=("${resume_flag[@]}")
fi

echo "arm=${arm} task=${task}"
echo "started_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
printf 'command='
printf '%q ' "${cmd[@]}"
printf '\n'
/usr/bin/time -v "${cmd[@]}"
echo "finished_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
