#!/bin/bash
#SBATCH --partition=rotskoff
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=24GB
#SBATCH --time=02:00:00
#SBATCH --job-name=dt_ladder
#SBATCH --output=/scratch/users/xsli/sarcomere/dt_ladder/dt_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/dt_ladder/dt_%j.err

set -euo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
unset PYTHONPATH
export OMP_NUM_THREADS=2 OMP_DYNAMIC=FALSE HDF5_USE_FILE_LOCKING=FALSE

ROOT=/scratch/users/xsli/sarcomere/dt_ladder
mkdir -p "${ROOT}"
EXE=/home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp/build_crossover_20260806/sarcomere

# Same 20 ms of PHYSICAL time at every dt, so observables are directly
# comparable and only the integrator resolution changes. The cap counters in
# the h5 are the stability instrument: they should stay ~0.
run () {
  DT=$1; N=$2; SAVE=$3; TAG=$4
  ${EXE} \
    --n_actins=160 --n_myosins=80 --Lx=3.86 --Ly=1 --Lz=1 \
    --dimension=3 --pbc_mask=111 \
    --initial_structure=myosin_x_noise \
    --actin_direction_noise=-1 --myosin_direction_noise=1 \
    --actin_length=1 --myosin_length=1.5 \
    --actin_radius=0.015 --myosin_radius=0.025 \
    --pre_phase_steps=0 --interaction_ramp_steps=0 \
    --physical_mobility=true --temperature_kelvin=310 \
    --thermal_noise_scale=1 --solvent_viscosity=0.001 \
    --actin_hydrodynamic_diameter=0.010 --myosin_hydrodynamic_diameter=0.035 \
    --max_actin_displacement=0.005 --max_myosin_displacement=0.005 \
    --max_actin_rotation=0.005 --max_myosin_rotation=0.005 \
    --sterics=true --actin_steric_k=2000 --myosin_steric_k=1000 --am_steric_k=1000 \
    --am_cutoff=0.10 --am_optimal=0.055 --am_surface_gap=0.015 \
    --k_am=0 --kappa_am=0 \
    --dynamic_am=true --am_offcenter_torque=true --am_spatial_heads=true \
    --am_heads_per_half=9 --am_sites_per_actin=182 --am_max_heads_per_actin=24 \
    --am_k_on=400 --am_k_weak_off=800 --am_k_powerstroke=800 --am_k_strong_off=640 \
    --am_unloaded_velocity=20 \
    --am_motor_k=1000 --am_power_stroke=0.005 --am_stall_force=4.5 \
    --am_tether_k=3000 --am_alignment_k=120 --am_head_transverse_k=100 \
    --am_strain_dependent_detachment=true --am_strain_detachment_theta=0.1 \
    --am_thin_activation=true --am_initialize_q_calcium_equilibrium=true \
    --directional=true --aa_local_multilink=true --aa_barbed_zone=0.20 \
    --k_aa=1000 --kappa_aa=100 --aa_perpendicular_k=50 --aa_status1_radial_k=75 \
    --k_on=2500 --k_off=1 --max_aa_bonds_per_actin=2 \
    --titin_mode=explicit_dynamic --titin_k=100 --titin_rest_length=0.665 \
    --titin_capture_distance=0.72 --titin_k_on=5000 --titin_k_off=0.01 \
    --titin_transition_theta=1.0 \
    --am_crossover_pull_through=true \
    --dt=${DT} --nsteps=${N} --save_every=${SAVE} --seed=20260806 \
    --filename=${ROOT}/${TAG}.h5 > ${ROOT}/${TAG}.log 2>&1 &
}

# 20 ms of physical time in every arm
run 5e-7  40000 2000 dt5e-7
run 1e-6  20000 1000 dt1e-6
run 2e-6  10000  500 dt2e-6
run 4e-6   5000  250 dt4e-6
run 8e-6   2500  125 dt8e-6
wait
echo "dt ladder complete"
ls -la ${ROOT}/ || true
