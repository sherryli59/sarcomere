#!/bin/bash
# Matched 200-ms A1600/M320 timing run for the principled sparse ACTN2-tip
# projection fix. Model inputs match job 40664866; the executable and duration
# differ. The archived launcher still supplies monitoring and validation.
#SBATCH --partition=rotskoff
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64GB
#SBATCH --time=1-00:00:00
#SBATCH --open-mode=append
#SBATCH --requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/A1600_M320_principled_projection_profile_200ms_20260825/logs/%x_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/A1600_M320_principled_projection_profile_200ms_20260825/logs/%x_%j.err

set -euo pipefail

density_case=${1:-}
case "${density_case}" in
    A1600_M320)
        n_actins=1600
        n_myosins=320
        ;;
    A1600_M640)
        n_actins=1600
        n_myosins=640
        ;;
    *)
        echo "usage: sbatch $0 {A1600_M320|A1600_M640}" >&2
        exit 2
        ;;
esac

base_script=/scratch/users/xsli/sarcomere/cellA_no_myomesin_Lx6_nonfatal_tip_audit_production_5s_20260821/source_snapshot/submit_cellA_no_myomesin_nonfatal_tip_audit_production_5s.sbatch
new_root=/scratch/users/xsli/sarcomere/A1600_M320_principled_projection_profile_200ms_20260825
profile_exe=${new_root}/bin/sarcomere_sparse_projection_profile_20260825
profile_exe_sha256=254f59812d1840314b1ea84036dda5d9591e0faf2a9ce66004ce9f75d22d7591
profile_nsteps=${PROFILE_NSTEPS:-20000}
actin_density=$(awk -v n="${n_actins}" 'BEGIN { printf "%.10g", n / 24.0 }')
myosin_density=$(awk -v n="${n_myosins}" 'BEGIN { printf "%.10g", n / 24.0 }')

test -r "${base_script}" || { echo "missing base script: ${base_script}" >&2; exit 1; }

render_script() {
    sed \
        -e "s#^run_kind=.*#run_kind=cellA_no_myomesin_Lx6_density_${density_case}_principled_projection_profile_200ms#" \
        -e "s#^root=\${run_root_base}\$#root=${new_root}#" \
        -e "s#^exe=.*#exe=${profile_exe}#" \
        -e "s#^expected_exe_sha256=.*#expected_exe_sha256=${profile_exe_sha256}#" \
        -e "s#^nsteps=.*#nsteps=${profile_nsteps}#" \
        -e "s#^n_actins=800\$#n_actins=${n_actins}#" \
        -e "s#^n_myosins=320\$#n_myosins=${n_myosins}#" \
        -e 's#^tag=A800_M320_cellA_Lx6_${arm}_seed${seed}_${duration_ms}ms$#tag=A${n_actins}_M${n_myosins}_cellA_Lx6_${arm}_seed${seed}_${duration_ms}ms#' \
        -e "s#^experiment=cell_A_no_myomesin_Lx6_axial_zdisc_tip_nonfatal_audit_production_5s#experiment=cell_A_no_myomesin_Lx6_density_${density_case}_principled_projection_profile_200ms#" \
        -e "s#^single_variable_change=.*#single_variable_change=projection_solver_and_runtime_profiling_vs_job_40664866#" \
        -e 's#^particle_counts_held_fixed=true#particle_counts_held_fixed=true#' \
        -e "s#^actin_density_per_um3=.*#actin_density_per_um3=${actin_density}#" \
        -e "s#^myosin_density_per_um3=.*#myosin_density_per_um3=${myosin_density}#" \
        -e "s#handle\[\"actin/center\"\].shape\[1\] == 800#handle[\"actin/center\"].shape[1] == ${n_actins}#" \
        -e "s#handle\[\"myosin/center\"\].shape\[1\] == 320#handle[\"myosin/center\"].shape[1] == ${n_myosins}#" \
        "${base_script}"
}

if [ "${DENSITY_DRY_RUN:-false}" = true ]; then
    render_script | bash -n
    echo "rendered_script_validation=PASS case=${density_case} actins=${n_actins} myosins=${n_myosins}"
    exit 0
fi

render_script | bash
