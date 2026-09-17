#!/bin/bash
# Matched density tests derived from the completed A800/M320, Lx=6 um,
# no-myomesin production run. The only model-input changes are particle counts.
# Submit as either:
#   sbatch ... A1600_M320
#   sbatch ... A1600_M640
#
# The archived production script is rendered in memory so these jobs use the
# exact frozen command, executable, monitoring, and post-run validation logic.
#SBATCH --partition=rotskoff
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64GB
#SBATCH --time=7-00:00:00
#SBATCH --open-mode=append
#SBATCH --requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/actin_density_A1600_no_myomesin_Lx6_production_5s_20260824/logs/%x_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/actin_density_A1600_no_myomesin_Lx6_production_5s_20260824/logs/%x_%j.err

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
new_root=/scratch/users/xsli/sarcomere/actin_density_A1600_no_myomesin_Lx6_production_5s_20260824
actin_density=$(awk -v n="${n_actins}" 'BEGIN { printf "%.10g", n / 24.0 }')
myosin_density=$(awk -v n="${n_myosins}" 'BEGIN { printf "%.10g", n / 24.0 }')

test -r "${base_script}" || { echo "missing base script: ${base_script}" >&2; exit 1; }

render_script() {
    sed \
        -e "s#^run_kind=.*#run_kind=cellA_no_myomesin_Lx6_density_${density_case}_production_5s#" \
        -e "s#^root=\${run_root_base}\$#root=${new_root}#" \
        -e "s#^n_actins=800\$#n_actins=${n_actins}#" \
        -e "s#^n_myosins=320\$#n_myosins=${n_myosins}#" \
        -e 's#^tag=A800_M320_cellA_Lx6_${arm}_seed${seed}_${duration_ms}ms$#tag=A${n_actins}_M${n_myosins}_cellA_Lx6_${arm}_seed${seed}_${duration_ms}ms#' \
        -e "s#^experiment=cell_A_no_myomesin_Lx6_axial_zdisc_tip_nonfatal_audit_production_5s#experiment=cell_A_no_myomesin_Lx6_density_${density_case}_production_5s#" \
        -e "s#^single_variable_change=.*#single_variable_change=particle_counts_A800_M320_to_${density_case}#" \
        -e 's#^particle_counts_held_fixed=true#particle_counts_held_fixed=false#' \
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
