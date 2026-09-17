#!/bin/bash
# Controlled one-link-per-pair A2400/M320 production run requested on
# 2026-08-27. Relative to the first links=1 run, restore the aggregate
# three-link pair stiffness in the single explicit link and compensate the
# empty-pair association hazard by tripling aa_encounter_k. The frozen
# executable uses the exact conservative spatial-incidence traversal.
#SBATCH --partition=rotskoff
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64GB
#SBATCH --time=7-00:00:00
#SBATCH --open-mode=append
#SBATCH --requeue
#SBATCH --job-name=A2400_M320_l1_kref1_kenc30
#SBATCH --output=/scratch/users/xsli/sarcomere/A2400_M320_links1_kref1_kenc30_spatial_5s_20260827/logs/%x_%j.out
#SBATCH --error=/scratch/users/xsli/sarcomere/A2400_M320_links1_kref1_kenc30_spatial_5s_20260827/logs/%x_%j.err

set -euo pipefail

base_script=/scratch/users/xsli/sarcomere/cellA_no_myomesin_Lx6_nonfatal_tip_audit_production_5s_20260821/source_snapshot/submit_cellA_no_myomesin_nonfatal_tip_audit_production_5s.sbatch
new_root=/scratch/users/xsli/sarcomere/A2400_M320_links1_kref1_kenc30_spatial_5s_20260827
frozen_exe=${new_root}/bin/sarcomere_spatial_default_20260827
frozen_exe_sha256=664fcb5d9ec0b311a55e22304f1a40f4cb3a32d1c00cc6cd1fb5795152f6494a
n_actins=2400
n_myosins=320
actin_density=100
myosin_density=13.33333333

test -r "${base_script}" || {
    echo "missing base script: ${base_script}" >&2
    exit 1
}
test -x "${frozen_exe}" || {
    echo "missing frozen executable: ${frozen_exe}" >&2
    exit 1
}

render_script() {
    sed \
        -e 's#^run_kind=.*#run_kind=cellA_no_myomesin_Lx6_A2400_M320_links1_kref1_kenc30_spatial_production_5s#' \
        -e "s#^root=\${run_root_base}\$#root=${new_root}#" \
        -e "s#^exe=.*#exe=${frozen_exe}#" \
        -e "s#^expected_exe_sha256=.*#expected_exe_sha256=${frozen_exe_sha256}#" \
        -e "s#^n_actins=800\$#n_actins=${n_actins}#" \
        -e "s#^n_myosins=320\$#n_myosins=${n_myosins}#" \
        -e 's#^tag=A800_M320_cellA_Lx6_${arm}_seed${seed}_${duration_ms}ms$#tag=A2400_M320_links1_kref1_kenc30_spatial_cellA_Lx6_${arm}_seed${seed}_${duration_ms}ms#' \
        -e 's#^experiment=cell_A_no_myomesin_Lx6_axial_zdisc_tip_nonfatal_audit_production_5s#experiment=cell_A_no_myomesin_Lx6_A2400_M320_links1_kref1_kenc30_spatial_production_5s#' \
        -e 's#^single_variable_change=.*#single_variable_change=A2400_links1_stiffness_reference_3_to_1_and_encounter_k_10_to_30#' \
        -e 's#^particle_counts_held_fixed=true#particle_counts_held_fixed=false#' \
        -e "s#^actin_density_per_um3=.*#actin_density_per_um3=${actin_density}#" \
        -e "s#^myosin_density_per_um3=.*#myosin_density_per_um3=${myosin_density}#" \
        -e 's#^max_aa_bonds_per_actin=4$#max_aa_bonds_per_actin=4\naa_prebundle_links_per_pair=1\naa_stiffness_reference_links=1\naa_encounter_k=30\naa_zdisc_tip_projection_traversal=spatial_incidence_exact_conservative_broadphase#' \
        -e 's#--aa_prebundle_links_per_pair=3#--aa_prebundle_links_per_pair=1#' \
        -e 's#--aa_stiffness_reference_links=3#--aa_stiffness_reference_links=1#' \
        -e 's#--aa_encounter_k=10#--aa_encounter_k=30#' \
        -e "s#handle\[\"actin/center\"\].shape\[1\] == 800#handle[\"actin/center\"].shape[1] == ${n_actins}#" \
        -e "s#handle\[\"myosin/center\"\].shape\[1\] == 320#handle[\"myosin/center\"].shape[1] == ${n_myosins}#" \
        -e 's#        "neighbor_valency": int(meta.get("max_aa_bonds_per_actin", -1)) == 4,#        "links_per_pair": int(meta.get("aa_prebundle_links_per_pair", -1)) == 1,\n        "stiffness_reference_links": int(meta.get("aa_stiffness_reference_links", -1)) == 1,\n        "encounter_k": math.isclose(float(meta.get("aa_encounter_k", -1)), 30.0),\n        "neighbor_valency": int(meta.get("max_aa_bonds_per_actin", -1)) == 4,#' \
        -e '/expected_solver = "mobility_weighted_projected_gauss_seidel_g_ge_0"/a\    expected_traversal = "spatial_incidence_exact_conservative_broadphase"' \
        -e '/if text_attr("aa_zdisc_tip_nonpenetration") != expected_solver:/i\    if text_attr("aa_zdisc_tip_projection_traversal") != expected_traversal:\n        raise SystemExit("wrong ACTN2 tip projection traversal")' \
        "${base_script}"
}

if [ "${DENSITY_DRY_RUN:-false}" = true ]; then
    rendered=$(mktemp)
    trap 'rm -f "${rendered}"' EXIT
    render_script > "${rendered}"
    bash -n "${rendered}"
    grep -q -- '--aa_prebundle_links_per_pair=1' "${rendered}"
    grep -q -- '--aa_stiffness_reference_links=1' "${rendered}"
    grep -q -- '--aa_encounter_k=30' "${rendered}"
    grep -q -- 'shape\[1\] == 2400' "${rendered}"
    grep -q -- 'aa_zdisc_tip_projection_traversal' "${rendered}"
    echo "rendered_script_validation=PASS actins=${n_actins} myosins=${n_myosins} links_per_pair=1 stiffness_reference_links=1 encounter_k=30 traversal=spatial_incidence"
    exit 0
fi

render_script | bash
