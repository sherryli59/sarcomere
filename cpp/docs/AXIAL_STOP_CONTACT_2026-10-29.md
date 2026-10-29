# Axial-stop connector force (2026-10-29)

Work copy of the active source (`zdisclink_20261020/cpp`, copied 2026-10-29). The active and durable source trees are unchanged.
`axial_stop.patch` holds the full change against that copy.

## What changed

The per-step ACTN2/myosin-tip force now uses `compute_zdisc_axial_stop_contact`. It is the C++ form of the standalone axial-stop law in `analysis/standalone_actn2_axial_contact.py`.

    u(x) = g(rho) Phi(max(0, H - sqrt(z^2 + eps^2))),   H = L/2 + aa_zdisc_tip_onset
    force on a connector point = axial part of -grad u only (no radial component)

- Phi is quadratic up to the saturation depth, then linear (`aa_zdisc_tip_stop_stiffness` 1000 pN/um, `aa_zdisc_tip_stop_force_max` 30 pN).
- g is 1 up to `aa_zdisc_tip_stop_gate_radius` (130 nm), then a C1 smoothstep to 0 at `aa_zdisc_tip_stop_gate_cutoff` (140 nm). Axial smoothing uses `aa_zdisc_tip_stop_epsilon` (5 nm).
- Energy and endpoint forces are normalized integrals over the whole connector, with endpoint shares (1-t) and t.
- The integral uses 16-node Gauss-Legendre on the analytically clipped support, split at the Phi saturation depth and at the gate radius.
- The myosin gets the opposite axial force and the exact reaction torque. Actin torques use the material attachment sites.
- All periodic images are summed.
- Every distinct actin pair contributes. Redundant links on one pair share one budget: the link with the greatest energy summed over both tips, with ties going to the lower link id. This replaces "one strongest connector per tip".

This force is exactly conservative along the axis at fixed radial positions. It is not the full gradient of u: a point loop that crosses the gate edge does net work Phi*(g_in - g_out).

Unchanged: `zdisc_myosin_tip_gap` and `compute_zdisc_myosin_tip_contact` still serve the exclusion audit and diffusion gate, and the hard projection, with the `aa_zdisc_tip_transverse_*` geometry. The barbed-tip end-stop still uses `aa_zdisc_tip_stiffness` and `aa_zdisc_tip_force_max`.

Diagnostics:
- `aa_zdisc_tip_force` (actin) is now the sum of |F_i| and |F_j|.
- `aa_zdisc_tip_force_{plus,minus}` (myosin) is the total axial force over the selected pairs.
- `gap/rho/gate_{plus,minus}` describe the deepest contacting quadrature node.

New-file and resume metadata record `aa_zdisc_tip_contact_force_law = full_connector_integral_axial_stop_pair_max_energy_v1` and the five `aa_zdisc_tip_stop_*` values.

## Validation

Release build with assertions in all tests: 46/46 CTest tests passed (Slurm 50175945). Binary SHA256: see `build/binary.sha256`.

- `zdisc_axial_stop_contact`: 44 fixtures generated from the standalone Python law with 4096 nodes.
  - 16-node C++ max error is 6.5e-10 pN um in energy and 2.3e-7 pN in force/torque. That is the reference's own accuracy.
  - 1274 random contacts are axial-only to 1e-12, with equal and opposite reaction and balanced angular momentum.
  - Axial finite differences (rigid, one endpoint, myosin) agree to 3.8e-8 pN.
  - Also covered: gate, onset and saturation edges; per-half selection; loop work (0 inside the gate, Phi*(g_in-g_out) across it); periodic translation and a non-nearest image.
- `zdisc_axial_stop_pair_budget`: the full Sarcomere force pass for a single link, redundant links in both orders (the max-energy link acts alone), and two distinct pairs (both add).

## Not done

- No simulator run has used this force yet.
- Resuming a checkpoint written before this change is allowed with the default `aa_zdisc_tip_stop_*` values. Such a resume switches the law abruptly; the plan's contact ramp is not implemented.
- Hard projection defaults to true and is not changed here.
- The durable repository (`$GROUP_HOME/.../sarcomere/cpp`) is not patched.
