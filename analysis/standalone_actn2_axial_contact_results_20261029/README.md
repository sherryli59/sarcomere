# Axial-stop ACTN2 contact (follow-up to standalone_actn2_contact_results_20261028)

The weighted-energy law (alpha=3) was meant to push the myosin back along its axis but also pushed sideways: on a 40 nm outward connector with 20 nm burial it gives 19.9 pN axial and 30.3 pN radial. When one connector was lost, the leftover radial pushes moved the myosin sideways (escape_instance_1: 102 nm sideways, 1067 nm past the plane).

This run replaces it with an axial stop:

    U(x) = g(rho) * Phi(depth(z)),   depth = max(0, H - sqrt(z^2 + eps^2)),   H = 800 nm
    force on connector point = -(dU/dz) * u        (myosin gets the equal and opposite axial reaction)
    g(rho) = 1 for rho <= 130 nm, C1 smoothstep to 0 at 140 nm

The radial position only decides whether contact applies (through g). The radial gradient term -Phi g' rho_hat is dropped, so contact never pushes either body sideways. Phi, k=1000 pN/um and the 30 pN slope are unchanged. Each connector has a 30 pN axial capacity, and alpha no longer appears in the law.

## Tradeoff

The force is the exact axial part of -grad U. At fixed radial positions it is conservative, and U is the work needed to push the connector out along the axis. It is not the full gradient of any energy. A point loop that stays inside the gate does zero work. A loop that crosses the gate edge does Phi*(g_in - g_out): the measured 399.67 pN nm matches the analytic value to 1e-12. In this rig all connectors stay inside the gate, so the work balance still closes: the residual is 7.9 pN nm at dt=1e-7, compared with 8.4 pN nm for the weighted law. Centre motion changes by 0.011 nm under dt/quadrature refinement. For weak birth, U can be used as the penalty and read as the axial work to eject the connector.

## Results (same rig, integrator, hazard law and seed as 20261028; only the contact kernel is swapped)

The weighted rows reproduce the 20261028 numbers exactly: 79.7% / 75% / 47.

| 40 pN drive | Weighted alpha=3 | Axial stop |
|---|---|---|
| 4 connectors, fixed | holds, -35 nm | holds, -35 nm |
| 3 connectors (one removed) | passes 1531 nm, **102 nm sideways** | holds, -30 nm, 0 sideways |
| 2 connectors | passes 1540 nm, 60 nm sideways | holds, -20 nm, 0 sideways |
| 1 connector (capacity 30 < 40 pN) | passes, 50 nm sideways | passes axially (405 nm at 2 ms), 0 sideways |
| 64-realization failure: lose >= 1 | 79.7% | 51.6% |
| passed plane | 75.0% (47 with live links) | 9.4% (4 with live links) |
| max myosin sideways shift | 102 nm | 0 nm (exact) |

Losing connectors now only lowers the axial capacity, which is 30 pN per live connector. Failure is less frequent because the contact hazard depends on endpoint force magnitude, and the radial component no longer adds to that magnitude. In realization 1 (the saved escape instance), no connector is lost and the myosin stays at -35 nm.

The integrator checks pass. Lateral force is below 4e-15 pN, the axial finite-difference error is 1.6e-8 pN, and the torque balance is 2e-15 pN um.

## Limits

The rig omits body rotations. An axial force applied off-axis still exerts a torque on the myosin, so under asymmetric loss a production myosin would be bent or tilted rather than pushed sideways. Other forces (motors, Brownian noise, network springs) can still move it laterally. Parameters are illustrative.

Files: `standalone_actn2_axial_contact.py` (imports the 20261028 script and reuses its integrator unchanged), `summary.json`, `trajectories.npz`, `axial_vs_weighted.png`. To reproduce, run on a compute node from `analysis/`: `python standalone_actn2_axial_contact.py --output <new empty dir>`. The run took 80 s on one CPU (Slurm job 50105633).

## Movies (analysis/animate_sarcomere_3d.py)

`movies/` holds side-by-side movies (weighted energy left, axial stop right) for the saved escape run 1 and for the run with C1 removed. Each has transverse (X–Y), end-on and isometric views, plus `*_side_by_side_final.png`. Each single-law, single-view movie is also kept.

- `export_rig_h5.py` writes the rig output into `h5/` in the renderer's HDF5 format. The rig simulates only connector endpoints and the myosin. Each endpoint is drawn on a 1 µm actin, 100 nm from its barbed end (the +x actin at 80 nm, the -x partner at 120 nm, matching the standalone weak-site geometry). These actins translate rigidly with the endpoint and are display-only. Actin force arrows are the contact force on that endpoint; the myosin arrow is the contact reaction.
- `render_rig_with_pipeline.py` runs the renderer unchanged except for a camera zoom (end-on view 3×). The renderer scales force arrows per frame, so arrow lengths are not comparable between panels or frames.
- Regenerate with `sbatch render_axial_contact_rig.sbatch <this directory>` (GPU partition, about 5 min).
EOF' 2>&1 | grep -v "^Automated"