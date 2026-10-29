# Interactive myosin / connector energy

Open `connector_energy_3d.html` in a modern browser. The file is standalone:
Plotly 3.1.0 and compressed numerical meshes are embedded; no server, external
requests, or installation is needed to view it. A browser with WebGL and
DecompressionStream is required.

This visualizes the **revised proposed energy**, not the old 100–130 nm taper,
a calibrated parameter set, a trajectory, or an implemented production law.
The simulator source and production trajectories were not changed.

**Revision note:** this file displays the isotropic alpha=1 predecessor.
The implementation plan now proposes alpha=3 weighting of radial escape to
extend axial capture across the actin shell. That weighted candidate has
not been regenerated here; the 30 pN whole-force bound stated below applies
to this alpha=1 visualization. At alpha=3 with the same potential slope,
the bounds would be 30 pN axial and 90 pN radial/net.

Parameters: myosin physical length 1500 nm and radius 50 nm; mechanical
cylinder radius 130 nm and half-length 800 nm (50 nm approach); smoothing
5 nm; k=1000 pN/um; force bound 30 pN per independent connector. Physical
steric-core energy is not included. Exact formation geometry uses the
physical 750 nm half-length, separately from this energy.

Modes: point probe; 40 nm connector along transverse Y, axial X, or diagonal
X–Y. Finite-connector energy is the normalized integral along the connector.
Each shell is an equal-energy surface in **connector-center coordinates** at
fixed orientation, not a shell of contact points on a particular connector.

Controls: rotate/zoom/pan; whole/tip/cross-section views; four energy shells;
surface opacity; physical myosin and boundary layers; Z=0 energy slice;
example connectors and proportional net-force arrows; selected endpoint
forces and myosin reaction; editable selected connector center; arrow scale.
Every example is an independent hypothetical connector, not a simultaneous
network. Select examples through the table or marker click. Hover for values.

Energy is in pN nm; force is pN; coordinates nm. Arrow length defaults to
4 nm/pN and is never normalized. The myosin reaction is the opposite of the
selected connector's net force; associated torque is not drawn.

Numerics: marching tetrahedra, nonuniform axial grid and 7.5 nm transverse
spacing. Fields use 32-point Gauss integration and readouts 64-point. The
surfaces are interpolated visual approximations. Analytic point gradients,
finite-connector translation and endpoint gradients, force bounds, and
64/128-point convergence were checked. JavaScript values were independently
compared with Python references for all eight examples and four modes.
See validation.json and js_validation.json. Actual headless browser rendering
could not be checked: Chromium requires newer glibc than the compute node
provides. See browser_validation.json. No system changes were attempted.
The standalone HTML should be downloaded and opened in a local modern browser;
the app cannot open remote filesystem file URLs directly in its local browser.

Sources are two directories above this README: build_interactive_energy.py,
interactive_energy_template.html, check_interactive_energy.cjs,
check_energy_browser.cjs, and build_energy_viewer.sh. Generation and checks
run inside Slurm. Plotly source is from
https://cdn.plot.ly/plotly-3.1.0.min.js; its embedded license header is retained.

The side force is intentional for this candidate: contact seeks the nearest
exit from a finite cylinder. Whether this holds a Z-disc under motor loading
remains a separate simulation/calibration question.
