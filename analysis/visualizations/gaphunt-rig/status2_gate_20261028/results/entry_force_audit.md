# Soft repulsion audit of newly mature entry clips

Input: `r27zt0npg130_seed20260805.h5`, original production continuation.
Machine-readable measurements: `entry_force_audit.json`; analysis: `../audit_entry_forces.py`.

The status 1-to-2 gate and the soft tip force use different geometry. The gate
checks the entire connector against the finite 130 nm cylinder. The tip force
uses only the connector midpoint, with a radial taper from 100 to 130 nm.

| Link / myosin | First saved excluded time | Connector minimum radius in axial slab | Midpoint radius | Reconstructed candidate tip force |
|---|---:|---:|---:|---:|
| 168764 / 32 minus | 2.110 s | 129.801 nm | 131.995 nm | 0 pN |
| 171917 / 255 plus | 2.195 s | 128.872 nm | 154.852 nm | 0 pN |

These forces are reconstructed at the saved pose. Recorded force diagnostics
were evaluated before the final integrator motion within that step. They are
not exactly simultaneous, and the trajectory does not store the selected link
ID. The recorded selected patch radii at those frames are 101.064 and
102.016 nm, respectively, substantially different from the tracked midpoints.
Thus the recorded 13.888 and 19.741 pN tip forces cannot be read as forces on
the tracked connectors.

The kernel applies `force * outward_axis` and no radial force. Its radial gate
changes the magnitude of an axial force; it does not push a connector outward
radially. It therefore has no direct radial restoring component against the
lateral entry in the second clip. The separate connector/body contact force
has an outer range of only 50 + 5 + 15 = 70 nm, leaving a 70–130 nm region
without that radial contact force.

Only the strongest candidate patch per myosin tip receives the tip force.
For example, at frame 423 the tracked first connector requests approximately
8.219 pN at its saved pose, while the recorded selected patch receives
16.915 pN and has a very different radius/gap. A 20 pN tip budget is a maximum
for one selected connector, split into 10 pN per actin at saturation; it is
not a guaranteed force on every excluded connector.

The sum of all forces is capped afterward. These mature actins have a
26.373 pN cap, not the 7.127 pN immature-actin cap. The two actins in the
second clip hit the cap on every step of the sampled intervals near entry.
Opposing forces and this cap can further limit the motion, but do not explain
the earlier midpoint detection gap. The saved data cannot isolate the
post-cap share of an individual contact contribution.

Relevant active source:
- `cpp/src/interaction.cpp:1207`: midpoint construction.
- `cpp/src/interaction.cpp:1227`: midpoint radial cutoff.
- `cpp/src/interaction.cpp:1246`: purely axial patch force.
- `cpp/src/sarcomere_actn2.cpp:2328`: strongest patch selection.
- `cpp/src/sarcomere.cpp:716`: sum-force cap.

Conclusion: the implemented force law is not a continuous full-connector
130 nm exclusion mechanism. Matching that intended region would require
full-connector contact detection, a restoring force opposing both lateral
and axial entry, and a force allocation that does not leave other overlapping
connectors unforced. Increasing the existing force maximum alone does not
close these gaps. No model implementation was changed for this audit.
