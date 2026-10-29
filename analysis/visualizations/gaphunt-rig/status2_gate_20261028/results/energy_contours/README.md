# Proposed energy contours

These plots illustrate the proposed harmonic A-band contact energy, not the
currently implemented force law or measured trajectory data. Stiffness
1000 pN/um is illustrative and uncalibrated. The separate physical-core
steric energy is omitted; the 50 nm myosin body is shown for reference.

- A: point-energy density around a 1.5 um myosin.
- B: the same density near the plus tip; gap is positive outside the tip.
- C: normalized energy integral over a 40 nm straight radial connector at
  fixed axial gap. The vertical coordinate is its midpoint radius.
- D: radial energy sections at the physical tip plane.

Color bars are logarithmic. Energy below 1 pN nm is blank in contour panels.
White arrows indicate normalized negative energy gradients, not magnitude.
For the radial connector, the W integral and its radial derivative are
evaluated analytically. Finite-difference checks verify plotted forces.

At gap zero and radius 115 nm, the point force is 25 pN axial and 62.5 pN
radial. At midpoint radius 140 nm, a 40 nm radial connector has energy
28.935 pN nm because its near endpoint is at 120 nm; a point at 140 nm has
zero energy. These values use the illustrative stiffness above.

Source: `../../plot_proposed_energy.py`.
Final rendering job: 50023349. PNG, PDF, SVG and parameter JSON are saved
alongside this README. No simulation implementation was changed.
