# Local opposing-actin pointed-end contact sandbox

Date: 2026-08-05

## Interaction

The sandbox adds a default-off unilateral contact between the material
pointed-end caps of antiparallel actins that occupy the same 3-D lattice lane.
It uses a shifted WCA force, applies equal and opposite forces and torques, and
retains the signed axial ordering of the two caps after an attempted crossing.
It does not inspect ACTN2 status, A-M occupancy, or a prescribed overlap.

Parameters used below:

- decorated cap diameter: 0.012 um
- WCA energy scale: 0.00428 pN um (approximately kBT at 310 K)
- force ceiling sweep: 0.5, 1, 2, 5, 10, 20, and 40 pN
- thermal noise: zero for the mechanism-isolation runs

## Results

1. Free actin, 100 ms: no pointed-end contacts; free translocation passes.
2. Normal balanced cell, 100 ms: no pointed-end contacts. The contact law is
   inert at the desired 0.34--0.42 um A-M overlap.
3. True short periodic cell (box x = 3.12 um), 50 ms:
   - contact disabled: mean A-M overlap 0.890 -> 0.922 um
   - contact enabled, 10 pN ceiling: 0.890 -> 0.784 um
   - the contact carries the opposing motor load, but remains much deeper than
     the target finite overlap.
4. Existing `balanced_buried` initializer, 100 ms: zero pointed-end contacts;
   A-M overlap remains approximately 0.99 um on the displaced central halves.

Changing the force ceiling changes the transient and peak load but not the
geometric set point. The contact begins only when opposing thin-filament tips
meet. A segmented bending model would change compliance after contact, not
move contact onset to 0.35--0.45 um A-M overlap.

## Interpretation

The existing `balanced_buried` initializer translates the two thick filaments
by 0.55 um while leaving both Z-pair positions and all opposing-actin pointed
ends unchanged. It is therefore a thick-filament off-centering challenge, not
a compressed/double-overlap challenge. A local actin-actin compression law
must correctly return zero force in that state. Its physical restoring force
would instead come from tension in the opposite half of a bilateral Z-to-M
link such as explicit titin.

The pointed-end contact remains default-off and is not suitable for production
or for selecting the target finite A-M overlap.

## Verification

The complete serial CTest suite passes: 28/28.
