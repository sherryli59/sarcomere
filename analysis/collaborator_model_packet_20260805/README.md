# Collaborator model packet — 2026-08-05

This packet intentionally excludes every later production trajectory with collapsed ACTN2 retention.

## Model brief

- [Current coarse-grained model with equations](sarcomere-model-brief.pdf)

## Interactive 3-D

- [Final 200,000-step snapshots (simplified)](final-snapshots-3d.html)
- [Successful 1280-actin production trajectory explorer](successful-production-3d.html)

The simplified final-snapshot view is the recommended collaborator visual. Myosin has one color; actin polarity is light versus dark blue; status 1 is dashed and can be hidden. Both views use the two successful 100 ms model-19 runs (`PRODALIGN120`, seeds 2080201 and 2080202) and show only ACTN2-bound actins plus myosins attached to those actins.

## Minimal animations

- [minimal-free-actin-walkoff.mp4](animations/minimal-free-actin-walkoff.mp4) — Free actin on one myosin: weak/strong heads, material binding sites, connector strain, transverse tether, force, torque, and walk-off.
- [minimal-one-sided-crosslink.mp4](animations/minimal-one-sided-crosslink.mp4) — One motor-bound actin crosslinked to an unsupported actin: the pair co-moves while ACTN2 remains turnover-capable.
- [minimal-force-balanced.mp4](animations/minimal-force-balanced.mp4) — Two myosins and four actins from the registered start: endogenous opposing motor load and a stalled A-A register.
- [minimal-annealing-start.mp4](animations/minimal-annealing-start.mp4) — The same balanced cell started with excess A-A overlap and reduced A-M overlap: sliding toward the preferred A-A register before stalling.

## Production data provenance

- `/scratch/users/xsli/sarcomere/prod_titin_serial_20260802/data/PRODALIGN120_A1280_M640_h35_100ms_seed2080201.h5`
- `/scratch/users/xsli/sarcomere/prod_titin_serial_20260802/data/PRODALIGN120_A1280_M640_h35_100ms_seed2080202.h5`

The 3-D explorer is evidence for A-A maturation and nematic alignment. It is not presented as proof that finite A-M placement is solved.
