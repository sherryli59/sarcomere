# Titin-like balanced production run

Date: 2026-07-31  
Seed: 1001  
Release duration: 1.000 s after 50 ms fixed-lattice chemistry pre-equilibration  
Integrator: production `sarcomere`, 3-D full translational and rotational Brownian dynamics, `dt=1 us`  
Output cadence: 5 ms (201 committed frames)

## Promoted sandbox behavior

- Balanced 4-actin/2-myosin periodic topology, with no fixed filaments.
- All 35 weak/strong head slots per myosin half, spatial crown anchors, connector strain, off-center force and torque.
- Three turnover-capable local ACTN2 sites per A-A pair.
- Every attached ACTN2 state carries the preferred 0.10-um A-A register force; status 1/2 still controls lifetime.
- Reciprocal, always-present barbed-end-to-owning-myosin-center titin-like springs (`k=220 pN/um`, `L0=1.125 um`).
- Physical anisotropic rod mobility, translational diffusion, rotational diffusion, and thermal torque at 310 K.

## Production-density result: 141 heads per half

| release time | mean A-A (um) | mean A-M (um) | A0-A1 3-D (um) | A2-A3 3-D (um) | live ACTN2 | weak + strong heads |
|---:|---:|---:|---:|---:|---:|---:|
| 0 ms | 0.1600 | 0.4950 | 0.1678 | 0.1678 | 4 | 32 + 37 |
| 5 ms | 0.0915 | 0.5290 | 0.1000 | 0.0992 | 6 | 46 + 30 |
| 100 ms | 0.0873 | 0.5337 | 0.0957 | 0.0977 | 6 | 41 + 30 |
| 500 ms | 0.0811 | 0.5372 | 0.0908 | 0.0932 | 6 | 36 + 39 |
| 1000 ms | 0.0868 | 0.5331 | 0.0955 | 0.0951 | 6 | 47 + 25 |

Mean A-M overlap during the terminal 500 ms was 0.5357 um, with fitted drift
of -0.00070 um/s, compared with 1.63 um/s during the first 20 ms. Mean head
occupancy over the terminal 500 ms was 73.9. Two empty ACTN2 sites rebound
after release and matured; no mature link detached during this 1-s seed.

## Sandbox-matched port validation: 35 heads per half

| release time | mean A-A (um) | mean A-M (um) | A0-A1 3-D (um) | A2-A3 3-D (um) | live ACTN2 | weak + strong heads |
|---:|---:|---:|---:|---:|---:|---:|
| 0 ms | 0.1600 | 0.4950 | 0.1678 | 0.1678 | 5 | 17 + 16 |
| 5 ms | 0.0979 | 0.5249 | 0.1017 | 0.1047 | 6 | 9 + 9 |
| 100 ms | 0.0917 | 0.5316 | 0.0997 | 0.1019 | 6 | 13 + 20 |
| 500 ms | 0.0887 | 0.5330 | 0.0984 | 0.1045 | 6 | 13 + 15 |
| 1000 ms | 0.0899 | 0.5313 | 0.0981 | 0.1050 | 6 | 11 + 16 |

Mean A-M overlap during the terminal 500 ms was 0.5325 um. Its fitted drift
was 0.00131 um/s, compared with 1.29 um/s during the first 20 ms. Thus the
trajectory shows rapid slide followed by a finite-overlap plateau rather than
continued walking.

Two mature ACTN2 molecules detached after release and both pairs recovered:

- A2-A3: detach 26.682 ms, rebind 26.693 ms, remature 26.702 ms.
- A0-A1: detach 365.534 ms, rebind 453.150 ms, remature 453.151 ms.

All saved geometry and orientations are finite. No actin or myosin force,
torque, or rotation cap fired. The final resume checkpoint is at absolute KMC
step 1,050,000 in both trajectories. A copied checkpoint from the 35-head
trajectory was also resumed from step 70,000 to 75,000 successfully.
