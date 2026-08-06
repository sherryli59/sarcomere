# Mechanism intervention runbook

This runbook preserves the requested order: establish the status-2-only
baseline first, then run history switching, then match constant-lifetime
kinetics.

## 1. Build and verify status-2-only stabilization

```bash
cmake -S cpp -B cpp/build_release -DBUILD_TESTING=ON -DBUILD_BENCHMARK=OFF
cmake --build cpp/build_release --target sarcomere catch_bond_kinetics_test -j 8
ctest --test-dir cpp/build_release --output-on-failure -R catch_bond_kinetics
```

With `base_lifetime=0`, status 1 now survives its formation step and is removed
at its next KMC decision. Only status 2 receives `directional_base_lifetime`
and `lifetime_coeff * f_load`.

Run the paired medium/high status-2-only controls before the switching arms:

```bash
sbatch --array=0-5 cpp/history_switch_status2.sbatch
```

The defaults are A2000/M1000 in 12x2x2, three replicas, medium lambda
0.00316227766, high lambda 0.1, and 400000 total steps. Environment variables
prefixed `HISTORY_` override these values.

## 2. History-switch test

After the control arms finish, submit the switching histories:

```bash
sbatch --array=6-11 cpp/history_switch_status2.sbatch
```

Arms 6-8 assemble at medium lambda to step 200000 and then switch high. Arms
9-11 run high to 100000, medium to 300000, then high to 400000. Each transition
uses the existing resume checkpoint, so positions, bonds, bond ages, and RNG
state continue without reinitialization.

Extract a jointly selected M/Z order over all axial modes and summarize each
kinetic phase:

```bash
python analysis/extract_sarcomeric_order.py \
  /scratch/users/xsli/sarcomere/data_history_status2 \
  --output analysis/history_status2_order.csv \
  --stride 10 --box 12 2 2 --yz-bins 2 2 --save-every 100 \
  --phi3d-n-max-x 30 --phi3d-n-max-perp 5 --shared-all-nx

python analysis/analyze_history_switch.py \
  analysis/history_status2_order.csv \
  --output-dir analysis/history_status2_summary
```

The four reported contrasts test failure to assemble at high lambda, retention
after medium-lambda assembly, rescue during the medium pulse, and retention
after restoring high lambda.

## 3. Match constant status-2 lifetime controls

Run a short two-dimensional calibration grid. Constant lifetime is supplied
through `directional_base_lifetime`, not `base_lifetime`, so status 1 remains
short-lived in both the catch and constant arms.

```bash
sbatch cpp/status2_lifetime_match.sbatch
```

The default A2000/M1000 grid has two replicas of one catch reference plus
`tau_const in {0.001,0.002,0.003,0.004,0.006}` crossed with
`k_on in {500,1000,2000}`. It is intentionally a calibration run of 150000
steps. Rank candidates using the tail CB2 occupancy and the mean of binding and
unbinding event rates:

```bash
python analysis/match_lifetime_controls.py \
  /scratch/users/xsli/sarcomere/data_status2_lifetime_match \
  --output-dir analysis/status2_lifetime_match \
  --reference-pattern catch_lt0p003162
```

The best one or two candidates should then be rerun to 400000 steps with three
or more replicas and compared using the same sarcomeric order as the catch
reference. Matching both occupancy and turnover prevents interpreting a simple
change in bound fraction as evidence for load selectivity.
