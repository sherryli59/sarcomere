# Production configuration

- `prodargs.sh`: production argv (hph27 campaign, burial block rule-only, titin and myomesin off, Z-disc stop gated on ACTN2). Copied from `burialzone_20261020/prodargs.sh`.
- The compiled defaults in `src/main.cpp` equal these values, apart from the per-run seed, output file and step counts. Shell variables resolve to HPH=27, KILL=0.050, RATE=1, SLAB=false and CUTOFF=0.
- `sync_defaults.py` checks that the compiled defaults match `prodargs.sh` and exits 1 when they differ. `--apply` rewrites the defaults.
- `run_production.sh EXE OUT.h5 SEED [NSTEPS] [SAVE_EVERY] [flags...]` runs with those defaults. For example, the 54-head variant: `run_production.sh ./build/sarcomere out.h5 20260805 20000 200 --am_heads_per_half=54`.

After changing a production flag, update `prodargs.sh` and run `python sync_defaults.py --apply`.

## Slim checkpoints (resume from the last frame only)

A resume reads only the last committed frame and the `/resume` snapshot, so a
trajectory can be cut to about 87 MB (from 7.6 GB at 2.0 s) before it is copied
to another machine. Each run that resumes from a slim file also writes about 1.5 GB
instead of about 9 GB, because it no longer carries the full history.

```bash
python cpp/production/slim_checkpoint.py SOURCE.h5 SLIM.h5
```

The script keeps:
- the last frame of every per-frame dataset;
- that frame's ACTN2 link rows, with offset starts rebased to 0;
- `/meta` and `/resume` unchanged;
- the event logs, with zero rows.

Dataset dtypes, chunking and maxshape are preserved. It then checks every kept
value against the source and exits non-zero on a mismatch. It refuses any
dataset it does not recognize (for example `/myomesin`, written when `k_mm > 0`)
rather than guessing.

On Sherlock, run it on a compute node, then copy the result from the Mac.
The 5080 has no route to Sherlock.

```bash
srun -p dev -c 1 --mem=4G -t 00:15:00 $GROUP_HOME/sherry/mamba_root/envs/sarcomere/bin/python cpp/production/slim_checkpoint.py $SCRATCH/sarcomere/RUN/data/X.h5 $SCRATCH/sarcomere/RUN/data/X_slim.h5
```

```bash
scp -3 sherlock:/scratch/users/xsli/sarcomere/RUN/data/X_slim.h5 5080:runs/
```

Resume from a slim file exactly as from the full one: copy it, then run with
`--resume=true --filename=COPY.h5`. In the output, the cut frame is frame 0,
not its original index. Find it by step, for example
`np.flatnonzero(f["state/current_step"][:, 0] == 200000)[0]`.

These scripts hard-code frame 400 and need that change before they can read
slim-based runs:
- `analysis/visualizations/gaphunt-rig/axialstop_20261029/run_resume.py`
  (`n == 401` and the `[400]` start-pose checks);
- `status2_gate_20261028/analyze_resume.py`;
- `zcheck.py`;
- `resab_an.py`.

`axialstop_20261029/analyze_ab.py` and `status2_gate_20261028/connector_metrics.py`
already find the start frame by step. History from before the cut (earlier frames
and event rows) is not in a slim file. Anything that needs it, such as the birth
steps of links formed before the cut, has to read the original.

`validate_slim_checkpoint.py FULL.h5 SLIM.h5 ARGV BINARY WORKDIR` resumes both
files side by side on one thread, which is deterministic. It then requires
every appended frame, link row, event row and the final `/resume` snapshot to
match bit for bit. 2026-10-29, `hph27m24_seed20260805.h5`, 200 steps saved
every 100: identical across 299 datasets, including 458 new link events and
7886 head transitions, with both the axial-stop binary (5e3bf3a6) and the
previous-law binary (c28e75eb).
