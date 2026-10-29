# Production configuration

- `prodargs.sh`: production argv (hph27 campaign, burial block rule-only, titin and myomesin off, Z-disc stop gated on ACTN2). Copied from `burialzone_20261020/prodargs.sh`.
- The compiled defaults in `src/main.cpp` equal these values, apart from the per-run seed, output file and step counts. Shell variables resolve to HPH=27, KILL=0.050, RATE=1, SLAB=false and CUTOFF=0.
- `sync_defaults.py` checks that the compiled defaults match `prodargs.sh` and exits 1 when they differ. `--apply` rewrites the defaults.
- `run_production.sh EXE OUT.h5 SEED [NSTEPS] [SAVE_EVERY] [flags...]` runs with those defaults. For example, the 54-head variant: `run_production.sh ./build/sarcomere out.h5 20260805 20000 200 --am_heads_per_half=54`.

After changing a production flag, update `prodargs.sh` and run `python sync_defaults.py --apply`.
