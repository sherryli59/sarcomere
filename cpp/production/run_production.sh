#!/bin/bash
# Production run with the compiled defaults, which equal the production
# argv in prodargs.sh (checked by sync_defaults.py).  Only per-run values are
# passed.  Under Slurm, wrap with srun.
# usage: run_production.sh EXE OUTPUT.h5 SEED [NSTEPS] [SAVE_EVERY] [extra flags...]
set -euo pipefail
EXE=${1:?simulator executable}; OUT=${2:?output .h5}; SEED=${3:?seed}
NSTEPS=${4:-20000}; SAVE_EVERY=${5:-200}
shift $(( $# < 5 ? $# : 5 ))
exec "$EXE" --seed="$SEED" --filename="$OUT" --nsteps="$NSTEPS" \
     --save_every="$SAVE_EVERY" "$@"
