#!/bin/bash
#SBATCH --partition=rotskoff
#SBATCH --account=rotskoff
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=8GB
#SBATCH --time=24:00:00
#SBATCH --array=1
#SBATCH --job-name=actn2_local_prod5s
#SBATCH --open-mode=append
#SBATCH --requeue
#SBATCH --output=/scratch/users/xsli/sarcomere/actn2_local_capture_production_A800_M320_5s_dt10us_20260819/logs/%x_%A_%a.out
#SBATCH --error=/scratch/users/xsli/sarcomere/actn2_local_capture_production_A800_M320_5s_dt10us_20260819/logs/%x_%A_%a.err

set -euo pipefail

production_root=/scratch/users/xsli/sarcomere/actn2_local_capture_production_A800_M320_5s_dt10us_20260819
wrapper=/scratch/users/xsli/sarcomere/actn2_local_capture_pilot_A800_M320_100ms_dt10us_20260819/bin/run_actn2_local_capture_pilot_binary_20260819.sh
base_launcher=/scratch/users/xsli/sarcomere/cellA_sparse_myomesin_A800_M320_5s_dt10us_20260818/source_snapshot_nobarrier/cell_a_sparse_myomesin_5s_20260818.sbatch
run_kind=thick50_currentdist_valency3_actn2_localcapture_pool360_production5s
wrapper_hash=$(sha256sum "${wrapper}" | awk '{print $1}')

export CELL_A_EXE="${wrapper}"
export CELL_A_EXE_SHA256="${wrapper_hash}"
export CELL_A_RUN_KIND="${run_kind}"
export CELL_A_NSTEPS=500000
export CELL_A_SAVE_EVERY=500
export CELL_A_SNAPSHOT=/scratch/users/xsli/sarcomere/cellA_sparse_myomesin_A800_M320_5s_dt10us_20260818/source_snapshot_nobarrier

bash "${base_launcher}"
