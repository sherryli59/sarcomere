#!/bin/bash
#SBATCH -p normal,hns,owners
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --job-name=build_sarcomere
#SBATCH -o /scratch/users/xsli/sarcomere/logs0624_myo1800/build_%j.out

set -e
eval "$(micromamba shell hook --shell bash)"
micromamba activate sarcomere
cmake --build /home/groups/rotskoff/sherry/sarcomere/sarcomere/cpp/build -j 8
echo "BUILD_EXIT_OK"
