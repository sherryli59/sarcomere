export OMP_NUM_THREADS=12
rm langevin_fix.out
DIMENSION=${1:-3}
numactl --interleave=all nohup build/sarcomere  --n_actins=2000 --n_myosins=200 --filename=../data/1000actin.h5 --Lx=5 --Ly=5 --Lz=5 --dimension=${DIMENSION} --nsteps=100000 --save_every=200 --n_fixed_myosins=0 --base_lifetime=0.0001 --resume > langevin_fix.out &
wait
python ../analysis/visualize_traj.py --filename=../data/1000actin.h5 --Lx=5 --Ly=5 --Lz=5 --frame_dir=1000actin --myosin_radius=0.15
