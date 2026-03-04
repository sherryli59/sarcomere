export OMP_NUM_THREADS=12
rm langevin_fix.out
DIMENSION=${1:-3}
numactl --interleave=all nohup build/sarcomere --directional=false --n_actins=1800 --n_myosins=1000 --filename=../data/1800actin.h5 --Lx=12 --Ly=2 --Lz=2 --initial_structure=myosin_x_noise --dimension=${DIMENSION} --nsteps=100000 --save_every=200 --n_fixed_myosins=0 --base_lifetime=0.0001 > langevin_fix.out &
wait
python ../analysis/visualize_traj.py --filename=../data/1800actin.h5 --Lx=12 --Ly=2 --Lz=2 --frame_dir=1800actin --myosin_radius=0.06