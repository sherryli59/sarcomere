export OMP_NUM_THREADS=12
rm langevin_fix.out
numactl --interleave=all nohup build/sarcomere  --n_actins=2000 --n_myosins=250 --filename=../data/1000actin.h5 --initial_structure=partial --Lx=5 --Ly=5 --Lz=5 --nsteps=100000 --save_every=200 --n_fixed_myosins=0 --resume > langevin_fix.out &
wait
python ../analysis/visualize_traj.py --filename=../data/1000actin.h5 --print_frame=-1 --Lx=5 --Ly=5 --Lz=5 --frame_dir=1000actin --myosin_radius=0.15
