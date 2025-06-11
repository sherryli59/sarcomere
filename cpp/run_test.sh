export OMP_NUM_THREADS=32
rm langevin_fix.out
numactl --interleave=all nohup build/sarcomere  --n_actins=1400 --n_myosins=100 --filename=data/debug.h5 --initial_structure=partial --Lx=4.5 --Ly=4.5 --Lz=4.5 --nsteps=100000 --save_every=200 --n_fixed_myosins=0 > langevin_fix.out &
wait
python visualize_traj3d.py --filename=data/debug.h5 --print_frame=-1 --Lx=4.5 --Ly=4.5 --Lz=4.5