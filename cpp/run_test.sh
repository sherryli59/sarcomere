export OMP_NUM_THREADS=12
rm langevin_fix.out
numactl --interleave=all nohup build/sarcomere  --n_actins=500 --n_myosins=100 --filename=../data/debug.h5 --initial_structure=partial --Lx=5 --Ly=4 --Lz=5 --nsteps=100000 --save_every=200 --n_fixed_myosins=0 > langevin_fix.out &
wait
python ../analysis/visualize_traj.py --filename=../data/debug.h5 --print_frame=-1 --Lx=5 --Ly=4 --Lz=5 --frame_dir=new_data
