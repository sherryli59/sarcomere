export OMP_NUM_THREADS=12
rm langevin_fix.out
numactl --interleave=all nohup build/sarcomere  --n_actins=1500 --n_myosins=150 --filename=../data/debug.h5 --initial_structure=partial --Lx=5 --Ly=4 --Lz=4 --nsteps=110000 --save_every=200 --n_fixed_myosins=0 --resume > langevin_fix.out &
wait
python ../analysis/visualize_traj.py --filename=../data/debug.h5 --print_frame=-1 --Lx=5 --Ly=4 --Lz=4 --frame_dir=new_data
