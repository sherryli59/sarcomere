rm langevin_fix.out
nohup build/sarcomere  --n_actins=1400 --n_myosins=100 --filename=data/langevin_fix2.h5 --initial_structure=partial --Lx=6 --Ly=4 --Lz=3.2 --nsteps=10000 --save_every=200 --n_fixed_myosins=1 > langevin_fix.out &
wait
python visualize_traj3d.py --filename=data/langevin_fix2.h5 --print_frame=-1 --Lx=6 --Ly=4 --Lz=3.2