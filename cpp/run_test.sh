rm langevin_fix.out
nohup build/sarcomere  --n_actins=800 --n_myosins=36 --filename=data/langevin_fix2.h5 --initial_structure=partial --Lx=12 --Ly=5 --nsteps=30000 --save_every=200 --fix_myosin=5 > langevin_fix.out &
wait
python visualize_traj.py --filename="data/langevin_fix2.h5" --frame_dir="0206" --Lx=12 --Ly=8 --myosin_radius=0.2 --actin_length=1.0 --myosin_length=1.5


