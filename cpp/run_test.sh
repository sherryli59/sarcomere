rm langevin_fix.out
nohup build/sarcomere  --n_actins=1200 --n_myosins=36 --filename=data/langevin_fix2.h5 --initial_structure=partial --Lx=12 --Ly=5 --Lz=3 --nsteps=5 --save_every=200 --n_fixed_myosins=5 > langevin_fix.out &
