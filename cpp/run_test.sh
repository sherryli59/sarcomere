rm langevin_fix.out
nohup build/sarcomere  --n_actins=1400 --n_myosins=40 --filename=data/langevin_fix2.h5 --initial_structure=partial --Lx=7 --Ly=5 --Lz=3 --nsteps=100000 --save_every=200 --n_fixed_myosins=5 > langevin_fix.out &
