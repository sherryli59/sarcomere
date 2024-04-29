import actin_gpu as actin 
import matplotlib.pyplot as plt
import torch
import numpy as np
if __name__ == "__main__":

    # n_filaments = 100
    # n_bundles = 10
    # n_crosslinks = 50
    # boxL = 30.

    # n_filaments = 10
    # n_bundles = 2
    # n_crosslinks = 5
    # boxL = 10.

    n_filaments = 20
    n_bundles = 4
    n_crosslinks = 10
    boxL = 15
    #fix random seed
    torch.manual_seed(42)
    sarcomere = actin.SarcomereModel(n_filaments, n_bundles, n_crosslinks, boxL, boxL)
    traj = np.load("frames_medium/trajectory.npy",allow_pickle=True).item()
    # #access the 'myosin' key in the dictionary traj
    sarcomere.load_traj(traj)
    #sarcomere.run_simulation(1e-3, 1., 100000, make_animation=True, frame_dir="frames_small")
    sarcomere.run_simulation(1e-3, 1., 200000, make_animation=True, frame_dir="frames_medium")
    #sarcomere.run_simulation(1e-3, 1., 200000, make_animation=True, frame_dir="frames")
