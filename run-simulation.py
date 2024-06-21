import actin 
import matplotlib.pyplot as plt
import torch
import numpy as np
import argparse

def parse():
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--n_filaments", type=int, default=100)
    argparser.add_argument("--n_bundles", type=int, default=10)
    argparser.add_argument("--n_crosslinks", type=int, default=50)
    argparser.add_argument("--boxL", type=float, default=30.)
    argparser.add_argument("--boxW", type=float, default=30.)
    argparser.add_argument("--n_steps", type=int, default=100000)
    argparser.add_argument("--dt", type=float, default=1e-3)
    argparser.add_argument("--D", type=float, default=1.)
    argparser.add_argument("--frame_dir", type=str, default="frames")
    argparser.add_argument("--seed", type=int, default=42)
    argparser.add_argument("--resume", action="store_true")
    return argparser.parse_args()

if __name__ == "__main__":
    args = parse()
    n_filaments = args.n_filaments
    n_bundles = args.n_bundles
    n_crosslinks = args.n_crosslinks
    n_steps = args.n_steps
    boxL = args.boxL
    boxW = args.boxW
    seed = args.seed   
    dt = args.dt     
    D = args.D
    frame_dir = args.frame_dir
    torch.manual_seed(seed)
    sarcomere = actin.SarcomereModel(n_filaments, n_bundles, n_crosslinks,boxW, boxL)
    if args.resume:
        traj = np.load("{}/trajectory.npy".format(frame_dir),allow_pickle=True).item()
        sarcomere.load_traj(traj)

    sarcomere.run_simulation(dt, D, n_steps, make_animation=True, frame_dir=frame_dir)
