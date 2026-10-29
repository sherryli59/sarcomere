"""MP4 of Cluster E (2 myosins, 8 actins), control arm, 2.0 -> 2.4 s every 5 ms."""
import os, sys, h5py
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import build_cluster_traj as B, render_mp4 as R
f0 = h5py.File("%s/data/r27ctrl_seed20260805.h5" % B.CAMP, "r", locking=False)
clusters, adj = B.find_clusters(f0, f0["actin/center"].shape[1], f0["myosin/center"].shape[1])
chosen = B.pick(clusters); c = chosen[4]
assert (len(c["myosins"]), len(c["actins"])) == (2, 8), (len(c["myosins"]), len(c["actins"]))
run, raw = B.build("r27ctrl", "control", c, "E", adj, None, return_raw=True)
out = os.path.join(HERE, "videos", "cluster_E_control.mp4"); os.makedirs(os.path.dirname(out), exist_ok=True)
R.render(raw, out, fps=12)
print("wrote", out, os.path.getsize(out) / 1e6, "MB, frames", len(raw), "myosins", c["myosins"], "actins", c["actins"])
