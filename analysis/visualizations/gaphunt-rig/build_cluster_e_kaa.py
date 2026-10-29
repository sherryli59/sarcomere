"""Cluster E only (2 myosins, 8 actins; membership fixed at the shared 2.0 s state), control vs
kappa_aa 3 -> 30, from 2.0 s to the latest COMMITTED frame of the running kappa_aa arm."""
import os, sys, json, h5py
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import build_cluster_traj as B
g = h5py.File("%s/data/r27kaa30_seed20260805.h5" % B.CAMP, "r", locking=False)
B.LAST = int(g["meta/committed_frames"][0]) - 1; g.close()
f0 = h5py.File("%s/data/r27ctrl_seed20260805.h5" % B.CAMP, "r", locking=False)
clusters, adj = B.find_clusters(f0, f0["actin/center"].shape[1], f0["myosin/center"].shape[1])
c = B.pick(clusters)[4]; assert (len(c["myosins"]), len(c["actins"])) == (2, 8)
runs = [B.build(arm, lab, c, "E", adj, None) for arm, lab in (("r27ctrl", "control"), ("r27kaa30", "kappa_aa 3 -> 30"))]
for r in runs:
    r["setupHtml"] = r["setupHtml"].replace("followed from 2.0 to 2.4 s", "followed from 2.0 s to %.3f s (latest committed frame of the running kappa_aa arm)" % (B.LAST * 0.005))
tpl = open(os.path.join(HERE, "cluster-traj-template.html"), encoding="ascii").read()
out = tpl.replace("/*__RUNS__*/", json.dumps(runs, separators=(",", ":"))).replace("__TITLE__", "Cluster E: control vs kappa_aa 30 (running)")
i0 = out.index("<script"); out = out[:i0] + "".join(ch if ord(ch) < 128 else "\\u%04x" % ord(ch) for ch in out[i0:])
assert out.isascii()
dst = os.path.join(HERE, "cluster-e-kaa30.html"); open(dst, "w", encoding="ascii").write(out)
print("wrote %s (%.2f MB), frames 400-%d (2.000-%.3f s)" % (dst, os.path.getsize(dst) / 1e6, B.LAST, B.LAST * 0.005))
