"""Is the head-count effect really sign-flipped by the alignment potential?

Paired per-seed interaction:
    I(seed) = [metric(hph27,k=30) - metric(hph9,k=30)] - [metric(hph27,k=0) - metric(hph9,k=0)]
with an SE over the 4 shared seeds. Reported for S (near-tautological under a
nematic coupling) and for ACTN2 links / eligible pairs (the stronger measures).
All at a matched 80 ms.
"""
import os, glob, collections
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import numpy as np, h5py
FULL = "/scratch/users/xsli/sarcomere/barbed_capsule_normal_full_20261002/data"
LAD  = "/scratch/users/xsli/sarcomere/hph_ladder_20261005/data"
AL   = "/scratch/users/xsli/sarcomere/amalign_20261005/data"
CELL = {"A":   f"{FULL}/*kaa3_seed*_5000ms_capsule_normal.h5",
        "A9":  f"{AL}/*kaa3_seed*_5000ms_align9.h5",
        "C":   f"{LAD}/*kaa3_seed*_5000ms_hph27.h5",
        "C27": f"{AL}/*kaa3_seed*_5000ms_align27.h5"}
def read(pat):
    out = {}
    for p in sorted(glob.glob(pat)):
        seed = p.split("seed")[1][:8]
        with h5py.File(p, "r", locking=False) as f:
            a = f["meta"].attrs
            cut = float(a["aa_preferred_centerline"]) + float(a["aa_capture_margin"])
            maxdot = float(a["aa_association_max_direction_dot"]); L = float(a["actin_length"])
            box = np.array([a["box_x"], a["box_y"], a["box_z"]], float)
            cf = int(np.asarray(f["meta/committed_frames"]).reshape(-1)[-1])
            cf = min(cf, f["actin/center"].shape[0], f["state/current_step"].shape[0])
            st = np.asarray(f["state/current_step"][:cf]).reshape(-1)
            k = int(np.searchsorted(st, 80*100, "right")) - 1
            if k < 1: continue
            u = np.asarray(f["actin/direction"][k], float); u /= np.linalg.norm(u, axis=1, keepdims=True)
            Q = 1.5*np.einsum("ia,ib->ab", u, u)/u.shape[0] - 0.5*np.eye(3)
            S = float(np.linalg.eigvalsh(Q).max())
            links = float(np.asarray(f["state/aa_pairs_current_count"][k]).reshape(-1)[-1])
            ac = np.asarray(f["actin/center"][k], float); barbed = ac - 0.5*L*u
            d = barbed[:, None, :] - barbed[None, :, :]; d -= box*np.rint(d/box)
            dist = np.linalg.norm(d, axis=2); np.fill_diagonal(dist, 1e9)
            dot = np.einsum("ik,jk->ij", u, u)
            elig = int(((dist < cut) & (dot < maxdot)).sum())//2
            out[seed] = dict(S=S, links=links, elig=float(elig))
    return out
D = {k: read(v) for k, v in CELL.items()}
seeds = sorted(set.intersection(*[set(v) for v in D.values()]))
print("shared seeds: %s\n" % ",".join(seeds))
for key, nm in (("S", "nematic S"), ("links", "ACTN2 links"), ("elig", "eligible pairs")):
    e0 = np.array([D["C"][s][key]  - D["A"][s][key]  for s in seeds])   # head effect, k=0
    e30 = np.array([D["C27"][s][key] - D["A9"][s][key] for s in seeds]) # head effect, k=30
    I = e30 - e0
    se = lambda v: v.std(ddof=1)/np.sqrt(v.size)
    print("%-15s head effect (hph27 - hph9)" % nm)
    print("    at k=0 : %+9.3f  (SE %7.3f, t=%6.2f)  %s" % (e0.mean(), se(e0), e0.mean()/se(e0),
          "HURTS" if e0.mean() < 0 else "helps"))
    print("    at k=30: %+9.3f  (SE %7.3f, t=%6.2f)  %s" % (e30.mean(), se(e30), e30.mean()/se(e30),
          "HURTS" if e30.mean() < 0 else "helps"))
    print("    INTERACTION %+9.3f  (SE %7.3f, t=%6.2f)   sign flipped? %s\n"
          % (I.mean(), se(I), I.mean()/se(I), "YES" if e0.mean()*e30.mean() < 0 else "NO - same sign"))
