"""Why don't ACTN2-linked actins align with each other?  r27ctrl, 2.0 -> 2.4 s every 25 ms.
For every linked actin i (and separately the 8 Cluster E actins):
  * partners j (live ACTN2 links, status 1 or 2); partner mean direction n_p (Q-average)
  * angle to n_p, median angle to individual partners, partner spread (median angle of partners to n_p)
  * expected kappa_aa drive toward n_p from the per-link law (k = kappa_aa / reference_links,
    drive_i = sum_j k (u_i.u_j) u_j, projected on the tangent toward n_p) vs recorded tb_actn2_align
  * every recorded source projected on the same tangent (kT)
  * link attachment position: distance from the barbed end along each actin
  * partners' strong heads: fraction of partners with >=1 strong head; angle partner->its myosin;
    angle between i's primary myosin and the partner's myosin"""
import os, sys, json, collections, numpy as np, h5py
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import build_cluster_traj as B
KT = 1.380649e-23 * 310 / 1e-18
SRC = ["actn2_align", "actn2_lever", "actn2_perp", "actn2_axial", "actn2_contact", "aa_lateral", "aa_steric", "am_head", "am_steric"]
U = lambda d: d / np.linalg.norm(d, axis=-1, keepdims=True)
ang = lambda a, b: np.degrees(np.arccos(np.clip(np.abs(np.sum(a * b, -1)), 0, 1)))
ARM = sys.argv[1] if len(sys.argv) > 1 else "r27ctrl"; T0, T1 = (int(sys.argv[2]), int(sys.argv[3])) if len(sys.argv) > 3 else (400, 480)
f = h5py.File("%s/data/%s_seed20260805.h5" % (B.CAMP, ARM), "r", locking=False); a = dict(f["meta"].attrs)
nA = f["actin/center"].shape[1]; nM = f["myosin/center"].shape[1]; La = float(a["actin_length"])
k_link = float(a.get("kappa_aa", 3.0)) / float(a.get("aa_stiffness_reference_links", 3))
f0 = h5py.File("%s/data/r27ctrl_seed20260805.h5" % B.CAMP, "r", locking=False)
clusters, adj = B.find_clusters(f0, nA, nM); E = B.pick(clusters)[4]; EA = set(E["actins"])
k_link = float(a.get("kappa_aa", 3.0)) / float(a.get("aa_stiffness_reference_links", 3))
if ARM == "r27kaa30": k_link = 30.0 / 3.0   # meta@ keeps the original kappa_aa on a resume
rows = []
for t in range(T0, T1 + 1, 5):
    u = U(f["actin/direction"][t]); m = U(f["myosin/direction"][t])
    off = f["catch_bond/local_link_frame_offsets"][t]; o, c = int(off[1]), int(off[2])
    ids = f["catch_bond/local_link_identity"][o:o + c]; cs = f["catch_bond/local_link_coordinates"][o:o + c]
    part = collections.defaultdict(list)
    for r, q in zip(ids, cs):
        i, j = int(r[2]), int(r[3]); part[i].append((j, float(q[0]))); part[j].append((i, float(q[1])))
    hs = f["state/am_head_state_current"][t]; hb = f["state/am_head_bound_actin_current"][t].astype(int)
    mi, hi = np.nonzero(hs == 2); ai = hb[mi, hi]
    strong_my = collections.defaultdict(collections.Counter)
    for x, y in zip(ai, mi): strong_my[int(x)][int(y)] += 1
    tb = {s: np.stack([f["actin/tb_%s_%s" % (s, k)][t, :, 0] for k in "xyz"], 1) for s in SRC}
    for i, lst in part.items():
        J = np.array([j for j, _ in lst]); sj = np.array([s for _, s in lst])
        uj = u[J] * np.sign(u[J] @ u[i])[:, None]
        Qp = (uj[:, :, None] * uj[:, None, :]).mean(0); npd = np.linalg.eigh(Qp)[1][:, -1]; npd *= np.sign(npd @ u[i])
        tg = npd - u[i] * (npd @ u[i]); nt = np.linalg.norm(tg)
        if nt < 1e-9: continue
        tg /= nt
        exp_drive = float(np.sum([k_link * (u[i] @ u[j]) * (u[j] @ tg) for j in J])) / KT
        r = dict(t=t, i=i, inE=i in EA, nlinks=len(J), th_mean=float(ang(u[i], npd)), th_each=float(np.median(ang(u[i], uj))),
                 spread=float(np.median(ang(uj, npd[None]))), attach_from_barbed=float(np.median(sj + 0.5 * La)),
                 exp_align=exp_drive)
        for s in SRC: r[s] = float(tb[s][i] @ tg) / KT
        # partners' strong heads
        ps = [j for j in J if strong_my[int(j)]]
        r["partners_strong_frac"] = len(ps) / len(J)
        prim_i = max(strong_my[i], key=strong_my[i].get) if strong_my[i] else None
        r["self_strong"] = bool(strong_my[i])
        pm_ang = [float(ang(u[j], m[max(strong_my[int(j)], key=strong_my[int(j)].get)])) for j in ps]
        r["partner_to_its_myosin"] = float(np.median(pm_ang)) if pm_ang else np.nan
        if prim_i is not None and ps:
            r["myo_disagree"] = float(np.median([ang(m[prim_i], m[max(strong_my[int(j)], key=strong_my[int(j)].get)]) for j in ps]))
            r["same_myosin_frac"] = float(np.mean([prim_i in strong_my[int(j)] for j in ps]))
        else:
            r["myo_disagree"] = np.nan; r["same_myosin_frac"] = np.nan
        rows.append(r)
R = {k: np.array([r[k] for r in rows], float) for k in rows[0]}
def summ(sel, lab):
    s = sel
    print("== %s: %d actin-frames" % (lab, s.sum()))
    print("   ACTN2 links per actin %.1f; link attachment %.0f nm from the barbed end (median)" % (np.median(R["nlinks"][s]), 1000 * np.median(R["attach_from_barbed"][s])))
    print("   angle to partner MEAN %.1f deg | median angle to individual partners %.1f deg | partner spread around their mean %.1f deg"
          % (np.median(R["th_mean"][s]), np.median(R["th_each"][s]), np.median(R["spread"][s])))
    print("   kappa_aa drive toward partner mean: expected from per-link law %.0f kT, recorded tb_actn2_align %.0f kT (corr %.3f)"
          % (np.mean(R["exp_align"][s]), np.mean(R["actn2_align"][s]), np.corrcoef(R["exp_align"][s], R["actn2_align"][s])[0, 1]))
    ok = s & (np.abs(R["aa_steric"]) < 5000)
    print("   recorded drive toward partner mean, mean (sd) kT: " + ", ".join("%s %+.0f (%.0f)" % (k, R[k][ok].mean(), R[k][ok].std()) for k in SRC))
    print("   partners with >=1 strong head %.0f%%; actin itself strong-bound %.0f%%; partner to its own myosin %.1f deg; "
          "partner's myosin vs this actin's myosin %.1f deg; partner on the SAME myosin %.0f%%"
          % (100 * np.mean(R["partners_strong_frac"][s]), 100 * np.mean(R["self_strong"][s]), np.nanmedian(R["partner_to_its_myosin"][s]),
             np.nanmedian(R["myo_disagree"][s]), 100 * np.nanmean(R["same_myosin_frac"][s])))
    for lo, hi in ((0, 10), (10, 20), (20, 40), (40, 90)):
        b = s & (R["th_mean"] >= lo) & (R["th_mean"] < hi) & (np.abs(R["aa_steric"]) < 5000)
        if b.sum() < 10: continue
        print("     angle to partner mean %2d-%2d deg (n %5d): align %+5.0f lever %+6.0f axial %+6.0f lateral %+6.0f a-a steric %+6.0f heads %+5.0f"
              % (lo, hi, b.sum(), R["actn2_align"][b].mean(), R["actn2_lever"][b].mean(), R["actn2_axial"][b].mean(), R["aa_lateral"][b].mean(), R["aa_steric"][b].mean(), R["am_head"][b].mean()))
summ(np.ones(len(rows), bool), "all linked actins")
summ(R["inE"] > 0, "Cluster E actins")
json.dump({"k_link_pN_um": k_link, "rows": rows}, open("/scratch/users/xsli/sarcomere/gaphunt_20261023/torque_out/partner_an_%s_%d_%d.json" % (ARM, T0, T1), "w"))
