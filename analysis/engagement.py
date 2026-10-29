"""Engagement ratios: who is actually participating in BOTH networks?

  - actins with any ACTN2 link / with a status-2 (mature) link
  - actins with a strong myosin head
  - actins with BOTH  <- the doubly-engaged population a sarcomere needs
  - myosins with a strong head
  - myosins COUPLED: a strong head on an actin that itself carries ACTN2
    (and the stricter version: on an actin carrying a status-2 link)

Link identity layout is the VISUALIZATION_PIPELINE.md contract (actin IDs in
columns 2 and 3, status in column 4; offsets rows are (step, offset, count)).
The attrs are not persisted in these files, so it is ASSERTED below.

Reported at a matched 600 ms (all five arms have passed it) and at 5000 ms for
the arms that are complete.
"""
import os, glob, collections
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import numpy as np, h5py
D = "/scratch/users/xsli/sarcomere"
ARMS = [("A   hph9/spk3  k=0 ", f"{D}/barbed_capsule_normal_full_20261002/data/*kaa3_seed*_5000ms_capsule_normal.h5"),
        ("B   hph9/spk1  k=0 ", f"{D}/hph_ladder_20261005/data/*kaa3_seed*_5000ms_spokes1.h5"),
        ("C   hph27/spk3 k=0 ", f"{D}/hph_ladder_20261005/data/*kaa3_seed*_5000ms_hph27.h5"),
        ("A9  hph9/spk3  k=30", f"{D}/amalign_20261005/data/*kaa3_seed*_5000ms_align9.h5"),
        ("C27 hph27/spk3 k=30", f"{D}/amalign_20261005/data/*kaa3_seed*_5000ms_align27.h5")]

def engagement(p, t_ms):
    with h5py.File(p, "r", locking=False) as f:
        cf = int(np.asarray(f["meta/committed_frames"]).reshape(-1)[-1])
        cf = min(cf, f["actin/center"].shape[0], f["state/current_step"].shape[0],
                 f["state/am_head_state_current"].shape[0])
        st = np.asarray(f["state/current_step"][:cf]).reshape(-1)
        k = int(np.searchsorted(st, t_ms*100, "right")) - 1
        if k < 1: return None
        # GUARD: searchsorted CLAMPS to the last frame when the run has not
        # reached t_ms, which silently mislabels a 1350 ms frame as 5000 ms.
        if float(st[k]) < t_ms*100 * 0.999:
            return None
        nA = f["actin/center"].shape[1]; nM = f["myosin/center"].shape[1]
        ident = f["catch_bond/local_link_identity"]
        off = np.asarray(f["catch_bond/local_link_frame_offsets"])
        j = int(np.searchsorted(off[:, 0], st[k], "left"))
        if j >= off.shape[0] or int(off[j, 0]) != int(st[k]): return None
        s0 = int(off[j, 1]); n = int(off[j, 2])
        has_aa = np.zeros(nA, bool); has_aa2 = np.zeros(nA, bool)
        if n > 0:
            blk = np.asarray(ident[s0:s0+n])
            i1 = blk[:, 2].astype(int); i2 = blk[:, 3].astype(int); stt = blk[:, 4].astype(int)
            assert i1.max() < nA and i2.max() < nA and i1.min() >= 0, "actin IDs not in cols 2,3"
            assert set(np.unique(stt).tolist()) <= {0, 1, 2}, "col 4 is not status: %s" % np.unique(stt)
            has_aa[i1] = True; has_aa[i2] = True
            m2 = stt == 2
            has_aa2[i1[m2]] = True; has_aa2[i2[m2]] = True
        hs = np.asarray(f["state/am_head_state_current"][k])
        ha = np.asarray(f["state/am_head_bound_actin_current"][k]).astype(int)
        strong = (hs == 2) & (ha >= 0)
        a_strong = np.zeros(nA, bool); m_strong = np.zeros(nM, bool)
        m_coupled = np.zeros(nM, bool); m_coupled2 = np.zeros(nM, bool)
        mi, hh = np.nonzero(strong)
        for m_, h_ in zip(mi.tolist(), hh.tolist()):
            ai = int(ha[m_][h_])
            a_strong[ai] = True; m_strong[m_] = True
            if has_aa[ai]: m_coupled[m_] = True
            if has_aa2[ai]: m_coupled2[m_] = True
        return dict(t=float(st[k])*1e-5*1e3,
                    a_aa=has_aa.mean(), a_aa2=has_aa2.mean(), a_strong=a_strong.mean(),
                    a_both=(has_aa & a_strong).mean(), a_both2=(has_aa2 & a_strong).mean(),
                    m_strong=m_strong.mean(), m_coup=m_coupled.mean(), m_coup2=m_coupled2.mean())

for T in (600, 1300, 5000):
    print("\n################ matched t = %d ms ################" % T)
    print("  %-20s %5s | %6s %6s %6s | %6s %6s | %6s %6s %6s"
          % ("arm", "n", "aa", "aa2", "strg", "both", "both2", "myStr", "coup", "coup2"))
    print("  %s" % ("-"*96))
    for tag, pat in ARMS:
        rows = [r for r in (engagement(p, T) for p in sorted(glob.glob(pat))) if r]
        if not rows: print("  %-20s   --  NOT REACHED %d ms yet" % (tag, T)); continue
        g = lambda k: float(np.mean([r[k] for r in rows]))
        print("  %-20s %5d | %5.1f%% %5.1f%% %5.1f%% | %5.1f%% %5.1f%% | %5.1f%% %5.1f%% %5.1f%%   [t=%.0f ms]"
              % (tag, len(rows), 100*g("a_aa"), 100*g("a_aa2"), 100*g("a_strong"),
                 100*g("a_both"), 100*g("a_both2"),
                 100*g("m_strong"), 100*g("m_coup"), 100*g("m_coup2"), g("t")))
print("""
  aa    = actins with any ACTN2 link        both  = actins with ACTN2 AND a strong head
  aa2   = actins with a status-2 link       both2 = actins with status-2 AND a strong head
  strg  = actins with a strong head         myStr = myosins with any strong head
  coup  = myosins with a strong head on an ACTN2-carrying actin
  coup2 = same, requiring the actin to carry a STATUS-2 link""")
