"""Which link of the chain fails in production?  A-A -> A-M -> M-M.

Part A (final frame): ACTN2-linked actin pairs' angle; A-M pair angle vs crowns
and vs partner count; how much an actin's partner myosins disagree; M-M angle
for myosins sharing an actin vs random myosin pairs.
Part B (last 1 s, every frame): follow every actin-myosin pair; for each
continuous episode with >= 3 crowns lasting >= 50 ms, does its angle fall?
"""
import h5py, numpy as np, json, os, sys, collections
os.environ.setdefault("HDF5_USE_FILE_LOCKING","FALSE")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hph_an import anchor_table
CAMP="/scratch/users/xsli/sarcomere/gaphunt_20261023/data"
def ang(a,b): return np.degrees(np.arccos(np.clip(np.abs(np.sum(a*b,axis=-1)),0,1)))
def pairs_at(f,t,crn):
    hs=f["state/am_head_state_current"][t]; hb=f["state/am_head_bound_actin_current"][t].astype(int)
    mi,hi=np.nonzero((hs==1)|(hs==2)); ai=hb[mi,hi]
    d=collections.defaultdict(set); n=collections.Counter()
    for a,m,h in zip(ai,mi,hi): d[(a,m)].add(crn[h]); n[(a,m)]+=1
    return {k:(len(v),n[k]) for k,v in d.items()}
def unitdir(f,grp,t):
    d=f[grp+"/direction"][t]; return d/np.linalg.norm(d,axis=1,keepdims=True)
def med(x): return float(np.median(x)) if len(x) else float("nan")
out={}
for tag,path in (("hph81",CAMP+"/hph81m24_seed20260805.h5"),("hph27",CAMP+"/hph27m24_seed20260805.h5")):
    f=h5py.File(path,"r",locking=False); a=dict(f["meta"].attrs)
    n=int(min(f["meta/committed_frames"][0],f["actin/center"].shape[0])); t=n-1
    _,crn=anchor_table(int(a["am_heads_per_half"]),float(a["myosin_length"]))
    ua=unitdir(f,"actin",t); um=unitdir(f,"myosin",t)
    P=pairs_at(f,t,crn)
    # A-A: ACTN2-linked actin pairs
    off=f["catch_bond/local_link_frame_offsets"][t]; o,c=int(off[1]),int(off[2])
    ids=f["catch_bond/local_link_identity"][o:o+c]
    aa=ang(ua[ids[:,2].astype(int)],ua[ids[:,3].astype(int)]); st=ids[:,4].astype(int)
    rng=np.random.default_rng(1); ri=rng.integers(0,ua.shape[0],(20000,2)); aa_rand=ang(ua[ri[:,0]],ua[ri[:,1]])
    # A-M by crowns, and by partner count
    partners=collections.defaultdict(list)
    for (ai_,mi_),(cr,nh) in P.items(): partners[ai_].append(mi_)
    rows=[(ai_,mi_,cr,len(partners[ai_]),float(ang(ua[ai_],um[mi_]))) for (ai_,mi_),(cr,nh) in P.items()]
    R=np.array(rows,float)
    def am_bin(mask): return med(R[mask,4]), int(mask.sum())
    by_crown={lab:am_bin(m) for lab,m in (("1",R[:,2]==1),("2",R[:,2]==2),("3",R[:,2]==3),("4-5",(R[:,2]>=4)&(R[:,2]<=5)),("6+",R[:,2]>=6))}
    solo=R[:,3]==1
    solo_by_crown={lab:am_bin(solo&m) for lab,m in (("1-2",R[:,2]<=2),("3+",R[:,2]>=3))}
    multi_by_crown={lab:am_bin((~solo)&m) for lab,m in (("1-2",R[:,2]<=2),("3+",R[:,2]>=3))}
    # partner disagreement: for actins with >=2 partners, median pairwise angle among its myosins
    dis=[]; best=[]; worst=[]
    for ai_,ms in partners.items():
        if len(ms)<2: continue
        v=um[ms]; pa=[float(ang(v[i],v[j])) for i in range(len(ms)) for j in range(i+1,len(ms))]
        dis.append(np.median(pa)); aa_=ang(ua[ai_][None,:],v); best.append(aa_.min()); worst.append(aa_.max())
    # M-M: myosins sharing an actin vs random myosin pairs
    sh=[float(ang(um[ms[i]],um[ms[j]])) for ms in partners.values() for i in range(len(ms)) for j in range(i+1,len(ms))]
    rj=rng.integers(0,um.shape[0],(20000,2)); mm_rand=ang(um[rj[:,0]],um[rj[:,1]])
    out[tag]=dict(t_ms=t*5.0,
        AA=dict(status1=med(aa[st==1]),n1=int((st==1).sum()),status2=med(aa[st==2]),n2=int((st==2).sum()),random=med(aa_rand)),
        AM_by_crown=by_crown, AM_single_partner=solo_by_crown, AM_multi_partner=multi_by_crown,
        partner_disagreement=dict(median=med(dis),n=len(dis),best=med(best),worst=med(worst)),
        MM=dict(sharing_an_actin=med(sh),n=len(sh),random=med(mm_rand)))
    # Part B: episodes with >=3 crowns, last 1 s, every frame (hph81 only -- slowest part)
    if tag=="hph81":
        T0=max(0,n-200); hist=collections.defaultdict(list)
        for tt in range(T0,n):
            uat=unitdir(f,"actin",tt); umt=unitdir(f,"myosin",tt); Pt=pairs_at(f,tt,crn)
            for (ai_,mi_),(cr,nh) in Pt.items():
                hist[(ai_,mi_)].append((tt,cr,float(ang(uat[ai_],umt[mi_]))))
        eps=[]
        for k,seq in hist.items():
            run=[]
            for (tt,cr,an) in seq:
                if cr>=3 and (not run or tt==run[-1][0]+1): run.append((tt,cr,an))
                else:
                    if len(run)>=10: eps.append(run)
                    run=[(tt,cr,an)] if cr>=3 else []
            if len(run)>=10: eps.append(run)
        st_=np.array([e[0][2] for e in eps]); en=np.array([e[-1][2] for e in eps]); mn=np.array([min(x[2] for x in e) for e in eps])
        ln=np.array([len(e)*5.0 for e in eps])
        hi_=st_>=20
        out[tag]["episodes"]=dict(n=len(eps),length_ms_median=med(ln),start=med(st_),end=med(en),min=med(mn),
            frac_fell_5deg=float(((st_-en)>5).mean()) if len(eps) else float("nan"),
            started_above_20=dict(n=int(hi_.sum()),start=med(st_[hi_]),end=med(en[hi_]),min=med(mn[hi_]),
                frac_fell_5deg=float(((st_[hi_]-en[hi_])>5).mean()) if hi_.any() else float("nan")))
json.dump(out,open("chain_an.json","w"),indent=1)
for tag,r in out.items():
    print("==== %s, t = %.0f ms"%(tag,r["t_ms"]))
    A=r["AA"]; print("A-A  ACTN2-linked actin pairs: status-1 %.1f deg (n %d), status-2 %.1f deg (n %d); random actin pair %.1f"%(A["status1"],A["n1"],A["status2"],A["n2"],A["random"]))
    print("A-M  angle by crowns on the pair: "+", ".join("%s: %.1f (n %d)"%(k,v[0],v[1]) for k,v in r["AM_by_crown"].items()))
    print("     single-partner actins: "+", ".join("%s crowns %.1f (n %d)"%(k,v[0],v[1]) for k,v in r["AM_single_partner"].items())
          +" | multi-partner: "+", ".join("%s crowns %.1f (n %d)"%(k,v[0],v[1]) for k,v in r["AM_multi_partner"].items()))
    D=r["partner_disagreement"]; print("     actins with >=2 partners (n %d): partner myosins disagree by %.1f deg; actin to best partner %.1f, to worst %.1f"%(D["n"],D["median"],D["best"],D["worst"]))
    M=r["MM"]; print("M-M  myosins sharing an actin %.1f deg (n %d) vs random myosin pairs %.1f"%(M["sharing_an_actin"],M["n"],M["random"]))
    if "episodes" in r:
        E=r["episodes"]; H=E["started_above_20"]
        print("TIME pair episodes with >=3 crowns for >=50 ms (last 1 s): n %d, median length %.0f ms, angle start %.1f -> end %.1f (min %.1f), %.0f%% fell >5 deg"
              %(E["n"],E["length_ms_median"],E["start"],E["end"],E["min"],100*E["frac_fell_5deg"]))
        print("     of those that STARTED above 20 deg: n %d, %.1f -> %.1f (min %.1f), %.0f%% fell >5 deg"%(H["n"],H["start"],H["end"],H["min"],100*H["frac_fell_5deg"]))
