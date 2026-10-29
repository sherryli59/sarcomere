"""Angle vs AXIAL SPAN of a pair's bound heads, hph27 vs hph81.
Span = max - min anchor_axial among the pair's bound heads (one myosin half
holds them, so this is the lever arm of the aligning couple).  If angle vs
span is one curve across head counts, span is the variable and crown count was
only its proxy.  Also the time test by span: pairs holding span >= 200 nm for
>= 50 ms -- does the angle fall?"""
import h5py, numpy as np, json, os, sys, collections
os.environ.setdefault("HDF5_USE_FILE_LOCKING","FALSE")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hph_an import anchor_table
CAMP="/scratch/users/xsli/sarcomere/gaphunt_20261023/data"
def ang(a,b): return np.degrees(np.arccos(np.clip(np.abs(np.sum(a*b,axis=-1)),0,1)))
def pair_spans(f,t,anc):
    hs=f["state/am_head_state_current"][t]; hb=f["state/am_head_bound_actin_current"][t].astype(int)
    mi,hi=np.nonzero((hs==1)|(hs==2)); ai=hb[mi,hi]
    lo=collections.defaultdict(lambda: 9.0); hi_=collections.defaultdict(lambda: -9.0); n=collections.Counter()
    for a,m,h in zip(ai,mi,hi):
        k=(a,m); z=anc[h]; lo[k]=min(lo[k],z); hi_[k]=max(hi_[k],z); n[k]+=1
    return {k:(1000*(hi_[k]-lo[k]),n[k]) for k in n}
BINS=[(0,1,"0 (one crown)"),(1,100,"1-100"),(100,200,"100-200"),(200,300,"200-300"),(300,450,"300-450"),(450,2000,">450")]
out={}
for tag,path in (("hph27",CAMP+"/hph27m24_seed20260805.h5"),("hph81",CAMP+"/hph81m24_seed20260805.h5")):
    f=h5py.File(path,"r",locking=False); a=dict(f["meta"].attrs)
    n=int(min(f["meta/committed_frames"][0],f["actin/center"].shape[0])); t=n-1
    anc,_=anchor_table(int(a["am_heads_per_half"]),float(a["myosin_length"]))
    ua=f["actin/direction"][t]; ua=ua/np.linalg.norm(ua,axis=1,keepdims=True)
    um=f["myosin/direction"][t]; um=um/np.linalg.norm(um,axis=1,keepdims=True)
    P=pair_spans(f,t,anc)
    sp=np.array([v[0] for v in P.values()]); an=np.array([float(ang(ua[k[0]],um[k[1]])) for k in P])
    out[tag]=dict(span_median=float(np.median(sp)),bins=[(lab,float(np.median(an[(sp>=lo)&(sp<hi)])) if ((sp>=lo)&(sp<hi)).sum()>=20 else None,int(((sp>=lo)&(sp<hi)).sum())) for lo,hi,lab in BINS])
    # time test by span, last 1 s
    T0=max(0,n-200); hist=collections.defaultdict(list)
    for tt in range(T0,n):
        uat=f["actin/direction"][tt]; uat=uat/np.linalg.norm(uat,axis=1,keepdims=True)
        umt=f["myosin/direction"][tt]; umt=umt/np.linalg.norm(umt,axis=1,keepdims=True)
        for k,(s,nh) in pair_spans(f,tt,anc).items(): hist[k].append((tt,s,float(ang(uat[k[0]],umt[k[1]]))))
    res={}
    for thr in (200,300):
        eps=[]
        for k,seq in hist.items():
            run=[]
            for (tt,s,a_) in seq:
                if s>=thr and (not run or tt==run[-1][0]+1): run.append((tt,s,a_))
                else:
                    if len(run)>=10: eps.append(run)
                    run=[(tt,s,a_)] if s>=thr else []
            if len(run)>=10: eps.append(run)
        if not eps: res[thr]=None; continue
        st=np.array([e[0][2] for e in eps]); en=np.array([e[-1][2] for e in eps]); ln=np.array([5.0*len(e) for e in eps])
        hi_=st>=20
        res[thr]=dict(n=len(eps),length=float(np.median(ln)),start=float(np.median(st)),end=float(np.median(en)),
                      fell=float(((st-en)>5).mean()),n20=int(hi_.sum()),start20=float(np.median(st[hi_])) if hi_.any() else None,
                      end20=float(np.median(en[hi_])) if hi_.any() else None,fell20=float(((st[hi_]-en[hi_])>5).mean()) if hi_.any() else None)
    out[tag]["episodes"]=res
json.dump(out,open("span_an.json","w"),indent=1)
print("median angle (deg) to the bound myosin by the pair's axial head span (nm); n pairs in brackets")
print("%-16s"%"span"+"".join("%22s"%t for t in out))
for i,(lo,hi,lab) in enumerate(BINS):
    print("%-16s"%lab+"".join("%22s"%(("%.1f (%d)"%(out[t]["bins"][i][1],out[t]["bins"][i][2])) if out[t]["bins"][i][1] is not None else "- (%d)"%out[t]["bins"][i][2]) for t in out))
print("%-16s"%"median span"+"".join("%22s"%("%.0f nm"%out[t]["span_median"]) for t in out))
for t in out:
    for thr,r in out[t]["episodes"].items():
        if r: print("%s episodes holding span >= %d nm for >= 50 ms: n %d, median %.0f ms, angle %.1f -> %.1f, %.0f%% fell >5 deg | started >20 deg: n %d, %s -> %s, %s fell"
                    %(t,thr,r["n"],r["length"],r["start"],r["end"],100*r["fell"],r["n20"],
                      "%.1f"%r["start20"] if r["start20"] else "-","%.1f"%r["end20"] if r["end20"] else "-",
                      "%.0f%%"%(100*r["fell20"]) if r["fell20"] is not None else "-"))
