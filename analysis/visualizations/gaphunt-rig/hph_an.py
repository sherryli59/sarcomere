"""Head-count production runs (hph 27/54/81, sliding force fixed by msat=24).

Per sampled frame: actin & myosin nematic S; EXACT heads per actin-myosin
pair and distinct crowns per pair (from state/am_head_bound_actin_current);
partners per bound actin; median angle between an actin and the myosin it is
bound to; ACTN2 links and status-2 per actin; and the per-pair total motor
force, reconstructed head by head as compute_am_two_spring_head does it:
z = (c_a + clamp(target,+-La/2) u_a - c_m - anchor_axial u_m) . u_a,
strain = stroke - z, strong & strain>0 -> F = min(k*ms*strain, stall*ms),
ms = 1/min(bound heads on that pair, am_motor_saturation_heads).
"""
import h5py, numpy as np, math, json, os, sys
os.environ.setdefault("HDF5_USE_FILE_LOCKING","FALSE")
CAMP="/scratch/users/xsli/sarcomere/gaphunt_20261023/data"
def nem(d):
    u=d/np.maximum(np.linalg.norm(d,axis=1,keepdims=True),1e-12)
    Q=1.5*np.einsum("ni,nj->ij",u,u)/len(u)-0.5*np.eye(3); return float(np.linalg.eigvalsh(Q)[-1])
def anchor_table(hph, Lm, spokes=3, bare=0.10):
    bear=0.5*(1-bare); cc=max(1,(hph+spokes-1)//spokes); out=np.empty(2*hph); crown=np.empty(2*hph,int)
    for half in (0,1):
        sg=-1.0 if half==0 else 1.0
        for h in range(hph):
            out[half*hph+h]=sg*(0.5*bare*Lm+((h//spokes+0.5)/cc)*bear*Lm); crown[half*hph+h]=half*1000+h//spokes
    return out,crown
def frame_stats(f,t,a,anc,crn):
    box=np.array([a["box_x"],a["box_y"],a["box_z"]],float)
    La=float(a["actin_length"]); kmot=float(a["am_motor_k"]); stall=float(a["am_stall_force"])
    msat=int(a["am_motor_saturation_heads"])
    ac=f["actin/center"][t]; ad=f["actin/direction"][t]; ad=ad/np.linalg.norm(ad,axis=1,keepdims=True)
    mc=f["myosin/center"][t]; md=f["myosin/direction"][t]; md=md/np.linalg.norm(md,axis=1,keepdims=True)
    hs=f["state/am_head_state_current"][t]; hb=f["state/am_head_bound_actin_current"][t].astype(int)
    tg=f["state/am_head_target_axial_current"][t]; sp=f["state/am_head_stroke_progress_current"][t]
    bound=(hs==1)|(hs==2)
    mi,hi=np.nonzero(bound); ai=hb[mi,hi]
    if len(mi)==0:                      # production starts with no heads bound
        lk=f["actin/aa_local_link_count"][t][:,0]; s2=f["actin/aa_local_status2_count"][t][:,0]
        nan=float("nan")
        return dict(S_actin=nem(ad),S_myosin=nem(md),pairs=0,heads_per_pair=nan,heads_per_pair_median=nan,
            crowns_per_pair=nan,frac_pairs_3crowns=nan,partners_per_bound_actin=nan,pair_angle_median=nan,
            links_per_actin=float(lk.mean()),status2_per_actin=float(s2.mean()),strong_per_actin=0.0,
            max_heads_per_pair=0,pair_force_mean=0.0,pair_force_p99=0.0,pair_force_max=0.0)
    key=ai.astype(np.int64)*100000+mi
    uk,inv,cnt=np.unique(key,return_inverse=True,return_counts=True)
    pa=(uk//100000).astype(int); pm=(uk%100000).astype(int)
    crowns=np.array([len(set(crn[hi[inv==k]])) for k in range(len(uk))])
    # pair angle
    pang=np.degrees(np.arccos(np.clip(np.abs(np.sum(ad[pa]*md[pm],axis=1)),0,1)))
    # motor force per pair
    strong=hs[mi,hi]==2
    A=mc[mi]+anc[hi][:,None]*md[mi]
    d=ac[ai]-A; d=d-box*np.round(d/box)
    P=d+np.clip(tg[mi,hi],-La/2,La/2)[:,None]*ad[ai]
    z=np.sum(P*ad[ai],axis=1); strain=sp[mi,hi]-z
    ms=1.0/np.minimum(cnt[inv],msat)
    F=np.where(strong&(strain>0),np.minimum(kmot*ms*strain,stall*ms),0.0)
    Fpair=np.bincount(inv,weights=F,minlength=len(uk))
    partners=np.bincount(pa,minlength=ac.shape[0]); bound_actins=partners>0
    lk=f["actin/aa_local_link_count"][t][:,0]; s2=f["actin/aa_local_status2_count"][t][:,0]
    return dict(S_actin=nem(ad),S_myosin=nem(md),pairs=int(len(uk)),heads_per_pair=float(cnt.mean()),
        heads_per_pair_median=float(np.median(cnt)),crowns_per_pair=float(crowns.mean()),
        frac_pairs_3crowns=float((crowns>=3).mean()),partners_per_bound_actin=float(partners[bound_actins].mean()),
        pair_angle_median=float(np.median(pang)),links_per_actin=float(lk.mean()),status2_per_actin=float(s2.mean()),
        strong_per_actin=float((hs==2).sum()/ac.shape[0]),max_heads_per_pair=int(cnt.max()),
        pair_force_mean=float(Fpair[Fpair>0].mean()) if (Fpair>0).any() else 0.0,
        pair_force_p99=float(np.percentile(Fpair,99)),pair_force_max=float(Fpair.max()))
if __name__=="__main__":
    runs=[a.split(":") for a in sys.argv[1].split(",")]          # hph:seed,...
    every=int(sys.argv[2]) if len(sys.argv)>2 else 20
    out={}
    for hph,seed in runs:
        p="%s/hph%sm24_seed%s.h5"%(CAMP,hph,seed)
        f=h5py.File(p,"r",locking=False); a=dict(f["meta"].attrs)
        n=int(min(f["meta/committed_frames"][0],f["actin/center"].shape[0]))
        anc,crn=anchor_table(int(a["am_heads_per_half"]),float(a["myosin_length"]),int(float(a.get("am_spokes_per_crown",3))))
        fm=float(f["state/current_step"][1,0]-f["state/current_step"][0,0])*float(a["dt"])*1e3
        ts=list(range(0,n,every)); ts+= [n-1] if ts[-1]!=n-1 else []
        rows=[]
        for t in ts:
            r=frame_stats(f,t,a,anc,crn); r["t_ms"]=t*fm; rows.append(r)
        out["hph%s_%s"%(hph,seed)]=dict(hph=int(hph),seed=int(seed),frames=n,rows=rows)
        last=rows[-1]
        print("hph %s seed %s  t=%.0f ms: S actin %.3f myosin %.3f | heads/pair %.2f crowns/pair %.2f (>=3 crowns %.0f%%) partners %.2f | pair angle %.1f | links %.2f status2 %.2f | pair force mean %.2f p99 %.2f max %.2f pN | max heads/pair %d"
              %(hph,seed,last["t_ms"],last["S_actin"],last["S_myosin"],last["heads_per_pair"],last["crowns_per_pair"],
                100*last["frac_pairs_3crowns"],last["partners_per_bound_actin"],last["pair_angle_median"],last["links_per_actin"],
                last["status2_per_actin"],last["pair_force_mean"],last["pair_force_p99"],last["pair_force_max"],last["max_heads_per_pair"]))
    json.dump(out,open(sys.argv[3] if len(sys.argv)>3 else "hph_an.json","w"))
