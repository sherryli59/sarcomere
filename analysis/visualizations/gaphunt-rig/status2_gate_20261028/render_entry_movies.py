"""Extract and animate newly mature links entering a myosin exclusion cylinder.

Uses the original production trajectory, not the diagnostic continuation.
Every animation pose is an actual saved frame (5 ms sampling); repeated frames
slow playback, with no coordinate interpolation.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="status2-movies-mpl-"))
import h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from audit_new_mature import blocked_connectors, unit, wrap

ARM = "r27zt0npg130"
R = .130
ACTIN_COLORS = ("#168d82", "#dd8154")
CLEAR, BLOCKED = "#7655ae", "#ce3349"


def geometry(p0, p1, center, axis, box, length):
    d = wrap(p1-p0, box)
    mid = wrap(p0+.5*d-center, box)
    p = mid-.5*d
    z0, dz = float(p@axis), float(d@axis)
    z1 = z0+dz
    clearance = max(min(z0,z1)-length/2, -length/2-max(z0,z1))
    rp, rd = p-z0*axis, d-dz*axis
    vv = rd@rd
    tfull = np.clip(-rp@rd/vv,0,1) if vv>0 else 0.0
    radial_all = float(np.linalg.norm(rp+tfull*rd))
    if dz == 0:
        lo,hi = (0.,1.) if abs(z0)<=length/2 else (1.,0.)
    else:
        a,b=(-length/2-z0)/dz,(length/2-z0)/dz
        lo,hi=max(0.,min(a,b)),min(1.,max(a,b))
    if lo>hi:
        return dict(blocked=False,axial_overlap=False,rho=None,
                    rho_all=radial_all,clearance=clearance,depth=None,point=None)
    t=np.clip(-rp@rd/vv,lo,hi) if vv>0 else lo
    rho=float(np.linalg.norm(rp+t*rd))
    return dict(blocked=bool(rho<=R),axial_overlap=True,rho=rho,
                rho_all=radial_all,clearance=clearance,
                depth=float(length/2-abs(z0+t*dz)),point=p+t*d)


def load_trajectory(path):
    with h5py.File(path,"r",locking=False) as f:
        attrs=dict(f["meta"].attrs)
        box=np.array([attrs["box_x"],attrs["box_y"],attrs["box_z"]])
        ec={k:i for i,k in enumerate(attrs["aa_local_link_events_v3_columns"].decode().split(","))}
        ev=f["catch_bond/local_link_events_v3"][:]
        birth=ev[(ev[:,ec["old_status"]]==0)&(ev[:,ec["new_status"]]==1)&(ev[:,ec["step"]]>200000)]
        born={int(r[ec["link_id"]]):int(r[ec["step"]]) for r in birth}
        mats=ev[(ev[:,ec["old_status"]]==1)&(ev[:,ec["new_status"]]==2)&(ev[:,ec["step"]]>200000)]
        mat={int(r[ec["link_id"]]):int(r[ec["step"]]) for r in mats}
        death={int(r[ec["link_id"]]):dict(step=int(r[ec["step"]]),reason=int(r[ec["detachment_reason"]]))
               for r in ev[(ev[:,ec["new_status"]]==0)&(ev[:,ec["step"]]>200000)]}
        nf=int(f["meta/committed_frames"][0])
        frames={}
        off=f["catch_bond/local_link_frame_offsets"]
        for frame in range(400,nf):
            start,count=map(int,off[frame,1:3])
            ident=f["catch_bond/local_link_identity"][start:start+count]
            sites=f["catch_bond/local_link_coordinates"][start:start+count]
            ca=f["actin/center"][frame]
            ua=unit(f["actin/direction"][frame])
            p0=ca[ident[:,2]]+sites[:,0,None]*ua[ident[:,2]]
            p1=ca[ident[:,3]]+sites[:,1,None]*ua[ident[:,3]]
            frames[frame]=dict(frame=frame,step=int(f["state/current_step"][frame,0]),
                ca=ca,ua=ua,mc=f["myosin/center"][frame],
                mu=unit(f["myosin/direction"][frame]),ident=ident,sites=sites,p0=p0,p1=p1,
                index={int(lid):k for k,lid in enumerate(ident[:,1])})
    return dict(attrs=attrs,box=box,frames=frames,birth=born,mat=mat,death=death)


def link_geometry(data, frame, lid, m):
    fr=data["frames"][frame]
    k=fr["index"].get(lid)
    if k is None:return None
    result=geometry(fr["p0"][k],fr["p1"][k],fr["mc"][m],fr["mu"][m],
                    data["box"],float(data["attrs"]["myosin_length"]))
    result.update(status=int(fr["ident"][k,4]),index=k,frame=frame)
    return result


def select_cases(data):
    first_clear={}
    entries=[]
    entered=set()
    lm=float(data["attrs"]["myosin_length"])
    for frame,fr in data["frames"].items():
        if frame==400:continue
        q=np.flatnonzero((fr["ident"][:,4]==2)&np.isin(fr["ident"][:,1],list(data["mat"])))
        hit,depth,rho,which,ax=blocked_connectors(
            fr["p0"][q],fr["p1"][q],fr["mc"],fr["mu"],data["box"],lm)
        for j,k in enumerate(q):
            lid=int(fr["ident"][k,1])
            if lid not in first_clear:
                first_clear[lid]=(frame,not bool(hit[j]))
            elif hit[j] and first_clear[lid][1] and lid not in entered:
                entered.add(lid)
                # Require both initial weak formation and maturation to occur
                # after the resume, to eliminate any inherited-cohort ambiguity.
                if lid in data["birth"]:
                    entries.append(dict(lid=lid,m=int(which[j]),entry=frame,first=first_clear[lid][0]))
        if frame%20==0:print("SELECTION_SCAN",frame,flush=True)
    candidates=[]
    for c in entries:
        lid,m,entry=c["lid"],c["m"],c["entry"]
        previous=link_geometry(data,entry-1,lid,m)
        if previous is None or previous["status"]!=2 or previous["blocked"]:continue
        before=[]
        for frame in range(entry-1,c["first"]-1,-1):
            g=link_geometry(data,frame,lid,m)
            if g is None or g["status"]!=2 or g["blocked"]:break
            before.append(frame)
        after=[]
        for frame in range(entry,min(entry+18,max(data["frames"]))+1):
            g=link_geometry(data,frame,lid,m)
            if g is None or g["status"]!=2:break
            after.append((frame,g))
        streak=0
        for _,g in after:
            if not g["blocked"]:break
            streak+=1
        bad=[(frame,g) for frame,g in after if g["blocked"]]
        if len(before)<2 or len(after)<5 or streak<2 or not bad:continue
        deepest=max(bad,key=lambda item:item[1]["depth"])
        if deepest[1]["depth"]<.030:continue
        kind="lateral" if previous["axial_overlap"] else "axial"
        # Favor a visible approach, persistent entry, and a sizeable depth.
        c.update(kind=kind,clear_frames=len(before),mature_frames_after=len(after),
                 blocked_streak=streak,deepest_frame=deepest[0],
                 max_depth_nm=1000*deepest[1]["depth"],
                 min_radius_nm=1000*min(g["rho"] for _,g in bad),
                 score=1000*deepest[1]["depth"]+12*streak+4*min(len(before),8))
        candidates.append(c)
    candidates.sort(key=lambda c:c["score"],reverse=True)
    chosen=[]
    for kind in ("axial","lateral"):
        available=[c for c in candidates if c["kind"]==kind and c["lid"] not in {p["lid"] for p in chosen}]
        fresh=[c for c in available if (500*c["entry"]-data["mat"][c["lid"]])*.01<=100]
        if fresh:available=fresh
        if available:chosen.append(available[0])
    if len(chosen)<2:
        chosen.extend(c for c in candidates if c["lid"] not in {p["lid"] for p in chosen})
        chosen=chosen[:2]
    assert len(chosen)==2,"Could not isolate two sufficiently sampled entries"
    for i,c in enumerate(chosen,1):
        c["case"]=i
        c["start"]=max(400,c["first"]-1,c["entry"]-10)
        c["end"]=min(max(data["frames"]),c["entry"]+14)
        # End the clip shortly after an observed detachment, if one occurs.
        d=data["death"].get(c["lid"])
        if d:c["end"]=min(c["end"],int(np.ceil(d["step"]/500))+1)
        c["mature_step"]=data["mat"][c["lid"]]
        c["birth_step"]=data["birth"][c["lid"]]
        c["death"]=d
        fr=data["frames"][c["entry"]];k=fr["index"][c["lid"]]
        c["actins"]=fr["ident"][k,2:4].astype(int).tolist()
        midpoint=fr["p0"][k]+.5*wrap(fr["p1"][k]-fr["p0"][k],data["box"])
        axial=float(wrap(midpoint-fr["mc"][c["m"]],data["box"])@fr["mu"][c["m"]])
        c["end_sign"]=1 if axial>=0 else -1
        print("SELECTED",json.dumps(c),flush=True)
    return chosen,candidates


def tube(ax, xlo, xhi, radius, color, alpha, outline=False):
    theta=np.linspace(0,2*np.pi,40)
    # A local crop through the shaft must not look like another physical end.
    for x in ((xlo,xhi) if xlo<=-1500 else (xhi,)):
        ax.plot(np.full_like(theta,x),radius*np.cos(theta),radius*np.sin(theta),
                color=color,lw=1.0,alpha=.55 if outline else .35)
    if outline:
        for theta0 in np.linspace(0,2*np.pi,8,endpoint=False):
            ax.plot([xlo,xhi],[radius*np.cos(theta0)]*2,[radius*np.sin(theta0)]*2,
                    color=color,lw=.6,alpha=.22)
    else:
        x,t=np.meshgrid(np.linspace(xlo,xhi,2),theta)
        ax.plot_surface(x,radius*np.cos(t),radius*np.sin(t),color=color,
                        alpha=alpha,shade=False,linewidth=0,antialiased=False)


def prepare_series(data,c):
    rows=[]
    old_e1=None
    first_basis_frame=data["frames"][c["entry"]]
    k=first_basis_frame["index"][c["lid"]]
    u0=c["end_sign"]*first_basis_frame["mu"][c["m"]]
    mid=first_basis_frame["p0"][k]+.5*wrap(first_basis_frame["p1"][k]-first_basis_frame["p0"][k],data["box"])
    initial=wrap(mid-first_basis_frame["mc"][c["m"]],data["box"])
    initial-=float(initial@u0)*u0
    initial/=np.linalg.norm(initial)
    for frame in range(c["start"],c["end"]+1):
        fr=data["frames"][frame]
        u=c["end_sign"]*fr["mu"][c["m"]]
        e1=initial.copy() if old_e1 is None else old_e1.copy()
        e1-=float(e1@u)*u;e1/=np.linalg.norm(e1)
        e2=np.cross(u,e1);old_e1=e1
        basis=np.stack([u,e1,e2],axis=1)
        row=dict(frame=frame,step=fr["step"],time_s=fr["step"]*1e-5,
                 age_ms=(fr["step"]-c["mature_step"])*.01,exists=False,status=0,
                 basis=basis)
        k=fr["index"].get(c["lid"])
        if k is None:
            # Track the actual actins even before this molecule is born or
            # after it detaches; no connector is drawn for an absent link.
            tip=.5*float(data["attrs"]["myosin_length"])*u
            segments=[]
            half=.5*float(data["attrs"]["actin_length"])
            for a in c["actins"]:
                center=wrap(fr["ca"][a]-fr["mc"][c["m"]],data["box"])
                endpoints=np.stack([center-half*fr["ua"][a],center+half*fr["ua"][a]])
                segments.append(1000*((endpoints-tip)@basis))
            row["actins"]=segments
            rows.append(row);continue
        d=wrap(fr["p1"][k]-fr["p0"][k],data["box"])
        mid=wrap(fr["p0"][k]+.5*d-fr["mc"][c["m"]],data["box"])
        p0,p1=mid-.5*d,mid+.5*d
        si,sj=fr["sites"][k]
        i,j=c["actins"]
        ci,cj=p0-si*fr["ua"][i],p1-sj*fr["ua"][j]
        tip=.5*float(data["attrs"]["myosin_length"])*u
        transform=lambda x:1000*((np.asarray(x)-tip)@basis)
        half=.5*float(data["attrs"]["actin_length"])
        ai=np.stack([ci-half*fr["ua"][i],ci+half*fr["ua"][i]])
        aj=np.stack([cj-half*fr["ua"][j],cj+half*fr["ua"][j]])
        g=link_geometry(data,frame,c["lid"],c["m"])
        row.update(exists=True,status=int(fr["ident"][k,4]),connector=transform([p0,p1]),
                   actins=[transform(ai),transform(aj)],geometry=g,
                   closest=transform(g["point"]) if g["point"] is not None else None)
        rows.append(row)
    return rows


def clip_segment(segment, bounds):
    p,d=segment[0],segment[1]-segment[0]
    lo,hi=0.,1.
    for k,(lower,upper) in enumerate(bounds):
        if d[k]==0:
            if not lower<=p[k]<=upper:return None
        else:
            a,b=(lower-p[k])/d[k],(upper-p[k])/d[k]
            lo=max(lo,min(a,b));hi=min(hi,max(a,b))
            if lo>hi:return None
    return np.stack([p+lo*d,p+hi*d])


def draw(data,c,series,index,fig):
    fig.clear()
    current=series[index]
    available=[r for r in series if r["exists"]]
    connector_pts=np.concatenate([r["connector"] for r in available])
    xlo=min(-230.,float(connector_pts[:,0].min())-90)
    xhi=max(170.,float(connector_pts[:,0].max())+90)
    xlo=max(-1500.,xlo)
    radial=max(200.,float(np.linalg.norm(connector_pts[:,1:],axis=1).max())+60)
    radius_display=min(350.,radial)
    gs=fig.add_gridspec(2,2,height_ratios=[1.8,1],width_ratios=[1.2,1],
                       left=.05,right=.97,bottom=.075,top=.86,hspace=.34,wspace=.23)
    ax3=fig.add_subplot(gs[0,0],projection="3d")
    ax3.set_proj_type("ortho")
    tube(ax3,xlo,0,50,"#74808b",.32)
    tube(ax3,xlo,0,130,"#6b96b8",.055)
    tube(ax3,xlo,0,130,"#6b96b8",0,outline=True)
    ax3.scatter(0,0,0,s=50,color="#39434d",marker="D")
    ax3.set_xlim(xlo,xhi);ax3.set_ylim(-radius_display,radius_display);ax3.set_zlim(-radius_display,radius_display)
    ax3.set_box_aspect([xhi-xlo,2*radius_display,2*radius_display])
    ax3.view_init(elev=23,azim=-62)
    ax3.set_xlabel("axial distance from tip (nm)",labelpad=8)
    ax3.set_ylabel("transverse y (nm)",labelpad=6)
    ax3.set_zlabel("transverse z (nm)",labelpad=5)
    ax3.set_title("3D local crop · camera follows myosin",loc="left",fontsize=11)
    ax3.tick_params(labelsize=8)
    ax=fig.add_subplot(gs[0,1])
    ax.fill_between([xlo,0],0,130,color="#dbeaf4",label="130 nm exclusion region")
    ax.fill_between([xlo,0],0,50,color="#9ca6b0",alpha=.65,label="50 nm myosin body")
    ax.axvline(0,color="#536575",lw=1,ls="--")
    ax.axhline(130,color="#6b96b8",lw=1,ls="--")
    ax.text(xlo+8,136,"130 nm",fontsize=9,color="#527895")
    ax.text(7,8,"tip plane",fontsize=9,color="#536575")
    ax.set_xlim(xlo,xhi);ax.set_ylim(0,radius_display)
    ax.set_aspect("equal",adjustable="box")
    ax.set_xlabel("axial distance from selected tip (nm)")
    ax.set_ylabel("distance from myosin axis, ρ (nm)")
    ax.set_title("Axial–radial geometry · full connector",loc="left",fontsize=11)
    g=current.get("geometry")
    blocked=bool(g and g["blocked"])
    link_color=BLOCKED if blocked else CLEAR
    if "actins" in current:
        for segment,color in zip(current["actins"],ACTIN_COLORS):
            bounds=[(xlo,xhi),(-radius_display,radius_display),(-radius_display,radius_display)]
            clipped=clip_segment(segment,bounds)
            if clipped is not None:ax3.plot(*clipped.T,color=color,lw=3.5)
            if all(lower<=segment[0,k]<=upper for k,(lower,upper) in enumerate(bounds)):
                ax3.scatter(*segment[0],color=color,s=32,marker="o")
            s=np.linspace(0,1,100)[:,None]
            points=segment[0]+s*(segment[1]-segment[0])
            ax.plot(points[:,0],np.linalg.norm(points[:,1:],axis=1),color=color,lw=2.4,alpha=.9)
            ax.scatter(segment[0,0],np.linalg.norm(segment[0,1:]),color=color,s=25)
    if current["exists"]:
        link=current["connector"]
        ax3.plot(*link.T,color=link_color,lw=5)
        ax3.scatter(*link.T,color=link_color,s=25)
        points=link[0]+np.linspace(0,1,80)[:,None]*(link[1]-link[0])
        ax.plot(points[:,0],np.linalg.norm(points[:,1:],axis=1),color=link_color,lw=4)
        if current["closest"] is not None:
            p=current["closest"]
            ax3.scatter(*p,color=link_color,s=55,marker="X")
            ax.scatter(p[0],np.linalg.norm(p[1:]),color=link_color,s=65,marker="X",zorder=10)
    handles=[Line2D([],[],color=color,lw=3,label=f"actin {a}") for color,a in zip(ACTIN_COLORS,c["actins"])]
    handles += [Line2D([],[],color=CLEAR,lw=3,label="ACTN2 clear"),
                Line2D([],[],color=BLOCKED,lw=3,label="ACTN2 excluded"),
                Line2D([],[],color="#74808b",lw=6,label=f"myosin {c['m']}")]
    fig.legend(handles=handles,loc="upper center",bbox_to_anchor=(.50,.945),ncol=5,frameon=False,fontsize=10)
    ages=np.array([r["age_ms"] for r in series])
    rho=np.array([1000*r["geometry"]["rho"] if r["exists"] and r["geometry"]["rho"] is not None else np.nan for r in series])
    clearance=np.array([1000*r["geometry"]["clearance"] if r["exists"] else np.nan for r in series])
    for col,values,threshold,label,color in (
        (0,rho,130,"minimum radius within axial overlap (nm)","#527895"),
        (1,clearance,0,"axial clearance of connector (nm)","#7655ae"),
    ):
        tr=fig.add_subplot(gs[1,col])
        if col==0:tr.axhspan(0,130,color="#fcebed",alpha=.6)
        tr.axhline(threshold,color="#ce3349",ls="--",lw=1,label=f"{threshold} nm boundary")
        tr.plot(ages,values,color=color,lw=1.5,alpha=.45,marker="o",ms=3)
        tr.plot(ages[:index+1],values[:index+1],color=color,lw=2.2,marker="o",ms=4)
        tr.axvline(0,color="#343a40",ls=":",lw=1)
        tr.axvline(current["age_ms"],color=link_color,lw=1.5)
        if np.isfinite(values[index]):tr.scatter(current["age_ms"],values[index],color=link_color,s=55,zorder=10)
        tr.set_xlim(ages.min()-2,ages.max()+2)
        tr.set_xlabel("time since 1→2 maturation event (ms)")
        tr.set_ylabel(label,fontsize=9)
        tr.grid(alpha=.17)
        if col==0:
            finite=values[np.isfinite(values)]
            if len(finite):tr.set_ylim(max(0,float(finite.min())-20),max(160,float(finite.max())+15))
            tr.set_title("Radius is defined only where axial overlap exists",fontsize=10,loc="left")
        else:
            tr.set_title("Negative clearance = axial overlap",fontsize=10,loc="left")
    state="LINK ABSENT" if not current["exists"] else (
        "EXCLUDED · STILL STATUS 2" if blocked and current["status"]==2 else
        "CLEAR · STATUS 2" if current["status"]==2 else "STATUS 1")
    rho_text="no axial overlap" if not g or g["rho"] is None else f"ρmin = {1000*g['rho']:.1f} nm"
    title_kind="axial entry" if c["kind"]=="axial" else "lateral entry"
    fig.suptitle(
        f"Case {c['case']} · {title_kind} after maturation · ACTN2 {c['lid']} / myosin {c['m']}\n"
        f"t = {current['time_s']:.3f} s   |   mature age {current['age_ms']:.1f} ms   |   {state}   |   {rho_text}",
        fontsize=15,fontweight="bold",y=.995,color=BLOCKED if blocked else "#273440")
    fig.text(.5,.017,"Original production trajectory r27zt0npg130 · poses every 5 ms · no interpolated motion · dots: barbed tips · shafts clipped to local view",
             ha="center",fontsize=9,color="#536575")


def serialize(row):
    if isinstance(row,dict):return {k:serialize(v) for k,v in row.items()}
    if isinstance(row,(list,tuple)):return [serialize(v) for v in row]
    if isinstance(row,np.ndarray):return row.tolist()
    if isinstance(row,np.generic):return row.item()
    return row


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):raise RuntimeError("Render in Slurm")
    out=args.output;out.mkdir(parents=True,exist_ok=True)
    source=Path(os.environ["SCRATCH"])/f"sarcomere/gaphunt_20261023/data/{ARM}_seed20260805.h5"
    data=load_trajectory(source)
    cases,candidates=select_cases(data)
    (out/"selection.json").write_text(json.dumps(dict(selected=cases,candidates=candidates),indent=2)+"\n")
    products=[]
    matplotlib.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"axes.spines.top":False,"axes.spines.right":False})
    for case in cases:
        series=prepare_series(data,case)
        case["poster_frame"]=max(
            (r for r in series if r["exists"] and r["geometry"]["blocked"]),
            key=lambda r:r["geometry"]["depth"])["frame"]
        stem=f"case{case['case']}_{case['kind']}_link{case['lid']}_myosin{case['m']}"
        movie=out/(stem+".mp4")
        fig=plt.figure(figsize=(16,9),dpi=100,facecolor="white")
        writer=FFMpegWriter(fps=12,codec="libx264",metadata={"title":stem,"comment":str(source)},
                           extra_args=["-pix_fmt","yuv420p","-crf","19","-movflags","+faststart","-threads","1"])
        with writer.saving(fig,str(movie),100):
            for k,row in enumerate(series):
                draw(data,case,series,k,fig)
                holds=3+(9 if k in (0,len(series)-1) else 0)+(6 if row["frame"]==case["entry"] else 0)
                for _ in range(holds):writer.grab_frame()
                if row["frame"]==case["poster_frame"]:
                    fig.savefig(out/(stem+".png"),dpi=130,facecolor="white")
                if k%5==0:print("RENDER",stem,k,len(series),flush=True)
        plt.close(fig)
        manifest=dict(
            input=str(source),source_size=source.stat().st_size,source_mtime_ns=source.stat().st_mtime_ns,
            arm=ARM,case=case,slurm_job_id=os.environ["SLURM_JOB_ID"],
            source_binary_sha256=data["attrs"].get("ab_binary_sha256",b"").decode(),
            render_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            description="Newly born and newly matured link first observed clear globally, then entering a fixed myosin's finite 130-nm cylinder while status 2.",
            display="Comoving myosin tip frame, continuous transverse basis, nm; actual 5-ms saved poses with repeated frames for slow playback.",
            sampling_caveat="Maturation time is from events; gate-time position is not stored. Motion between saved poses is not reconstructed.",
            series=serialize(series),video=str(movie),poster=str(out/(stem+".png")))
        (out/(stem+".json")).write_text(json.dumps(manifest,indent=2)+"\n")
        # Assert the selected trajectory shows at least two clear mature
        # poses followed by excluded mature poses for this same myosin.
        before=[r for r in series if r["exists"] and r["status"]==2 and r["frame"]<case["entry"]]
        after=[r for r in series if r["exists"] and r["status"]==2 and r["frame"]>=case["entry"]]
        assert sum(not r["geometry"]["blocked"] for r in before)>=2
        assert sum(r["geometry"]["blocked"] for r in after)>=2
        products.append(manifest)
        print("VIDEO_DONE",movie,flush=True)
    joined=out/"newly_mature_entries.mp4"
    # Concatenate the two equal-size/equal-codec streams without reencoding.
    concat=out/"concat.txt"
    concat.write_text("".join("file '"+Path(r["video"]).name+"'\n" for r in products))
    subprocess.run(["ffmpeg","-nostdin","-y","-v","error","-f","concat","-safe","0",
                    "-i",str(concat),"-c","copy","-movflags","+faststart",str(joined)],check=True)
    (out/"manifest.json").write_text(json.dumps(dict(input=str(source),videos=[r["video"] for r in products],
        combined=str(joined),cases=cases,sampling_ms=5),indent=2)+"\n")
    print("ALL_MOVIES_COMPLETE",out,flush=True)


if __name__=="__main__":main()
