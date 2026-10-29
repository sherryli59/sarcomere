"""Plot the proposed soft energy, not a calibrated/implemented force law."""
import json
import os
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle

K=1000.0                 # pN/um, illustrative
HALF=.750; APPROACH=.050 # um
FLAT=.100; OUTER=.130
ELL=.040                # illustrative radial connector length, um


def weight(r):
    q=np.clip((np.abs(r)-FLAT)/(OUTER-FLAT),0,1)
    return (1-q)**2*(1+2*q)


def weight_prime(r):
    a=np.abs(r);q=np.clip((a-FLAT)/(OUTER-FLAT),0,1)
    return np.where((a>FLAT)&(a<OUTER),6*q*(q-1)/(OUTER-FLAT)*np.sign(r),0)


def primitive(r):
    """Integral from 0 to r of W(abs(s)), including signed r."""
    a=np.abs(r);q=np.clip((a-FLAT)/(OUTER-FLAT),0,1)
    f=np.where(a<=FLAT,a,FLAT+(OUTER-FLAT)*(q-q**3+.5*q**4))
    return np.sign(r)*f


def averaged_weight(r):
    return np.clip((primitive(r+.5*ELL)-primitive(r-.5*ELL))/ELL,0,1)


def averaged_prime(r):
    return (weight(r+.5*ELL)-weight(r-.5*ELL))/ELL


def energy(g,r,segment=False):
    delta=np.maximum(APPROACH-g,0)
    return .5*K*delta**2*(averaged_weight(r) if segment else weight(r))


def force(g,r,segment=False):
    delta=np.maximum(APPROACH-g,0)
    w=averaged_weight(r) if segment else weight(r)
    wp=averaged_prime(r) if segment else weight_prime(r)
    return K*delta*w,-.5*K*delta**2*wp


def validate():
    assert np.isclose(energy(0,.08),1.25)
    assert np.isclose(energy(0,.115),.625)
    assert energy(.05,.08)==0
    assert energy(0,.13)==0
    assert energy(0,.14,True)>0
    assert np.isclose(1000*energy(0,.14,True),28.935185185,rtol=1e-8)
    assert np.allclose(force(0,.115),[25,62.5])
    for segment in (False,True):
        for g,r in ((-.1,.115),(0,.12),(.02,.105),(0,.14)):
            e=1e-7
            fz,fr=force(g,r,segment)
            numz=-(energy(g+e,r,segment)-energy(g-e,r,segment))/(2*e)
            numr=-(energy(g,r+e,segment)-energy(g,r-e,segment))/(2*e)
            assert np.allclose([fz,fr],[numz,numr],rtol=1e-6,atol=1e-5)


def contours(ax,x,y,e,maximum):
    levels=np.geomspace(1,maximum,25)
    fill=ax.contourf(x,y,np.ma.masked_less(e,1),levels=levels,
        cmap="magma",norm=LogNorm(vmin=1,vmax=maximum),extend="min")
    selected=[v for v in (10,100,1000,10000,100000) if v<maximum]
    lines=ax.contour(x,y,e,levels=selected,colors="#eeeeee",linewidths=.65,alpha=.8)
    ax.clabel(lines,fmt=lambda v:f"{v:g}",fontsize=7,inline=True)
    ax.set_facecolor("#f5f6f8")
    return fill


def zoom(ax,segment=False):
    g=np.linspace(-150,100,601);r=np.linspace(0,170,501)
    gx,ry=np.meshgrid(g,r)
    u=1000*energy(gx*.001,ry*.001,segment)
    fill=contours(ax,g,r,u,20000)
    ax.axvline(0,color="#29b6d1",lw=1.5)
    ax.axvline(50,color="#218b8b",lw=1.2,ls="--")
    ax.axhline(100,color="#777777",lw=.9,ls=":")
    ax.axhline(130,color="#29b6d1",lw=1.2,ls="--")
    if segment:ax.axhline(150,color="#52a37c",lw=.9,ls=":")
    ax.add_patch(Rectangle((-150,0),150,50,fill=False,edgecolor="#7d8794",lw=.8,hatch="///"))
    ax.text(-143,23,"myosin core",color="#a3afbc",fontsize=8)
    xx,yy=np.meshgrid(np.arange(-135,45,25),np.arange(15,156,20))
    fz,fr=force(xx*.001,yy*.001,segment)
    norm=np.hypot(fz,fr)
    valid=(1000*energy(xx*.001,yy*.001,segment)>15)&~((xx<0)&(yy<50))&(norm>1e-10)
    ax.quiver(xx[valid],yy[valid],10*fz[valid]/norm[valid],10*fr[valid]/norm[valid],
        angles="xy",scale_units="xy",scale=1,color="white",width=.004,
        headwidth=3.2,headlength=4,alpha=.85)
    ax.text(2,164,"tip plane",color="#087f97",fontsize=8,rotation=90,va="top")
    ax.text(52,164,"force onset",color="#127d7d",fontsize=8,rotation=90,va="top")
    ax.set(xlim=(-150,100),ylim=(0,170),xlabel="axial gap g from plus tip (nm)",
        ylabel="connector midpoint radius (nm)" if segment else "point radius ρ (nm)")
    ax.set_aspect("equal",adjustable="box")
    return fill


def main():
    assert os.environ.get("SLURM_JOB_ID"),"Render in Slurm"
    validate()
    out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,
                         "axes.titlesize":11,"axes.titleweight":"bold"})
    fig,axes=plt.subplots(2,2,figsize=(13.2,9.8),layout="constrained")
    fig.suptitle("Proposed A-band contact energy",fontsize=19,fontweight="bold",y=.992)
    fig.get_layout_engine().set(hspace=.10,wspace=.08,rect=(0,.075,1,.825))
    fig.text(.5,.950,"Illustrative k = 1000 pN/µm  ·  myosin half-length = 750 nm  ·  approach = 50 nm  ·  taper = 100–130 nm",
             ha="center",fontsize=10,color="#475569")
    ax=axes[0,0]
    z=np.linspace(-900,900,1001);r=np.linspace(0,170,401)
    zz,rr=np.meshgrid(z,r)
    delta=np.maximum(HALF+APPROACH-np.abs(zz*.001),0)
    u=1000*.5*K*delta**2*weight(rr*.001)
    full=contours(ax,z,r,u,320000)
    for x in (-750,750):ax.axvline(x,color="#29b6d1",lw=1.2)
    for x in (-800,800):ax.axvline(x,color="#218b8b",lw=.9,ls="--")
    ax.axhline(100,color="#777777",lw=.9,ls=":")
    ax.axhline(130,color="#29b6d1",lw=1.2,ls="--")
    ax.add_patch(Rectangle((-750,0),1500,50,fill=False,edgecolor="#9ca3af",lw=.8,hatch="///"))
    ax.set(xlabel="axial position z from myosin center (nm)",ylabel="point radius ρ (nm)",
           title="A  Point-energy density around the full myosin",xlim=(-900,900),ylim=(0,170))
    cb=fig.colorbar(full,ax=ax,pad=.02,fraction=.045,ticks=[1,100,10000,100000])
    cb.set_label("energy (pN nm; logarithmic color)")
    b=zoom(axes[0,1])
    axes[0,1].set_title("B  Near the plus tip: point-energy density")
    c=zoom(axes[1,0],True)
    axes[1,0].set_title("C  Full-connector energy: 40 nm radial segment")
    axes[1,0].plot([0],[140],"o",color="#52a37c",mec="white",ms=6,zorder=10)
    axes[1,0].annotate("midpoint 140 nm\nnear endpoint 120 nm",(0,140),(-138,151),
        fontsize=8,color="#17623e",arrowprops=dict(arrowstyle="->",color="#17623e"))
    for fill,ax in ((b,axes[0,1]),(c,axes[1,0])):
        cb=fig.colorbar(fill,ax=ax,pad=.02,fraction=.045,ticks=[1,10,100,1000,10000])
        cb.set_label("energy (pN nm; logarithmic color)")
    ax=axes[1,1]
    radius=np.linspace(0,170,1001)
    ax.plot(radius,1000*energy(0,radius*.001),lw=2.4,label="point-energy density")
    ax.plot(radius,1000*energy(0,radius*.001,True),lw=2.4,color="#269063",label="40 nm full connector")
    for x in (100,130,150):ax.axvline(x,color="#9ca3af",lw=.8,ls=":" if x!=130 else "--")
    ax.scatter([140],[1000*energy(0,.140,True)],color="#269063",zorder=5)
    ax.annotate("Full connector at midpoint 140 nm:\nU = 28.9 pN nm; point U = 0",(140,28.935),(48,470),
        fontsize=9,arrowprops=dict(arrowstyle="->",color="#269063"),
        bbox=dict(boxstyle="round,pad=.4",fc="white",ec="#cbd5e1",alpha=.94))
    ax.set(title="D  Energy at the physical tip plane (g = 0)",xlabel="point / connector midpoint radius (nm)",
           ylabel="energy (pN nm)",xlim=(0,170),ylim=(-35,1370))
    ax.grid(alpha=.18);ax.legend(loc="upper right",frameon=False,fontsize=9)
    fig.text(.5,.027,"A/B show u = ½ k W(ρ) [h + a − |z|]₊²; C integrates u over a 40 nm radial connector (normalized by length).\n"
             "White arrows: normalized −∇U direction, not force magnitude; U < 1 pN nm is blank. Separate core steric energy is omitted. Proposal, not trajectory data.",
             ha="center",va="center",fontsize=9,color="#475569")
    for suffix in ("png","pdf","svg"):
        fig.savefig(out/("proposed_energy_contours."+suffix),dpi=210,facecolor="white")
    plt.close(fig)
    metadata=dict(stiffness_pN_per_um=K,half_length_nm=1000*HALF,approach_nm=1000*APPROACH,
        radial_flat_nm=1000*FLAT,radial_cutoff_nm=1000*OUTER,example_connector_length_nm=1000*ELL,
        soft_contact_only=True,calibrated=False,implemented=False,
        arrow_normalization="direction only; equal arrow lengths",
        examples=dict(point_tip_r115_energy_pN_nm=1000*energy(0,.115),
                      point_tip_r115_axial_radial_force_pN=list(force(0,.115)),
                      radial_connector_tip_midpoint140_energy_pN_nm=1000*energy(0,.14,True)),
        validation="Analytic W integral and finite-difference force checks passed",
        job=os.environ["SLURM_JOB_ID"])
    (out/"proposed_energy_contours.json").write_text(json.dumps(metadata,indent=2)+"\n")
    print("ENERGY_CONTOURS_DONE",out,flush=True)


if __name__=="__main__":main()
