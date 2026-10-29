"""Standalone 3D visualization of the proposed bounded connector energy.

Run numerical work in Slurm. Coordinates in um internally, display in nm.
No simulator or trajectory is modified. Vendor Plotly 3.1.0 separately.
"""
import argparse
import base64
import gzip
import json
from pathlib import Path
import platform
import numpy as np
from numpy.polynomial.legendre import leggauss
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

K, FC, EPS, R, H, HALF = 1000., 30., .005, .130, .800, .750
LEVELS = [10., 100., 500., 2000.]
COLORS = ['#37b9d6', '#41b899', '#e6ac45', '#df658c']
MODES = [
    ('point', 'Point probe', np.array([0., 0., 0.])),
    ('radial', '40 nm connector along Y', np.array([0., .040, 0.])),
    ('axial', '40 nm connector along X', np.array([.040, 0., 0.])),
    ('diagonal', '40 nm connector diagonal X–Y', np.array([.040, .040, 0.])/np.sqrt(2)),
]
EXAMPLES = [
    ('A', 'Approaching the + tip', [780., 40., 0.]),
    ('B', 'At the physical + tip', [750., 60., 0.]),
    ('C', 'Side/end transition', [750., 80., 0.]),
    ('D', 'Near the cylinder side', [500., 110., 20.]),
    ('E', 'Deep axial overlap', [0., 80., 0.]),
    ('F', 'Beyond the approach plane', [810., 80., 0.]),
    ('G', 'Midpoint outside 130 nm', [700., 140., 0.]),
    ('H', 'Approaching the − tip', [-780., 40., 30.]),
]


def point(x, force=False):
    x = np.asarray(x)
    rho = np.sqrt(x[..., 1]**2+x[..., 2]**2+EPS**2)
    za = np.sqrt(x[..., 0]**2+EPS**2)
    aa, bb = R-rho, H-za
    dd = np.sqrt((aa-bb)**2+EPS**2)
    p = np.maximum(.5*(aa+bb-dd), 0.)
    ps = FC/K
    e = np.where(p <= ps, .5*K*p*p, FC*(p-.5*ps))
    if not force:
        return e
    wa = .5*(1-(aa-bb)/dd)
    wb = 1-wa
    mag = np.minimum(K*p, FC)
    f = np.stack([mag*wb*x[..., 0]/za,
                  mag*wa*x[..., 1]/rho,
                  mag*wa*x[..., 2]/rho], axis=-1)
    return e, f


def connector(x, d, n=64, force=False):
    if not np.any(d):
        if force:
            e, f = point(x, True)
            return e, f, .5*f, .5*f
        return point(x)
    nodes, weights = leggauss(n)
    ts, weights = .5*(nodes+1), .5*weights
    e = np.zeros(np.asarray(x).shape[:-1])
    if force:
        fi, fj = np.zeros_like(x, dtype=float), np.zeros_like(x, dtype=float)
    # Iterate nodes, not points: bounded memory for the field mesh.
    for t, w in zip(ts, weights):
        xx = x+(t-.5)*d
        if force:
            ee, ff = point(xx, True)
            fi += w*(1-t)*ff
            fj += w*t*ff
        else:
            ee = point(xx)
        e += w*ee
    return (e, fi+fj, fi, fj) if force else e


def level_mesh(x, y, z, values, level):
    """Marching tetrahedra; mesh only cubes crossing an iso-energy level."""
    shape = values.shape
    idx = np.arange(values.size).reshape(shape)
    corners = np.stack([idx[:-1,:-1,:-1], idx[1:,:-1,:-1],
        idx[1:,1:,:-1], idx[:-1,1:,:-1], idx[:-1,:-1,1:],
        idx[1:,:-1,1:], idx[1:,1:,1:], idx[:-1,1:,1:]], axis=-1).reshape(-1,8)
    vals = values.ravel()
    cv = vals[corners]
    corners = corners[(cv.min(axis=1) < level)&(cv.max(axis=1) >= level)]
    coords = np.stack(np.meshgrid(x,y,z,indexing='ij'),axis=-1).reshape(-1,3)
    triangles = []
    for tet in ((0,1,2,6),(0,2,3,6),(0,3,7,6),(0,7,4,6),(0,4,5,6),(0,5,1,6)):
        ids = corners[:,tet]
        tv = vals[ids]
        code = ((tv >= level)*np.array([1,2,4,8])).sum(axis=1)
        for case in range(1,15):
            group = ids[code == case]
            if len(group) == 0:
                continue
            inside = [j for j in range(4) if case & (1<<j)]
            outside = [j for j in range(4) if not case & (1<<j)]
            def edge(a,b):
                ia, ib = group[:,a], group[:,b]
                q = ((level-vals[ia])/(vals[ib]-vals[ia]))[:,None]
                return coords[ia]+q*(coords[ib]-coords[ia])
            if len(inside) in (1,3):
                single, other = (inside[0],outside) if len(inside)==1 else (outside[0],inside)
                triangles.append(np.stack([edge(single,j) for j in other],axis=1))
            else:
                a,b = inside; c,d = outside
                ac,ad,bc,bd = edge(a,c),edge(a,d),edge(b,c),edge(b,d)
                triangles.extend([np.stack([ac,ad,bc],axis=1),np.stack([ad,bd,bc],axis=1)])
    tri = np.concatenate(triangles)
    verts, inverse = np.unique(np.round(1000*tri.reshape(-1,3),3),axis=0,return_inverse=True)
    faces = inverse.reshape(-1,3)
    return dict(vertices=verts.tolist(), faces=faces.tolist())


def validate():
    rng = np.random.default_rng(20261028)
    points = rng.uniform([-.86,-.17,-.17],[.86,.17,.17],(400,3))
    points = np.concatenate([points,np.array([p[2] for p in EXAMPLES])*.001])
    summary = {}
    for key, _, d in MODES:
        e, f, fi, fj = connector(points,d,64,True)
        numerical = np.zeros_like(f)
        step = 1.e-7
        for axis in range(3):
            dx = np.zeros(3); dx[axis] = step
            numerical[:,axis] = -(connector(points+dx,d)-connector(points-dx,d))/(2*step)
        err = float(np.max(np.abs(f-numerical)))
        assert err < 3.e-4, (key,err)
        assert np.all(np.linalg.norm(f,axis=1) <= FC+1.e-9)
        assert np.allclose(fi+fj,f)
        e128, f128, _, _ = connector(points,d,128,True)
        eq = float(np.max(np.abs(e-e128))*1000)
        fq = float(np.max(np.linalg.norm(f-f128,axis=1)))
        assert eq < .12 and fq < .12, (key,eq,fq)
        # Independent endpoint-coordinate finite differences.
        base = np.array([.750,.080,.015])
        if np.any(d):
            _, _, fa, fb = connector(base,d,64,True)
            for axis in range(3):
                dx=np.zeros(3);dx[axis]=step
                da = -(connector(base+dx/2,d-dx)-connector(base-dx/2,d+dx))/(2*step)
                db = -(connector(base+dx/2,d+dx)-connector(base-dx/2,d-dx))/(2*step)
                assert abs(da-fa[axis]) < 3.e-4
                assert abs(db-fb[axis]) < 3.e-4
        summary[key] = dict(max_gradient_error_pN=err,
            max_64_vs_128_energy_error_pNnm=eq,max_64_vs_128_force_error_pN=fq)
    assert point([.810,.080,0]) == 0
    assert connector(np.array([.700,.140,0]),MODES[1][2]) > 0
    return summary


def preview(path):
    fig=plt.figure(figsize=(14,6.5))
    ax=fig.add_axes([.01,.12,.49,.69],projection='3d')
    th=np.linspace(0,2*np.pi,65)
    xx,tt=np.meshgrid(np.linspace(-750,750,51),th)
    ax.plot_surface(xx,50*np.cos(tt),50*np.sin(tt),color='#8993a2',alpha=.45,linewidth=0)
    # Analytic point-energy isosurfaces of the smoothed minimum.
    for level,color in zip(LEVELS,COLORS):
        ee=level/1000
        p=np.sqrt(2*ee/K) if ee<=.5*FC*FC/K else ee/FC+.5*FC/K
        bmin=p+EPS**2/(4*(R-EPS-p))
        xmax=np.sqrt((H-bmin)**2-EPS**2)
        xs=xmax*np.cos(np.linspace(np.pi,0,125))
        b=H-np.sqrt(xs**2+EPS**2)
        a=p+EPS**2/(4*(b-p))
        rr=np.sqrt(np.maximum((R-a)**2-EPS**2,0))
        ax.plot_surface(1000*np.tile(xs,(len(th),1)),1000*np.outer(np.cos(th),rr),
            1000*np.outer(np.sin(th),rr),color=color,alpha=.15,linewidth=0)
    pts=np.array([e[2] for e in EXAMPLES]);es,fs=point(pts*.001,True)
    ax.scatter(*pts.T,c='#242b3a',s=17,depthshade=False)
    for (label,_,_),p,f in zip(EXAMPLES,pts,fs):
        ax.quiver(*p,*(4*f),color='#182334',arrow_length_ratio=.18,linewidth=1.5)
        ax.text(*p,label,fontsize=9)
    ax.set(xlabel='X: myosin axis (nm)',
        title='3D energy shells + force arrows\nArrow length: 4 nm per pN')
    ax.set_xticks([-750,0,750]);ax.set_yticks([-100,0,100]);ax.set_zticks([-100,0,100])
    ax.tick_params(labelsize=8,pad=0)
    ax.text2D(.04,.08,'Y / Z: transverse coordinates (nm)',transform=ax.transAxes,fontsize=9,color='#637389')
    ax.set_box_aspect([3.2,1,1]);ax.view_init(elev=24,azim=-65)
    bx=fig.add_axes([.55,.16,.36,.65])
    xs=np.linspace(550,850,401);ys=np.linspace(0,170,281)
    gx,gy=np.meshgrid(xs,ys);xyz=np.stack([gx,gy,np.zeros_like(gx)],axis=-1)*.001
    u,f=point(xyz,True);u*=1000
    fill=bx.contourf(xs,ys,np.ma.masked_less(u,1),levels=np.geomspace(1,3400,35),
        cmap='magma',norm=LogNorm(1,3400))
    cs=bx.contour(xs,ys,u,levels=LEVELS,colors=COLORS,linewidths=1)
    bx.clabel(cs,fmt=lambda v:f'{v:g} pN nm',fontsize=8)
    bx.axvline(750,color='#999',ls=':');bx.axvline(800,color='#258f9c',ls='--')
    bx.axhline(130,color='#258f9c',ls='--')
    for (label,_,_),p,ff in zip(EXAMPLES,pts,fs):
        if 550<=p[0]<=850 and p[2]==0:
            bx.scatter(p[0],p[1],s=22,c='white',edgecolor='#222',zorder=5)
            bx.quiver(p[0],p[1],4*ff[0],4*ff[1],angles='xy',scale_units='xy',scale=1,
                color='white',width=.006)
            bx.annotate(f'{label}: {np.linalg.norm(ff):.1f} pN',(p[0],p[1]),
                xytext=(4,7),textcoords='offset points',color='#197789',fontsize=8)
    bx.set(xlabel='X (nm)',ylabel='Y at Z = 0 (nm)',title='Point energy near the + tip',
        xlim=(550,850),ylim=(0,170),aspect='equal')
    cax=fig.add_axes([.93,.16,.014,.65])
    cb=fig.colorbar(fill,cax=cax,label='Energy (pN nm)',ticks=[1,10,100,1000])
    cb.ax.set_yticklabels(['1','10','100','1000'])
    fig.suptitle('Revised bounded-energy proposal · illustrative, not calibrated\n'
        'k = 1000 pN/µm · F* = 30 pN · smoothing = 5 nm · myosin radius = 50 nm',fontsize=12)
    fig.savefig(path,dpi=150)
    plt.close(fig)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--plotly',type=Path,required=True);args=ap.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    validation=validate();print(json.dumps(validation),flush=True)
    x=np.unique(np.r_[np.linspace(-.88,.88,111),np.linspace(-.83,-.65,37),np.linspace(.65,.83,37)])
    y=z=np.linspace(-.180,.180,49)
    xyz=np.stack(np.meshgrid(x,y,z,indexing='ij'),axis=-1)
    payload=dict(parameters=dict(k=K,forceBound=FC,epsilon=EPS,radius=R,
        halfLength=HALF,approach=H-HALF,kBT_pNnm=1.380649e-2*310),
        levels=LEVELS,colors=COLORS,examples=[dict(id=i,label=s,position=p) for i,s,p in EXAMPLES],modes=[])
    sx=np.linspace(-880,880,221);sy=np.linspace(-180,180,91)
    xx,yy=np.meshgrid(sx,sy);sxyz=np.stack([xx,yy,np.zeros_like(xx)],axis=-1)*.001
    for key,label,d in MODES:
        print(f'mesh {key}: {xyz.shape}',flush=True)
        values=connector(xyz,d,n=32)*1000
        meshes=[level_mesh(x,y,z,values,lev) for lev in LEVELS]
        print(f'{key} faces: {[len(m["faces"]) for m in meshes]}',flush=True)
        energy, force, fi, fj=connector(np.array([e[2] for e in EXAMPLES])*.001,d,n=64,force=True)
        payload['modes'].append(dict(key=key,label=label,displacement=(1000*d).tolist(),meshes=meshes,
            slice=dict(x=sx.tolist(),y=sy.tolist(),energy=np.round(connector(sxyz,d,32)*1000,4).tolist()),
            referenceEnergy=np.round(energy*1000,6).tolist(),referenceForce=force.tolist()))
    nodes,weights=leggauss(64)
    payload['quadrature']=dict(t=(.5*(nodes+1)).tolist(),w=(.5*weights).tolist())
    template=Path(__file__).with_name('interactive_energy_template.html').read_text()
    encoded=base64.b64encode(gzip.compress(json.dumps(payload,separators=(',',':')).encode(),compresslevel=6,mtime=0)).decode()
    html=template.replace('__PLOTLY_BUNDLE__',args.plotly.read_text()).replace('__ENERGY_DATA__',encoded)
    (args.output/'connector_energy_3d.html').write_text(html)
    # Compact payload needed for independent browser/Node checks, not all meshes.
    check={k:v for k,v in payload.items() if k!='modes'}
    check['modes']=[{k:v for k,v in m.items() if k not in ('meshes','slice')} for m in payload['modes']]
    (args.output/'connector_energy_examples.json').write_text(json.dumps(check,indent=2))
    (args.output/'validation.json').write_text(json.dumps(dict(python=platform.python_version(),numpy=np.__version__,checks=validation,
        implemented=False,calibrated=False,field_quadrature=32,probe_quadrature=64,
        soft_contact_only=True,mesh_note='Nonuniform X grid, 7.5 nm transverse spacing; interpolated iso-surfaces.'),indent=2))
    preview(args.output/'connector_energy_3d_preview.png')
    print(f'written {args.output}',flush=True)


if __name__=='__main__':
    main()
