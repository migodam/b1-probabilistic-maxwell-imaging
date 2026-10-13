"""One continuous coefficient prior; topology is derived after decoding.

The two compact bumps have disjoint supports. Inactive coefficients integrate
to genuine material atoms, not to extra topology experts. No thin layer is
removed: only exactly coincident interfaces are merged.
"""
from __future__ import annotations
import numpy as np

BOUNDS = np.array([[-3., 1.], [-3., 1.], [2.2, 2.8], [.04, .16]])
WIDTH = np.diff(BOUNDS, axis=1).ravel()
SCALES = np.array([.45, .18])
RMAX = SCALES*np.sqrt(1-1/np.sqrt(3.))
EVENT_LABELS = [f"central-{c}_annular-{s}" for c in range(4) for s in range(3)]

def prior_sample(rng, n):
    return rng.uniform(BOUNDS[:, 0], BOUNDS[:, 1], size=(n, 4))

def prior_logpdf(x):
    x=np.asarray(x); inside=((x>=BOUNDS[:,0]) & (x<=BOUNDS[:,1])).all(-1)
    return np.where(inside, -np.log(WIDTH).sum(), -np.inf)

def geometry(x):
    z=np.asarray(x)[..., :2]
    active=z < -1.
    a=np.where(active, -z, 1.)
    return SCALES*np.sqrt(np.maximum(1-1/np.sqrt(a), 0.))

def from_geometry(r, h, er=2.5, ei=.1):
    g=np.array([r,h],float)
    if np.any(g<0) or np.any(g>RMAX+1e-14): raise ValueError('geometry outside declared support')
    z=np.where(g>0, -(1-(g/SCALES)**2)**-2, 0.)
    # Exact physical endpoints should encode exact coefficient endpoints; the
    # inverse formula can round -3 to -3.0000000000000004. Not a sample clip.
    z=np.where(g==RMAX,-3.,z)
    return np.r_[z,er,ei]

def event_ids(x):
    g=geometry(x)
    c=np.where(g[...,0]==0,0,np.where(g[...,0]<=.1,1,np.where(g[...,0]<=.2,2,3)))
    s=np.where(g[...,1]==0,0,np.where(g[...,1]<=.06,1,2))
    return c*3+s

def prior_event_masses():
    def mass(t,a): return ((1-(t/a)**2)**-2-1)/4
    a=mass(.1,.45); b=mass(.2,.45); d=mass(.06,.18)
    return np.outer([.5,a,b-a,.5-b],[.5,d,.5-d]).ravel()

def event_intervals():
    """Coefficient intervals in each physical event; None denotes inactive atom.

    Inactive coefficients remain uniform on [-1,1] conditional on material.
    Their probability .5 is retained explicitly in all integrals.
    """
    def z(t,a): return -(1-(t/a)**2)**-2
    return ([None,(z(.1,.45),-1.),(z(.2,.45),z(.1,.45)),(-3.,z(.2,.45))],
            [None,(z(.06,.18),-1.),(-3.,z(.06,.18))])

def void_volume(x):
    g=geometry(x); r,h=g[...,0],g[...,1]
    return r**3+np.where(h>0,(.65+h)**3-(.65-h)**3,0.)

def decode(x):
    x=np.asarray(x,float)
    if x.shape!=(4,) or not np.isfinite(prior_logpdf(x)): raise ValueError('invalid material latent')
    r,h=geometry(x); voids=[]
    if r>0: voids.append((0.,float(r)))
    if h>0: voids.append((float(.65-h),float(.65+h)))
    edges=sorted({v for interval in voids for v in interval if v>0}|{1.})
    eps=[]; lo=0.
    for hi in edges:
        mid=(lo+hi)/2
        eps.append(1.+0j if any(a<=mid<b for a,b in voids) else complex(x[2],x[3]))
        lo=hi
    return dict(radii=np.array(edges),eps_layers=np.array(eps),void_intervals=voids,
                central_radius=float(r),annular_halfwidth=float(h),event_id=int(event_ids(x)))

def material_distance(x,y):
    a,b=decode(x),decode(y); edges=sorted(set(a['radii'])|set(b['radii'])|{0.})
    value=0.
    for lo,hi in zip(edges[:-1],edges[1:]):
        mid=(lo+hi)/2
        ea=a['eps_layers'][np.searchsorted(a['radii'],mid)]
        eb=b['eps_layers'][np.searchsorted(b['radii'],mid)]
        value+=(hi**3-lo**3)*abs(ea-eb)**2
    return float(np.sqrt(value)/abs(1.5+.1j))

def image(x, n=151):
    d=decode(x); axis=np.linspace(-1.05,1.05,n); X,Y=np.meshgrid(axis,axis); r=np.hypot(X,Y)
    eps=np.ones(r.shape,complex); inside=r<=1.
    eps[inside]=complex(x[2],x[3])
    for a,b in d['void_intervals']: eps[inside & (r>=a) & (r<b)]=1.
    return axis,eps
