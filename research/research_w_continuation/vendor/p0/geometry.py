"""Frozen total-material worlds; generation is separate from search input."""
import numpy as np

def antennas():
    tx=np.c_[np.linspace(-1.4,1.4,4),np.full(4,-.3)]
    rx=np.c_[np.linspace(-1.5,1.5,12),np.full(12,-.3)]
    return tx,rx

def in_query(points):
    x,y=np.asarray(points).T
    return (abs(x)<=1.2)&(y>=1.6)&(y<=2.8)

def wall_values(points,seed):
    """Shared epsilon_infinity and conductivity, not per-frequency fitting."""
    x,y=np.asarray(points).T
    rng=np.random.default_rng(seed);p=rng.uniform(0,2*np.pi,5)
    low=.64+.09*x+.055*np.sin(2.4*x+p[0])
    thick=.48+.085*np.cos(2.1*x+p[1])
    support=(abs(x)<1.48+.08*np.sin(4*y+p[2]))&(y>low)&(y<low+thick)
    angle=np.arctan2(y-.92,x+.20);rad=np.hypot((x+.20)/.25,(y-.92)/.15)
    cavity=rad<1+.18*np.sin(3*angle+p[3])
    occupied=support&~cavity
    er=1+occupied*(3.1+.65*np.cos(2.3*x+p[4])+.25*np.sin(5*y))
    conductivity=occupied*(.018+.009*(1+np.sin(3*x+2*y)))
    return er,conductivity

def material_grid(backend,seed,frequency=None,subsamples=4,probe=None):
    # Deterministic area averaging of one common continuous geometry. This is
    # a coefficient discretization, not a continuous-PDE error certificate.
    f=backend.freq if frequency is None else frequency
    er=np.zeros(backend.N);sigma=np.zeros(backend.N)
    offsets=(np.arange(subsamples)+.5)/subsamples-.5
    for a in offsets:
        for b in offsets:
            pts=backend.points+backend.h*np.array([a,b])
            re,loss=wall_values(pts,seed)
            re=np.where(in_query(pts),1.,re);loss=np.where(in_query(pts),0.,loss)
            if probe is not None:
                center,radius,contrast=probe
                target=np.linalg.norm(pts-np.asarray(center),axis=1)<=radius
                re=np.where(target,1+contrast,re);loss=np.where(target,.02,loss)
            er+=re;sigma+=loss
    return er/subsamples**2+1j*sigma/(subsamples**2*f)

def probes():
    return [(tuple((x,y)),.12,c) for y in [1.85,2.15,2.45]
            for x in [-.9,-.3,.3,.9] for c in [3.,7.]]

def target_centroid(points,epsilon):
    mask=in_query(points);weights=np.maximum(np.asarray(epsilon).real-1,0)*mask
    mass=weights.sum()
    return np.full(2,np.nan) if mass==0 else weights@np.asarray(points)/mass
