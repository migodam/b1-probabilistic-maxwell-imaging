"""Preregistered synthetic labels; import/use only in offline evaluator."""
from __future__ import annotations
import numpy as np

def backgrounds(points,config):
    weak=complex(*config['background_weak'])
    out=np.full(len(points),weak,complex)
    radius=np.linalg.norm(points,axis=1)
    ring=out.copy()
    ring[(radius>=config['ring_radius_inner'])&(radius<=config['ring_radius_outer'])]=complex(*config['background_ring'])
    return {'weak':out,'ring':ring}

def rectangle(x,y,limits):
    lo,hi,bottom,top=limits
    return (x>=lo)&(x<=hi)&(y>=bottom)&(y<=top)

def material_pattern(points,parent,config):
    """Coordinates/masks are frozen before data generation, no image-dependent choice."""
    x,y=np.asarray(points)[:,:2].T
    c=config['parent_masks'][str(parent)]
    if parent==2401:
        cx,cy=c['elliptical_center']; r=np.hypot(x-cx,(y-cy)/c['y_axis_ratio'])
        p=((r>=c['ring'][0])&(r<=c['ring'][1])).astype(float)
        p+=c['inner_disk_weight']*(r<=c['inner_disk_radius'])
    elif parent==2402:
        cx,cy=c['elliptical_center'];r=np.hypot((x-cx)/c['x_axis_ratio'],y-cy)
        p=((r>=c['ring1'][0])&(r<=c['ring1'][1])).astype(float)
        p+=c['ring2_weight']*((r>=c['ring2'][0])&(r<=c['ring2'][1]))
    else:
        p=np.zeros(len(x),float)
        for lim in c.get('rectangles',[c.get('rectangle')]):
            if lim is not None: p=np.maximum(p,rectangle(x,y,lim))
        cx,cy=c['disk_center']
        p=np.maximum(p,c['disk_weight']*(np.hypot(x-cx,y-cy)<=c['disk_radius']))
    phase=complex(*config['perturbation_complex_phase'])
    return p*phase/abs(phase)

def noisy_observation(clean,yb,noise_db,rng):
    if noise_db is None: return clean.copy(),0.,np.zeros_like(clean)
    sigma=float(np.sqrt(np.mean(abs(clean-yb)**2))*10**(-noise_db/20))
    noise=sigma/np.sqrt(2)*(rng.standard_normal(clean.shape)+1j*rng.standard_normal(clean.shape))
    return clean+noise,sigma,noise
