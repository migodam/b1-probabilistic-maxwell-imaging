"""Risk actions, never posterior distributions or decoder-legality claims."""
import numpy as np
from .risk_models import load_predictor
from .residual_risk import load_residual_predictor

def passive_box(chi,config,roi):
    """Projection onto a convex rectangle containing every original P1 field.

    For every legal truth this cannot increase squared material error. The
    projected mean need not equal one decoded material or fit y by itself.
    """
    a=np.asarray(chi).copy();host=complex(*config['host_chi'])
    a[...,roi]=np.clip(a[...,roi].real,host.real,host.real+config['inclusion_increment_real_interval'][1])+1j*np.clip(a[...,roi].imag,host.imag,host.imag+config['inclusion_increment_imag_interval'][1])
    outside=np.ones(576,bool);outside[np.asarray(roi,int)]=False;a[...,outside]=host
    return a

def empirical_blend_selection(anchor,neural,truth):
    """Fixed candidate grid; validation selection, not an unbiased risk claim."""
    rows=[]
    for alpha in np.linspace(0,1,21):
        a=anchor+alpha*(neural-anchor)
        rows.append(dict(alpha=float(alpha),complex_MSE=float(np.mean(abs(a-truth)**2))))
    return min(rows,key=lambda x:(x['complex_MSE'],x['alpha'])),rows

def prior_material_mean(config,roi):
    # Zero-mean symmetric Gaussian f and h(s)+h(-s)=1 give E[fraction]=1/2.
    host=complex(*config['host_chi']);a=np.full(576,host,dtype=complex)
    delta=complex(np.mean(config['inclusion_increment_real_interval']),np.mean(config['inclusion_increment_imag_interval']))
    a[np.asarray(roi,int)]+=delta/2
    return a
