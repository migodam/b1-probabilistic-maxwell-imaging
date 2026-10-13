"""Gaussian real packing with a proper, bounded, shared calibration gain.

The truncated-normal marginal is an untruncated Gaussian marginal multiplied
by the posterior/prior interval probabilities. The determinant is retained.
"""
from __future__ import annotations
import numpy as np
from scipy.special import log_ndtr

def log_cdf_difference(a,b):
    a,b=np.broadcast_arrays(a,b)
    if np.any(b<=a): raise ValueError('normal interval must be increasing')
    # Reflect positive tails before subtracting to avoid subtracting two ones.
    high=np.where(a>0,log_ndtr(-a),log_ndtr(b))
    low=np.where(a>0,log_ndtr(-b),log_ndtr(a))
    return high+np.log(-np.expm1(low-high))

def conditional_gain(F,y,sigma,sd=.15):
    F=np.asarray(F); sigma=np.asarray(sigma); a=F/sigma; yy=np.asarray(y)/sigma
    v=1/(np.sum(a*a,axis=-1)+1/sd**2)
    m=v*(np.sum(a*yy,axis=-1)+1/sd**2)
    return m,v

def log_likelihood(F,y,sigma,gain_sd=.15,gain_bounds=(.7,1.3)):
    F=np.asarray(F); sigma=np.broadcast_to(np.asarray(sigma), (F.shape[-1],))
    if np.any(sigma<=0): raise ValueError('positive fixed real noise required')
    r=(np.asarray(y)-F)/sigma; f=F/sigma
    ff=np.sum(f*f,axis=-1); rf=np.sum(r*f,axis=-1)
    den=1+gain_sd**2*ff
    q=np.sum(r*r,axis=-1)-gain_sd**2*rf*rf/den
    logdet=2*np.log(sigma).sum()+np.log(den)
    out=-.5*(q+logdet+F.shape[-1]*np.log(2*np.pi))
    if gain_sd:
        m,v=conditional_gain(F,y,sigma,gain_sd); s=np.sqrt(v); lo,hi=gain_bounds
        correction=log_cdf_difference((lo-m)/s,(hi-m)/s)-log_cdf_difference((lo-1)/gain_sd,(hi-1)/gain_sd)
        out=out+correction
    if not np.all(np.isfinite(out)): raise FloatingPointError('nonfinite declared likelihood')
    return out

def context_loglik(x,c,sd=.01):
    from .material import void_volume
    return -.5*((void_volume(x)-c)/sd)**2-np.log(sd*np.sqrt(2*np.pi))
