"""Ordinary smoothing/shrink actions on a physics material point; no neural model."""
import numpy as np
from scipy.ndimage import gaussian_filter
from .fusion import prior_material_mean


def apply_action(chi,cfg,roi,sigma,blend,prior_shrink):
    """Filter spatial axes independently within each parent, then shrink to prior mean."""
    sigma,blend,prior_shrink=float(sigma),float(blend),float(prior_shrink)
    if sigma<0 or not 0<=blend<=1 or not 0<=prior_shrink<=1:
        raise ValueError('nonnegative sigma and blend/shrink in[0,1] required')
    source=np.asarray(chi)
    if source.shape[-1]!=576 or not np.iscomplexobj(source) or not np.isfinite(source).all():
        raise ValueError('finite complex material with final dimension576 required')
    grid=source.reshape(-1,24,24)
    # Zero sigma on the parent axis prevents any foreign-parent information.
    if sigma==0:blur=source.copy()
    else:
        blur=(gaussian_filter(grid.real,(0,sigma,sigma),mode='reflect',truncate=4)+
              1j*gaussian_filter(grid.imag,(0,sigma,sigma),mode='reflect',truncate=4)).reshape(source.shape)
    smooth=source+blend*(blur-source)
    result=(1-prior_shrink)*smooth+prior_shrink*prior_material_mean(cfg,roi)
    outside=np.ones(576,dtype=bool);outside[np.asarray(roi,dtype=int)]=False
    result[...,outside]=complex(*cfg['host_chi'])
    return result
