"""Ordinary final-stage strong-mode release using mature native CIE code.

Only the final-stage current dictionary changes: [VD,Q]. Strong coefficients
start at zero through new frequency keys; old Fourier coefficients carry over
unchanged. No historical source file is edited. Not thread-safe by design:
the experiment is sequential CPU1 and restores the forwarding hook afterward.
"""
import numpy as np
from research.research_fbg_r3_experiments.src import cie_tsom_baseline as baseline

def run_strong_release(model,Y,m0,config):
    cfg=dict(config);last_mf=int(cfg['stages'][-1]['mf'])
    if any(int(s['mf'])==last_mf for s in cfg['stages'][:-1]):
        raise ValueError('Final release mask must be distinct from prior stages')
    original=baseline.fourier_apcs
    def dictionary(shape,mf,VD,cost=None):
        Q,info=original(shape,mf,VD,cost)
        if int(mf)!=last_mf:
            return Q,dict(info,released_DPCS_rank=0,current_basis_columns=Q.shape[1],
                          independent_current_spatial_rank=info['apcs_rank'])
        # New names are absent from old_coefficients, hence exact zero seed.
        keys=[('DPCS_RELEASE',i)for i in range(VD.shape[1])]
        info=dict(info,frequencies=keys+info['frequencies'],released_DPCS_rank=VD.shape[1],
                  current_basis_columns=VD.shape[1]+Q.shape[1],
                  independent_current_spatial_rank=VD.shape[1]+info['apcs_rank'])
        if cost is not None:cost['released_DPCS_columns']+=VD.shape[1]
        return np.column_stack([VD,Q]),info
    baseline.fourier_apcs=dictionary
    try:
        result=baseline.run_cie(model,Y,m0,cfg)
        result['control']='ordinary strong-mode release final stage only; same data term/readout/PR'
        return result
    finally:baseline.fourier_apcs=original
