#!/usr/bin/env python3
"""One public-state cost comparison: old wide QR vs exact physical small Gram.

This is a single kernel timing, not a reconstruction or a speedup theorem.
No ground truth, optimization or derivative finite differences are used.
"""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
import argparse,json,time
from pathlib import Path
import numpy as np
from scipy.linalg import qr,lstsq
from src.physics import Model

def measured(fn):
    t,c=time.perf_counter(),time.process_time();out=fn()
    return out,{'wall_seconds':time.perf_counter()-t,'cpu_seconds':time.process_time()-c}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',required=True);a=ap.parse_args();dest=Path(a.output)
    if dest.exists():raise FileExistsError(str(dest))
    root=Path(__file__).resolve().parent;ops=dict(np.load(root/'public_data/operators.npz',allow_pickle=False))
    case=next(iter(sorted((root/'public_data/dev').glob('*.npz'))));Y=np.load(case,allow_pickle=False)['Y']
    G,q,L=ops['G'],ops['E0'],ops['S'];n,s=q.shape;r=L.shape[0];m=np.ones(n)
    def old_offline():
        u,_=qr(q,mode='full',check_finite=False);Z=u[:,s:].conj().T
        return Z,Z@G,np.linalg.pinv(q)
    (Z,ZG,qdag),offline=measured(old_offline)
    def target_setup():
        W0=Y@qdag;return W0,L+W0@G
    (W0,V0),target=measured(target_setup)
    def old_gradient():
        A=Z-ZG*m[None,:];B=W0-V0*m[None,:]
        coef,_,rank,_=lstsq(A.T,-B.T,cond=None,lapack_driver='gelsy',check_finite=False)
        W=W0+coef.T@Z;V=L+W@G;E=W-V*m[None,:]
        return float(np.vdot(E,E).real),-2*np.real(np.sum(E.conj()*V,axis=0)),int(rank)
    (old_value,old_grad,rank),oldtime=measured(old_gradient)
    model=Model(ops)
    def new_gradient():
        state=model.forward(m,label='PUBLIC_SCALING_BENCHMARK')
        return model.smallgram(state,Y,gradient=True)
    (value,g,details),newtime=measured(new_gradient)
    out={'scope':'Single public-material1 kernel timing; not repeated precision validation or an imaging comparison',
         'public_case':str(case.relative_to(root)),'truth_access':'NONE','N':n,'sources':s,'receivers':r,
         'old_fixed_setup':offline,'old_per_target_setup':target,
         'old_wide_QR_objective_gradient':dict(oldtime,QR_shape=[n,n-s],LS_RHS=r,rank=rank,objective=old_value,gradient_norm=float(np.linalg.norm(old_grad))),
         'small_Gram_objective_gradient':dict(newtime,objective=value,gradient_norm=float(np.linalg.norm(g)),Gram=details,counts=dict(model.counts)),
         'wall_ratio_old_kernel_over_new_kernel':oldtime['wall_seconds']/newtime['wall_seconds'],
         'limitations':['One fixed public state, one call per representation, CPU machine also running two imaging workers',
                       'Includes old QR target setup separately; new LU/source/adjoint charged',
                       'Not a universal wall advantage over Full-J or ordinary data adjoint; same physical zeros only in invertible domain']}
    dest.parent.mkdir(parents=True,exist_ok=True);dest.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))

if __name__=='__main__':main()
