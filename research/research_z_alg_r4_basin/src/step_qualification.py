"""One substantive experiment: cost and conservatism of a resolvent certificate."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key]='1'
import json,sys,time
from pathlib import Path
import numpy as np
from scipy.linalg import lu_solve
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT.parents[1]))
from research.research_z_alg_r3_digits.src.physics import Model
HISTORY=ROOT.parent/'research_z_alg_r3_digits'

def main():
    begun,cpu=time.perf_counter(),time.process_time()
    model=Model(dict(np.load(HISTORY/'public_data/operators.npz',allow_pickle=False)))
    rows=[]
    cases=[('DEV7','dev_v2','dev_2a5c221387b1'),('holdout7','holdout','obs_c3604cf3555d'),
           ('holdout9','holdout','obs_42c48e795766'),('holdout8','holdout','obs_a63df7f20566')]
    with threadpool_limits(limits=1):
        for label,split,case in cases:
            target=np.load(HISTORY/'public_data'/('dev' if split=='dev_v2' else split)/(case+'.npz'))['Y']
            old=json.loads((HISTORY/'results'/split/'fullj_native'/(case+'.json')).read_text())
            material=np.asarray(old['top1']['m'])
            state=model.forward(material)
            jac=model.jacobian(state)
            residual=state['y']-target
            gradient=jac.T@model.pack(residual)
            direction=-gradient/np.sqrt(np.mean(gradient**2))
            with model.operation('DENSE_CERTIFICATE_PREPARATION',rhs=model.n+model.nr,
                                 certificate_inverse_RHS=model.n,certificate_receiver_RHS=model.nr):
                GM=lu_solve(model.factor,model.G.T,trans=1,check_finite=False).T
                SM=lu_solve(model.factor,model.S.T,trans=1,check_finite=False).T
                rownorm2=np.sum(abs(GM)**2,axis=1)
                snorm=float(np.linalg.norm(SM,2))
            rnorm=float(np.linalg.norm(residual))
            for radius in (.001,.01,.1,1.):
                trial=np.clip(material+radius*direction,.5,20.)
                h=trial-material
                beta=float(np.sqrt(np.dot(h*h,rownorm2)))
                tangent=(jac@h).reshape(model.ns,model.nr,2)
                tangent=(tangent[:,:,0]+1j*tangent[:,:,1]).T
                error_bound=snorm*beta/(1-beta)*np.linalg.norm(h[:,None]*state['Q']) if beta<1 else None
                qualified=bool(beta<1 and np.linalg.norm(residual+tangent)+error_bound<rnorm)
                nextstate=model.forward(trial)
                actual_error=float(np.linalg.norm(nextstate['y']-state['y']-tangent))
                rows.append(dict(case=label,requested_rms=radius,actual_rms=float(np.sqrt(np.mean(h*h))),
                    beta_F=beta,invertibility_sufficient_condition=beta<1,descent_qualified=qualified,
                    rigorous_absolute_error_bound=float(error_bound) if error_bound is not None else None,
                    actual_remainder=actual_error,actual_data_descent=bool(np.linalg.norm(nextstate['y']-target)<rnorm),
                    before_relative_residual=rnorm/np.linalg.norm(target),
                    after_relative_residual=float(np.linalg.norm(nextstate['y']-target)/np.linalg.norm(target))))
    out=dict(status='OFFLINE_DEVELOPMENT_ONLINE_COMPUTABLE_QUANTITIES',
       purpose='Measure whether a mathematically sufficient resolvent bound usefully qualifies finite steps; no formula precision verification',
       rows=rows,counts=dict(model.counts),wall_seconds=time.perf_counter()-begun,cpu_seconds=time.process_time()-cpu,
       conclusion='The Frobenius bound is sufficient, never necessary. Rejection does not prove singularity or failure.')
    (ROOT/'diagnostics/operator_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'qualified':sum(r['descent_qualified'] for r in rows),'actual_descents':sum(r['actual_data_descent'] for r in rows),
                      'rows':len(rows),'counts':out['counts'],'wall_seconds':out['wall_seconds']}))

if __name__=='__main__':main()
