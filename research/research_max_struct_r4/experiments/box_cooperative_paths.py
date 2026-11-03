"""Refinement: solve a coupled box quadratic instead of clipping its eigenvector.

Uses the same already-paid full physical Hessian. Two deterministic starts;
this is a nonconvex quadratic heuristic, not a global QP certificate.
"""
from pathlib import Path
import sys,time,json,traceback
import numpy as np
from scipy import linalg as la
from scipy.optimize import minimize
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
sys.path.insert(0,str(REPO))
from research.research_fbg_r3_experiments.src.physics import MaxwellModel
from research.research_max_struct_r2.experiments.run_blind import dump,ALLOW,metric
from research.research_max_struct_r4.experiments.cooperative_paths import refine
from research.research_z_alg_r4_basin.src.rational_rays import scan_ray

def run():
    tick,cpu=time.perf_counter(),time.process_time()
    out=ROOT/'results/box_cooperative_paths';out.mkdir(parents=True,exist_ok=True)
    path=out/'RECEIPT.json'
    if path.exists():raise FileExistsError('Preserve development result')
    cfg=json.loads((ROOT/'configs/cooperative_paths.json').read_text());model=None
    receipt=dict(run_id='box_cooperative_paths',status='RUNNING',stage='DEV_STRUCTURAL_REFINEMENT',truth_access='NONE',
        reason='Clipping the negative eigenvector destroyed its interference cancellation and made curvature positive',
        refinement='Use the signed exact Hessian in a full bounded quadratic; keep exact original finite-ray acceptance',
        quadratic_starts='two box endpoints determined by +/- least full Hessian eigenvector',
        quadratic_maxiter=500,ray_config=cfg,configuration='frozen in this source before execution',
        no_global_QP_or_global_ray_certificate=True)
    dump(path,receipt)
    try:
        with np.load(ROOT/'results/cooperative_paths/full_geometry.npz') as f:
            m=f['m'];g=f['gradient'];H=f['H']
        with np.load(REPO/f"research/research_max_struct_r2/data_public/{cfg['case_id']}.npz") as f:
            public={k:f[k] for k in f.files if k in ALLOW}
        model=MaxwellModel(public);Y=public['Y'];state=model.state(m);base=model.physical_score(m,Y)
        eig,v=la.eigh(H,subset_by_index=[0,0]);seeds=[];records=[]
        for sign in [-1,1]:
            initial=np.where(sign*v[:,0]>0,20.-m,.02-m)
            history=[]
            def quadratic(d):
                h=H@d;value=float(g@d+.5*d@h);grad=g+h
                history.append(value)
                return 1e6*value,1e6*grad
            fit=minimize(quadratic,initial,jac=True,method='L-BFGS-B',bounds=list(zip(.02-m,20.-m)),
                options=dict(maxiter=500,maxfun=2000,ftol=1e-12,gtol=.01,maxls=30))
            d=fit.x
            rec=dict(sign=sign,qp_stop=str(fit.message),qp_nit=int(fit.nit),qp_nfev=int(fit.nfev),
                     quadratic_value=float(fit.fun/1e6),qp_history=history,gradient_dot=float(g@d),
                     curvature=float(d@H@d),active_moved=int(np.count_nonzero((m<=.020001)&(d>1e-6))))
            proposed,rec=scan_ray(model._physics,state,d,Y,float(np.vdot(Y,Y).real),model.reg,cfg,'box_qp_'+str(sign),rec)
            records.append(rec);seeds+=proposed
            print(json.dumps({k:rec[k] for k in ('sign','qp_stop','quadratic_value','gradient_dot','curvature','active_moved')}
                             |dict(best_ray_Phi=min(x['raw_score'] for x in proposed))),flush=True)
        dump(out/'RAYS.json',records)
        best=min(seeds,key=lambda x:x['raw_score'])
        seed=best['m'] if best['raw_score']<base else m.copy()
        raw=metric(model,Y,seed)
        if raw['objective']<base:final,corr=refine(model,Y,seed)
        else:final=m.copy();corr=dict(stop='NO_ACTUAL_RAY_IMPROVEMENT',cost={})
        dump(out/'CORRECTOR.json',corr)
        finalrow=min([metric(model,Y,x) for x in (m,seed,final)],key=lambda x:x['objective'])
        st=model.state(finalrow['m'])
        np.savez_compressed(out/'prediction.npz',m=finalrow['m'],Y_prediction=st['Y'])
        receipt.update(status='SEALED_DEVELOPMENT_PREDICTION',anchor_Phi=base,best_ray_Phi=best['raw_score'],
                       final_Phi=finalrow['objective'],physical_residual=finalrow['relative_residual'],
                       material_distance=float(la.norm(finalrow['m']-m)),corrector_stop=corr['stop'])
    except Exception as e:
        receipt.update(status='FAILED',error=repr(e),traceback=traceback.format_exc());raise
    finally:
        receipt.update(physical_cost=dict(model.cost) if model else {},wall_seconds=time.perf_counter()-tick,
                       cpu_seconds=time.process_time()-cpu,shared_geometry='already charged cooperative_paths receipt; no second Hessian charge')
        dump(path,receipt)
    print(json.dumps(receipt))

if __name__=='__main__':
    with threadpool_limits(limits=1):run()
