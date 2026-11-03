"""Finite cooperative paths using the existing exact R4-BASIN ray backend."""
from pathlib import Path
import sys,time,json,traceback
from collections import Counter
import numpy as np
from scipy import linalg as la
from scipy.optimize import minimize
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
sys.path.insert(0,str(REPO))
from research.research_fbg_r3_experiments.src.physics import MaxwellModel
from research.research_max_struct_r2.experiments.run_blind import dump,ALLOW,metric
from research.research_max_struct_r4.src.cooperative_geometry import curvature
from research.research_z_alg_r4_basin.src.rational_rays import scan_ray,feasible_interval


def refine(model,Y,m,maximum=1000):
    before=Counter(model.cost);history=[];scale=1e6;start=time.perf_counter()
    best=[m.copy(),model.physical_score(m,Y)]
    def objective(x):
        value,g=model.data_value_gradient(x,Y)
        history.append([float(value),float(la.norm(g,np.inf))])
        if value<best[1]:best[:]=[x.copy(),float(value)]
        return scale*value,scale*g
    result=minimize(objective,m,method='L-BFGS-B',jac=True,bounds=[(.02,20.)]*model.n,
        options=dict(maxiter=maximum,maxfun=4*maximum,ftol=1e-12,gtol=.01,maxls=30))
    return best[0],dict(stop=str(result.message),iterations=int(result.nit),evaluations=int(result.nfev),
        history=history,cost=dict(Counter(model.cost)-before),wall_seconds=time.perf_counter()-start)


def run():
    started,cpu=time.perf_counter(),time.process_time()
    out=ROOT/'results/cooperative_paths';out.mkdir(parents=True,exist_ok=True)
    path=out/'RECEIPT.json'
    if path.exists():raise FileExistsError('Preserve completed development outcome')
    cfg=json.loads((ROOT/'configs/cooperative_paths.json').read_text())
    model=None
    receipt=dict(run_id='cooperative_full_paths',stage='DEV_CALIBRATION',status='RUNNING',truth_access='NONE',
                 configuration=cfg,limitations='finite sampled scalar search is not global ray certification')
    dump(path,receipt)
    try:
        with np.load(REPO/f"research/research_max_struct_r2/data_public/{cfg['case_id']}.npz") as f:
            public={k:f[k] for k in f.files if k in ALLOW}
        m=np.load(REPO/'research/research_max_struct_r3/results/development/scaled_original_restart/prediction.npz')['m']
        model=MaxwellModel(public);Y=public['Y'];value,gradient=model.data_value_gradient(m,Y)
        geom=curvature(model,m,Y,np.arange(model.n));state=geom['state']
        np.savez_compressed(out/'full_geometry.npz',m=m,gradient=gradient,H=geom['H'],GN=geom['GN'],
                            feedback=geom['feedback'],regularizer=geom['regularizer'])
        seeds={};rays=[];summaries={};cost_before=Counter(model.cost)
        for arm,H in [('cooperative',geom['H']),('weak_gn',geom['GN']+geom['regularizer'])]:
            with model._physics.operation('FULL_DIRECTION_EIGEN',symmetric_eigenproblems=1,eigen_dimension=model.n):
                vals,vecs=la.eigh(H,subset_by_index=[0,0],check_finite=False)
            proposals=[]
            for sign in [-1,1]:
                p=sign*vecs[:,0].copy()
                p[(m<=.02+1e-6)&(p<0)]=0
                p[(m>=20.-1e-6)&(p>0)]=0
                p/=np.sqrt(np.mean(p*p))
                lo,hi=feasible_interval(m,p,.02,20.)
                predicted_linear=float(gradient@p);predicted_quadratic=float(p@geom['H']@p)
                record=dict(selected_eigenvalue=float(vals[0]),direction_source=arm,sign=sign,
                    gradient_dot=predicted_linear,exact_curvature=predicted_quadratic,
                    feedback_curvature=float(p@geom['feedback']@p),GN_curvature=float(p@geom['GN']@p),
                    quadratic_crossing=(2*predicted_linear/-predicted_quadratic if predicted_quadratic<0 else None))
                found,record=scan_ray(model._physics,state,p,Y,float(np.vdot(Y,Y).real),model.reg,
                                      cfg,arm+str(sign),record)
                proposals+=found;rays.append(record)
                print(json.dumps(dict(arm=arm,sign=sign,eigenvalue=float(vals[0]),interval=[lo,hi],
                    curvature=predicted_quadratic,best_ray_Phi=min(x['raw_score'] for x in found),anchor_Phi=value)),flush=True)
            best=min(proposals,key=lambda x:x['raw_score'])
            seeds[arm]=best['m'] if best['raw_score']<value else m.copy()
            summaries[arm]=dict(best_ray_Phi=best['raw_score'],anchor_Phi=float(value),
                               ray_material_distance=float(la.norm(seeds[arm]-m)),ray_label=best['ray'])
        dump(out/'RAYS.json',rays)
        receipt['ray_cost']=dict(Counter(model.cost)-cost_before)
        for arm,seed in seeds.items():
            raw=metric(model,Y,seed)
            if raw['objective']<value:
                final,corr=refine(model,Y,seed)
            else:
                final=m.copy();corr=dict(stop='NO_PHYSICAL_RAY_IMPROVEMENT',cost={})
            rows=[metric(model,Y,x) for x in (m,seed,final)];best=min(rows,key=lambda x:x['objective'])
            state2=model.state(best['m'])
            np.savez_compressed(out/f'{arm}_prediction.npz',m=best['m'],Y_prediction=state2['Y'])
            dump(out/f'{arm}_CORRECTOR.json',corr)
            summaries[arm].update(final_Phi=best['objective'],physical_residual=best['relative_residual'],
                                  material_distance=float(la.norm(best['m']-m)),corrector_stop=corr['stop'])
        receipt.update(status='SEALED_DEVELOPMENT_PREDICTIONS',arms=summaries,
                       candidate_selection='same original Phi, original prior, includes unchanged incumbent')
    except Exception as e:
        receipt.update(status='FAILED',error=repr(e),traceback=traceback.format_exc());raise
    finally:
        receipt.update(physical_cost=dict(model.cost) if model else {},
                       wall_seconds=time.perf_counter()-started,cpu_seconds=time.process_time()-cpu)
        dump(path,receipt)
    print(json.dumps(receipt))

if __name__=='__main__':
    with threadpool_limits(limits=1):run()
