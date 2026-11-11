"""R11 exposed bad endpoints: identical physical starts, no truth access."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key]='1'
from pathlib import Path
import sys,json,time,argparse,platform
import numpy as np
from scipy.optimize import minimize
from scipy.ndimage import gaussian_filter
from scipy.stats import chi2
ROOT=Path(__file__).resolve().parents[1]
REPO=ROOT.parents[1]
sys.path.insert(0,str(REPO))
if os.environ.get('MAXWELL_REUSE_REPO'):
    sys.path.insert(0,os.environ['MAXWELL_REUSE_REPO'])
from research.research_fbg_r3_experiments.src.physics import MaxwellModel
from research.research_max_struct_r12.src.field_current import DenseFieldCurrentProjectors,PhysicalAudit,run_projection

def arrays(path):
    with np.load(path,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}

class BudgetEnd(Exception):pass

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--case',required=True)
    ap.add_argument('--method',choices=['direct','multistart','alternating_projection','douglas_rachford'],required=True)
    ap.add_argument('--budget',type=float,default=240);ap.add_argument('--steps',type=int,default=24000)
    a=ap.parse_args();out=ROOT/'RESULTS/basin_v1'/f'{a.case}__{a.method}'
    out.mkdir(parents=True,exist_ok=False)
    public=arrays(ROOT/'data_public/basin/operators.npz')
    public.update(arrays(ROOT/'data_public/basin'/f'{a.case}.npz'))
    m0=arrays(ROOT/'data_public/basin'/f'{a.case}_initial.npz')['m']
    Y=public['Y'];sigma=float(public['sigma']);delta=sigma*np.sqrt(chi2.ppf(.99,2*Y.size)/2)
    config=dict(case=a.case,method=a.method,budget_cpu=a.budget,proposal_fraction=.6,scale=6.,omega=1.,
        steps=a.steps,audit_every=50,box=[.02,20.],noise_complex_rms=sigma,raw_delta=delta,
        prior='inherited smooth=1e-4 tv=1e-4 epsilon=.01',selection='minimum original Maxwell Phi including inherited prior',
        role='EXPOSED_DEVELOPMENT',truth_access='NONE',start='R11 saved bad material, recomputed complete physical fields',
        restart_protocol='include bad endpoint; cyclic sigma-pixels2 smoothed Gaussian perturbations RMS1,4,10, seeds120120 onward; clipped same box',
        hardware=platform.platform(),python=sys.version,baseline='frozen MaxwellModel exact prior and forward')
    (out/'CONFIG.json').write_text(json.dumps(config,indent=2))
    cpu=time.process_time();wall=time.perf_counter();model=MaxwellModel(public);history=[];mats=[]
    best=[float('inf'),m0.copy()];before={};count=0
    def audit(m):
        nonlocal before,count
        state=model.state(m);value=model.physical_score(m,Y)
        if value<best[0]:best[:]=[value,m.copy()]
        counts=dict(model.cost);inc={k:float(v-before.get(k,0)) for k,v in counts.items()};before=counts
        metrics=dict(relative_data_residual=float(np.linalg.norm(state['Y']-Y)/np.linalg.norm(Y)),
            noise_feasible=bool(np.linalg.norm(state['Y']-Y)<=delta),prior=float(model.reg.evaluate(m)[0]),
            lower_active=int(np.sum(m<=.02+1e-8)),upper_active=int(np.sum(m>=20-1e-8)))
        history.append(dict(candidate=count,cpu=time.process_time()-cpu,wall=time.perf_counter()-wall,
            Phi=value,incumbent=best[0],**metrics,cost=counts));mats.append(m.copy());count+=1
        return PhysicalAudit(value,state['j'].copy(),state['Q'].copy(),metrics,inc)
    def optimize(start,deadline,maxiter=3000):
        last=[start.copy()];status='not_started'
        def fun(m):
            if time.process_time()-cpu>=deadline:raise BudgetEnd()
            return model.data_value_gradient(m,Y)
        def callback(m):
            last[0]=m.copy();audit(m)
        try:
            sol=minimize(fun,start,jac=True,method='L-BFGS-B',bounds=[(.02,20.)]*model.n,
                callback=callback,options=dict(maxiter=maxiter,ftol=1e-14,gtol=1e-10,maxls=30))
            last[0]=sol.x.copy();status=str(sol.message)
        except BudgetEnd:status='cpu_budget'
        audit(last[0]);return status
    audit(m0);proposal_stop=[];projcost=None
    if a.method in ('alternating_projection','douglas_rachford'):
        pr=DenseFieldCurrentProjectors(model.G,model.E,model.S,Y,delta=delta,scale=6.,lower=.02,upper=20.)
        run=run_projection(a.method,pr,m0,audit,steps=a.steps,audit_every=50,
            cpu_budget_seconds=max(0.,.6*a.budget-(time.process_time()-cpu)))
        run.save(out/'projection');proposal_stop.append(run.stop_reason)
        projcost=dict(setup_cpu=run.setup_cpu_seconds,setup_wall=run.setup_wall_seconds,
            setup_and_total=dict(pr.costs),iteration_cost=run.costs)
    else:
        proposal_stop.append(optimize(m0,.6*a.budget))
        if a.method=='multistart':
            restart=0
            while time.process_time()-cpu<.6*a.budget:
                rng=np.random.default_rng(120120+restart)
                v=gaussian_filter(rng.normal(size=model.shape),2.).ravel();v/=np.sqrt(np.mean(v*v))
                start=np.clip(m0+[1.,4.,10.][restart%3]*v,.02,20.)
                deadline=min(.6*a.budget,time.process_time()-cpu+24.)
                proposal_stop.append(optimize(start,deadline,1000));restart+=1
    candidate=best[1].copy();proposal_cpu=time.process_time()-cpu
    # Unused proposal allowance is not converted into a privileged corrector budget.
    stop=optimize(candidate,proposal_cpu+.4*a.budget)
    final=audit(best[1]);grad=model.data_value_gradient(best[1],Y)[1]
    np.savez_compressed(out/'prediction.npz',m=best[1],J=final.J,Q=final.Q,candidates=np.asarray(mats))
    (out/'history.json').write_text(json.dumps(history,indent=2))
    receipt=dict(status='COMPLETE',method=a.method,case=a.case,objective=best[0],
        **final.metrics,projected_gradient_inf=float(np.max(abs(best[1]-np.clip(best[1]-grad,.02,20.)))),
        proposal_termination=proposal_stop,corrector_termination=stop,projection_cost=projcost,
        physics_cost=dict(model.cost),proposal_cpu=proposal_cpu,cpu_seconds=time.process_time()-cpu,
        wall_seconds=time.perf_counter()-wall,truth_access='NONE',scoring_status='SEALED_NOT_SCORED')
    (out/'RECEIPT.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt),flush=True)

if __name__=='__main__':main()
