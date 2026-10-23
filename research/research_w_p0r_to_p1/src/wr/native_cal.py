"""Round 2: affine material coordinates; PHYSICAL EXISTENCE, not online.

The material/noise contract and the original allowance are unchanged. The
search chart removes coefficient saturation and allocates local wall freedom.
"""
import sys,json,time
import numpy as np
from scipy.optimize import minimize
from .native_world import LinearWorld
from .physical_witness import mesh_for
from p0.geometry import probes
from .physics import Forward,pack,packj
from .accounting import ROOT,ScopedBook,write,BudgetStop
from .optimization import trust_region

def bounded_task(j,k,metric,rms):
    """Sequential linear task maximization; not a global certificate."""
    u,_,_=np.linalg.svd(k,full_matrices=False)
    trials=[]
    for index in range(3):
        direction=u[:,index];x=np.zeros(j.shape[1])
        for iteration in range(3):
            c=k.T@direction
            fit=minimize(lambda q:-float(c@q),x,jac=lambda q:-c,
                method='SLSQP',bounds=[(-1.,1.)]*len(x),constraints=[
                dict(type='ineq',fun=lambda q:.8**2-float((j@q)@(j@q)),
                     jac=lambda q:-2*j.T@(j@q)),
                dict(type='ineq',fun=lambda q:rms**2-float(q@metric@q),
                     jac=lambda q:-2*metric@q)],
                options=dict(maxiter=160,ftol=1e-11))
            # Numerical subproblem failure cannot authorize an infeasible edit.
            x=np.clip(fit.x,-1.,1.)
            factor=min(1.,.8/max(np.linalg.norm(j@x),1e-30),
                       rms/max(np.sqrt(max(0.,x@metric@x)),1e-30))
            x*=factor
            kx=k@x
            trials.append(dict(start=index,iteration=iteration,success=bool(fit.success),
                message=str(fit.message),predicted_cal=float(np.linalg.norm(j@x)),
                predicted_task=float(np.linalg.norm(kx)),material_RMS=float(np.sqrt(max(0,x@metric@x))),x=x.copy()))
            if np.linalg.norm(kx)>0:direction=kx/np.linalg.norm(kx)
    best=max(trials,key=lambda r:r['predicted_task'])
    return best['x'],trials

def main(seed,band='narrow'):
    scope=f'cal_{seed}_{band}';book=ScopedBook('B2',scope)
    out=ROOT/f'results/B2/native_cal_{seed}_{band}';out.mkdir(parents=True,exist_ok=False)
    freq=json.loads((ROOT/'configs/ORIGINAL_P0.json').read_text())['frequencies'][band]
    p=probes()[8];model=LinearWorld(seed);tm=LinearWorld(seed,p);zero=np.zeros(model.n)
    folder=str(out.relative_to(ROOT/'results'))
    fw=Forward(model,mesh_for(seed,p),freq,book,folder+'/cal',phase='B2_native_cal',identity='PHYSICAL_EXISTENCE_SEARCH',h=.145,pml_strength=6,retain_current=True)
    ft=Forward(tm,mesh_for(seed,p),freq,book,folder+'/task',phase='B2_native_task',identity='PHYSICAL_EXISTENCE_SEARCH',h=.145,pml_strength=6)
    d0,j0=fw.evaluate(zero,jac=True);t0,k0=ft.evaluate(zero,jac=True)
    j,k=packj(j0),packj(k0)
    grid=np.stack(np.meshgrid(np.linspace(-2,2,201),np.linspace(0,3,151)),axis=-1).reshape(-1,2)
    dm=model.derivative(zero,grid,1.);material0=model.coefficient(zero,1.)(grid)
    metric=np.real(dm.conj().T@dm)/len(grid)
    target,linear=bounded_task(j,k,metric,.04)
    write(out/'LINEAR_SEARCH.json',linear)
    wanted=float(np.sqrt(target@metric@target))
    # The matched random/ordinary controls get the identical actual radius.
    rng=np.random.default_rng(2026102315+seed)
    rand=rng.normal(size=model.n);rand*=wanted/max(np.sqrt(rand@metric@rand),1e-30)
    rand*=min(1.,1/max(abs(rand)))
    vals,vec=np.linalg.eigh(metric);active=vals>vals.max()*1e-8
    transform=vec[:,active]/np.sqrt(vals[active])[None,:]
    _,directions=np.linalg.eigh(transform.T@(j.T@j)@transform)
    weak=None
    for col in range(directions.shape[1]):
        candidate=transform@directions[:,col]*wanted
        if max(abs(candidate))<=1.:
            weak=candidate;break
    np.savez_compressed(out/'BASE.npz',J=j0,K=k0,cal=d0,task=t0,material_metric=metric,grid=grid,target=target,material0=material0)
    # Actual full-physics derivative regression for the changed chart, paid.
    fd=[];v=target/max(np.linalg.norm(target),1e-30)
    for h in [1e-3,1e-4,1e-5]:
        a,_=fw.evaluate(h*v);b,_=fw.evaluate(-h*v)
        err=float(np.linalg.norm((a-b).ravel()/(2*h)-j0@v)/max(np.linalg.norm(j0@v),1e-30))
        fd.append(dict(step=h,relative_error=err));write(out/'DERIVATIVES.json',fd)
    if min(r['relative_error'] for r in fd)>1e-5:raise RuntimeError('native chart Full-J mismatch')
    rows=[]
    for method,proposal in [('TASK_AWARE',target),('RANDOM_MATCHED',rand),('FULL_J_WEAK',weak)]:
        if proposal is None:
            rows.append(dict(method=method,status='NOT_RUN_NO_BOX_FEASIBLE_MATCHED_DIRECTION',wanted_RMS=wanted));write(out/'RESULTS.json',rows);continue
        start_rhs=book.totals()['B2']['rhs'];start_cpu=time.process_time();start_wall=time.perf_counter()
        raw,_=fw.evaluate(proposal);rawt,_=ft.evaluate(proposal)
        delta=dm@proposal;norm=float(np.vdot(delta,delta).real)
        pg=np.real(delta.conj()@dm)/norm
        history=[]
        def ev(x,jac=False):
            d,jx=fw.evaluate(x,jac=jac);progress=float(pg@x)
            return np.r_[pack(d-d0),100*(progress-1)],None if jx is None else np.vstack([packj(jx),100*pg])
        def callback(r):history.append(r);write(out/f'{method}_ITERATIONS.json',history)
        fit=trust_region(ev,proposal,np.full(model.n,-1.),np.full(model.n,1.),
             max_iter=8,max_calls=20,radius=.25,max_radius=.75,l2_radius=False,
             callback=callback,objective_floor=.05,fatal_exceptions=(Exception,))
        d,_=fw.evaluate(fit.x);td,_=ft.evaluate(fit.x)
        field=model.coefficient(fit.x,1.)(grid)
        dc=float(np.linalg.norm(d-d0)/.01);dt=float(np.linalg.norm(td-t0)/.01)
        row=dict(method=method,seed=seed,band=band,chart='108 affine material hats',
            raw_material_RMS=float(np.sqrt(proposal@metric@proposal)),wanted_RMS=wanted,
            material_RMS=float(np.sqrt(fit.x@metric@fit.x)),signed_progress=float(pg@fit.x),
            raw_d_cal=float(np.linalg.norm(raw-d0)/.01),raw_d_T=float(np.linalg.norm(rawt-t0)/.01),
            d_cal=dc,d_T=dt,cal_threshold=bool(dc<=1),task_threshold=bool(dt>=3),
            stop=str(fit.message),iterations=fit.nit,gradient=float(fit.projected_gradient),
            rhs=book.totals()['B2']['rhs']-start_rhs,cpu_seconds=time.process_time()-start_cpu,
            wall_seconds=time.perf_counter()-start_wall,status='NOT_YET_INDEPENDENTLY_VALIDATED',
            identity='ORACLE_CONSTRUCTED_PHYSICAL_WITNESS',posterior_probability=None)
        np.savez_compressed(out/f'{method}_FINAL.npz',x=fit.x,proposal_x=proposal,cal=d,task=td,
              cal_base=d0,task_base=t0,material_grid=grid,material=field,reference_material=material0)
        rows.append(row);write(out/'RESULTS.json',rows);print(json.dumps(row),flush=True)
    write(out/'COMPLETE.json',dict(totals=book.totals(),no_new_allowance=True))

if __name__=='__main__':main(int(sys.argv[1]),*sys.argv[2:])
