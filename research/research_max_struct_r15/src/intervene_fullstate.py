"""Fixed public-data full-state interventions; no truth argument or lookup."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):os.environ[k]='1'
import argparse,json,sys,time,traceback
from pathlib import Path

class EndBudget(BaseException):pass

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('root','legacy-repo','source-checkpoint','public','geometry'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--case-id',required=True);p.add_argument('--method',required=True)
    p.add_argument('--band',choices=('base','U11','U01','U00'),default='base');p.add_argument('--sign',type=int,choices=(-1,1),default=1);p.add_argument('--cpu',type=int)
    a=p.parse_args();start=time.perf_counter();root=a.root.resolve()
    import psutil
    process=psutil.Process()
    if a.cpu is not None:process.cpu_affinity([a.cpu])
    if os.name=='nt':process.nice(psutil.NORMAL_PRIORITY_CLASS)
    sys.path[:0]=[str(root/'src'),str(root/'runtime'),str(a.legacy_repo),str(a.legacy_repo/'vendor')]
    import numpy as np
    from twofold_solver import TwofoldG19
    from research.research_fbg_r3_experiments.src.physics import MaxwellModel
    ident=a.case_id+'__'+a.method+'__'+a.band+(str(a.sign) if a.band!='base' else '')
    out=root/'RESULTS'/'fullstate_interventions_v1'/ident;out.mkdir(parents=True,exist_ok=False)
    def write(name,value):
        def encode(x):
            if isinstance(x,np.generic):return x.item()
            if isinstance(x,np.ndarray):return x.tolist()
            raise TypeError(type(x).__name__)
        (out/name).write_text(json.dumps(value,default=encode,indent=2,allow_nan=False)+'\n')
    with np.load(a.public,allow_pickle=False) as f:public={k:f[k].copy() for k in f.files}
    model=MaxwellModel(public);Y=public['Y'];norm2=float(np.vdot(Y,Y).real);delta=float(public['delta'])
    rows=[];materials=[];best=None;solver=None;status='COMPLETE';error=None;cap=60.
    def check():
        if time.process_time()>=cap:raise EndBudget()
    def audit(m,metadata):
        nonlocal best
        check();state=model.state(m);res=state['Y']-Y;D=.5*float(np.vdot(res,res).real)/norm2;R=float(model.reg.evaluate(m)[0]);t=time.process_time()
        row=dict(candidate=len(rows),D=D,R=R,Phi=D+R,relative_residual=float(np.linalg.norm(res)/np.sqrt(norm2)),residual_over_delta=float(np.linalg.norm(res)/delta),
            CPU_seconds=t,eligible=bool(t<=cap and np.isfinite(m).all() and m.min()>=float(public['lower'])-1e-12 and m.max()<=float(public['upper'])+1e-12),metadata=metadata,
            physics_cost=dict(model.cost))
        rows.append(row);materials.append(np.asarray(m).copy())
        if row['eligible'] and (best is None or row['Phi']<best['row']['Phi']):best=dict(m=np.array(m,copy=True),row=row)
        check()
    native=[];perturbation={};updates=64 if a.band=='base' else 16
    try:
        cfg=dict(regularizer=lambda m:model.reg.evaluate(m)[:2],aplus_regularizer=model.reg,budget_check=check,r15_data_delta=delta)
        solver=TwofoldG19.from_checkpoint(model,Y,cfg,a.source_checkpoint,out=out,audit_callback=audit,method_override=a.method)
        # The inherited counters are retained in the state; report their initial
        # values separately so a fork never makes old operations newly paid.
        inherited_aux=dict(solver.cost);inherited_projection={} if solver.projector is None else dict(solver.projector.costs)
        if a.band!='base':
            with np.load(a.geometry,allow_pickle=False) as g:U=g[a.band].copy()
            rng=np.random.default_rng(20261116);h=U@rng.standard_normal(U.shape[1]);h=np.asarray(h.real);h/=max(np.linalg.norm(h),1e-300)
            before=solver.m.copy();requested=a.sign*1e-4*max(np.linalg.norm(before),1.)*h
            solver.m=np.clip(before+requested,solver.lo,solver.hi);realized=solver.m-before
            np.savez_compressed(out/'PERTURBATION.npz',requested=requested,realized=realized,unperturbed_material=before)
            perturbation=dict(band=a.band,sign=a.sign,seed=20261116,requested_norm=float(np.linalg.norm(requested)),realized_norm=float(np.linalg.norm(realized)),
                clipped_coordinates=int(np.count_nonzero(requested!=realized)),interpretation='Material-only finite perturbation; all current/PR/reflection memory retained')
        solver.checkpoint(out/'FORK_INITIAL.npz');audit(solver.m,dict(source='complete-state fork initial',perturbation=perturbation))
        for local in range(1,updates+1):
            check()
            try:row=solver.step()
            except StopIteration:status='SCHEDULE_COMPLETE';break
            scalar={k:v for k,v in row.items() if not isinstance(v,np.ndarray)};scalar.update(local_update=local,CPU_seconds=time.process_time());native.append(scalar)
            if local in (1,16,64):
                np.savez_compressed(out/f'STATE_{local:03d}.npz',m=solver.m,J=solver.J,Q=solver.Q,A=solver.A,
                    stage=np.asarray(solver.stage),iteration=np.asarray(solver.iteration),beta=np.asarray([solver.beta0,solver.beta1]),
                    memory_X=solver.memory_X,memory_Q=solver.memory_Q)
        audit(solver.m,dict(source='complete-state fork final',completed_updates=len(native)))
        solver.checkpoint(out/'FORK_FINAL.npz')
    except EndBudget:status='WHOLE_PROCESS_CPU'
    except Exception as exc:status='FAILED';error=repr(exc);traceback.print_exc()
    finally:
        if best is not None:np.savez_compressed(out/'primary.npz',m=best['m']);write('primary.json',best['row'])
        np.savez_compressed(out/'CANDIDATES.npz',m=np.array(materials) if materials else np.empty((0,model.n)))
        write('HISTORY.json',rows);write('NATIVE_HISTORY.json',native)
        write('RECEIPT.json',dict(status=status,error=error,case_id=a.case_id,method=a.method,band=a.band,sign=a.sign,source_checkpoint=str(a.source_checkpoint),
            requested_updates=updates,completed_updates=len(native),CPU_cap=cap,whole_process_CPU_seconds=time.process_time(),wall_seconds=time.perf_counter()-start,
            physics_cost=dict(model.cost),auxiliary_cumulative={} if solver is None else dict(solver.cost),
            auxiliary_inherited=locals().get('inherited_aux',{}),fixed_projection_cumulative={} if solver is None or solver.projector is None else dict(solver.projector.costs),
            fixed_projection_inherited=locals().get('inherited_projection',{}),corrector=False,truth_access='NONE',perturbation=perturbation))
    return 1 if status=='FAILED' else 0

if __name__=='__main__':sys.exit(main())
