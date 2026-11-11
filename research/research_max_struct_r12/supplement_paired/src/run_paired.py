"""Paid, public-only adapter over frozen Maxwell/G19/AP/DR algorithms."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
import sys,time,json,argparse,traceback
from pathlib import Path
from collections import Counter

class PhaseEnd(BaseException):pass

def dump(path,obj):
    import numpy as np
    def cv(x):
        if isinstance(x,np.ndarray):return x.tolist()
        if isinstance(x,np.generic):return x.item()
        if isinstance(x,Path):return str(x)
        raise TypeError(type(x).__name__)
    Path(path).write_text(json.dumps(obj,indent=2,default=cv,allow_nan=False)+'\n')

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,required=True);ap.add_argument('--legacy-repo',type=Path,required=True)
    ap.add_argument('--task-id',required=True);ap.add_argument('--cpu',type=int,required=True)
    args=ap.parse_args();root=args.root.resolve();wall=time.perf_counter()
    import psutil
    proc=psutil.Process();proc.cpu_affinity([args.cpu])
    if os.name=='nt':proc.nice(psutil.NORMAL_PRIORITY_CLASS)
    sys.path.insert(0,str(args.legacy_repo));sys.path.insert(0,str(root/'runtime'))
    import numpy as np
    from scipy.stats import chi2
    from research.research_fbg_r3_experiments.src.physics import MaxwellModel
    from research.research_max_struct_r8.src import benchmark_solver as solver
    from research.research_max_struct_r11.src.material_readouts import prepare_full,full_material_readout
    from research.research_max_struct_r10.src.twod_adapter import noise_ball_rank
    from field_current import DenseFieldCurrentProjectors,PhysicalAudit,run_projection
    tasks=json.loads((root/'configs/TASKS.json').read_text())['tasks']
    task=next(t for t in tasks if t['task_id']==args.task_id)
    cfg=json.loads((root/'configs/FROZEN_PROTOCOL.json').read_text());limits=cfg['budgets'][task['cohort']]
    cap=float(limits['whole_child']);deadline=float(limits['proposal']) if task['method']!='direct' else cap
    out=root/'RESULTS'/args.task_id;out.mkdir(parents=True,exist_ok=False)
    source=root/task['input_file']
    with np.load(source,allow_pickle=False) as z:public={k:z[k].copy() for k in z.files}
    m0=public.pop('m_init');Y=public['Y'];sigma=float(public['sigma'])
    delta=float(sigma*np.sqrt(chi2.ppf(.99,2*Y.size)/2));model=MaxwellModel(public)
    # Explicitly retain inherited material prior; it is not reweighted for rescue.
    phase='initial';history=[];materials=[];best=None;secondary=None;phase_best={};projection_history=[]
    snapshots={};physics_before={};g19_cost={};projection_cost=None;projection_defect=None
    timings={'imports_and_input_CPU_seconds':time.process_time()};terminations={}
    def check():
        if time.process_time()>=deadline:raise PhaseEnd()
    def record(m,state,value=None,gradient=None):
        nonlocal best,secondary
        residual=state['Y']-Y;R=float(model.reg.evaluate(m)[0]);D=.5*float(np.vdot(residual,residual).real/np.vdot(Y,Y).real)
        val=float(D+R if value is None else value)
        row=dict(m=np.asarray(m).copy(),j=state['j'].copy(),Q=state['Q'].copy(),Y_prediction=state['Y'].copy(),objective=val)
        if gradient is not None:row['gradient']=np.asarray(gradient).copy()
        t=time.process_time();eligible=t<=deadline and t<=cap
        meta=dict(candidate=len(history),CPU_seconds=t,wall_seconds=time.perf_counter()-wall,phase=phase,Phi=val,D=D,R=R,
            relative_residual=float(np.linalg.norm(residual)/np.linalg.norm(Y)),residual_norm=float(np.linalg.norm(residual)),
            residual_over_delta=float(np.linalg.norm(residual)/delta),noise_feasible=bool(np.linalg.norm(residual)<=delta),
            legal=bool(np.min(m)>=.02 and np.max(m)<=20 and np.isfinite(m).all()),eligible=eligible,
            deadline=deadline,physics_cost=dict(model.cost),lower_active=int(np.sum(m<=.02+1e-8)),upper_active=int(np.sum(m>=20-1e-8)))
        row['metadata']=meta
        history.append(meta);materials.append(row['m'])
        if eligible and (best is None or val<best['objective']):best=row
        if eligible and (phase not in phase_best or val<phase_best[phase]['objective']):phase_best[phase]=row
        if eligible and meta['legal'] and meta['noise_feasible'] and (secondary is None or (R,val)<(secondary['metadata']['R'],secondary['objective'])):secondary=row
        return row
    class Observer(solver.PaidStateObserver):
        def record(self,m,state,value=None,gradient=None):
            row=record(m,state,value,gradient)
            self.best=best or row
            check()
            return row
        def state(self,*a,**kw):check();return super().state(*a,**kw)
        def physical_score(self,*a,**kw):check();return super().physical_score(*a,**kw)
        def data_value_gradient(self,*a,**kw):check();return super().data_value_gradient(*a,**kw)
    observer=Observer(model,Y,lambda name,row:None,gradient_checkpoints=())
    def physical(m):
        nonlocal physics_before
        check();state=model.state(m);row=record(m,state)
        now=dict(model.cost);cost={k:now.get(k,0)-physics_before.get(k,0) for k in set(now)|set(physics_before)};physics_before=now
        check()
        return PhysicalAudit(row['objective'],row['j'],row['Q'],metrics=row['metadata'],costs=cost)
    def checkpoint(item):
        projection_history.append(dict(iteration=item.iteration,CPU_seconds=time.process_time(),shadow_gap=item.gap,
            Phi=item.objective,projector_cost=dict(item.costs)))
        check()
    status='COMPLETE';error=None;raw=None;corrected=None
    try:
        check();initial=record(m0,model.state(m0));snapshots['initial']=initial;check();initial_cpu=time.process_time()
        if task['method']=='direct':
            phase='direct';observer.phase='corrector';deadline=cap
            try:_,terminations['direct']=solver.refine_observed(observer,Y,m0,maximum=6000,accepted_checkpoints=())
            except PhaseEnd:terminations['direct']='WHOLE_CHILD_CPU_CAP'
            raw=initial;corrected=best
        else:
            phase='proposal';observer.phase='auxiliary'
            if task['method'] in ('ap','dr'):
                projector=DenseFieldCurrentProjectors(model.G,model.E,model.S,Y,delta=delta,scale=6.,lower=.02,upper=20.,grouped_ports=True)
                try:
                    run=run_projection(task['method'],projector,m0,physical,steps=3000,audit_every=50,omega=1.,
                        cpu_budget_seconds=max(0,deadline-time.process_time()),checkpoint_callback=checkpoint)
                    terminations['proposal']=run.stop_reason
                except PhaseEnd:terminations['proposal']='PROPOSAL_CPU_CAP'
                projection_cost=dict(projector.costs)
            else:
                decision=noise_ball_rank(model.S,Y,sigma,out/'RECEIVER_DPCS.npz');dump(out/'DPCS_RULE.json',decision)
                common=dict(host=.1,material_box=[.02,20.],dpcs_rank=decision['rank'])
                numerical=solver.arm_config('a_plus',common,Y,model.E,[dict(mf=12,maxiter=1500)])
                numerical.update(fft_geometry=dict(coords=public['coords'],radius=float(public.get('cell_radius',public['radius'])),k=float(public['k'])),
                    regularizer=lambda m:model.reg.evaluate(m)[:2],aplus_regularizer=model.reg,
                    physical_sigma=float(np.linalg.norm(Y)),receipt_path=str(out/'G19_EXECUTION.json'),
                    budget_check=check,initialize_physical_current=True,cost_callback=lambda receipt:g19_cost.update(receipt))
                numerical=prepare_full(model,numerical);numerical['material_readout']=full_material_readout
                # Observe the existing first-stage least-squares projection; do not
                # repeat its SVD or change its output/affine admissible current family.
                old_lstsq=np.linalg.lstsq
                def observed_lstsq(A,b,*a,**kw):
                    nonlocal projection_defect
                    result=old_lstsq(A,b,*a,**kw)
                    if projection_defect is None and b.shape==initial['j'].shape:
                        remainder=A@result[0]-b;J=initial['j']+remainder;Q=model.E+model.G@J
                        projection_defect=dict(relative_current_loss=float(np.linalg.norm(remainder)/np.linalg.norm(initial['j'])),
                            constitutive_defect=float(np.linalg.norm(J-m0[:,None]*Q)/max(np.linalg.norm(J),1e-300)),
                            relative_port_change=float(np.linalg.norm(model.S@remainder)/np.linalg.norm(Y)),rank=int(result[2]),
                            singular_values=result[3].tolist(),extra_diagnostic_G_RHS=model.ns,extra_diagnostic_S_RHS=model.ns,
                            semantics='Inherited initialize_physical_current LS into unchanged original affine DPCS/Fourier space')
                    return result
                np.linalg.lstsq=observed_lstsq
                try:
                    result=solver.run_cie_fft(observer,Y,m0,numerical);terminations['proposal']=result['status']
                except PhaseEnd:
                    terminations['proposal']='PROPOSAL_CPU_CAP';g19_cost['status']='PROPOSAL_CPU_PREFIX_HANDOFF'
                finally:np.linalg.lstsq=old_lstsq
            raw=best;snapshots['raw']=raw
            phase='corrector';observer.phase='corrector';deadline=min(cap,time.process_time()+limits['corrector'])
            try:_,terminations['corrector']=solver.refine_observed(observer,Y,raw['m'],maximum=3000,accepted_checkpoints=())
            except PhaseEnd:terminations['corrector']='CORRECTOR_CPU_CAP'
            corrected=best
    except PhaseEnd:
        status='CPU_BUDGET_REACHED';terminations['outer']='WHOLE_CHILD_OR_SETUP_CAP'
    except Exception:
        status='FAILED';error=traceback.format_exc()
    snapshots.update(raw=raw,corrected=corrected,primary=best,secondary=secondary)
    # Final output copies only already paid physical candidates, never extra FW.
    for name,row in snapshots.items():
        if row is not None:
            np.savez_compressed(out/(name+'.npz'),**{k:row[k] for k in ('m','j','Q','Y_prediction')})
            dump(out/(name+'.json'),row['metadata'])
    candidates=np.asarray(materials);np.savez_compressed(out/'CANDIDATES.npz',m=candidates)
    dump(out/'HISTORY.json',history);dump(out/'PROJECTION_HISTORY.json',projection_history)
    prefixes=[]
    for point in cfg['prefixes'][task['cohort']]:
        eligible=[r for r in history if r['eligible'] and r['CPU_seconds']<=point]
        selected=min(eligible,key=lambda r:r['Phi']) if eligible else None
        prefixes.append(dict(CPU_cap=point,candidate=selected['candidate'] if selected else None))
    dump(out/'PREFIXES.json',prefixes)
    receipt=dict(status=status,task=task,case_id=task['case_id'],method=task['method'],truth_access='NONE',
        whole_child_CPU_seconds=time.process_time(),wall_seconds=time.perf_counter()-wall,CPU_cap=cap,phase_limits=limits,
        timer_scope='Absolute process CPU includes Python imports, data/model setup, observations and serialization; external controller adds full exit timing',
        eligible_states=sum(r['eligible'] for r in history),excluded_states=sum(not r['eligible'] for r in history),
        primary_candidate=best['metadata']['candidate'] if best else None,secondary_candidate=secondary['metadata']['candidate'] if secondary else None,
        raw_candidate=raw['metadata']['candidate'] if raw else None,corrected_candidate=corrected['metadata']['candidate'] if corrected else None,
        terminations=terminations,error=error,physics_cost=dict(model.cost),projection_cost=projection_cost,G19_execution=g19_cost,
        initial_current_projection=projection_defect,delta_CN99=delta,sigma_complex_RMS=sigma,imports_and_input=timings,
        affinity=proc.cpu_affinity(),priority=proc.nice(),code_version=cfg['status'],scoring_status='SEALED_NOT_SCORED')
    dump(out/'RECEIPT.json',receipt);print(json.dumps({k:receipt[k] for k in ('status','case_id','method','whole_child_CPU_seconds','eligible_states','error')}),flush=True)
    return 1 if status=='FAILED' or best is None else 0
if __name__=='__main__':raise SystemExit(main())
