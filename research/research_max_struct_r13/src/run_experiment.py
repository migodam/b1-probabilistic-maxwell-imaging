"""Public-input R13 runner: inherited Maxwell physics and paid material search."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import argparse, json, sys, time, traceback
from pathlib import Path


class BudgetEnd(BaseException):
    pass


def dump(path, value):
    import numpy as np
    def encode(x):
        if isinstance(x, np.ndarray): return x.tolist()
        if isinstance(x, np.generic): return x.item()
        if isinstance(x, Path): return str(x)
        raise TypeError(type(x).__name__)
    Path(path).write_text(json.dumps(value, indent=2, default=encode, allow_nan=False)+'\n')


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--legacy-repo',type=Path,required=True)
    parser.add_argument('--manifest',type=Path,required=True)
    parser.add_argument('--task-id',required=True)
    parser.add_argument('--cpu',type=int)
    args=parser.parse_args(); root=args.root.resolve(); wall=time.perf_counter()
    import psutil
    proc=psutil.Process()
    if args.cpu is not None and hasattr(proc, 'cpu_affinity'): proc.cpu_affinity([args.cpu])
    if os.name=='nt': proc.nice(psutil.NORMAL_PRIORITY_CLASS)
    sys.path[:0]=[str(args.legacy_repo),str(root/'runtime'),str(root/'src')]
    import numpy as np
    from scipy.stats import chi2
    from scipy.optimize import minimize
    from research.research_fbg_r3_experiments.src.physics import MaxwellModel
    from research.research_max_struct_r8.src import benchmark_solver as solver
    from research.research_max_struct_r10.src.twod_adapter import noise_ball_rank
    from field_current import DenseFieldCurrentProjectors
    from material_events import norm_pair, event_start, finite_threshold, angular_branch_margin
    manifest=json.loads(args.manifest.read_text())
    task=next(t for t in manifest['tasks'] if t['task_id']==args.task_id)
    cfg=manifest['protocol']; cap=float(task.get('CPU_cap',cfg['CPU_cap']))
    proposal_fraction=float(cfg.get('proposal_fraction',.6))
    deadline=cap if task['method']=='direct' else cap*proposal_fraction
    out=root/'RESULTS'/task['campaign']/task['task_id'];out.mkdir(parents=True,exist_ok=False)
    with np.load(root/task['input_file'],allow_pickle=False) as f: public={k:f[k].copy() for k in f.files}
    m0=public.pop('m_init'); Y=public['Y']; sigma=float(public['sigma'])
    components=int(public.get('components',1)); lower=float(public.get('lower',.02)); upper=float(public.get('upper',20.))
    kind=str(public.get('kind','synthetic'))
    if kind=='fresnel':
        from vector_adapter import SourceSpecificPublicModel
        model=SourceSpecificPublicModel(public)
    elif components==3:
        from vector_adapter import VectorPublicModel
        model=VectorPublicModel(public)
    else:
        model=MaxwellModel(public)
    delta=float(public.get('delta',sigma*np.sqrt(chi2.ppf(.99,2*Y.size)/2)))
    normY2=float(np.vdot(Y,Y).real)
    history=[]; materials=[]; best=None; secondary=None; phase_best={}; phase='initial'
    projection_history=[]; event_history=[]; auxiliary_cost={}; auxiliary_receipt={}
    P=None; snapshots={}; stop={}; error=None; status='COMPLETE'; final_states={}

    def check():
        if time.process_time() >= deadline: raise BudgetEnd()

    def record(m,state,value=None,gradient=None,extra=None):
        nonlocal best,secondary
        rr=state['Y']-Y
        D=model.data_objective(rr,Y) if hasattr(model,'data_objective') else .5*float(np.vdot(rr,rr).real)/normY2
        R=float(model.reg.evaluate(m)[0]); t=time.process_time()
        legal=bool(np.isfinite(m).all() and np.min(m)>=lower-1e-12 and np.max(m)<=upper+1e-12)
        meta=dict(candidate=len(history),CPU_seconds=t,wall_seconds=time.perf_counter()-wall,
                  phase=phase,Phi=D+R,D=D,R=R,relative_residual=float(np.linalg.norm(rr)/np.sqrt(normY2)),
                  residual_over_delta=float(np.linalg.norm(rr)/delta),noise_feasible=bool(np.linalg.norm(rr)<=delta),
                  legal=legal,eligible=bool(t<=deadline and t<=cap),deadline=deadline,
                  lower_active=int(np.sum(m<=lower+1e-8)),upper_active=int(np.sum(m>=upper-1e-8)),
                  physics_cost=dict(model.cost),auxiliary=extra)
        row=dict(m=np.asarray(m).copy(),j=state['j'].copy(),Q=state['Q'].copy(),Y_prediction=state['Y'].copy(),metadata=meta)
        history.append(meta);materials.append(row['m'])
        if meta['eligible'] and legal:
            if best is None or meta['Phi']<best['metadata']['Phi']: best=row
            if phase not in phase_best or meta['Phi']<phase_best[phase]['metadata']['Phi']: phase_best[phase]=row
            if meta['noise_feasible'] and (secondary is None or (R,D+R)<(secondary['metadata']['R'],secondary['metadata']['Phi'])): secondary=row
        return row

    def physical(m,extra=None):
        check();state=model.state(m);row=record(m,state,extra=extra);check();return row

    def g19_audit(m,state,value,metadata):
        record(m,state,value,extra=metadata);check()

    def correct(start):
        def fun(m):
            check();value,grad=model.data_value_gradient(m,Y)
            record(m,model.state(m),value,grad);check()
            return value*1e6,grad*1e6
        result=minimize(fun,start,method='L-BFGS-B',jac=True,bounds=[(lower,upper)]*model.n,
                        options=dict(maxiter=6000 if task['method']=='direct' else 3000,
                                     maxfun=24000 if task['method']=='direct' else 12000,
                                     ftol=1e-12,gtol=.01,maxls=30))
        return dict(message=str(result.message),iterations=int(result.nit),evaluations=int(result.nfev))

    try:
        initial=physical(m0);snapshots['initial']=initial
        if task['method']=='direct':
            phase='direct'
            try: stop['direct']=correct(m0)
            except BudgetEnd: stop['direct']='WHOLE_PROCESS_CPU'
            snapshots['raw']=initial;snapshots['corrected']=best
        else:
            phase='proposal'
            try:
                method=task['method']
                if (kind=='fresnel' or components==3) and method in ('cie','dual_cie','global_dual','a_plus','g19','g19_full_apcs','g19_full_current'):
                    from vector_adapter import run_vector_auxiliary,run_fresnel_auxiliary
                    adapter=run_fresnel_auxiliary if kind=='fresnel' else run_vector_auxiliary
                    options=dict(audit_every=cfg['audit_every'])
                    if kind=='fresnel':options['cost_callback']=auxiliary_receipt.update
                    result=adapter(model,Y,m0,method,check,g19_audit,**options)
                    auxiliary_receipt.update(result)
                    stop['proposal']=result.get('status')
                elif method in ('cie','dual_cie','global_dual','a_plus','g19','g19_full_apcs','g19_full_current'):
                    from full_current import run_full_current
                    if kind=='classic':
                        from physics_inputs import classic_auxiliary_config
                        dec=dict(rank=int(public['baseline_dpcs_rank']),rule='inherited_paper_rank_not_CN99_selection')
                        numerical=classic_auxiliary_config(public,method,model,Y,out)
                    else:
                        dec=noise_ball_rank(model.S,Y,sigma,out/'RECEIVER_DPCS.npz')
                        common=dict(host=.1,material_box=[lower,upper],dpcs_rank=dec['rank'])
                        name={'cie':'cie_single','dual_cie':'cie_dual','global_dual':'global_dual'}.get(method,'a_plus')
                        numerical=solver.arm_config(name,common,Y,model.E,[dict(mf=12,maxiter=1500)])
                    dump(out/'DPCS_RULE.json',dec)
                    numerical.update(fft_geometry=dict(coords=public['coords'],radius=float(public.get('cell_radius',public['radius'])),k=float(public['k'])),
                        regularizer=lambda m:model.reg.evaluate(m)[:2],aplus_regularizer=model.reg,
                        physical_sigma=float(np.linalg.norm(Y)),receipt_path=str(out/'AUXILIARY_EXECUTION.json'),
                        cost_callback=lambda r:auxiliary_receipt.update(r),
                        # R12 G19 explicitly warms physical current. Historical
                        # CIE/Dual/Global and original R7 A+ leave this false.
                        initialize_physical_current=method.startswith('g19'))
                    if method=='a_plus':
                        from research.research_max_struct_r7.src.aplus_readout import mobius_material_readout
                        numerical['material_readout']=mobius_material_readout
                    mode='full_apcs' if method=='g19_full_apcs' else 'full_current' if method=='g19_full_current' else 'restricted'
                    result=run_full_current(model,Y,m0,numerical,mode=mode,budget_check=check,
                        audit_callback=g19_audit,audit_every=cfg['audit_every'],g19=method.startswith('g19'))
                    stop['proposal']=result.get('status')
                else:
                    scale=float(cfg.get('scale',6.))
                    P=DenseFieldCurrentProjectors(model.G,model.E,getattr(model,'projector_S',model.S),Y,delta,scale=scale,
                        lower=lower,upper=upper,components=components,grouped_ports=True)
                    check(); x,q=initial['j']/scale,initial['Q'].copy()
                    rng=np.random.default_rng(int(task.get('seed',13013)))
                    prep=int(cfg.get('AP_preparation',128)) if method.startswith(('event_','normal_','random_','ap_dr','reset_','online_')) else 0
                    for k in range(prep):
                        check();cx,cq,cm=P.material(x,q)
                        if method.startswith('online_'):
                            ax,aq=P.affine_ball(2*cx-x,2*cq-q)
                            x,q=x+ax-cx,q+aq-cq
                        else:
                            x,q=P.affine_ball(cx,cq)
                        if (k+1)%cfg['audit_every']==0: physical(cm,dict(preparation_iteration=k+1))
                    if prep:
                        if method.startswith('online_'):
                            # Form a single A point from the CURRENT DR shadow,
                            # rather than pay a separate AP convergence phase.
                            cx,cq,_=P.material(x,q);x,q=P.affine_ball(cx,cq)
                        mode='random_integer' if method=='random_integer' else 'integer' if 'integer' in method else 'double' if method.startswith('normal_') else 'random' if method.startswith('random_') or method=='online_random' else 'fractional'
                        if method=='ap_dr':
                            event_history.append(dict(status='AP_PREPARATION_THEN_ORDINARY_DR',preparation=prep))
                        else:
                            x,q,event,th=event_start(P,x,q,mode,rng)
                            event['preparation']=prep;event['method']=method
                            event_history.append(event)
                            np.savez_compressed(out/'EVENT_GEOMETRY.npz',**{k:v for k,v in th.items() if isinstance(v,np.ndarray)})
                            if method.startswith('reset_'):
                                _,_,cm=P.material(x,q);row=physical(cm,dict(event_memory_reset=True))
                                x,q=row['j']/scale,row['Q'].copy()
                    prev_x=None;prev_q=None;prev_f=None;prev_r=None
                    for k in range(int(cfg['steps'])):
                        check();cx,cq,cm=P.material(x,q)
                        if method in ('ap','anderson'):
                            ax,aq=P.affine_ball(cx,cq)
                            if method=='anderson':
                                # Depth-one Anderson on the AP fixed-point map.
                                # One projection per step; coefficients clipped to
                                # [-1,1] to bound the affine extrapolation norm.
                                r=np.r_[ (ax-x).ravel(),(aq-q).ravel() ]
                                if prev_r is not None:
                                    dr=r-prev_r; den=float(np.vdot(dr,dr).real)
                                    alpha=float(np.clip(np.vdot(dr,r).real/den,-1.,1.)) if den>0 else 0.
                                    nx=ax-alpha*(ax-prev_x);nq=aq-alpha*(aq-prev_q)
                                else: nx,nq=ax,aq
                                prev_x,prev_q,prev_r=ax.copy(),aq.copy(),r.copy()
                                x,q=nx,nq
                            else: x,q=ax,aq
                        else:
                            ax,aq=P.affine_ball(2*cx-x,2*cq-q)
                            x,q=x+ax-cx,q+aq-cq
                        if (k+1)%cfg['audit_every']==0 or k+1==cfg['steps']:
                            gap=norm_pair(ax-cx,aq-cq)
                            physical(cm,dict(iteration=k+1,shadow_gap=gap))
                            projection_history.append(dict(iteration=k+1,CPU_seconds=time.process_time(),shadow_gap=gap,cost=dict(P.costs)))
                    final_states=dict(X=x,Q=q)
                    stop['proposal']='FIXED_STEPS_COMPLETE'
            except BudgetEnd: stop['proposal']='PROPOSAL_PROCESS_CPU'
            snapshots['raw']=best
            if best is not None:
                phase='corrector'
                deadline=min(cap,time.process_time()+cap*(1-proposal_fraction))
                try: stop['corrector']=correct(best['m'])
                except BudgetEnd: stop['corrector']='CORRECTOR_PROCESS_CPU'
                snapshots['corrected']=phase_best.get('corrector',best)
    except BudgetEnd:
        status='CPU_CAP';stop['outer']='IMPORT_SETUP_OR_INITIAL_CAP'
    except Exception:
        status='FAILED';error=traceback.format_exc()
    snapshots['primary']=best;snapshots['secondary']=secondary
    for name,row in snapshots.items():
        if row is not None:
            np.savez_compressed(out/(name+'.npz'),**{k:row[k] for k in ('m','j','Q','Y_prediction')})
            dump(out/(name+'.json'),row['metadata'])
    np.savez_compressed(out/'CANDIDATES.npz',m=np.asarray(materials))
    if final_states: np.savez_compressed(out/'FINAL_AUXILIARY_STATE.npz',**final_states)
    dump(out/'HISTORY.json',history);dump(out/'PROJECTION_HISTORY.json',projection_history)
    dump(out/'EVENTS.json',event_history)
    prefixes=[]
    for point in cfg.get('prefixes',[15,30,60,90,120]):
        rows=[h for h in history if h['eligible'] and h['CPU_seconds']<=point]
        selected=min(rows,key=lambda h:h['Phi']) if rows else None
        prefixes.append(dict(CPU_cap=point,candidate=selected['candidate'] if selected else None))
    dump(out/'PREFIXES.json',prefixes)
    receipt=dict(status=status,task=task,protocol=cfg,truth_access='NONE',error=error,termination=stop,
        whole_process_CPU_seconds=time.process_time(),wall_seconds=time.perf_counter()-wall,
        physical_cost=dict(model.cost),projector_cost=dict(P.costs) if P is not None else None,
        auxiliary_receipt=auxiliary_receipt,raw_candidate=snapshots.get('raw',{}).get('metadata',{}).get('candidate'),
        primary_candidate=best['metadata']['candidate'] if best else None,
        secondary_candidate=secondary['metadata']['candidate'] if secondary else None,
        eligible_candidates=sum(h['eligible'] for h in history),overrun_candidates=sum(not h['eligible'] for h in history),
        components=components,material_DOF=model.n,state_real_DOF=4*model.G.shape[0]*model.ns,
        delta=delta,sigma_complex_RMS=sigma,
        noise_qualification=str(public.get('sigma_role','iid circular complex RMS; CN99 raw ball')),
        prior=dict(name=getattr(model.reg,'name','inherited_scalar_regularizer'),
            smooth=getattr(model.reg,'smooth',None),tv=getattr(model.reg,'tv',None),epsilon=getattr(model.reg,'epsilon',None)),
        affinity=proc.cpu_affinity() if hasattr(proc,'cpu_affinity') else None,
        endpoint_RSS_bytes=proc.memory_info().rss,
        memory_scope='end-of-process resident snapshot, not measured peak',
        timer_scope='absolute process CPU including imports/setup/proposal/audits/corrector/export; external exit CPU separately')
    dump(out/'RECEIPT.json',receipt)
    print(json.dumps(dict(task_id=args.task_id,status=status,CPU=receipt['whole_process_CPU_seconds'],candidates=len(history),error=error)),flush=True)
    return 1 if status=='FAILED' or best is None else 0


if __name__=='__main__': raise SystemExit(main())
