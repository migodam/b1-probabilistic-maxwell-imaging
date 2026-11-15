"""R15 public-input stateful G19/twofold experiment; truth is not an input.

Budgets are absolute whole-process CPU, including imports and preparation.
Completed over-budget operations are charged but cannot yield eligible candidates.
"""
import os
for _name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
              'VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_name]='1'
import argparse
import json
from pathlib import Path
import sys
import time
import traceback


STATEFUL_METHODS=('g19','g19_s19','twofold','twofold_hop','twofold_dr','twofold_ap','twofold_full_dr','twofold_reset')


class BudgetEnd(BaseException):
    pass


def encode(value):
    import numpy as np
    if isinstance(value,np.ndarray): return value.tolist()
    if isinstance(value,np.generic): return value.item()
    if isinstance(value,Path): return str(value)
    raise TypeError(type(value).__name__)


def dump(path,value):
    Path(path).write_text(json.dumps(value,default=encode,indent=2,allow_nan=False)+'\n',encoding='utf8')


def build_config(model,public,out):
    from research.research_max_struct_r8.src import benchmark_solver
    from research.research_max_struct_r10.src.twod_adapter import noise_ball_rank
    decision=noise_ball_rank(model.S,public['Y'],float(public['sigma']),out/'RECEIVER_DPCS.npz')
    cfg=benchmark_solver.arm_config('a_plus',dict(host=.1,material_box=[.02,20.],dpcs_rank=decision['rank']),
                                   public['Y'],model.E,[dict(mf=12,maxiter=1500)])
    cfg.update(fft_geometry=dict(coords=public['coords'],radius=float(public['cell_radius']),k=float(public['k'])),
               regularizer=lambda m:model.reg.evaluate(m)[:2],aplus_regularizer=model.reg,
               sourcewise_cg=False,aplus_mode='mobius',aplus_arc_candidates=True,
               aplus_inherited_gtol=True,aplus_sweeps=2,aplus_step_tol=1e-8)
    return cfg,decision


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--legacy-repo',type=Path,required=True)
    parser.add_argument('--manifest',type=Path,required=True)
    parser.add_argument('--task-id',required=True)
    parser.add_argument('--cpu',type=int)
    args=parser.parse_args();root=args.root.resolve();wall=time.perf_counter()
    import psutil
    proc=psutil.Process()
    if args.cpu is not None and hasattr(proc,'cpu_affinity'):proc.cpu_affinity([args.cpu])
    if os.name=='nt':proc.nice(psutil.NORMAL_PRIORITY_CLASS)
    sys.path[:0]=[str(root/'src'),str(root/'runtime'),str(args.legacy_repo),str(args.legacy_repo/'vendor')]
    import numpy as np
    from scipy.optimize import minimize
    from research.research_fbg_r3_experiments.src.physics import MaxwellModel
    from twofold_solver import TwofoldG19
    manifest=json.loads(args.manifest.read_text(encoding='utf8'))
    task=next(t for t in manifest['tasks'] if t['task_id']==args.task_id)
    protocol=manifest['protocol'];cap=float(task.get('CPU_cap',protocol['CPU_cap']))
    proposal_cap=cap*float(protocol.get('proposal_fraction',.6));deadline=proposal_cap
    out=root/'RESULTS'/task['campaign']/task['task_id'];out.mkdir(parents=True,exist_ok=False)
    with np.load(root/task['input_file'],allow_pickle=False) as f:public={k:f[k].copy() for k in f.files}
    model=MaxwellModel(public);Y=public['Y'];m0=public['m_init'].copy()
    normY=float(np.linalg.norm(Y));delta=float(public['delta']);lo=float(public['lower']);hi=float(public['upper'])
    history=[];pool=[];best=None;secondary=None;phase_best={};snapshots={};native=[]
    phase='initial';next_audit=2.;solver=None;status='COMPLETE';error=None;termination={}
    checkpoints={1,32,128,512,1500};memory_peak=int(proc.memory_info().rss);final_checkpoint_attempted=False
    checkpoint_errors=[]
    if task['method'] not in STATEFUL_METHODS:raise ValueError('Unregistered R15 stateful method')
    git_identity=manifest.get('git_identity')

    def check():
        if time.process_time()>=deadline:raise BudgetEnd()

    def record(m,state,extra=None):
        nonlocal best,secondary,memory_peak
        elapsed=time.process_time();res=state['Y']-Y
        D=.5*float(np.vdot(res,res).real)/(normY*normY);R=float(model.reg.evaluate(m)[0])
        legal=bool(np.isfinite(m).all() and m.min()>=lo-1e-12 and m.max()<=hi+1e-12)
        memory_peak=max(memory_peak,int(proc.memory_info().rss))
        row=dict(candidate=len(history),phase=phase,CPU_seconds=elapsed,wall_seconds=time.perf_counter()-wall,
                 D=D,R=R,Phi=D+R,relative_residual=float(np.linalg.norm(res)/normY),
                 residual_over_delta=float(np.linalg.norm(res)/delta),noise_feasible=bool(np.linalg.norm(res)<=delta),
                 legal=legal,eligible=bool(legal and elapsed<=deadline and elapsed<=cap),deadline=deadline,
                 lower_active=int(np.count_nonzero(m<=lo+1e-8)),upper_active=int(np.count_nonzero(m>=hi-1e-8)),
                 physics_cost=dict(model.cost),extra=extra)
        candidate=dict(m=np.array(m,copy=True),j=state['j'].copy(),Q=state['Q'].copy(),Y_prediction=state['Y'].copy(),metadata=row)
        history.append(row);pool.append(candidate['m'])
        if row['eligible']:
            if best is None or row['Phi']<best['metadata']['Phi']:best=candidate
            if phase not in phase_best or row['Phi']<phase_best[phase]['metadata']['Phi']:phase_best[phase]=candidate
            if row['noise_feasible'] and (secondary is None or (R,D+R)<(secondary['metadata']['R'],secondary['metadata']['Phi'])):secondary=candidate
        return candidate

    def physical(m,extra=None):
        check();state=model.state(m);candidate=record(m,state,extra);check();return candidate

    try:
        snapshots['initial']=physical(m0)
        config,decision=build_config(model,public,out)
        config['budget_check']=check
        config.update(protocol.get('r15',{}))
        config['r15_protocol']=dict(protocol.get('r15',{}))
        config['r15_data_delta']=delta
        config['stages']=[dict(stage) for stage in protocol['stages']]
        config['material_box']=[lo,hi]
        audit_interval=float(protocol.get('physical_audit_CPU_interval',2.))
        if audit_interval!=2.:raise ValueError('Frozen original2-process-CPU-second audit interval required')
        config_record={k:v for k,v in config.items() if not callable(v) and k not in ('aplus_regularizer',)}
        dump(out/'ACTUAL_CONFIG.json',dict(task=task,numerical=config_record,dpcs=decision,truth_access='NONE',
             runtime_contract='Parent-supplied TwofoldG19; same original physical record/corrector, periodic cycle audits for every stateful method',
             git_identity=git_identity))
        phase='proposal'
        def conditional_audit(m,metadata):
            return physical(m,dict(source='paid conditional DR/AP shadow or matched cycle opportunity',cycle_metadata=metadata))
        solver=TwofoldG19(model,Y,m0,config,task['method'],out,conditional_audit)
        try:
            while True:
                check();row=solver.step()
                arrays={k:v for k,v in row.items() if isinstance(v,np.ndarray)}
                scalar={k:v for k,v in row.items() if k not in arrays}
                scalar['CPU_seconds']=time.process_time();native.append(scalar)
                iteration=int(row.get('global_iteration',len(native)))
                if iteration in checkpoints or row.get('stage_end',False):
                    if iteration in checkpoints:solver.checkpoint(out/f'CHECKPOINT_{iteration:05d}.npz')
                    if row.get('stage_end',False):solver.checkpoint(out/f'CHECKPOINT_STAGE_{int(row.get("stage",-1)):02d}_{iteration:05d}.npz')
                if row.get('stage_end',False) or time.process_time()>=next_audit:
                    physical(solver.m,dict(iteration=iteration,stage=row.get('stage'),source='complete material update'))
                    next_audit=(int(time.process_time()/2)+1)*2.
                check()
        except StopIteration:termination['proposal']='SCHEDULE_COMPLETE'
        except BudgetEnd:termination['proposal']='WHOLE_PROCESS_CPU_PROPOSAL'
        snapshots['raw']=best  # Original C2 proposal incumbent includes initial.
        if solver is not None:
            # Export is charged and occurs before the corrector's shared cap.
            final_checkpoint_attempted=True
            solver.checkpoint(out/'CHECKPOINT_FINAL.npz')
            np.savez_compressed(out/'JOINT_FINAL.npz',m=solver.m,J=solver.J,Q=solver.Q)
        deadline=cap;phase='corrector'
        start=best['m'].copy() if best is not None else m0.copy()
        def objective(m):
            check();value,gradient=model.data_value_gradient(m,Y)
            record(m,model.state(m));check()
            return value*1e6,gradient*1e6
        try:
            opt=minimize(objective,start,jac=True,method='L-BFGS-B',bounds=[(lo,hi)]*model.n,
                         options=dict(maxiter=3000,maxfun=12000,ftol=1e-12,gtol=.01,maxls=30))
            termination['corrector']=dict(message=str(opt.message),iterations=int(opt.nit),evaluations=int(opt.nfev))
        except BudgetEnd:termination['corrector']='WHOLE_PROCESS_CPU'
        snapshots['corrected']=phase_best.get('corrector',best)
    except BudgetEnd:
        termination['setup']='WHOLE_PROCESS_CPU'
        status='INCOMPLETE_SETUP'
    except Exception as exc:
        status='FAILED';error=repr(exc);traceback.print_exc()
    finally:
        if solver is not None and not final_checkpoint_attempted:
            final_checkpoint_attempted=True
            try:solver.checkpoint(out/'CHECKPOINT_FINAL.npz')
            except Exception as checkpoint_exc:checkpoint_errors.append(repr(checkpoint_exc))
        snapshots['primary']=best;snapshots['secondary']=secondary
        for name,row in snapshots.items():
            if row is not None:
                np.savez_compressed(out/(name+'.npz'),**{k:v for k,v in row.items() if k!='metadata'})
                dump(out/(name+'.json'),row['metadata'])
        np.savez_compressed(out/'CANDIDATES.npz',m=np.array(pool) if pool else np.empty((0,model.n)))
        dump(out/'HISTORY.json',history);dump(out/'NATIVE_HISTORY.json',native)
        dump(out/'RECEIPT.json',dict(status=status,error=error,task=task,termination_reason=termination,
             CPU_seconds=time.process_time(),whole_process_CPU_seconds=time.process_time(),wall_seconds=time.perf_counter()-wall,CPU_cap=cap,proposal_CPU_cap=proposal_cap,
             physical_cost=dict(model.cost),auxiliary_cost=dict(solver.cost) if solver is not None else {},
             fixed_projection_cost=dict(solver.projector.costs) if solver is not None and solver.projector is not None else {},
             peak_sampled_RSS_bytes=memory_peak,RSS_scope='sampled at physical audits, not exact peak',
             candidates=len(history),eligible_candidates=sum(r['eligible'] for r in history),
             primary_candidate=best['metadata']['candidate'] if best is not None else None,
             secondary_candidate=secondary['metadata']['candidate'] if secondary is not None else None,
             truth_access='NONE',online_selection_rule='min eligible paid original Phi; secondary minR in originalCN99',
             physical_convention='exp(-i omega t); scalar TM disk pulse VIE; real chi=epsilon_r-1',
             material_dof=model.n,state_representation='inherited native Fourier APCS + fixed receiver DPCS',
             git_identity=git_identity,git_identity_null_reason=None if git_identity else 'Not supplied at manifest freeze',
             checkpoint_errors=checkpoint_errors,source_record_logic='Inherited R14 run_experiment record/corrector/export; C2 view fallbacks and CANDIDATES naming',
             cycle_snapshots='First4 before/after handled by parent TwofoldG19; fixed background cache once, paid same clock',
             unknown_fields=dict(currency='no metered bill',GPU_seconds='CPU workload',exact_memory_peak='not instrumented')))
    return 0 if status=='COMPLETE' else 1


if __name__=='__main__':sys.exit(main())
