"""Authorized public0000 integration check: fixed3/17 steps, no GT/main runs."""
from __future__ import annotations
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
import argparse,inspect,json,sys,time,traceback,importlib.util
from pathlib import Path
from copy import deepcopy
_IMPORT_WALL=time.perf_counter()


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--legacy-repo',type=Path,required=True);p.add_argument('--public',type=Path);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    root=a.root.resolve();W=a.legacy_repo.resolve();out=a.out.resolve()
    if out.exists():raise FileExistsError(out)
    out.mkdir(parents=True);sys.path[:0]=[str(root/'src'),str(root/'runtime'),str(W),str(W/'vendor')]
    import numpy as np
    from research.research_fbg_r3_experiments.src.physics import MaxwellModel
    from run_stateful import build_config
    from twofold_solver import TwofoldG19
    import frozen_full_current as original
    public_path=(a.public or root/'data_public/independent24/r14_new_0000.npz').resolve()
    if public_path.name!='r14_new_0000.npz' or any(x in ('sealed_truth','truth','private') for x in public_path.parts):raise ValueError('Only authorized public0000')
    with np.load(public_path,allow_pickle=False) as z:public={k:z[k].copy() for k in z.files}
    Y=public['Y'];m0=public['m_init'];executions=[];problems=[];results={};setup_config=None
    source=root.parent/'research_max_struct_r14/configs/independent24_v1.json';protocol=json.loads(source.read_text())['protocol']
    r15=dict(tauS=.1,tauG=.1,numerical_rank_rtol=1e-12,period=16,omega=1.,scale=6.,DRscale=6.,AP_preparation=0,onehop_resolvent_corrections=1)
    def write(path,obj):
        def encode(x):
            if isinstance(x,np.generic):return x.item()
            if isinstance(x,np.ndarray):return x.tolist()
            if isinstance(x,Path):return str(x)
            raise TypeError(type(x).__name__)
        with Path(path).open('x') as f:json.dump(obj,f,default=encode,indent=2,allow_nan=False);f.write('\n')
    def delta(after,before):return {k:after.get(k,0)-before.get(k,0) for k in set(after)|set(before) if after.get(k,0)!=before.get(k,0)}
    def metrics(x,y):
        if x is None or y is None:return dict(both_none=x is None and y is None,absolute=None,relative=None)
        x,y=np.asarray(x),np.asarray(y);d=float(np.linalg.norm(x-y));return dict(absolute=d,relative=d/max(float(np.linalg.norm(y)),1e-300),max_absolute=float(np.max(abs(x-y))) if x.size else 0.,array_equal=bool(np.array_equal(x,y)))
    def snapshot(s):return {k:None if getattr(s,k) is None else np.array(getattr(s,k),copy=True) for k in ['m','J','Q','A','old_grad','old_direction','memory_X','memory_Q']}
    def configuration(model,folder):
        model.state(m0) # Same paid initial physical warm-state opportunity.
        cfg,decision=build_config(model,public,folder)
        cfg.update(initialize_physical_current=True,stages=deepcopy(protocol['stages']),r15_protocol=deepcopy(r15),r15_data_delta=float(public['delta']),receipt_path=str(folder/'NATIVE_EXECUTION.json'))
        return cfg,decision
    def audit_for(model,items):
        def audit(m,metadata):
            state=model.state(m);res=state['Y']-Y;D=.5*float(np.vdot(res,res).real)/float(np.vdot(Y,Y).real);R=float(model.reg.evaluate(m)[0])
            items.append(dict(Phi=D+R,D=D,R=R,metadata=metadata,physical_cost=dict(model.cost),CPU_seconds=time.process_time()))
        return audit
    class Attempt:
        def __init__(self,name):self.name=name;self.folder=out/name;self.folder.mkdir();self.model=None;self.solver=None;self.extra={}
        def __enter__(self):self.cpu=time.process_time();self.wall=time.perf_counter();return self
        def __exit__(self,kind,error,tb):
            rec=dict(attempt=self.name,status='FAILED' if error else 'COMPLETE',error=repr(error) if error else None,CPU_seconds=time.process_time()-self.cpu,wall_seconds=time.perf_counter()-self.wall,
                physical_cost=dict(self.model.cost) if self.model else None,auxiliary_cost=dict(self.solver.cost) if self.solver else None,
                projector_cost=dict(self.solver.projector.costs) if self.solver and self.solver.projector else None,extra=self.extra,truth_access='NONE')
            if error:rec['traceback']=''.join(traceback.format_exception(kind,error,tb));problems.append(dict(attempt=self.name,error=repr(error)))
            executions.append(rec);write(self.folder/'RECEIPT.json',rec);return error is not None
    # Literal frozen baseline source budget-line contract, not another baseline.
    frozen=root/'runtime/frozen_C2_runner.py';equal=root/'runtime/equal_budget_C2_runner.py'
    literal_match=equal.read_text()==frozen.read_text().replace('deadline=min(cap,time.process_time()+cap*(1-proposal_fraction))','deadline=cap')
    if not literal_match:raise ValueError('C2 one-line budget provenance failed')
    native_rows=[];native_receipt={}
    with Attempt('native_g19_3steps') as t:
        t.model=MaxwellModel(public);cfg,decision=configuration(t.model,t.folder)
        cfg['max_iterations']=3;cfg['cost_callback']=native_receipt.update
        def observe(row):
            frame=inspect.currentframe().f_back
            while frame is not None and frame.f_code.co_name!='_run_cie':frame=frame.f_back
            if frame is None:raise RuntimeError('Original R8 callback loop frame missing')
            loc=frame.f_locals
            native_rows.append(dict(m=row['m'].copy(),J=row['J'].copy(),Q=loc['Q'].copy(),A=row['coefficients'].copy(),
                old_grad=loc['old_grad'].copy(),old_direction=loc['old_direction'].copy(),PR=np.array(loc['pr']),alpha=np.array(loc['alpha']),history=deepcopy(loc['history'][-1])))
            del frame
        cfg['iteration_callback']=observe
        original.run_full_current(t.model,Y,m0,cfg,mode='restricted',g19=True,audit_every=0)
        assert len(native_rows)==3
        t.extra.update(native_receipt=native_receipt,diagnostic_limit='max_iterations=3; only global step stop added, original betas/tolerances/readout unchanged',DPCS=decision)
    with Attempt('stateful_g19_3steps') as t:
        t.model=MaxwellModel(public);cfg,decision=configuration(t.model,t.folder);cfg['max_iterations']=3
        t.solver=TwofoldG19(t.model,Y,m0,cfg,'g19',t.folder,None);comparisons=[]
        for i in range(3):
            row=t.solver.step();state=snapshot(t.solver)
            comparisons.append(dict(step=i+1,arrays={k:metrics(state[k],native_rows[i][k]) for k in ['m','J','Q','A','old_grad','old_direction']},PR=metrics(row['PR'],native_rows[i]['PR']),alpha=metrics(row['alpha'],native_rows[i]['alpha'])))
        results['g19_3step_comparison']=comparisons;t.extra.update(comparison=results['g19_3step_comparison'],DPCS=decision)
    for method in ('twofold_hop','twofold_dr','twofold_ap','twofold_full_dr','twofold_reset'):
        with Attempt(method+'_17steps') as t:
            t.model=MaxwellModel(public);cfg,decision=configuration(t.model,t.folder);audits=[]
            t.solver=TwofoldG19(t.model,Y,m0,cfg,method,t.folder,audit_for(t.model,audits));seen=[];last=None
            for i in range(17):
                last=t.solver.step()
                seen.append(dict(iteration=t.solver.iterations,cycle_event=bool(last['cycle_event']),cycle_metadata=last.get('cycle_metadata'),
                    finite=bool(all(np.isfinite(v).all() for v in (t.solver.m,t.solver.J,t.solver.Q,t.solver.A))),
                    material_legal=bool(np.all(t.solver.m>=t.solver.lo-1e-12) and np.all(t.solver.m<=t.solver.hi+1e-12))))
            checkpoint=t.folder/'CHECKPOINT_17.npz';t.solver.checkpoint(checkpoint)
            cycles=[x for x in seen if x['cycle_event']]
            if not cycles:raise AssertionError('Original period16 did not produce a cycle')
            # Every arm checkpoint is continued one step. DR additionally restores
            # the actual pending-before-cycle checkpoint and replays the event.
            original_before=snapshot(t.solver);statecost=dict(t.solver.cost)
            restore_dir=t.folder/'restored_after17';restore_dir.mkdir();restore_audits=[]
            restored=TwofoldG19.from_checkpoint(t.model,Y,cfg,checkpoint,out=restore_dir,audit_callback=audit_for(t.model,restore_audits))
            restored_before=snapshot(restored);rowa=t.solver.step();rowb=restored.step()
            comparison=dict(before={k:metrics(restored_before[k],original_before[k]) for k in original_before},after={k:metrics(snapshot(restored)[k],snapshot(t.solver)[k]) for k in original_before},PR=metrics(rowb['PR'],rowa['PR']),alpha=metrics(rowb['alpha'],rowa['alpha']))
            replay=None
            if method=='twofold_dr':
                beforecycle=t.folder/'CYCLE_001_BEFORE.npz';replay_dir=t.folder/'restored_pending_cycle';replay_dir.mkdir();replay_audits=[]
                replay_solver=TwofoldG19.from_checkpoint(t.model,Y,cfg,beforecycle,out=replay_dir,audit_callback=audit_for(t.model,replay_audits))
                replay_row=replay_solver.step()
                # Compare to frozen AFTER checkpoint, avoiding another solve.
                aftercycle=t.folder/'CYCLE_001_AFTER.npz';after_dir=t.folder/'restored_after_cycle_reference';after_dir.mkdir()
                after_solver=TwofoldG19.from_checkpoint(t.model,Y,cfg,aftercycle,out=after_dir,audit_callback=None)
                replay=dict(cycle_event=replay_row['cycle_event'],comparison={k:metrics(snapshot(replay_solver)[k],snapshot(after_solver)[k]) for k in original_before},replay_auxiliary_cost=dict(replay_solver.cost),replay_projector_cost=dict(replay_solver.projector.costs),reference_auxiliary_cost=dict(after_solver.cost),physical_audits=replay_audits)
            t.extra.update(steps=seen,cycle_events=len(cycles),cycle_audits=audits,restore_one_step=comparison,pending_cycle_restore=replay,
                cost_at17=statecost,restore_auxiliary_cost=dict(restored.cost),continuation_native_cost_increment=delta(dict(t.solver.cost),statecost),
                checkpoint=str(checkpoint),DPCS=decision)
            results[method]=t.extra
    totals=dict(status='COMPLETE_WITH_FAILURES' if problems else 'COMPLETE',public=str(public_path),truth_access='NONE',main_experiment=False,attempts=len(executions),failures=problems,
        diagnostic_step_limits='Native/stateful G19 exactly3; periodic five methods17 plus one continuation each; DR pending-cycle replay once. Original period16/stages/betas/gtol unchanged; no scoring/corrector.',
        literal_C2_budget_only_source_match=literal_match,source_paths=dict(core=str(Path(sys.modules['dualmetric_core'].__file__)),original_full_current=original.__file__,solver=str(root/'src/twofold_solver.py')),
        recorded_attempt_CPU_seconds=sum(x['CPU_seconds'] for x in executions),recorded_attempt_wall_seconds=sum(x['wall_seconds'] for x in executions),whole_process_CPU_seconds=time.process_time(),wall_since_imports=time.perf_counter()-_IMPORT_WALL,
        costs_scope='Attempt CPU includes setup, checkpoints, all failed checks/continuations/replay; cumulative restored counters are nested historical counts, never sum them as fresh work. Model physical cost covers extra physical warm/audit calls; projector background counters separate. Native final stage audit after diagnostic stop is paid.',
        physical_LU_total=sum(x['physical_cost'].get('LU',0) for x in executions if x['physical_cost']),physical_actual_complex_RHS_total=sum(x['physical_cost'].get('actual_complex_RHS',0) for x in executions if x['physical_cost']),results=results)
    write(out/'RESULTS.json',totals)
    lines=['# R15 public integration check','',f"Status: {totals['status']}; attempts {len(executions)}, failures {len(problems)}. Only public r14_new_0000; no GT, main experiment, scoring or corrector.", '', 'Native G19 uses original frozen R13/R8 program with diagnostic max_iterations3; new stateful uses the same config and three step calls. Callback reads existing R8 frame arrays without replacing formulas. Physical warm state and original readout are preserved. Cost/trajectory equivalence is not assumed.','', 'Per-attempt receipts preserve every failure and full cumulative physical/auxiliary/projector counters. Restored historical counters are not fresh-work totals; attempt process CPU includes all continuation/replay and serialization work.','',f"Recorded attempt CPU {totals['recorded_attempt_CPU_seconds']:.6f}s; whole process CPU {totals['whole_process_CPU_seconds']:.6f}s; physical LU {totals['physical_LU_total']}, actual complex RHS {totals['physical_actual_complex_RHS_total']}."]
    if problems:lines += ['', 'Failures:']+['- '+x['attempt']+': '+x['error'] for x in problems]
    for name,value in results.items():
        if name=='g19_3step_comparison':
            lines+=['','## Native vs stateful G19']
            for v in value:lines.append('Step '+str(v['step'])+': '+', '.join(k+' rel='+str(d['relative']) for k,d in v['arrays'].items())+'; PR rel='+str(v['PR']['relative']))
        else:
            after=value['restore_one_step']['after'];lines+=['',f"## {name}",f"17 steps, cycle events={value['cycle_events']}; all finite={all(x['finite'] for x in value['steps'])}, all material legal={all(x['material_legal'] for x in value['steps'])}.",'Checkpoint next-step differences: '+', '.join(k+' abs='+str(v['absolute']) for k,v in after.items())+'.']
            if value['pending_cycle_restore']:lines.append('Pending-before-cycle restore differences: '+', '.join(k+' abs='+str(v['absolute']) for k,v in value['pending_cycle_restore']['comparison'].items())+'.')
    lines+=['','This checks software integration only. Parent owns interpretation and numerical scope. No parameter was retuned to force equality.']
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n');print(json.dumps({k:v for k,v in totals.items() if k!='results'}))
    return 0 if not problems else 1

if __name__=='__main__':sys.exit(main())
