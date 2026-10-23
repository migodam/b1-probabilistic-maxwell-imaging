"""Finite serial B2 tail after released P1A/P1B queues; no new budget.

Selection inspects completed cavity calibration distances only. It never reads
P1 truth/observations or target distances to nominate a material. Import is inert.
"""
from __future__ import annotations

import argparse
import fcntl
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
CHECKPOINT=ROOT/'B2_FINISH_QUEUE_CHECKPOINT.json'
FIRST='configs/FIRST_JOINT_NOMINATION.json'
SECOND='configs/SECOND_CAVITY_NOMINATION.json'
MAIN_FINISHED={'P1A_REFERENCE_FINISHED','P1A_REFERENCE_FINISHED_REVIEW_GATE_BEFORE_P1B'}
P1B_FINISHED={'NUMERICAL_WORK_FINISHED','P1B_NUMERICAL_WORK_FINISHED_PENDING_SCIENTIFIC_REVIEW',
              'NOT_RUN_P1A_REFERENCE_GATE_CLOSED'}
JOBS=(('b2-second-cavity-independent-validation','wr.endpoint_validation',SECOND,1080),
      ('b2-first-joint-task-argmax-grid-check','wr.argmax_validation',FIRST,40),
      ('b2-second-cavity-task-argmax-grid-check','wr.argmax_validation',SECOND,40))


def _read(path):return json.loads(path.read_text())


def _write(path,row):
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(row,indent=2,allow_nan=False)+'\n')
    temporary.replace(path)


def _queues():
    main=_read(ROOT/'QUEUE_CHECKPOINT.json') if (ROOT/'QUEUE_CHECKPOINT.json').exists() else {}
    p1b=_read(ROOT/'P1B_QUEUE_CHECKPOINT.json') if (ROOT/'P1B_QUEUE_CHECKPOINT.json').exists() else {}
    statuses=[str(main.get('status','')),str(p1b.get('status',''))]
    halted=any('HALTED' in value or value=='WAITING_FOR_PREDECESSOR_REVIEW' for value in statuses)
    released=main.get('status') in MAIN_FINISHED and p1b.get('status') in P1B_FINISHED
    return dict(main=main,P1B=p1b,halted=halted,released=released)


def _completed_validation(manifest):
    cfg=_read(ROOT/manifest)
    out=ROOT/'results/B2/validation'/cfg['id']
    for name in ('RESULTS.json','ALL_TASKS.json','MANIFEST.json','GRID_COMPARISON.json','joint8_PAIR.npz'):
        if not (out/name).exists():raise RuntimeError('completed endpoint artifact missing: '+str(out/name))
    if _read(out/'MANIFEST.json')!=cfg:raise RuntimeError('existing validation manifest differs')
    tasks=_read(out/'ALL_TASKS.json')
    if len(tasks)!=24 or {r['probe'] for r in tasks}!=set(range(24)):
        raise RuntimeError('endpoint validation has not completed all fixed24 tasks')
    if not all((out/f'probe{i}_PAIR.npz').exists() for i in range(24)):
        raise RuntimeError('endpoint validation raw pairs incomplete')
    gate=_read(out/'RESULTS.json')
    if not isinstance(gate.get('numerical_pass'),bool) or 'd_T_sup_24' not in gate:
        raise RuntimeError('endpoint final result marker missing')
    return out


def _nomination():
    rule_path=ROOT/'configs/SECOND_VALIDATION_NOMINATION_RULE.json'
    rule=_read(rule_path)
    if (rule.get('candidate_family')!='native_cavity_direct_2026102202_narrow'
            or rule.get('no_new_budget') is not True or rule.get('P1_prior_not_changed') is not True
            or rule.get('posterior_or_GT_selection') is not False or rule.get('max_validation_rhs')!=1120):
        raise RuntimeError('second nomination rule differs from registered bounded protocol')
    source=ROOT/'results/B2'/rule['candidate_family']
    if not (source/'COMPLETE.json').exists():return None,'CAVITY_TWO_METHODS_NOT_COMPLETE'
    complete=_read(source/'COMPLETE.json')
    if (complete.get('status')!='COMPLETED_TWO_NUMERICAL_METHODS' or complete.get('failure') is not None
            or set(complete.get('completed_methods',[]))!={'FULL_J','RANDOM_MULTISTART'}):
        return None,'CAVITY_TWO_METHODS_NOT_COMPLETE'
    # No target/GT file is opened. Read only frozen per-method calibration field.
    candidates=[]
    for method in ('FULL_J','RANDOM_MULTISTART'):
        path=source/f'{method}_RESULTS.json'
        latent=source/f'{method}_FINAL.npz'
        if not path.exists() or not latent.exists():return None,'CAVITY_METHOD_OUTPUT_NOT_COMPLETE'
        row=_read(path)
        if (row.get('status')!='NUMERICAL_METHOD_COMPLETED_NEEDS_INDEPENDENT_VALIDATION'
                or row.get('method')!=method or row.get('kind')!='cavity'
                or row.get('seed')!=2026102202 or row.get('band')!='narrow'):
            return None,'CAVITY_METHOD_OUTPUT_NOT_COMPLETE'
        distance=row.get('final_distance')
        if isinstance(distance,bool) or not isinstance(distance,(int,float)) or not math.isfinite(distance) or distance<0:
            raise RuntimeError('invalid completed cavity calibration distance')
        candidates.append(dict(method=method,d_cal=float(distance),x_file=str(latent.relative_to(ROOT))))
    selected=min(candidates,key=lambda row:row['d_cal'])
    if selected['d_cal']>1:return None,'NO_COMPLETED_CAVITY_CALIBRATION_CANDIDATE_LE_1'
    cfg=dict(id='cavity_seed2_narrow_bestcal',seed=2026102202,band='narrow',kind='cavity',
             chart='linear',x_file=selected['x_file'],cavities=[[[ -.85,.85],.06]],
             identity='ORACLE_CONSTRUCTED_PHYSICAL_WITNESS',selected_method=selected['method'],
             coarse_calibration_distance=selected['d_cal'],
             candidate_calibration_distances={r['method']:r['d_cal'] for r in candidates},
             nomination_rule='smaller completed FULL_J/RANDOM_MULTISTART coarse calibration distance <=1; fixed order breaks exact ties',
             registered_rule='configs/SECOND_VALIDATION_NOMINATION_RULE.json',
             selection_uses_only_calibration=True,GT_or_target_selection=False,
             gate_before_validation='PROVISIONAL_ONLY',controlled_P1_prior_unchanged=True,
             original_minimum_feature_status='UNRESOLVED')
    return cfg,None


def _job_state(label,module,manifest,proposed_nomination=None):
    directory=ROOT/'results/jobs'/label
    receipt=directory/'RECEIPT.json'
    cfg=(_read(ROOT/manifest) if (ROOT/manifest).exists()
         else proposed_nomination if manifest==SECOND else None)
    result_dir=(ROOT/'results/B2/validation'/cfg['id']) if cfg else None
    if receipt.exists():
        record=_read(receipt)
        if record.get('returncode')!=0 or record.get('status')!='COMPLETE':
            raise RuntimeError('failed receipt requires parent review, never automatic rerun: '+label)
        if record.get('command',[])[-3:]!=['-m',module,manifest]:
            raise RuntimeError('completed receipt command differs from requested job: '+label)
        if cfg is None:raise RuntimeError('receipt exists but nomination manifest missing')
        endpoint=_completed_validation(manifest)
        if module=='wr.argmax_validation':
            path=endpoint/'argmax_validation/TASK_ARGMAX_GRID_CHECK.json'
            if not path.exists() or not (endpoint/'argmax_validation/COMPLETE.json').exists():
                raise RuntimeError('argmax receipt complete but output incomplete')
            if _read(path).get('task_argmax_numerical_pass') is not True:
                raise RuntimeError('argmax numerical failure needs review')
        return 'COMPLETE'
    if directory.exists():raise RuntimeError('wrapper directory without complete receipt preserved for review: '+label)
    if result_dir is not None:
        output=result_dir if module=='wr.endpoint_validation' else result_dir/'argmax_validation'
        if output.exists():
            raise RuntimeError('prior output without this successful wrapper receipt; inspect, do not overwrite: '+str(output))
    return 'NOT_RUN'


def _budget(required):
    path=ROOT/'PHYSICS_LEDGER.jsonl'
    if path.exists():
        with path.open() as stream:
            fcntl.flock(stream,fcntl.LOCK_SH)
            rows=[json.loads(line) for line in stream if line.strip()]
    else:rows=[]
    reserves=[r for r in rows if r.get('stage')=='B2' and r.get('event')=='reserve']
    used=sum(r.get('rhs',0) for r in reserves if r.get('scope')=='validation')
    allocation=_read(ROOT/'configs/B2_ALLOCATION.json')['scopes']['validation']['rhs']
    authorization=_read(ROOT/'configs/AUTHORIZATION.json')['stages']['B2']
    stage_rhs=sum(r.get('rhs',0) for r in reserves)
    stage_lu=sum(r.get('factorizations',0) for r in reserves)
    cpu=sum(r.get('charged_cpu_seconds',0) for r in rows if r.get('stage')=='B2')
    if allocation!=2500:raise RuntimeError('registered validation scope cap changed; no implicit new allowance')
    result=dict(required_new_rhs=required,validation_rhs_used=used,validation_rhs_cap=allocation,
                validation_rhs_projected_max=used+required,stage_rhs_available=authorization['rhs']-stage_rhs,
                stage_LU_available=authorization['factorizations']-stage_lu,
                stage_CPU_available=authorization['cpu_seconds']-cpu,no_budget_transfer=True)
    if (used+required>allocation or required>result['stage_rhs_available']
            or math.ceil(required/4)>result['stage_LU_available'] or result['stage_CPU_available']<=0):
        raise RuntimeError('remaining existing B2/validation budget cannot fund whole remaining serial tail')
    return result


def main(max_wait_seconds=86400,poll_seconds=30):
    if not math.isfinite(max_wait_seconds) or max_wait_seconds<=0 or not 1<=poll_seconds<=60:
        raise ValueError('finite positive wait and 1..60 second polling required')
    # Only this finite coordinator is locked; other queues release via their
    # terminal checkpoints. No daemon, scheduler, reservation or budget is added.
    with (ROOT/'B2_FINISH_QUEUE.lock').open('a') as lock:
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:raise RuntimeError('another B2 tail coordinator owns the serial lock')
        deadline=time.monotonic()+max_wait_seconds
        while True:
            states=_queues()
            if states['halted']:
                result=dict(status='WAITING_FOR_PARENT_REVIEW',queues=states,no_new_physics=True)
                _write(CHECKPOINT,result);return result
            if states['released']:break
            if time.monotonic()>=deadline:
                result=dict(status='FINITE_WAIT_EXPIRED_NOT_RUN',queues=states,no_new_physics=True)
                _write(CHECKPOINT,result);return result
            _write(CHECKPOINT,dict(status='WAITING_FOR_SERIAL_P1A_P1B_RELEASE',queues=states,
                                   no_new_physics=True,deadline_monotonic=deadline))
            time.sleep(min(poll_seconds,max(.01,deadline-time.monotonic())))
        cfg,not_run=_nomination()
        if cfg is None:
            result=dict(status='NOT_RUN_'+not_run,no_new_physics=True,nomination_created=False)
            _write(CHECKPOINT,result);return result
        target=ROOT/SECOND
        if target.exists():
            if _read(target)!=cfg:raise RuntimeError('existing second nomination differs; never overwrite')
        # Determine unfinished jobs before creating a nomination; any previous
        # failed/partial wrapper output stops instead of spending/resetting again.
        states={label:_job_state(label,module,manifest,cfg) for label,module,manifest,_ in JOBS}
        required=sum(cap for label,_,_,cap in JOBS if states[label]!='COMPLETE')
        budget=_budget(required)
        if not target.exists():
            with target.open('x') as stream:
                json.dump(cfg,stream,indent=2,allow_nan=False);stream.write('\n');stream.flush();os.fsync(stream.fileno())
        completed=[label for label,_,_,_ in JOBS if states[label]=='COMPLETE']
        for label,module,manifest,cap in JOBS:
            if label in completed:continue
            queues=_queues()
            if queues['halted'] or not queues['released']:
                result=dict(status='WAITING_FOR_PARENT_REVIEW_SERIAL_RELEASE_CHANGED',queues=queues,completed=completed)
                _write(CHECKPOINT,result);return result
            remaining=sum(c for name,_,_,c in JOBS if name not in completed)
            budget=_budget(remaining)
            _write(CHECKPOINT,dict(status='RUNNING',label=label,completed=completed,budget_forecast=budget,
                                   per_job_rhs_max=cap,selection='coarse calibration only'))
            code=subprocess.run([str(Path(sys.executable).absolute()),str(ROOT/'tools/run_job.py'),
                                 'B2',label,str(Path(sys.executable).absolute()),'-m',module,manifest],
                                cwd=ROOT.parents[1]).returncode
            if code:
                result=dict(status='HALTED_REQUIRES_REVIEW',label=label,returncode=code,completed=completed,
                            failed_runs_preserved=True,no_automatic_retry=True)
                _write(CHECKPOINT,result);return result
            if _job_state(label,module,manifest)!='COMPLETE':raise RuntimeError('job missing successful receipt/output')
            completed.append(label)
        result=dict(status='B2_VALIDATION_TAIL_FINISHED_PENDING_PARENT_REVIEW',completed=completed,
                    budget_forecast=budget,nomination=SECOND,selection_uses_target_or_GT=False,
                    controlled_P1_prior_unchanged=True,no_new_budget=True)
        _write(CHECKPOINT,result);return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-wait-seconds',type=float,default=86400)
    parser.add_argument('--poll-seconds',type=float,default=30)
    args=parser.parse_args()
    try:
        result=main(args.max_wait_seconds,args.poll_seconds)
    except Exception as exc:
        _write(CHECKPOINT,dict(status='HALTED_REQUIRES_REVIEW',error_type=type(exc).__name__,
                               error=str(exc),no_automatic_retry=True))
        raise
    if result['status']!='B2_VALIDATION_TAIL_FINISHED_PENDING_PARENT_REVIEW':raise SystemExit(2)
