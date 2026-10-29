#!/usr/bin/env python3
"""Before-GT supplemental GN budget fill; immutable original source/results."""
import os
for _key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_key]='1'
import argparse,copy,json,shutil,time,traceback
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
from threadpoolctl import threadpool_limits
from src.solver import Run,load_public,BudgetStop,StageStop,json_safe
ROOT=Path(__file__).resolve().parent
CONTRACT=ROOT/'configs/BUDGET_FILL_CONTRACT.json'
SNAPSHOT=ROOT/'audit/BUDGET_FILL_SOURCE_SNAPSHOT'


def utc():return datetime.now(timezone.utc).isoformat()


def write_json(path,data):
    path=Path(path)
    if path.exists():raise FileExistsError('Preserving existing output '+str(path))
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(json_safe(data),indent=2,allow_nan=False)+'\n')


def freeze():
    if SNAPSHOT.exists():raise FileExistsError('Budget-fill snapshot already exists')
    SNAPSHOT.mkdir(parents=True)
    names=['budget_fill.py','configs/BUDGET_FILL_CONTRACT.json','src/solver.py','src/physics.py','src/multiscale.py','src/__init__.py']
    files=[]
    for name in names:
        target=SNAPSHOT/name;target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(ROOT/name,target);files.append(dict(path=name,bytes=target.stat().st_size))
    write_json(SNAPSHOT/'FREEZE.json',dict(frozen_utc=utc(),files=files,
        scope='Supplemental budget fairness control before any main GT scoring',
        after_public_Y_generation=True,fully_pre_generation_source_freeze=False,
        original_primary_source_unchanged=True,truth_access='NONE',hashes='NOT_COMPUTED'))


def check_freeze(contract):
    manifest=json.loads((SNAPSHOT/'FREEZE.json').read_text())
    for row in manifest['files']:
        if (ROOT/row['path']).read_bytes()!=(SNAPSHOT/row['path']).read_bytes():
            raise RuntimeError('Source changed after budget-fill freeze: '+row['path'])
    marker=ROOT/contract['original_prediction_freeze_marker']
    if not marker.exists():raise RuntimeError('All original predictions must freeze before continuation: '+str(marker))
    return manifest


def add_numeric(a,b):
    result=copy.deepcopy(a or {})
    for k,v in (b or {}).items():
        if isinstance(v,dict):result[k]=add_numeric(result.get(k,{}),v)
        elif isinstance(v,(int,float)) and not isinstance(v,bool):result[k]=result.get(k,0)+v
    return result


def run_fill(original_path,public_path,operators_path,output_path):
    contract=json.loads(CONTRACT.read_text());manifest=check_freeze(contract)
    original_path=Path(original_path);original=json.loads(original_path.read_text())
    method=original['method'];old_wall=float(original['wall_seconds']);remaining=max(0.,float(contract['total_numerical_seconds'])-old_wall)
    root=copy.deepcopy(original)
    initial=np.asarray(original['top1']['m'],float)
    if initial.shape!=(576,) or not np.all(np.isfinite(initial)):raise ValueError('Original online incumbent must have576finite material cells')
    counts0=original.get('accounting',{}).get('counts',{})
    inc=None;live=None;error=None;inc_wall=inc_cpu=0.;fill_status='SKIPPED_REMAINING_BELOW_ATOMIC_THRESHOLD'
    ledger=dict(original_input_path=str(original_path.resolve()),original_input_bytes=original_path.stat().st_size,
                original_wall_seconds=old_wall,remaining_seconds=remaining,
                source_freeze_utc=manifest['frozen_utc'],source_freeze_after_Y_generation=True,
                source_freeze_fully_pre_generation=False,initialization='original online top1 material; no material reset',
                reused_factorization=False,new_LU_charged=True,tail_method='reciprocal GN' if method=='fullj_reciprocal' else 'native GN',
                original_primary_unchanged=True,truth_access='NONE')
    if remaining>=float(contract['minimum_atomic_seconds']):
        begun,cpu=time.perf_counter(),time.process_time()
        try:
            cfg=dict(original['config']);cfg.update(wall_seconds=remaining,max_gn_iterations=int(contract['continuation_max_gn_iterations']))
            previous=[e for e in original.get('events',[]) if e.get('kind')=='GN_ITERATION' and e.get('stage')==24 and e.get('accepted')]
            if previous and np.isfinite(previous[-1].get('damping',np.nan)):
                cfg['initial_damping']=float(previous[-1]['damping'])
            ledger['continuation_initial_damping']=cfg['initial_damping']
            with threadpool_limits(limits=1):
                live=Run(load_public(public_path,operators_path),method,cfg)
                # Include public loading, new Run construction and all refactorization in the remaining cap.
                live.started,live.cpu_started=begun,cpu
                live.active_stage=24;live.stage_deadline=remaining;live.phase_deadline=remaining
                summary={};phase_status='OPTIMIZER_FINISHED'
                try:
                    _,summary=live.gn(initial,np.eye(live.model.n),remaining,
                                      reciprocal=method=='fullj_reciprocal',origin='budget_fill_original_tail',
                                      max_iterations=int(contract['continuation_max_gn_iterations']))
                    fill_status='CONTINUATION_CONVERGED_OR_NO_ACCEPTED_STEP'
                except (BudgetStop,StageStop):
                    phase_status='REMAINING_BUDGET_STOP';live.status='BUDGET_STOP';live.stop_reason='REMAINING_ORIGINAL_BUDGET'
                    fill_status='REMAINING_BUDGET_USED'
                except BaseException as exc:
                    error=dict(type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc())
                    phase_status='FAILED_PAID';live.status='FAILED_PAID_UNEXPECTED';live.stop_reason=type(exc).__name__
                    live.failures[type(exc).__name__]+=1;fill_status='FAILED_PAID_ORIGINAL_INCUMBENT_RETAINED'
                finally:
                    live.events.append(dict(kind='MULTISCALE_STAGE',stage=24,parameter_dimension=576,
                         fixed_physics_dimension=576,started_seconds=0.,deadline_seconds=remaining,
                         wall_seconds=live.elapsed(),status=phase_status,phase='BUDGET_FILL',
                         initial_material=initial,initial_coarse_material=initial,coarse_seed_counts=counts0,
                         counts=dict(live.model.counts),has_evaluated_candidate=bool(live.archive),
                         completed_gn_iterations=sum(e.get('kind')=='GN_ITERATION' for e in live.events),**summary))
                inc=live.finish()
        except BaseException as exc:
            error=dict(type=type(exc).__name__,message=str(exc),traceback=traceback.format_exc())
            fill_status='FAILED_INITIALIZATION_ORIGINAL_INCUMBENT_RETAINED'
            if live is not None:inc=live.finish()
        inc_wall,inc_cpu=time.perf_counter()-begun,time.process_time()-cpu
    inc_accounting=inc['accounting'] if inc is not None else dict(counts={},full_factorizations=0,full_complex_rhs=0)
    if error and live is None:inc_accounting['unavailable_counts_reason']='Run initialization failed before live counters were available'
    # All counters in incremental_accounting are new cost only. Root totals include original once.
    root['accounting']=add_numeric(original.get('accounting',{}),inc_accounting)
    root['wall_seconds']=old_wall+inc_wall
    root['cpu_seconds']=float(original.get('cpu_seconds',0.))+inc_cpu
    root['budget_seconds']=float(contract['total_numerical_seconds'])
    root['budget_overshoot_seconds']=max(0.,root['wall_seconds']-root['budget_seconds'])
    root['incremental_wall_seconds']=inc_wall;root['incremental_cpu_seconds']=inc_cpu
    root['incremental_accounting']=inc_accounting
    root['budget_fill_status']=fill_status;root['original_status']=original.get('status')
    root['budget_fill_error']=error
    root['phase_ledger']=[dict(phase='ORIGINAL_PRIMARY',input_path=ledger['original_input_path'],input_bytes=ledger['original_input_bytes'],
                              wall_seconds=old_wall,cpu_seconds=original.get('cpu_seconds'),cost_accounting='original root input; counted once'),
                          dict(phase='BUDGET_FILL',wall_seconds=inc_wall,cpu_seconds=inc_cpu,
                               cost_accounting='incremental_accounting only when aggregating original+supplement',**ledger)]
    root['budget_fill_contract']=contract
    root['original_root_input_path']=ledger['original_input_path'];root['original_root_input_bytes']=ledger['original_input_bytes']
    pool=copy.deepcopy(original.get('candidates',original.get('archive',[])))
    if inc is not None:
        for row in inc['candidates']:
            row=copy.deepcopy(row);row['discovery_seconds']=old_wall+row.get('discovery_seconds',0.)
            row['discovery_counts']=add_numeric(counts0,row.get('discovery_counts',{}));row['phase']='BUDGET_FILL';pool.append(row)
        events=[]
        for e in inc['events']:
            e=copy.deepcopy(e);e['phase']='BUDGET_FILL'
            for key in ('seconds','started_seconds','deadline_seconds'):
                if isinstance(e.get(key),(int,float)):e[key]+=old_wall
            events.append(e)
        root['events']=copy.deepcopy(original.get('events',[]))+[e for e in events if e.get('kind')!='MULTISCALE_STAGE']
        trace=copy.deepcopy(original.get('best_trace',[]))
        for item in inc.get('best_trace',[]):
            item=copy.deepcopy(item);item['seconds']=old_wall+item.get('seconds',0.)
            item['counts']=add_numeric(counts0,item.get('counts',{}));item['phase']='BUDGET_FILL'
            if not trace or item['objective']<trace[-1]['objective']:trace.append(item)
        root['best_trace']=trace
        # Preserve original stages and merge the new tail into the actual last576stage, never replace its prior cost.
        stages=copy.deepcopy(original.get('stages',[]))
        base=stages[-1] if stages else dict(kind='MULTISCALE_STAGE',stage=24,parameter_dimension=576,fixed_physics_dimension=576)
        extra=inc['stages'][-1] if inc['stages'] else {}
        base['counts']=add_numeric(base.get('counts',{}),inc_accounting.get('counts',{}))
        base['wall_seconds']=float(base.get('wall_seconds',0.))+inc_wall
        base['completed_gn_iterations']=int(base.get('completed_gn_iterations',0))+int(extra.get('completed_gn_iterations',0))
        base['has_evaluated_candidate']=bool(base.get('has_evaluated_candidate') or inc.get('candidates'))
        base['budget_fill_phase']=dict(status=fill_status,remaining_seconds=remaining,wall_seconds=inc_wall,counts=inc_accounting.get('counts',{}))
        if stages:stages[-1]=base
        else:stages=[base]
        root['stages']=stages
        root['events']=[e for e in root['events'] if not(e.get('kind')=='MULTISCALE_STAGE' and e.get('stage')==24)]+[copy.deepcopy(base)]
        root['failures']=add_numeric(original.get('failures',{}),inc.get('failures',{}))
        root['all_unknowns_released']=bool(original.get('all_unknowns_released') or inc.get('all_unknowns_released'))
    chosen=[]
    for row in sorted(pool,key=lambda x:x['objective']):
        if all(np.sqrt(np.mean((np.asarray(row['m'])-np.asarray(old['m']))**2))>float(original['config']['min_material_rms']) for old in chosen):
            chosen.append(row)
            if len(chosen)>=int(original['config']['archive_size']):break
    root['candidates']=chosen;root['archive']=chosen;root['top1']=chosen[0] if chosen else original['top1'];root['m']=root['top1']['m']
    root['truth_access']='NONE_PUBLIC_INPUT_ONLY_BEFORE_ANY_MAIN_GT_SCORE'
    root['supplemental_control']='BUDGET_FILL; separately frozen after publicY generation, before continuation/GT; original primary preserved'
    write_json(output_path,root)
    print(json.dumps(dict(output=str(output_path),method=method,budget_fill_status=fill_status,remaining_seconds=remaining,
                          incremental_wall_seconds=inc_wall,incremental_cpu_seconds=inc_cpu,wall_seconds=root['wall_seconds'])))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--freeze',action='store_true')
    p.add_argument('--original');p.add_argument('--public');p.add_argument('--operators');p.add_argument('--output')
    a=p.parse_args()
    if a.freeze:freeze();return
    if not all((a.original,a.public,a.operators,a.output)):p.error('Supply --original --public --operators --output')
    run_fill(a.original,a.public,a.operators,a.output)


if __name__=='__main__':main()
