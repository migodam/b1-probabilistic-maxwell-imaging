"""Public-only CPU CIE CLI; finite NOT_RUN plan by default, no scoring."""
import os
for _k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[_k]='1'
import argparse,json,sys,time,platform
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parents[1];sys.path.insert(0,str(HERE.parents[1]))
from research.research_gn_nn_preservation_r4.src.cie_reference import run_reference
from research.research_gn_nn_completion_r3.src.evaluate import save_result
from research.research_gn_nn_preservation_r4.experiments.run_test_partition import (load_contract,load_partition,require_host,online_records,complete_pair,run_fifo,canonical_costs,write_json,read_json,baseline744_complete,pair_host,DEADLINE,DEADLINE_TEXT)

def execute_one(task):
    tick=time.perf_counter()
    with np.load(HERE/'data/operators.npz',allow_pickle=False)as z:public={k:np.array(z[k])for k in ('G','E0','S')}
    with np.load(HERE/'data'/task['record']['public_path'],allow_pickle=False)as z:public['Y']=np.array(z['Y'])
    preparation=time.perf_counter()-tick
    receipt,arrays=run_reference(public,method=task['method'],noise=task['record']['noise'],preparation_seconds=preparation)
    receipt.update(obs_id=task['record']['obs_id'],parent_id=task['record']['parent_id'],split=task['record']['split'],run=task['run'],execution=dict(device='CPU',platform=platform.platform(),host_system=platform.system(),machine=platform.machine(),processor=platform.processor(),pid=os.getpid(),execution_partition=task.get('execution_partition'),global_online_index=task['global_online_index'],BLAS_threads=1,cold_case=True,public_initializer=.1,data_inputs='G/E0/S/Y/noise only; no truth, support or label',python_process_startup='excluded from per-case timer, included in queue elapsed'))
    return receipt,arrays

def main():
    p=argparse.ArgumentParser();p.add_argument('--split',choices=['val','test'],default='test');p.add_argument('--scope',choices=['figure13','all108'],default='figure13');p.add_argument('--obs-id',nargs='*',help='Explicit public VAL cost checks only');p.add_argument('--methods',nargs='+',choices=['cie_single','cie_dual'],default=['cie_dual']);p.add_argument('--workers',type=int,choices=[1,2,4],default=1);p.add_argument('--execute',action='store_true');p.add_argument('--test-methods-frozen',action='store_true');p.add_argument('--execution-partition',choices=['local','remote']);p.add_argument('--output',type=Path,default=HERE/'evaluation');p.add_argument('--worker-obs',help=argparse.SUPPRESS);a=p.parse_args()
    setup=time.perf_counter();online=online_records();global_indices={r['obs_id']:i for i,r in enumerate(online)}
    full_test=[r for r in online if r['split']=='test'];global_records=full_test
    if a.execution_partition:require_host(a.execution_partition,a.workers)
    elif a.workers>(4 if sys.platform=='win32' else 2):raise ValueError('At most2 local /4 Windows CPU workers')
    if a.split=='val':
        if not a.obs_id:raise ValueError('Explicit public VAL IDs required')
        records=[r for r in online if r['split']=='val' and r['obs_id'] in a.obs_id]
        if len(records)!=len(set(a.obs_id)):raise ValueError('Requested VAL observation missing')
        if a.execution_partition:raise ValueError('TEST host partition does not apply to VAL cost checks')
        global_records=records
    else:
        if a.obs_id:raise ValueError('TEST follows only frozen full108/figure13 plan and host allocation')
        if len(full_test)!=216 or len({r['parent_id']for r in full_test})!=108:raise ValueError('Full TEST denominator must remain216/108')
        records=full_test
        if a.scope=='figure13':
            selection=read_json(HERE/'audit/FIGURE_SELECTION.json');parents={r['parent_id']for r in selection['rows']}
            if len(parents)!=13:raise ValueError('Frozen figure subset requires13 mothers')
            records=[r for r in full_test if r['parent_id'] in parents]
        if a.execution_partition:
            allocation,_=load_partition();allowed=set(allocation[a.execution_partition]['obs_ids']);records=[r for r in records if r['obs_id'] in allowed]
        if a.execute:
            if not a.execution_partition or not a.test_methods_frozen:raise ValueError('TEST execution requires frozen methods and explicit assigned host')
            freeze,final,_,_=load_contract();registered={q['method'] for q in final['plans'] if q['kind']=='cie'}
            if any(m not in registered for m in a.methods):raise ValueError('Unregistered CIE TEST method')
            if not baseline744_complete():raise RuntimeError('Wait for all744 prior baseline pairs before CIE TEST CPU allocation')
    tasks=[dict(record=r,method=method,run=method+'_cold55',execution_partition=a.execution_partition,global_online_index=global_indices[r['obs_id']]) for r in records for method in a.methods]
    global_combinations=[dict(obs_id=r['obs_id'],parent_id=r['parent_id'],run=method+'_cold55',status='NOT_RUN',assigned_here=r['obs_id'] in {q['obs_id']for q in records})for r in global_records for method in a.methods]
    plan=dict(status='NOT_RUN' if not a.execute else 'EXECUTE_REQUESTED',split=a.split,scope=a.scope,methods=a.methods,parent_test_freeze_assertion=a.test_methods_frozen,execution_partition=a.execution_partition,workers=a.workers,BLAS_threads_per_worker=1,mothers=len({r['parent_id'] for r in records}),observations=len(records),tasks=len(tasks),global_mothers=108 if a.split=='test' else len({r['parent_id']for r in records}),global_observations=216 if a.split=='test' else len(records),global_plan_combinations=len(global_combinations),global_plan=global_combinations,observation_ids=[r['obs_id']for r in records],truth_access='NONE_PUBLIC_INPUT_ONLY',budget=dict(wall_seconds=55,max_LU=600,max_complex_RHS=36000),absolute_deadline_utc=DEADLINE_TEXT,comparison='cold CIE55s versus NN GN45s+post10s; actual costs separate',baseline_schedule='mf3/beta1 → mf6/beta.2 → mf10/beta.2; min50,max500/stage,coefficient_tol1e-8; laterdual beta_current1',prior='zero spatial prior; unchanged data-only adaptation',setup_wall_seconds=time.perf_counter()-setup,resource_coordination='same preassigned host per mother for every TEST method; own2local/4remote CPUs only')
    if not a.execute:print(json.dumps(plan,indent=2));return
    if a.worker_obs:
        matched=[t for t in tasks if t['record']['obs_id']==a.worker_obs]
        if len(matched)!=1:raise ValueError('Worker must select exactly one declared observation/method')
        t=matched[0];base=a.output/t['run']/a.worker_obs
        if complete_pair(base):
            if a.split=='test' and pair_host(base)!=a.execution_partition:raise RuntimeError('Existing CIE case host provenance differs; preserved')
            print(json.dumps(dict(status='REUSED_COMPLETE_PAIR',obs_id=a.worker_obs,new_physics_counts={})));return
        if base.with_suffix('.npz').exists() or base.with_suffix('.json').exists():raise RuntimeError('Preserve partial CIE case; interrupted cost unknown')
        if time.time()>=DEADLINE:print(json.dumps(dict(status='NOT_RUN',reason='ABSOLUTE_DEADLINE')));return
        receipt,arrays=execute_one(t);save_result(a.output/t['run'],a.worker_obs,receipt,arrays)
        print(json.dumps(dict(run=t['run'],obs_id=a.worker_obs,status=receipt['status'],wall_seconds=receipt['wall_seconds'],costcounts=receipt['costcounts'])),flush=True);return
    out=HERE/'delegated/cie';out.mkdir(parents=True,exist_ok=True);stamp=time.time_ns();write_json(out/f'PLAN_{stamp}.json',plan)
    pending=[];existing=[];held=[]
    for t in tasks:
        obs=t['record']['obs_id'];base=a.output/t['run']/obs
        if complete_pair(base):
            if a.split=='test' and pair_host(base)!=a.execution_partition:held.append(dict(run=t['run'],obs_id=obs,status='HOLD_EXISTING_HOST_PROVENANCE'));continue
            existing.append(dict(run=t['run'],obs_id=obs,status='EXISTING_PRESERVED'));continue
        if base.with_suffix('.npz').exists() or base.with_suffix('.json').exists():held.append(dict(run=t['run'],obs_id=obs,status='HOLD_PARTIAL_EXISTING_UNKNOWN_COST'));continue
        command=[sys.executable,str(Path(__file__).resolve()),'--split',a.split,'--scope',a.scope,'--methods',t['method'],'--workers','1','--execute','--worker-obs',obs,'--output',str(a.output)]
        if a.split=='test':command+=['--test-methods-frozen','--execution-partition',a.execution_partition]
        else:command+=['--obs-id',*a.obs_id]
        pending.append(dict(obs_id=obs,run=t['run'],command=command,log=str(out/f'CASE_{obs}_{t["method"]}_{stamp}.log')))
    began=time.perf_counter();results=run_fifo(pending,a.workers,out/f'QUEUE_{stamp}.json',label='CIE_'+str(a.execution_partition))
    entries=[]
    for t in tasks:
        base=a.output/t['run']/t['record']['obs_id']
        if complete_pair(base):
            r=read_json(base.with_suffix('.json'));entries.append(dict(run=t['run'],obs_id=t['record']['obs_id'],status=r['status'],selection_status=r.get('selection_status'),wall_seconds=r['wall_seconds'],costcounts=r['costcounts']))
    for row in global_combinations:
        base=a.output/row['run']/row['obs_id']
        if complete_pair(base):row['status']=read_json(base.with_suffix('.json')).get('status','COMPLETED')
    summary=dict(plan,status='FINITE_PUBLIC_QUEUE_CLOSED',queue_wall_seconds=time.perf_counter()-began,entries=entries,existing=existing,held=held,jobs=results,total_LU=sum(r['costcounts'].get('LU',0)for r in entries),total_complex_RHS=sum(r['costcounts'].get('actual_complex_RHS',0)for r in entries),cost_scope='canonical saved case costs once; includes reused cases, never add reused physics a second time',unknown_interruption_policy='UNKNOWN_NOT_ZERO; missing planned combinations stayNOT_RUN')
    write_json(out/f'SUMMARY_{stamp}.json',summary);print(json.dumps(summary),flush=True)
if __name__=='__main__':main()
