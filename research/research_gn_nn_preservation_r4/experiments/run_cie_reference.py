"""Public-only, CPU CIE baseline CLI. Default is a finite NOT_RUN plan."""
import os
for _k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[_k]='1'
import argparse,json,sys,time,platform
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
from multiprocessing import get_context
import numpy as np
HERE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(HERE.parents[1]))
from research.research_gn_nn_preservation_r4.src.cie_reference import run_reference
from research.research_gn_nn_completion_r3.src.evaluate import save_result

def execute_one(task):
    # Public allowlists only: no generic loader can accidentally import truth.
    tick=time.perf_counter()
    with np.load(HERE/'data/operators.npz',allow_pickle=False)as z:public={k:np.array(z[k])for k in ('G','E0','S')}
    with np.load(HERE/'data'/task['record']['public_path'],allow_pickle=False)as z:public['Y']=np.array(z['Y'])
    preparation=time.perf_counter()-tick
    receipt,arrays=run_reference(public,method=task['method'],noise=task['record']['noise'],preparation_seconds=preparation)
    receipt.update(obs_id=task['record']['obs_id'],parent_id=task['record']['parent_id'],split=task['record']['split'],run=task['run'],
        execution=dict(device='CPU',host_system=platform.system(),machine=platform.machine(),processor=platform.processor(),BLAS_threads=1,cold_case=True,public_initializer=.1,
            data_inputs='G/E0/S/Y/noise only; no truth, support or label',python_process_startup='excluded from per-case timer, included in queue elapsed'))
    return task,receipt,arrays

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--split',choices=['val','test'],default='test')
    p.add_argument('--scope',choices=['figure13','all108'],default='figure13')
    p.add_argument('--obs-id',nargs='*',help='For public VAL cost checks only; cannot select TEST after the frozen parent plan')
    p.add_argument('--methods',nargs='+',choices=['cie_single','cie_dual'],default=['cie_dual'])
    p.add_argument('--workers',type=int,choices=[1,2],default=1)
    p.add_argument('--execute',action='store_true')
    p.add_argument('--test-methods-frozen',action='store_true',help='Parent assertion that TEST method selection is frozen before this explicit execution')
    p.add_argument('--output',type=Path,default=HERE/'evaluation')
    a=p.parse_args();online=json.loads((HERE/'data/manifests/online.json').read_text())['records']
    if a.split=='val':
        if not a.obs_id:raise ValueError('Explicit public VAL observation ID(s) required for a cost check')
        records=[r for r in online if r['split']=='val'and r['obs_id']in a.obs_id]
        if len(records)!=len(set(a.obs_id)):raise ValueError('A requested VAL observation is missing or not VAL')
    else:
        if a.obs_id:raise ValueError('TEST identities must follow the preselected figure13 or full108 plan')
        records=[r for r in online if r['split']=='test']
        if a.scope=='figure13':
            selection=json.loads((HERE/'audit/FIGURE_SELECTION.json').read_text())
            parents={r['parent_id']for r in selection['rows']}
            if len(parents)!=13:raise ValueError('Frozen figure subset must have 13 mothers')
            records=[r for r in records if r['parent_id']in parents]
            if len(records)!=26:raise ValueError('Frozen subset must include both noise observations of all 13 mothers')
        elif len({r['parent_id']for r in records})!=108:raise ValueError('Full TEST plan is not 108 mothers')
    tasks=[dict(record=r,method=method,run=method+'_cold55')for r in records for method in a.methods]
    plan=dict(status='NOT_RUN'if not a.execute else'EXECUTE_REQUESTED',split=a.split,scope=a.scope,methods=a.methods,
        parent_test_freeze_assertion=a.test_methods_frozen,workers=a.workers,BLAS_threads_per_worker=1,
        mothers=len({r['parent_id']for r in records}),observations=len(records),tasks=len(tasks),
        observation_ids=[r['obs_id']for r in records],truth_access='NONE_PUBLIC_INPUT_ONLY',
        budget=dict(wall_seconds=55,max_LU=600,max_complex_RHS=36000),comparison='cold CIE55s versus common NN GN45s plus post10s; actual costs remain separate',
        baseline_schedule='mf3/beta1 → mf6/beta.2 → mf10/beta.2; min50, max500/stage, coefficient_tol1e-8; later dual beta_current1',
        prior='zero spatial prior, declared adaptation matching NN data-only correction',resource_coordination='parent coordinates CPU training; no remote i9 use')
    if not a.execute:print(json.dumps(plan,indent=2));return
    if a.split=='test'and not a.test_methods_frozen:raise ValueError('TEST execution requires the parent frozen-method assertion')
    out=HERE/'delegated/cie';out.mkdir(parents=True,exist_ok=True)
    stamp=time.time_ns();(out/f'PLAN_{stamp}.json').write_text(json.dumps(plan,indent=2))
    pending=[];existing=[]
    for task in tasks:
        dest=a.output/task['run'];base=dest/task['record']['obs_id']
        if base.with_suffix('.json').exists()or base.with_suffix('.npz').exists():existing.append(dict(run=task['run'],obs_id=task['record']['obs_id'],status='EXISTING_PRESERVED'))
        else:pending.append(task)
    started=time.perf_counter();entries=[]
    def save(result):
        task,receipt,arrays=result
        save_result(a.output/task['run'],task['record']['obs_id'],receipt,arrays)
        row=dict(run=task['run'],obs_id=task['record']['obs_id'],status=receipt['status'],selection_status=receipt['selection_status'],wall_seconds=receipt['wall_seconds'],costcounts=receipt['costcounts'],current_status=receipt['current_phase']['status'])
        entries.append(row);print(json.dumps(row),flush=True)
    if a.workers==1:
        for task in pending:save(execute_one(task))
    else:
        with ProcessPoolExecutor(max_workers=2,mp_context=get_context('spawn'))as pool:
            futures={pool.submit(execute_one,task):task for task in pending}
            for future in as_completed(futures):save(future.result())
    summary=dict(plan,status='FINITE_PUBLIC_QUEUE_CLOSED',queue_wall_seconds=time.perf_counter()-started,entries=entries,existing=existing,
        total_LU=sum(r['costcounts'].get('LU',0)for r in entries),total_complex_RHS=sum(r['costcounts'].get('actual_complex_RHS',0)for r in entries))
    (out/f'SUMMARY_{stamp}.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary),flush=True)

if __name__=='__main__':main()
