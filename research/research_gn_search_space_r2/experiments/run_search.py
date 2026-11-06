"""Finite task list, public-only prediction runner, independently sealable receipts."""
import argparse,json,os,sys,time,shutil
from concurrent.futures import ProcessPoolExecutor,as_completed
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from src.search import HIST,load_public,solve

def dump(path,value):path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
def tasks():
    rows=json.loads((HIST/'RESULTS_SUMMARY.json').read_text())['rows']
    rows=sorted([r for r in rows if r['backend']=='cpu'and r['method']=='native_trf'and r['amplitude']==.4],key=lambda r:(r['noise'],r['shape']))
    out=[]
    for r in rows:
        common={k:r[k]for k in ('observation_id','shape','noise','amplitude','task_id')};common['anchor_task_id']=common.pop('task_id')
        if r['noise']=='clean':
            for scale in (1.,.1,0.):out.append({**common,'experiment':'e2','arm':'continue_trf','reg_scale':scale,'task_id':f'e2_{r["observation_id"]}_reg{scale:g}'})
        for arm in ('continue_trf','haar_unfiltered','green_pullback','soft_jacobian'):out.append({**common,'experiment':'e3','arm':arm,'reg_scale':1.,'task_id':f'e3_{r["observation_id"]}_{arm}'})
    return out

def work(task):
    output=ROOT/'results'/task['experiment']/task['task_id'];start=time.time()
    if Path(str(output)+'.json').exists():return dict(task_id=task['task_id'],status='EXISTS')
    if task['experiment']=='e3'and task['arm']=='continue_trf'and task['noise']=='clean':
        source=ROOT/'results/e2'/f'e2_{task["observation_id"]}_reg1'
        if source.with_suffix('.json').exists():
            receipt=json.loads(source.with_suffix('.json').read_text());receipt.update(task=task,arrays=Path(str(output)+'.npz').name,reused_from=str(source.relative_to(ROOT)),incremental_compute_seconds=0.,reuse_semantics='Exact E2 full-regularization prediction and its charged original 45s budget; no second continuation run')
            shutil.copyfile(source.with_suffix('.npz'),Path(str(output)+'.npz'));dump(Path(str(output)+'.json'),receipt);return dict(task_id=task['task_id'],status='REUSED')
    try:
        with np.load(HIST/'predictions/endpoints.npz',allow_pickle=False)as a:
            i=list(a['task_id']).index(task['anchor_task_id']);anchor=a['m'][i].copy()
        result,arrays=solve(load_public(task['observation_id']),anchor,task)
        np.savez_compressed(Path(str(output)+'.npz'),**arrays)
        dump(Path(str(output)+'.json'),dict(task=task,status='OK'if result['status']!='FAILED'else'FAILED',result=result,arrays=Path(str(output)+'.npz').name,started_unix=start,closed_unix=time.time()))
        return dict(task_id=task['task_id'],status=result['status'],wall=result['wall_seconds'],top1=result['top1'])
    except Exception as exc:
        dump(Path(str(output)+'.json'),dict(task=task,status='FAILED',error=repr(exc),started_unix=start,closed_unix=time.time()));return dict(task_id=task['task_id'],status='FAILED',error=repr(exc))

def stationarity():
    from src.search import Search
    from threadpoolctl import threadpool_limits
    rows=[];arrays={}
    with threadpool_limits(limits=1):
        for task in tasks():
            if task['experiment']!='e2':continue
            path=ROOT/'results/e2'/(task['task_id']+'.npz')
            if not path.exists():continue
            with np.load(path,allow_pickle=False)as z:m=z['top1'].copy()
            run=Search(load_public(task['observation_id']),m,task);state,row=run.evaluate(m,'offline_stationarity_public_only')
            J=run.model.jacobian(state);gd=J.T@run.model.pack(state['y']-run.target)/run.scale2;gr=run.reg.evaluate(m)[1];g=gd+gr
            pg=g.copy();lower=m<=.02+1e-10;upper=m>=20.-1e-10;pg[lower&(g>0)]=0.;pg[upper&(g<0)]=0.
            v=np.ones(576);v[g<0]=20.-m[g<0];v[g>0]=m[g>0]-.02
            arrays[task['task_id']+'_gradient']=g;arrays[task['task_id']+'_data_gradient']=gd;arrays[task['task_id']+'_reg_gradient']=gr
            rows.append(dict(task=task,gradient_inf=float(np.linalg.norm(g,np.inf)),projected_gradient_inf=float(np.linalg.norm(pg,np.inf)),scipy_trf_scaled_optimality=float(np.linalg.norm(g*v,np.inf)),active_lower_count=int(lower.sum()),active_upper_count=int(upper.sum()),wall_seconds=run.elapsed(),costcounts=dict(run.model.counts)))
    dump(ROOT/'results/e2/ENDPOINT_STATIONARITY.json',dict(scope='Separate public-only diagnostic AFTER predictions; no prediction changed and costs charged separately',rows=rows,total_LU=sum(r['costcounts']['LU']for r in rows),total_complex_RHS=sum(r['costcounts']['actual_complex_RHS']for r in rows)))
    np.savez_compressed(ROOT/'results/e2/STATIONARITY_ARRAYS.npz',**arrays)

def main():
    p=argparse.ArgumentParser();p.add_argument('--stationarity',action='store_true');p.add_argument('--experiment',choices=['e2','e3','all'],default='all');p.add_argument('--workers',type=int,default=1);p.add_argument('--limit',type=int);p.add_argument('--arm');p.add_argument('--stopping-profile',choices=['frozen','tight'],default='frozen');args=p.parse_args()
    if args.stationarity:stationarity();return
    all_tasks=tasks();planned=[t for t in all_tasks if args.experiment=='all'or t['experiment']==args.experiment]
    if args.stopping_profile=='tight':
        planned=[{**t,'experiment':'e2','stopping_profile':'tight','task_id':t['task_id']+'_tight'}for t in planned if t['experiment']=='e2' and t['shape'] in ('digit2','digit7','ring','curved_bar')]
    if args.stopping_profile=='tight':dump(ROOT/'results/e2/STOPPING_MANIFEST.json',dict(tasks=planned,selection='4 prespecified shapes selected after observing early stopping and before new GT scores',profile=dict(ftol=1e-12,xtol=1e-12,gtol=1e-12),same_scientific_budget=True))
    if args.arm:planned=[t for t in planned if t['arm']==args.arm]
    for exp in ('e2','e3'):
        (ROOT/'results'/exp).mkdir(parents=True,exist_ok=True)
        dump(ROOT/'results'/exp/'MANIFEST.json',dict(schema='gn-search-r2-manifest-v1',tasks=[t for t in all_tasks if t['experiment']==exp],truth_in_solver=False,budget_per_task=dict(wall_seconds=45,max_LU=500,max_complex_RHS=30000)))
    planned=[t for t in planned if not (ROOT/'results'/t['experiment']/(t['task_id']+'.json')).exists()]
    if args.limit:planned=planned[:args.limit]
    with ProcessPoolExecutor(max_workers=args.workers,max_tasks_per_child=1)as executor:
        futures=[executor.submit(work,t)for t in planned]
        for f in as_completed(futures):
            row=f.result();print(json.dumps(row),flush=True)
            with (ROOT/'results'/args.experiment/'RUN_LOG.jsonl').open('a')as log:log.write(json.dumps(row)+'\n')
if __name__=='__main__':main()
