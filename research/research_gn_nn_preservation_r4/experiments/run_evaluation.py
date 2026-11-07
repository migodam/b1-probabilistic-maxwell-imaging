"""Public-only proposal verification; scoring is a separate explicit command."""
import argparse,json,os,sys,time,platform
from pathlib import Path
for _k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS'):
    os.environ[_k]='1'
STUDY=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(STUDY.parents[1]))
import numpy as np
from research.research_gn_nn_preservation_r4.src.evaluate import evaluate_arm,evaluate_control,evaluate_independent_baseline,save_result
from research.research_gn_nn_preservation_r4.experiments.run_test_partition import load_contract,require_host,DEADLINE,INDEPENDENT,pair_host

def read_npz(path):
    with np.load(path,allow_pickle=False) as z:return {k:z[k] for k in z.files}

def main():
    setup_started=time.perf_counter()
    p=argparse.ArgumentParser()
    p.add_argument('--runs',nargs='+',required=True)
    p.add_argument('--splits',nargs='+',default=['val'])
    p.add_argument('--partition-count',type=int,default=1,choices=[1,2,4])
    p.add_argument('--partition-index',type=int,default=0)
    p.add_argument('--execution-partition',choices=['local','remote'])
    p.add_argument('--obs-id',nargs='*')
    p.add_argument('--limit',type=int)
    a=p.parse_args()
    if not 0<=a.partition_index<a.partition_count:p.error('partition-index out of range')
    data=STUDY/'data'
    online=json.loads((data/'manifests/online.json').read_text())['records']
    global_indices={r['obs_id']:i for i,r in enumerate(online)}
    allowed_test=None
    if 'test' in a.splits:
        if not a.execution_partition:p.error('TEST needs the declared --execution-partition')
        freeze,final,allocation,_=load_contract();require_host(a.execution_partition,1)
        if any(run not in freeze['test_runs'] for run in a.runs):raise ValueError('Unregistered TEST run')
        allowed_test=set(allocation[a.execution_partition]['obs_ids'])
        if a.obs_id and any(obs not in allowed_test for obs in a.obs_id):raise ValueError('TEST observation assigned to another host')
    records=[r for r in online if r['split'] in a.splits and (not a.obs_id or r['obs_id'] in a.obs_id) and (r['split']!='test' or r['obs_id'] in allowed_test)]
    if a.limit:records=records[:a.limit]
    ops=read_npz(data/'operators.npz')
    setup_wall_seconds=time.perf_counter()-setup_started
    started=time.perf_counter(); entries=[]
    independent=INDEPENDENT
    for index,r in enumerate(records):
        if index%a.partition_count!=a.partition_index:continue
        # Rotate method order deterministically; no GT information enters it.
        order_index=global_indices[r['obs_id']] if r['split']=='test' else index
        runs=a.runs[order_index%len(a.runs):]+a.runs[:order_index%len(a.runs)]
        for run in runs:
            obs=r['obs_id'];dest=STUDY/'evaluation'/run
            if r['split']=='test' and time.time()>=DEADLINE:break
            if (dest/f'{obs}.json').exists():
                if r['split']=='test' and not (dest/f'{obs}.npz').exists():raise RuntimeError('Preserve incomplete existing TEST pair; interrupted costs unknown')
                if r['split']=='test' and pair_host(dest/obs)!=a.execution_partition:raise RuntimeError('Existing TEST case lacks matching assigned-host provenance; preserved')
                continue
            if r['split']=='test' and (dest/f'{obs}.npz').exists():raise RuntimeError('Preserve incomplete existing TEST NPZ; interrupted costs unknown')
            anchor_path=STUDY/'anchors'/f'{obs}.npz'
            pred=STUDY/'predictions'/run/f'{obs}.npz'
            if not anchor_path.exists() or (run!='continue_trf' and run not in independent and not pred.exists()):continue
            anchor=read_npz(anchor_path); public={**ops,**read_npz(data/r['public_path'])}
            noise='clean' if float(r['noise'])==0 else 'noisy'
            if run in independent:
                base_path=STUDY/'anchors'/f'{obs}{independent[run]}.npz'
                if not base_path.exists():continue
                base_receipt=json.loads(base_path.with_suffix('.json').read_text())
                receipt,arrays=evaluate_independent_baseline(public,anchor,read_npz(base_path),base_receipt,run,noise)
            elif run=='continue_trf':
                receipt,arrays=evaluate_control(public,anchor,control=run,noise=noise)
            else:
                info=json.loads(pred.with_suffix('.json').read_text())
                receipt,arrays=evaluate_arm(public,anchor,read_npz(pred),arm=info.get('arm',run),noise=noise,
                                           preparation_cost=info['preparation_cost'])
            receipt.update(obs_id=obs,parent_id=r['parent_id'],split=r['split'],run=run,
                           prediction_path=str(pred),anchor_path=str(anchor_path),
                           execution={'partition_index':a.partition_index,'partition_count':a.partition_count,
                                      'BLAS_threads':1,'platform':platform.platform(),'machine':platform.machine(),'processor':platform.processor(),'pid':os.getpid(),'execution_partition':a.execution_partition,'global_online_index':global_indices[obs],'method_rotation_index':order_index,'method_order':'complete online-manifest index rotation for TEST; legacy filtered-case index for VAL, public only'})
            save_result(dest,obs,receipt,arrays)
            item={'run':run,'obs_id':obs,'status':receipt['status'],'wall_seconds':receipt['wall_seconds'],
                  'selection_status':receipt['selection_status']}
            entries.append(item);print(json.dumps(item),flush=True)
    path=STUDY/'evaluation'/f'BATCH_{"_".join(a.splits)}_P{a.partition_index}_{time.time_ns()}.json'
    path.parent.mkdir(exist_ok=True)
    path.write_text(json.dumps({'truth_access':'NONE','wall_seconds':time.perf_counter()-started,
                               'setup_wall_seconds':setup_wall_seconds,'entries':entries,'runs':a.runs,'splits':a.splits,'execution_partition':a.execution_partition,'platform':platform.platform(),'machine':platform.machine(),'absolute_deadline_utc':'2026-11-07T10:44:27Z','interrupted_cost_policy':'unknown interrupted costs are not zero'},indent=2)+'\n')
    print(path,flush=True)

if __name__=='__main__':main()
