"""Public-only proposal verification; scoring is a separate explicit command."""
import argparse,json,os,sys,time
from pathlib import Path
for _k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS'):
    os.environ[_k]='1'
STUDY=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(STUDY.parents[1]))
import numpy as np
from research.research_gn_nn_preservation_r4.src.evaluate import evaluate_arm,evaluate_control,evaluate_independent_baseline,save_result

def read_npz(path):
    with np.load(path,allow_pickle=False) as z:return {k:z[k] for k in z.files}

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--runs',nargs='+',required=True)
    p.add_argument('--splits',nargs='+',default=['val'])
    p.add_argument('--partition-count',type=int,default=1,choices=[1,2])
    p.add_argument('--partition-index',type=int,default=0)
    p.add_argument('--obs-id',nargs='*')
    p.add_argument('--limit',type=int)
    a=p.parse_args()
    data=STUDY/'data'
    online=json.loads((data/'manifests/online.json').read_text())['records']
    records=[r for r in online if r['split'] in a.splits and (not a.obs_id or r['obs_id'] in a.obs_id)]
    if a.limit:records=records[:a.limit]
    ops=read_npz(data/'operators.npz')
    started=time.perf_counter(); entries=[]
    independent={'strict_gn45_anchor':'','weakreg45_anchor':'_weakreg1e-5','discrepancy45_anchor':'_discrepancy'}
    for index,r in enumerate(records):
        if index%a.partition_count!=a.partition_index:continue
        # Rotate method order deterministically; no GT information enters it.
        runs=a.runs[index%len(a.runs):]+a.runs[:index%len(a.runs)]
        for run in runs:
            obs=r['obs_id'];dest=STUDY/'evaluation'/run
            if (dest/f'{obs}.json').exists():continue
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
                                      'BLAS_threads':1,'method_order':'case-index rotation, public only'})
            save_result(dest,obs,receipt,arrays)
            item={'run':run,'obs_id':obs,'status':receipt['status'],'wall_seconds':receipt['wall_seconds'],
                  'selection_status':receipt['selection_status']}
            entries.append(item);print(json.dumps(item),flush=True)
    path=STUDY/'evaluation'/f'BATCH_{"_".join(a.splits)}_P{a.partition_index}_{time.time_ns()}.json'
    path.parent.mkdir(exist_ok=True)
    path.write_text(json.dumps({'truth_access':'NONE','wall_seconds':time.perf_counter()-started,
                               'entries':entries,'runs':a.runs,'splits':a.splits},indent=2)+'\n')
    print(path,flush=True)

if __name__=='__main__':main()
