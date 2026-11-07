"""Postseal scoring of exactly three contractedLagged DEV runs; no forward calls."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import sys,json,time,csv
from collections import Counter
import numpy as np
from scipy.stats import chi2
ROOT=Path(__file__).absolute().parents[1];sys.path.insert(0,str(ROOT.parents[1]))
from research.research_max_struct_r10.experiments.score_twod_comparison import material_metrics
OUT=ROOT/'delegated/diagnostics';R10=ROOT.parents[2]/'r10_worktree/research/research_max_struct_r10'

def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False,default=lambda v:v.item() if isinstance(v,np.generic) else v.tolist())+'\n')
def main():
    wall,cpu=time.perf_counter(),time.process_time();seal=json.loads((OUT/'CONTRACTED_ALL3_SEALED.json').read_text())
    if seal['status']!='ALL3_SEALED_BEFORE_OFFLINE_SCORING':raise RuntimeError('Three seals required')
    tasks=json.loads((ROOT/'configs/twod_dev_contracted_tasks.json').read_text())['tasks'];rows=[];checkpoints=[];cost=Counter()
    metas={r['case_id']:r for r in json.loads((R10/'configs/twod_generator.json').read_text())['cases']}
    for task in tasks:
      d=ROOT.parents[1]/task['output_dir'];rec=json.loads((d/'RECEIPT.json').read_text())
      if rec['status']!='SEALED_PUBLIC_PREDICTION' or rec['truth_access']!='NONE':raise RuntimeError('Invalid online seal')
      cid=task['case_id'];meta=metas[cid]
      with np.load(R10/'sealed_truth/twod'/f'{cid}.npz') as z:truth=z['m']
      with np.load(R10/'data_public/twod'/f'{cid}.npz') as z:Y=z['Y'];sigma=float(z['sigma'])
      amp=float(truth.max()-.1)
      with np.load(d/'R11_CAP_PREDICTION.npz') as z:
        metric=material_metrics(z['m'],truth,amp,meta['high_chi']);pred=z['Y_prediction']
      res=float(np.linalg.norm(pred-Y));cn=float(res/sigma);cn99=float(np.sqrt(chi2.ppf(.99,2*Y.size)/2))
      row=dict(case_id=cid,method='contractedLagged',status=rec['status'],experiment_scope='EXPOSED_DEV_NOT_CONFIRMATION',truth_access='OFFLINE_POST_ALL3_SEAL',
        whole_child_CPU_seconds=rec['whole_child_CPU_seconds'],CPU_scope='Absolute process_time sampled near child end; final external child rusage was not captured, serialization/exit tail unavailable',
        original_Phi=rec['cap_physical_objective'],relative_residual=float(res/np.linalg.norm(Y)),whitened_residual=cn,CN99_threshold=cn99,data_CN99_pass=bool(cn<=cn99),
        selected_stage=rec.get('auxiliary_selected_stage'),corrector_stop=rec.get('corrector',{}).get('stop'),corrector_iterations=rec.get('corrector',{}).get('iterations'),
        physical_cost=rec.get('physical_model_cost'),auxiliary_cost=rec.get('auxiliary_execution_receipt',{}).get('operation_counts'),**metric)
      rows.append(row);cost.update(rec.get('physical_model_cost',{}))
      cp=json.loads((d/'R11_CPU_CHECKPOINTS.json').read_text())
      with np.load(d/'R11_CHECKPOINT_MATERIALS.npz') as z:cm=z['m']
      for point in cp:
        if 'material_index' not in point:checkpoints.append(dict(case_id=cid,**point));continue
        mm=material_metrics(cm[point['material_index']],truth,amp,meta['high_chi'])
        checkpoints.append(dict(case_id=cid,**point,**mm))
    dump(OUT/'CONTRACTED_DEV_SCORES.json',rows);dump(OUT/'CONTRACTED_DEV_CHECKPOINT_SCORES.json',checkpoints)
    dump(OUT/'CONTRACTED_DEV_COST_LEDGER.json',dict(status='COMPLETE',executed_tasks=3,failed_tasks=0,retries=0,whole_child_CPU_seconds_sum=sum(r['whole_child_CPU_seconds'] for r in rows),physical_cost=dict(cost),
      physical_forward_during_scoring=0,physical_adjoint_during_scoring=0,truth_access='OFFLINE_AFTER_ALL3_SEALED',scoring_CPU_seconds=time.process_time()-cpu,scoring_wall_seconds=time.perf_counter()-wall,
      resource_note='Local Mac CPU, two lanes then one; BLAS threads=1 as inherited runner; no remote process changed',
      wiring_check='CONTRACTED_WIRING_CHECK.json; zero FW; process CPU not recorded for that short check',CPU_limit_note='480 CPU cap did not trigger; runs ended under mature schedule/stop criteria',CPU_precision_note='Near-end child process_time recorded; exact external exit CPU unavailable; no repeat executed to recover it'))
    cols=['case_id','whole_child_CPU_seconds','amplitude_NRMSE','fixed_support_IoU','material_success','original_Phi','data_CN99_pass','whitened_residual','corrector_stop','corrector_iterations']
    with (OUT/'CONTRACTED_DEV_SCORES.csv').open('w') as f:
      w=csv.DictWriter(f,fieldnames=cols,extrasaction='ignore');w.writeheader();w.writerows(rows)
    for r in rows:print({k:r[k] for k in cols})
if __name__=='__main__':main()
