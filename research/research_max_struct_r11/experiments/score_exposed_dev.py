"""Offline score of already exposed R10 objects; no inverse or online selection."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[k]='1'
from pathlib import Path
import argparse, csv, json, sys, time
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT.parents[1]))
from research.research_max_struct_r10.experiments.score_twod_comparison import material_metrics

def write_csv(path,rows):
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--historical-root',type=Path,required=True)
    ap.add_argument('--manifest',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args();started=time.process_time();wall=time.perf_counter()
    tasks=json.loads(a.manifest.read_text())['tasks'];rows=[];checkpoints=[]
    # Verify all requested public outcomes exist before this process reads any exposed GT.
    seals=[]
    for t in tasks:
        d=ROOT.parents[1]/t['output_dir'];r=json.loads((d/'RECEIPT.json').read_text())
        if r['status']!='SEALED_PUBLIC_PREDICTION':raise RuntimeError(t['task_id']+' is not sealed')
        with np.load(d/'R11_CAP_PREDICTION.npz',allow_pickle=False) as z: m=z['m'].copy()
        seals.append((t,d,r,m))
    a.output.mkdir(parents=True,exist_ok=True)
    (a.output/'SEALS_BEFORE_SCORE.json').write_text(json.dumps(dict(status='ALL_REQUESTED_PUBLIC_PREDICTIONS_SEALED',
        task_ids=[t['task_id'] for t in tasks],truth_access='NONE_YET',dataset='EXPOSED_R10_DEVELOPMENT'),indent=2)+'\n')
    for t,d,r,m in seals:
        case=t['case_id'];cohort='conditional_confirmation' if 'conditional' in case else 'twod'
        with np.load(a.historical_root/'sealed_truth'/cohort/(case+'.npz'),allow_pickle=False) as z: gt=z['m'].copy()
        amp=float(gt.max()-.1);high=float(gt.max())
        label=t.get('material_readout',{}).get('label',t['arm'])
        p=json.loads((d/'R11_CAP_PREDICTION.json').read_text())
        row=dict(case_id=case,method=label,task_id=t['task_id'],**material_metrics(m,gt,amp,high),
            original_Phi=p['objective'],relative_residual=p['physical_relative_residual'],
            corrector_stop=(r.get('corrector') or {}).get('stop'),
            corrector_iterations=(r.get('corrector') or {}).get('iterations'),
            truth_access='OFFLINE_EXPOSED_DEV_AFTER_ALL_REQUESTED_SEALS')
        rows.append(row)
        cp=json.loads((d/'R11_CPU_CHECKPOINTS.json').read_text())
        with np.load(d/'R11_CHECKPOINT_MATERIALS.npz',allow_pickle=False) as z: mm=z['m'].copy()
        for c in cp:
            entry=dict(case_id=case,method=label,**c)
            if c.get('material_index') is not None:entry.update(material_metrics(mm[c['material_index']],gt,amp,high))
            checkpoints.append(entry)
    write_csv(a.output/'SCORES.csv',rows);write_csv(a.output/'CPU_CHECKPOINT_SCORES.csv',checkpoints)
    (a.output/'RECEIPT.json').write_text(json.dumps(dict(status='COMPLETE',objects=len(rows),physical_forward=0,
        CPU_seconds=time.process_time()-started,wall_seconds=time.perf_counter()-wall,
        scope='EXPOSED_R10_DEVELOPMENT; no confirmation or equal hardware inference from cross-platform DEV'),indent=2)+'\n')
    print(json.dumps([{k:r[k] for k in ('case_id','method','amplitude_NRMSE','fixed_support_IoU','original_Phi')} for r in rows],indent=2))

if __name__=='__main__':main()
