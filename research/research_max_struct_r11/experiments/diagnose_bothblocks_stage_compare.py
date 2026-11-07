"""Saved-state stage/corrector comparison; complete-trajectory intervention, zero FW."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import sys,json,time,csv
import numpy as np
ROOT=Path(__file__).absolute().parents[1];sys.path.insert(0,str(ROOT.parents[1]))
from research.research_max_struct_r10.experiments.score_twod_comparison import material_metrics
OUT=ROOT/'delegated/diagnostics';R10=ROOT.parents[2]/'r10_worktree/research/research_max_struct_r10'

def dump(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False,default=lambda v:v.item() if isinstance(v,np.generic) else v.tolist())+'\n')
def main():
 wall,cpu=time.perf_counter(),time.process_time();rows=[];starts=[];meta={r['case_id']:r for r in json.loads((R10/'configs/twod_generator.json').read_text())['cases']}
 for cid in ('r10_twod_0001','r10_twod_0004','r10_twod_0022'):
  with np.load(R10/'sealed_truth/twod'/f'{cid}.npz') as z:gt=z['m']
  amp=float(gt.max()-.1)
  for method,d in [('R10_fullGreen19',R10/'RESULTS/twod'/f'{cid}__matrix_green19__noise_ball_99'),('R11_fullGreenBothBlocks',ROOT/'RESULTS/exposed_dev'/f'{cid}__fullGreenBothBlocks__480CPU')]:
   rec=json.loads((d/'RECEIPT.json').read_text());schedule=json.loads((d/'NUMERICAL_SCHEDULE.json').read_text())
   selected=None
   with np.load(d/'selected_before_corrector.npz') as z:selected=z['m'].copy()
   for name in ['stage_00_physical','stage_01_physical','stage_02_physical','stage_03_physical','selected_before_corrector']+[x.stem for x in sorted(d.glob('corrector_accepted_*.npz'))]+['prediction']:
    npz=d/(name+'.npz');jp=d/(name+'.json')
    if not npz.exists():rows.append(dict(case_id=cid,method=method,state=name,status='MISSING_SAVED_STATE'));continue
    info=json.loads(jp.read_text())
    with np.load(npz) as z:
     metrics=material_metrics(z['m'],gt,amp,meta[cid]['high_chi']);pg=float(np.max(abs(z['m']-np.clip(z['m']-z['gradient'],.02,20.)))) if 'gradient' in z.files else None
     step_from_selected=float(np.linalg.norm(z['m']-selected));m=z['m'].copy()
    stage=int(name[6:8]) if name.startswith('stage_') else rec.get('auxiliary_selected_stage') if name=='selected_before_corrector' else None
    rows.append(dict(case_id=cid,method=method,state=name,status='AVAILABLE',stage=stage,
      mf=schedule['stages'][stage]['mf'] if stage is not None and stage>=0 else None,
      original_Phi=info.get('objective',info.get('physical_objective')),amplitude_NRMSE=metrics['amplitude_NRMSE'],fixed_support_IoU=metrics['fixed_support_IoU'],material_success=metrics['material_success'],saved_PG_inf=pg,
      step_from_corrector_start=step_from_selected,selected_auxiliary_stage=rec.get('auxiliary_selected_stage'),corrector_stop=rec.get('corrector',{}).get('stop'),corrector_iterations=rec.get('corrector',{}).get('iterations'),
      state_CPU_seconds=info.get('cpu_seconds'),scope='Complete-trajectory intervention, not same-state counterfactual',truth_access='POSTSEAL_EXPOSED_GT_ONLY'))
   with np.load(d/'stage_01_joint.npz') as z:starts.append(dict(case_id=cid,method=method,stage1_m=z['m'].copy()))
 dump(OUT/'BOTHBLOCKS_STAGE_CORRECTOR_COMPARISON.json',rows)
 keys=list(dict.fromkeys(k for x in rows for k in x))
 with (OUT/'BOTHBLOCKS_STAGE_CORRECTOR_COMPARISON.csv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
 paired=[]
 for cid in ('r10_twod_0001','r10_twod_0004','r10_twod_0022'):
  for name in ['stage_00_physical','stage_01_physical','stage_02_physical','stage_03_physical','selected_before_corrector','prediction']:
   old=next(x for x in rows if x['case_id']==cid and x['method']=='R10_fullGreen19' and x['state']==name);new=next(x for x in rows if x['case_id']==cid and x['method']=='R11_fullGreenBothBlocks' and x['state']==name)
   paired.append(dict(case_id=cid,state=name,R10_amplitude_NRMSE=old['amplitude_NRMSE'],R11_amplitude_NRMSE=new['amplitude_NRMSE'],R10_fixed_IoU=old['fixed_support_IoU'],R11_fixed_IoU=new['fixed_support_IoU'],R10_Phi=old['original_Phi'],R11_Phi=new['original_Phi'],R10_selected_stage=old['selected_auxiliary_stage'],R11_selected_stage=new['selected_auxiliary_stage'],scope='Separate full alternating trajectories; identical stage index does not imply identical J or m'))
 dump(OUT/'BOTHBLOCKS_PAIRED_STAGE_SUMMARY.json',paired)
 dump(OUT/'BOTHBLOCKS_STAGE_COMPARISON_RECEIPT.json',dict(status='COMPLETE',new_forward=0,new_adjoint=0,new_optimization=0,states=len(rows),CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,truth_access='POSTSEAL_EXPOSED_GT'))
 for x in paired:
  if x['case_id']=='r10_twod_0022':print(x)
if __name__=='__main__':main()
