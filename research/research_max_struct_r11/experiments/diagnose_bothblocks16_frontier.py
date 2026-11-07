"""Postseal saved-incumbent CPU/LU/RHS material frontier, no new FW/optimization."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import sys,json,time,csv
from collections import defaultdict
import numpy as np
ROOT=Path(__file__).absolute().parents[1];REPO=ROOT.parents[1];OUT=ROOT/'delegated/diagnostics/bothblocks_confirmation';sys.path.insert(0,str(REPO))
from research.research_max_struct_r10.experiments.score_twod_comparison import material_metrics

def main():
 wall,cpu=time.perf_counter(),time.process_time();common=json.loads((OUT/'COMMON_PHYSICS_CAP_SELECTION_BEFORE_GT.json').read_text());seal=json.loads((OUT/'ALL80_TERMINAL_PREDICTION_SEAL.json').read_text())
 if common['status']!='SEALED_PUBLIC_COST_AND_PHI_BEFORE_GT' or len(seal['entries'])!=80 or not (OUT/'SCORING_RECEIPT.json').exists():raise RuntimeError('Main root-authorized scoring plus all80/common-Phi seals required')
 config=json.loads((ROOT/'configs/twod_bothblocks16_generator.json').read_text());meta={x['case_id']:x for x in config['cases']};truth={}
 for cid in meta:
  with np.load(ROOT/'sealed_truth/bothblocks16'/f'{cid}.npz') as z:truth[cid]=z['m'].copy()
 commons=[]
 with np.load(OUT/'COMMON_PHYSICS_CAP_SEALED_MATERIALS.npz') as z:cm=z['m']
 for x in common['rows']:
  r=dict(x);cid=x['case_id'];r.update(family=meta[cid]['family'],regime=meta[cid]['regime'],material_success=False)
  if x['material_index'] is not None:r.update(**material_metrics(cm[x['material_index']],truth[cid],float(truth[cid].max()-.1),meta[cid]['high_chi']))
  commons.append(r)
 (OUT/'COMMON_PHYSICS_CAP_SCORES.json').write_text(json.dumps(commons,indent=2)+'\n')
 with (OUT/'COMMON_PHYSICS_CAP_SCORES.csv').open('w') as f:
  keys=list(dict.fromkeys(k for r in commons for k in r));w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(commons)
 fields=['case_id','method','point_kind','source_index','CPU_seconds','LU','complex_RHS','original_Phi','amplitude_NRMSE','fixed_support_IoU','material_success','CPU_eligibility_scope'];count=0;terminals=[]
 with (OUT/'PAID_INCUMBENT_FRONTIER.csv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
  for e in seal['entries']:
   t=e['task'];cid=t['case_id'];d=REPO/t['output_dir'];p=d/'R11_INCUMBENT_TRAJECTORY.npz';gt=truth[cid];amp=float(gt.max()-.1)
   if not p.exists():continue
   with np.load(p) as z:
    m=z['m'];nr=np.sqrt(np.mean((m-gt)**2,axis=1))/amp;support=m>1.;ts=gt>1.;union=np.sum(support|ts,axis=1);inter=np.sum(support&ts,axis=1);iou=np.divide(inter,union,out=np.ones_like(inter,dtype=float),where=union>0);legal=(m.min(axis=1)>=.02-1e-12)&(m.max(axis=1)<=20+1e-12)
    for i in range(len(m)):
     row=dict(case_id=cid,method=t['method'],point_kind='PAID_PHI_INCUMBENT',source_index=i,CPU_seconds=float(z['CPU_seconds'][i]),LU=int(z['LU'][i]),complex_RHS=int(z['complex_RHS'][i]),original_Phi=float(z['original_Phi'][i]),amplitude_NRMSE=float(nr[i]),fixed_support_IoU=float(iou[i]),material_success=bool(nr[i]<=.15 and iou[i]>=.7 and legal[i]),CPU_eligibility_scope='Original whole-child CPU<=480; exact saved state, no interpolation');w.writerow(row);count+=1
   rec=json.loads((d/'RECEIPT.json').read_text());pc=rec.get('physical_model_cost',{});j=e['terminal_job'];end=dict(case_id=cid,method=t['method'],point_kind='ACTUAL_TERMINAL_SPEND_SEPARATE',CPU_seconds=j['child_process_CPU_seconds'],LU=pc.get('LU'),complex_RHS=pc.get('actual_complex_RHS'),original_Phi=e['cap_original_Phi'],CPU_eligibility_scope='Actual entire child spend including paid overruns; prediction remains best eligible cap state')
   if e['cap_prediction_available']:
    with np.load(e['prediction_path']) as z:mt=material_metrics(z['m'],gt,amp,meta[cid]['high_chi'])
    end.update(amplitude_NRMSE=mt['amplitude_NRMSE'],fixed_support_IoU=mt['fixed_support_IoU'],material_success=mt['material_success'])
   terminals.append(end)
 (OUT/'ACTUAL_TERMINAL_SPEND_POINTS.json').write_text(json.dumps(terminals,indent=2)+'\n')
 with (OUT/'ACTUAL_TERMINAL_SPEND_POINTS.csv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(terminals)
 summary=[]
 for method in sorted({r['method'] for r in commons}):
  rr=[r for r in commons if r['method']==method];v=[r for r in rr if r.get('amplitude_NRMSE') is not None];summary.append(dict(method=method,objects=16,available=len(v),material_passes=sum(r['material_success'] for r in rr),mean_amplitude_NRMSE=float(np.mean([r['amplitude_NRMSE'] for r in v])) if v else None,median_amplitude_NRMSE=float(np.median([r['amplitude_NRMSE'] for r in v])) if v else None,scope='Secondary per-object shared LU cap; does not replace frozen primary CPU-cap endpoint'))
 (OUT/'COMMON_PHYSICS_CAP_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
 (OUT/'FRONTIER_SCORING_RECEIPT.json').write_text(json.dumps(dict(status='COMPLETE',paid_frontier_rows=count,terminal_actual_points=len(terminals),common_cap_states=len(commons),new_forward=0,new_adjoint=0,new_optimization=0,CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,truth_access='After root-authorized main scoring and all80/common selections sealed'),indent=2)+'\n')
 print(json.dumps(dict(frontier_rows=count,terminal_points=len(terminals),common_summary=summary)))
if __name__=='__main__':main()
