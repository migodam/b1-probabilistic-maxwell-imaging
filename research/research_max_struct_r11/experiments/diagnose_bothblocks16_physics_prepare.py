"""Public-cost/common-LU selection seal, after80 terminal and before any GT scoring."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import json,time
import numpy as np
ROOT=Path(__file__).absolute().parents[1];REPO=ROOT.parents[1];OUT=ROOT/'delegated/diagnostics/bothblocks_confirmation'

def main():
 wall,cpu=time.perf_counter(),time.process_time();seal=json.loads((OUT/'ALL80_TERMINAL_PREDICTION_SEAL.json').read_text())
 if len(seal['entries'])!=80 or seal['status']!='ALL80_TERMINAL_BEFORE_ANY_CONFIRMATION_GT_SCORING':raise RuntimeError('All80 terminal required')
 rows=[];mats=[]
 for cid in sorted({e['task']['case_id'] for e in seal['entries']}):
  entries=[e for e in seal['entries'] if e['task']['case_id']==cid];attained={}
  for e in entries:
   d=REPO/e['task']['output_dir'];p=d/'R11_PAID_PHYSICAL_STATES.json';paid=json.loads(p.read_text()) if p.exists() else []
   eligible=[x for x in paid if x.get('eligible_at_cap') and x['CPU_seconds']<=480]
   attained[e['task']['method']]=max((x['physical_model_cost'].get('LU',0) for x in eligible),default=0)
  cap=min(attained.values()) if len(attained)==5 else 0
  for e in entries:
   t=e['task'];d=REPO/t['output_dir'];p=d/'R11_INCUMBENT_TRAJECTORY.npz'
   row=dict(case_id=cid,method=t['method'],common_physical_LU_cap=cap,arm_attained_eligible_LU=attained[t['method']],all5_attained_eligible_LU=attained,selection_rule='Min originalPhi among stored CPU-eligible paid incumbents with LU<=min(all5 eligible attained maxLU)',truth_access='NONE_SELECTION',material_index=None)
   if p.exists() and cap>0:
    with np.load(p) as z:
     allowed=np.flatnonzero((z['LU']<=cap)&(z['CPU_seconds']<=480))
     if len(allowed):
      i=int(allowed[np.argmin(z['original_Phi'][allowed])]);row.update(status='SEALED_COMMON_LU_PAID_STATE',source=str(p),source_index=i,original_Phi=float(z['original_Phi'][i]),selected_CPU_seconds=float(z['CPU_seconds'][i]),selected_LU=int(z['LU'][i]),selected_complex_RHS=int(z['complex_RHS'][i]),material_index=len(mats));mats.append(z['m'][i].copy())
     else:row['status']='NO_ELIGIBLE_STORED_INCUMBENT'
   else:row['status']='NO_ALL5_ELIGIBLE_COMMON_CAP'
   rows.append(row)
 np.savez_compressed(OUT/'COMMON_PHYSICS_CAP_SEALED_MATERIALS.npz',m=np.asarray(mats))
 (OUT/'COMMON_PHYSICS_CAP_SELECTION_BEFORE_GT.json').write_text(json.dumps(dict(status='SEALED_PUBLIC_COST_AND_PHI_BEFORE_GT',rows=rows,definition='Per-object min of five arms max eligible paid physical LU; actual terminal/overrun spend separate',CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,new_forward=0,truth_access='NONE'),indent=2)+'\n')
 print(json.dumps(dict(status='SEALED_BEFORE_GT',states=len(mats),objects=16)))
if __name__=='__main__':main()
