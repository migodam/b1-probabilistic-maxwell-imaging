"""Read-only sealed-public-prediction replay comparison; no GT or rerun."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
import argparse,json,time
from pathlib import Path
import numpy as np

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--reference',type=Path,required=True);ap.add_argument('--replay',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();t,c=time.perf_counter(),time.process_time()
 with np.load(a.reference,allow_pickle=False) as x,np.load(a.replay,allow_pickle=False) as y:
  rows={}
  for key in ('m','j','Q','Y_prediction'):
   d=np.asarray(y[key])-np.asarray(x[key]);den=max(float(np.linalg.norm(x[key])),1e-300)
   rows[key]=dict(shape=list(x[key].shape),max_abs_difference=float(np.max(abs(d))),relative_norm_difference=float(np.linalg.norm(d))/den,exact_array_equal=bool(np.array_equal(x[key],y[key])))
 receipt=dict(status='EXACT_REPLAY' if all(r['exact_array_equal'] for r in rows.values()) else 'NUMERICAL_DIFFERENCE_REQUIRES_ROOT_REVIEW',reference=str(a.reference),replay=str(a.replay),fields=rows,truth_access='NONE',parameter_selection='NONE',physics_runs=0,CPU_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-t)
 a.output.write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))
if __name__=='__main__':main()
