"""Paid post-seal prediction; public config and selected material only, no truth."""
from pathlib import Path
import os,sys,json,argparse,time
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='1'
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from gsda.fem_model import FEMModel
from gsda.accounting import BudgetBook,dump
from gsda.guard import OnlineGuard

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--methods',nargs='+',default=['B0_FULL_FWI','B4_GSDA']);ap.add_argument('--cases',nargs='+',default=[f'c{i:03d}' for i in range(1,10)]);a=ap.parse_args()
 book=BudgetBook(job='heldout_predictions');guard=OnlineGuard(ROOT).install();status='FAILED';rows=[]
 try:
  c=json.loads((ROOT/'inputs/online/c001/config.json').read_text());c['frequencies']=c['heldout_frequencies'];model=FEMModel(c,book,'R4_HELDOUT_POSTSEAL')
  for method in a.methods:
   for case in a.cases:
    d=ROOT/'results'/method/case
    if not (d/'ONLINE_SEALED.json').exists():continue
    if (d/'heldout_prediction.npz').exists():continue
    with np.load(d/'endpoint.npz') as f:z=f['z'].copy()
    before=book.snapshot(pid=os.getpid());cpu=time.process_time();wall=time.perf_counter();pred=model.predict(z);after=book.snapshot(pid=os.getpid())
    np.savez_compressed(d/'heldout_prediction.npz',data=pred.data,z=z,frequencies=c['frequencies'])
    row=dict(case=case,method=method,role='post-seal heldout prediction, not online selection or new information',frequencies=c['frequencies'],full_fe=after['full_fe']-before['full_fe'],rhs=after['rhs']-before['rhs'],cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,diagnostics=pred.diagnostics)
    dump(d/'HELDOUT_RECEIPT.json',row);rows.append(row)
  status='COMPLETED'
 except BaseException as e:book.failure('HELDOUT','predict',e);raise
 finally:
  record=dict(rows=rows,receipt=book.receipt(status=status),truth_access=False)
  dump(ROOT/'results/heldout_predictions'/('HELDOUT_PREDICTIONS_RECEIPT_'+str(os.getpid())+'.json'),record)
  dump(ROOT/'results/heldout_predictions/HELDOUT_PREDICTIONS_RECEIPT.json',record)
  guard.save(ROOT/'results/heldout_predictions'/('HELDOUT_ACCESS_SUMMARY_'+str(os.getpid())+'.json'))
  guard.save(ROOT/'results/heldout_predictions/HELDOUT_ACCESS_SUMMARY.json')
if __name__=='__main__':main()
