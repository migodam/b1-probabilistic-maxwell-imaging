"""Reconstruct from measurements under an explicit, separate deployment budget.

Never resumes or resets the closed R4 research clock. No truth/digit ID input.
"""
from pathlib import Path
import os,sys,json,argparse,time,datetime
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='1'
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'vendor'))
from gsda.accounting import BudgetBook,dump
from gsda.guard import OnlineGuard
from gsda.fem_model import FEMModel

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--measurement',type=Path,required=True);ap.add_argument('--config',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--options',type=Path,required=True);ap.add_argument('--method',choices=['B0','B1','B2','B3','B4'],default='B0');ap.add_argument('--deployment-budget-seconds',type=int,required=True);a=ap.parse_args()
 if not 1<=a.deployment_budget_seconds<=10800:raise ValueError('explicit deployment budget must be 1..10800 s')
 out=a.output.resolve()
 if out.exists() and any(out.iterdir()):raise FileExistsError('new deployment requires empty output')
 out.mkdir(parents=True,exist_ok=True);acct=out/'accounting';acct.mkdir();start=time.time();utc=lambda t:datetime.datetime.fromtimestamp(t,datetime.timezone.utc).isoformat()
 dump(acct/'RUN_CLOCK.json',dict(start_utc=utc(start),deadline_utc=utc(start+a.deployment_budget_seconds),no_new_heavy_after_utc=utc(start+a.deployment_budget_seconds),single_job_rss_bytes=4294967296,mode='EXPLICIT_NEW_DEPLOYMENT_NOT_R4_CAMPAIGN_EXTENSION'))
 book=BudgetBook(root=acct,job='deployment_'+a.method);guard=OnlineGuard(ROOT).install();status='FAILED'
 try:
  c=json.loads(a.config.read_text());opts=json.loads(a.options.read_text())
  with np.load(a.measurement,allow_pickle=False) as f:
   if f.files!=['y']:raise ValueError('online measurement contains only y')
   y=f['y'].copy()
  model=FEMModel(c,book,'DEPLOY_'+a.method)
  if a.method=='B0':
   from gsda.fwi import invert
   result=invert(model,y,out,opts.get('fwi_options',opts))
  else:
   from gsda.currents import invert
   co=dict(opts.get('csi_options',{}))
   if a.method in ('B3','B4'):co['spectral_preconditioner']=True
   result=invert(model,y,out,co,opts.get('fwi_options'),absorption=a.method in ('B2','B3'),risk_controlled=a.method=='B4')
  from gsda.artifacts import save_outputs
  save_outputs(out,c)
  status=result['status'];dump(out/'reconstruction.json',result)
 finally:
  guard.save(out/'FILE_ACCESS_SUMMARY.json');dump(out/'runtime_cost.json',book.receipt(status=status))
if __name__=='__main__':main()
