from pathlib import Path
import os,sys,argparse,json,time
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='1'
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from gsda.accounting import BudgetBook,dump
from gsda.fem_model import FEMModel
from gsda.currents import CurrentGreen
from gsda.spectral import build_modes,compare_state

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--cases',nargs='+',default=['c010','c011','c012']);ap.add_argument('--rank',type=int,default=48);ap.add_argument('--frequency',type=float,default=1.);a=ap.parse_args()
 out=ROOT/'results/spectral';out.mkdir(exist_ok=True,parents=True);book=BudgetBook(job='spectral_pilot');status='FAILED';rows=[]
 try:
  c=json.loads((ROOT/'configs/base.json').read_text());c['frequencies']=[a.frequency]
  model=FEMModel(c,book,'R4_GD_SPECTRAL');g=CurrentGreen(model,0)
  modes=build_modes(g,a.rank);np.savez_compressed(out/'GD_MODES.npz',**{k:v for k,v in modes.items() if isinstance(v,np.ndarray)})
  for case in a.cases:
   with np.load(ROOT/'inputs/private'/case/'truth.npz') as f:z=f['z']
   row,arrays=compare_state(g,z,modes);row['case']=case;rows.append(row);np.savez_compressed(out/f'{case}.npz',z=z,**arrays);dump(out/'SUMMARY.json',dict(status='ACTUAL_FULL_WAVE_DIAGNOSTIC',rows=rows,mode_cpu=modes['cpu_seconds'],mode_wall=modes['elapsed_seconds'],not_online_truth_optimization=True))
  status='COMPLETED'
 except BaseException as e:book.failure('SPECTRAL','pilot',e);raise
 finally:dump(out/'RECEIPT.json',book.receipt(status=status))
if __name__=='__main__':main()
