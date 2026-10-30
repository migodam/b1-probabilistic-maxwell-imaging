"""Further public FEM level for pixel-world responses; independent of inversion."""
from pathlib import Path
import os,sys,argparse,json,time
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='1'
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'vendor'))
from gsda.accounting import BudgetBook,dump
from gsda.fem_model import FEMModel

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--cases',nargs='+',default=['c000','c006','c010','c011','c012']);a=ap.parse_args()
 out=ROOT/'results/independent_validation';out.mkdir(parents=True,exist_ok=True);book=BudgetBook(job='independent_fine');status='FAILED';rows=[]
 try:
  c=json.loads((ROOT/'configs/base.json').read_text());c.update(h=.12,targetstep=.05,wall_refine=2)
  model=FEMModel(c,book,'R4_FINE_VALIDATION')
  for case in a.cases:
   with np.load(ROOT/'inputs/private'/case/'truth.npz') as f:z=f['z'].copy();coarse=f['clean'].copy()
   pred=model.predict(z);err=float(np.linalg.norm(pred.data-coarse)/.01)
   arrays=dict(z=z,coarse_data=coarse,fine_data=pred.data)
   row=dict(case=case,truth_white_mesh_distance=err,mesh=dict(h=.12,targetstep=.05,wall_refine=2),scope='empirical further public FEM refinement; no continuum error bound',original_W_tolerance=.1,within_original_tolerance=bool(err<=.1),diagnostics=pred.diagnostics)
   ep=ROOT/'results/B0_FULL_FWI'/case/'endpoint.npz'
   if ep.exists() and (ep.parent/'ONLINE_SEALED.json').exists():
    with np.load(ep) as f:zh=f['z'];dc=f['data']
    pp=model.predict(zh);row['endpoint_white_mesh_distance']=float(np.linalg.norm(pp.data-dc)/.01);arrays.update(estimated_z=zh,estimated_coarse_data=dc,estimated_fine_data=pp.data)
   np.savez_compressed(out/f'{case}.npz',**arrays);rows.append(row);dump(out/'SUMMARY.json',dict(status='ACTUAL_FURTHER_FEM_LEVEL',rows=rows,config=c))
  status='COMPLETED'
 except BaseException as e:book.failure('VALIDATION','fine',e);raise
 finally:dump(out/'RECEIPT.json',book.receipt(status=status))
if __name__=='__main__':main()
