"""Post-seal diagnostic only: does the wrong estimate outscore the legal truth?"""
from pathlib import Path
import os,sys,json,time,numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from gsda.fwi import FWIOptions,regularization
from gsda.accounting import BudgetBook,dump

def main():
 book=BudgetBook(job='offline_truth_objective');out=ROOT/'evaluation';out.mkdir(exist_ok=True);rows=[];options=FWIOptions(**json.loads((ROOT/'configs/main.json').read_text())['fwi_options'])
 for ep in sorted((ROOT/'results').glob('*/c*/endpoint.npz')):
  d=ep.parent;seal=d/'ONLINE_SEALED.json'
  if not seal.exists():continue
  with np.load(ep) as f:z=f['z'];pred=f['data']
  with np.load(ROOT/'inputs/private'/d.name/'truth.npz') as f:zt=f['z'];clean=f['clean']
  with np.load(ROOT/'inputs/online'/d.name/'observation.npz') as f:y=f['y']
  if pred.shape!=y.shape:continue
  selected_options=options
  if d.parent.name=='B0_TV0_ABLATION':selected_options=FWIOptions(**json.loads((ROOT/'configs/fwi_tv0_c006.json').read_text())['fwi_options'])
  prior,_,terms=regularization(zt,6,(7,7),(2,3),selected_options);priorhat,_,termshat=regularization(z,6,(7,7),(2,3),selected_options)
  td=.5*float(np.sum(abs(clean-y)**2))/.01**2;pd=.5*float(np.sum(abs(pred-y)**2))/.01**2
  rows.append(dict(method_dir=d.parent.name,case=d.name,truth_in_exact_same_material_class=True,truth_data_objective=td,truth_prior=prior,truth_prior_components=terms,truth_full_objective=td+prior,estimated_data_objective=pd,estimated_prior=priorhat,estimated_prior_components=termshat,estimated_full_objective=pd+priorhat,estimated_minus_truth_full_objective=pd+priorhat-td-prior,estimated_outscores_truth=bool(pd+priorhat<td+prior),clean_mean_white_distance=float(np.linalg.norm(pred-clean)/.01),not_noise_normalized=True,not_exact_fiber_or_global_optimum=True))
 dump(out/'TRUTH_OBJECTIVE_OFFLINE.json',dict(rows=rows,new_physics=0,online_used_truth=False,limits='lower criterion does not prove exact nonidentifiability or global optimum; no refitting or branch reselection'))
 book.receipt(status='COMPLETED',new_physics=0,scoring_only=True)
 print('postseal diagnostics',len(rows))
if __name__=='__main__':main()
