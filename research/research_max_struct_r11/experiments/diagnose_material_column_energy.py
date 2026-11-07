"""Ten saved-state material Jacobian column energies; zero physical solves, no GT."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
from collections import Counter
import sys,json,csv,time,traceback
import numpy as np
from scipy.stats import spearmanr
ROOT=Path(__file__).absolute().parents[1];sys.path.insert(0,str(ROOT.parents[1]))
from research.research_max_struct_r10.src.lowrank_denominator import LowRankDenominator
OUT=ROOT/'delegated/diagnostics/material_column_energy'
def dump(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def quant(x):
 x=np.asarray(x);v=x[np.isfinite(x)]
 return dict(zip(('min','q10','median','q90','max'),map(float,np.quantile(v,[0,.1,.5,.9,1])))) if len(v) else None
def ratio(a,b):return np.divide(a,b,out=np.full_like(a,np.nan),where=b!=0)
def corr(a,b):
 good=np.isfinite(a)&np.isfinite(b)
 return float(spearmanr(a[good],b[good]).statistic) if good.sum()>1 and np.ptp(a[good])>0 and np.ptp(b[good])>0 else None
def main():
 wall=time.perf_counter();protocol=json.loads((OUT/'FROZEN_DIAGNOSTIC.json').read_text());rows=[];cost=Counter()
 receipt=dict(status='RUNNING',GT_access='NONE',physical_LU=0,physical_FW=0,physical_adjoint=0,new_inverse=0,classification='FIXED_STATE_EXPOSED_OFFLINE_DIAGNOSTIC_NO_GT')
 try:
  for st in protocol['states']:
   count=Counter();beg=time.process_time();sid=st['case_id']+'_s'+str(st['stage'])
   with np.load(st['state'],allow_pickle=False) as z:m,J,Q=[z[k].copy() for k in ('m','J','Q')]
   with np.load(st['basis'],allow_pickle=False) as z:U,B=[z[k].copy() for k in ('U','B')]
   with np.load(st['source_metric'],allow_pickle=False) as z:weights=z['state_weights'].copy()
   with np.load(st['public'][0],allow_pickle=False) as z:E=z['E'].copy()
   # Hydrate the exact saved production operator; no fresh B=PGPU multiplication needed.
   op=LowRankDenominator.__new__(LowRankDenominator);op.U=U;op.B=B;op.projection='galerkin';op.n,op.rank=U.shape;op.cost=count
   beta=float(st['beta']);f=op.factor(m,beta);w=op.solve(f,J-m[:,None]*Q);qeff=Q+beta*w-op.action(w)
   columns=op.solve(f,np.eye(len(m),dtype=complex));inv_energy=np.sum(np.abs(columns)**2,axis=0)
   qeff0=(Q+beta*J)/(1+beta*m[:,None]);inv0_energy=1/(1+beta*m)**2
   def energy(x):
    count['source_metric_diagonal_products']+=1
    return np.sum(np.abs(x)**2*weights[None,:],axis=1)
   qe=energy(qeff);qe0=energy(qeff0);ee=energy(E);de=energy(E+beta*J)
   gn=inv_energy*qe;gn0=inv0_energy*qe0
   qr,ar,gr=ratio(qe,qe0),ratio(inv_energy,inv0_energy),ratio(gn,gn0)
   distance=np.minimum(m-.02,20-m);normalized_distance=distance/(20-.02)
   lower=m<=.02+1e-12;upper=m>=20-1e-12;interior=~(lower|upper)
   data=dict(m=m,beta=np.array(beta),state_weights=weights,Qeff=qeff,Qeff_ordinary=qeff0,Qeff_energy=qe,Qeff_ordinary_energy=qe0,inverse_column_norm2=inv_energy,inverse_column_norm2_ordinary=inv0_energy,GN_diagonal_without_common_beta2rho=gn,GN_diagonal_ordinary_without_common_beta2rho=gn0,Qeff_energy_ratio=qr,inverse_column_norm2_ratio=ar,total_column_energy_ratio=gr,distance_to_box=distance,normalized_distance_to_box=normalized_distance,lower_active=lower,upper_active=upper,interior=interior,excitation_degeneracy_distance_per_row=ratio(np.sqrt(de),np.sqrt(ee)))
   np.savez_compressed(OUT/(sid+'.npz'),**data)
   groups={name:dict(count=int(mask.sum()),Qeff_energy_ratio=quant(qr[mask]),inverse_column_norm2_ratio=quant(ar[mask]),total_column_energy_ratio=quant(gr[mask])) for name,mask in [('exact_lower_bound',lower),('exact_upper_bound',upper),('interior',interior)]}
   row=dict(state_id=sid,case_id=st['case_id'],stage=st['stage'],beta=beta,rank=op.rank,material_quantiles=quant(m),distance_to_box_quantiles=quant(distance),Qeff_energy_ratio=quant(qr),inverse_column_norm2_ratio=quant(ar),total_column_energy_ratio=quant(gr),Qeff_energy_Green=quant(qe),Qeff_energy_ordinary=quant(qe0),GN_diagonal_Green=quant(gn),GN_diagonal_ordinary=quant(gn0),zero_ordinary_Qeff_columns=int(np.count_nonzero(qe0==0)),zero_Green_Qeff_columns=int(np.count_nonzero(qe==0)),spearman_distance_Qeff_ratio=corr(distance,qr),spearman_distance_inverse_ratio=corr(distance,ar),spearman_distance_total_ratio=corr(distance,gr),bound_groups=groups,excitation_degeneracy_distance_weighted=float(np.sqrt(de.sum()/ee.sum())),excitation_degeneracy_distance_row_quantiles=quant(ratio(np.sqrt(de),np.sqrt(ee))),auxiliary_cost=dict(count),state_CPU_seconds=time.process_time()-beg,raw_array=str(OUT/(sid+'.npz')))
   rows.append(row);cost.update(count)
  dump(OUT/'ROWS.json',rows)
  flat=[]
  for r in rows:
   for component in ('Qeff_energy_ratio','inverse_column_norm2_ratio','total_column_energy_ratio'):
    flat.append(dict(state_id=r['state_id'],component=component,**r[component],spearman_distance_to_box=r['spearman_distance_'+{'Qeff_energy_ratio':'Qeff_ratio','inverse_column_norm2_ratio':'inverse_ratio','total_column_energy_ratio':'total_ratio'}[component]],lower_active=r['bound_groups']['exact_lower_bound']['count'],upper_active=r['bound_groups']['exact_upper_bound']['count'],excitation_degeneracy_distance_weighted=r['excitation_degeneracy_distance_weighted']))
  with (OUT/'QUANTILES.csv').open('w') as f:
   cw=csv.DictWriter(f,fieldnames=list(flat[0]));cw.writeheader();cw.writerows(flat)
  receipt['status']='COMPLETE'
 except Exception:
  receipt.update(status='FAILED',traceback=traceback.format_exc());raise
 finally:
  receipt.update(states_completed=len(rows),auxiliary_cost=dict(cost),CPU_seconds_near_write=time.process_time(),wall_seconds=time.perf_counter()-wall,physical_counts=dict(LU=0,FULL_FORWARD=0,DATA_ADJOINT=0,actual_complex_RHS=0),prior_scope='No prior Hessian included: auxiliary residual Jacobian GN only; common rho*beta^2 omitted equally in ratios',failure_attempts=0 if receipt['status']=='COMPLETE' else 1)
  dump(OUT/'RECEIPT.json',receipt)
 print(json.dumps(dict(status=receipt['status'],states=len(rows),auxiliary_cost=dict(cost))))
if __name__=='__main__':main()
