"""Zero-FW split of scalar contraction and nonlocal feedback metric derivatives."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import sys,json,time
from collections import Counter
import numpy as np
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).absolute().parents[1];sys.path.insert(0,str(ROOT.parents[1]))
from research.research_max_struct_r11.src.diagnostic_core import cosine,projected_gradient
OUT=ROOT/'delegated/diagnostics'

def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False,default=lambda v:v.item() if isinstance(v,np.generic) else v.tolist())+'\n')

def stats(g,gphi,m,lo,hi):
    pg=projected_gradient(m,g,lo,hi)
    return dict(norm=float(np.linalg.norm(g)),inf=float(np.linalg.norm(g,np.inf)),
      cosine_with_Phi_gradient=cosine(g,gphi),unit_descent_linear_work_Phi=float(-np.dot(g,gphi)),
      projected_gradient_norm=float(np.linalg.norm(pg)),projected_unit_descent_linear_work_Phi=float(-np.dot(pg,gphi)),
      projected_descent_cos_negative_Phi_gradient=cosine(pg,gphi))

def main():
    wall,cpu=time.perf_counter(),time.process_time();inv=json.loads((OUT/'STATE_INVENTORY.json').read_text());rows=[]
    with threadpool_limits(limits=1):
      for obj in inv['objects']:
        with np.load(Path(obj['folder'])/'DENOMINATOR_BASIS.npz') as z:U=z['U'];B=z['B']
        for state in obj['stages']:
          sid=f"{obj['case_id']}_s{state['stage']}"
          with np.load(OUT/f'{sid}_START.npz') as z:m=z['m'];J=z['J'];Q=z['Q'];w=z['w'];a=z['a'];gphi=z['physical_gradient'];gc=z['g_const'];gp=z['g_prior'];gm=z['g_metric']
          beta=state['metadata']['beta_material'];kw=B@(U.conj().T@w)
          gb=-beta**3*np.real(np.sum(a.conj()*w,axis=1))
          gk=beta**2*np.real(np.sum(a.conj()*kw,axis=1))
          # T=A D_d^-1, so freezing T gives dw/dm=-A^-1 (Q+beta J)/(1+beta m).
          # a_T=T^-*(w Os)=(1+beta*m)*a, giving this mapped-m derivative.
          contracted=-beta**2*np.real(np.sum(a.conj()*(Q+beta*J)/(1+beta*m[:,None]),axis=1))+gp
          lo,hi=np.full(m.size,.02),np.full(m.size,20.)
          row=dict(state_id=sid,beta=beta,truth_access='NO_NEW_TRUTH_ACCESS; no candidate selection',
            reconstruction_relative=float(np.linalg.norm(gb+gk-gm)/max(np.linalg.norm(gm),1e-300)),
            beta_to_K_norm_ratio=float(np.linalg.norm(gb)/max(np.linalg.norm(gk),1e-300)),
            cosine_beta_K=cosine(gb,gk),cosine_const_beta=cosine(gc,gb),cosine_const_K=cosine(gc,gk),
            gradients={k:stats(g,gphi,m,lo,hi) for k,g in dict(g_beta=gb,g_K=gk,g_const=gc,g_metric=gm,g_contracted_lagged=contracted).items()})
          np.savez_compressed(OUT/f'{sid}_METRIC_SPLIT.npz',g_beta=gb,g_K=gk,g_contracted_lagged=contracted)
          rows.append(row)
    dump(OUT/'METRIC_SPLIT.json',rows)
    dump(OUT/'METRIC_SPLIT_RECEIPT.json',dict(status='COMPLETE',states=len(rows),physical_forward=0,physical_adjoint=0,truth_access='NONE_NEW',
      operation_cost=dict(saved_lowrank_Green_actions=10,saved_lowrank_Green_complex_RHS=120),wall_seconds=time.perf_counter()-wall,cpu_seconds=time.process_time()-cpu,
      mapped_m_derivative='g_contracted_lagged=-rho beta^2 Re sum conj(a) (Q+beta J)/(1+beta m)+gprior; a=A^-*(w Omega); only T frozen'))
    for r in rows:
      print(r['state_id'],'beta/K %.4g'%r['beta_to_K_norm_ratio'],'cos_const_beta %.4g'%r['cosine_const_beta'],'cos_const_K %.4g'%r['cosine_const_K'],'norm beta/K %.5g/%.5g'%(r['gradients']['g_beta']['norm'],r['gradients']['g_K']['norm']),'contracted projected Phi cos %.4g'%r['gradients']['g_contracted_lagged']['projected_descent_cos_negative_Phi_gradient'])
if __name__=='__main__':main()
