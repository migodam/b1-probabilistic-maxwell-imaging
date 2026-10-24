"""First fresh-anchor derivative/adjoint check, charged to Stage I."""
from pathlib import Path
import sys,json
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3r4.core import OUT,Book,Target,save,write
import numpy as np
def main():
 b=Book('I','fresh_anchor_derivative_adjoint',max_wall=300);status='COMPLETE';error=None
 try:
  t=Target(b,keep_fields=True);parent='R4DEV_00000';y=t.observation(parent)
  base=t.load_cache(OUT/f'results/contexts/cache/{parent}.npz','R4 paid dev context');J=t.jacobian(base)
  rng=np.random.default_rng(20261025031);v=np.r_[rng.normal(size=64),0.,0.];v/=np.linalg.norm(v)
  rows=[]
  for h in [.01,.001,.0001]:
   plus=t.evaluate(base.q+h*v,f'{parent}_FD_plus_{h}');minus=t.evaluate(base.q-h*v,f'{parent}_FD_minus_{h}')
   finite=t.white(plus.state.data-minus.state.data)/(2*h);exact=J@v
   rows.append(dict(h=h,relative_error=float(np.linalg.norm(finite-exact)/np.linalg.norm(exact)),
       plus=plus.key,minus=minus.key));del plus,minus
  # An actual source-adjoint action: the cached J is NOT allowed to answer it.
  with b.physical('source_adjoint_audit_8'):
   adj=t.backend().solve_anchor(t.backend().S.conj().T@(base.state.data-y),base.state,adjoint=True,role='source_adjoint_audit')
  hchi=np.sum(np.conj(adj)*base.state.field,axis=1)
  gradient=2*np.real(t.prior.derivative(base.q).T@hchi)/t.sigma**2
  gj=J.T@t.white(base.state.data-y)
  gerror=float(np.linalg.norm(gradient-gj)/max(np.linalg.norm(gj),1e-30))
  save(b.jobdir/'DERIVATIVES.npz',v=v,Jv=J@v,source_gradient=gradient,receiver_gradient=gj)
  if rows[-1]['relative_error']>1e-5 or gerror>1e-9:raise RuntimeError('paid derivative contract failed')
  write(OUT/'AUDIT/FRESH_PHYSICS_DERIVATIVE.json',dict(status='PASS',parent=parent,FD=rows,adjoint_gradient_relative_error=gerror,
    GT_used=False,all_new_forward_and_adjoint_paid=True,new_RHS_expected=56,new_states_expected=6,original_backend_unchanged=True))
 except BaseException as e:status='FAILED';error=repr(e);raise
 finally:b.close(status,error)
if __name__=='__main__':main()
