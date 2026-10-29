"""Local exact full physics with efficient trial LU/J; no truth access."""
import sys
from pathlib import Path
import numpy as np
from scipy.linalg import lu_factor,lu_solve
P=Path(__file__).resolve().parents[1];sys.path.insert(0,str(P/'src'))
from physics import Model
class FastModel(Model):
 def forward(self,m,label='trial',jac=False,independent=False):
  if independent:return super().forward(m,label=label,jac=jac,independent=True)
  m=np.asarray(m,dtype=float)
  if m.shape!=(self.n,) or not np.all(np.isfinite(m)):raise ValueError('finite real material shape required')
  self.state=self.factor=None
  with self.operation('FULL_FORWARD',rhs=self.ns,LU=1,new_states=1,Green_matvec_columns=self.ns,label=label):
   A=np.eye(self.n,dtype=complex)-m[:,None]*self.G;b=m[:,None]*self.E0;fac=lu_factor(A,check_finite=False);j=lu_solve(fac,b,check_finite=False);E=self.E0+self.G@j;y=self.S@j;s={'m':m.copy(),'j':j,'E':E,'y':y,'L':A,'state_relative_residual':float(np.linalg.norm(A@j-b)/max(np.linalg.norm(b),1e-300)),'resolvent_sigma_min':None}
  self.state,self.factor=s,fac
  if jac:self.jacobian(s)
  return s
 def jacobian(self,state=None):
  state=self.state if state is None else state
  if 'J' in state:return state['J']
  if state is not self.state or self.factor is None:raise RuntimeError('J requires retained matching LU')
  if self.n<=self.nr:
   with self.operation('FULL_MATERIAL_J_TANGENT',rhs=self.n,tangent_RHS=self.n,material_J_builds=1):
    inverse_columns=lu_solve(self.factor,np.eye(self.n),check_finite=False);transfer=self.S@inverse_columns;J=np.einsum('ri,is->sri',transfer,state['E']).reshape(self.ns*self.nr,self.n);state['Jcomplex']=J;state['J']=self.pack_jac(J)
   return state['J']
  return super().jacobian(state)
