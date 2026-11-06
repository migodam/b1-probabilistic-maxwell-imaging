"""Source-dependent receiver extension of established CIE/Global/A readouts.

Each source has its original SVD DPCS and native projected IDFT APCS. PR+, exact
quadratic line step, scalar/onehop material readouts and physical candidate
selection follow R6. No absent receiver is turned into a zero observation.
Only diagonal source metrics are supported; source mixing is not claimed.
"""
from collections import Counter
import time
import json
from pathlib import Path
import numpy as np
from scipy.linalg import svd
from research.research_max_struct_r2.src.cie_global import _baseline_readout
from research.research_max_struct_r4.src.cie_onehop import OneHopAbsorption,onehop_material_readout
from research.research_max_struct_r7.src.aplus_readout import mobius_material_readout


def run(model,Y,m0,arm,config):
 start=time.process_time();wall=time.perf_counter();cost=Counter();G,E,S=model.G,model.E,np.moveaxis(model.Ss,2,0).copy()
 n,ns=E.shape;m=m0.copy();reg=config['regularizer'];lo,hi=(np.full(n,x) for x in config['material_box'])
 wd=1/np.sum(abs(Y)**2,axis=0);ws=1/np.sum(abs(E)**2,axis=0);Os=np.diag(ws);Ls=np.diag(np.sqrt(ws));rho=1.
 sourcewise=arm in ('cie_single','cie_dual');blocks=OneHopAbsorption(G,model.shape,cost) if arm in ('a_material','a_plus') else None
 JD=np.zeros_like(E);V=[];ranks=[]
 for s in range(ns):
  u,q,vh=svd(S[s],full_matrices=False,check_finite=False);rank=int(np.count_nonzero(q>config['dpcs_rtol']*q[0]));v=vh[:rank].conj().T
  JD[:,s]=v@((u[:,:rank].conj().T@Y[:,s])/q[:rank]);V.append(v);ranks.append(rank);cost['receiver_SVD']+=1
 receive=lambda J:np.matmul(S,J.T[:,:,None])[:,:,0].T
 apply=lambda U,A:np.matmul(U,A.T[:,:,None])[:,:,0].T
 adjoint=lambda U,R:np.matmul(U.conj().transpose(0,2,1),R.T[:,:,None])[:,:,0].T
 GJD=G@JD;SJD=receive(JD);coeff={};history=[];candidates=[];stage_states=[]
 def validate(z,stage):
  value=model.physical_score(z,Y);candidates.append(dict(m=z.copy(),objective=value,stage=stage));cost['physical_validations']+=1
 validate(m,-1)
 for stageidx,st in enumerate(config['stages']):
  mf=st['mf'];bm=st['beta_material'];bj=st['beta_current'];shape=model.shape
  axes=[np.unique(np.r_[np.arange(mf),np.arange(z-mf,z)]) for z in shape];freq=np.array(np.meshgrid(*axes,indexing='ij')).reshape(2,-1).T
  grid=np.indices(shape).reshape(2,-1).T;F=np.exp(2j*np.pi*((grid/np.array(shape))@freq.T))/np.sqrt(n)
  U=np.stack([F-v@(v.conj().T@F) for v in V],axis=0)
  # This is the native projected Fourier operator. We omit only the inherited
  # reconstruction via an SVD which removes roundoff nulls then returns same
  # nonorthogonal coefficient coordinates; do NOT orthonormalize U.
  A=np.array([coeff.get(tuple(k),np.zeros(ns,complex)) for k in freq])
  GU=np.stack([G@U[s] for s in range(ns)],axis=0)
  SU=np.stack([S[s]@U[s] for s in range(ns)],axis=0)
  cost['G_columns']+=U.shape[2]*ns;cost['S_columns']+=U.shape[2]*ns
  Vcurrent=GU+bj*U
  oldg=oldd=None
  for it in range(st['maxiter']):
   r=bj*m/(1+bj*m);T=bj*JD-r[:,None]*(E+GJD+bj*JD)
   UA=apply(U,A);VA=apply(Vcurrent,A)
   D=SJD-Y+apply(SU,A);R=T+bj*UA-r[:,None]*VA
   grad=adjoint(SU,D)*wd+(bj*adjoint(U,R)-adjoint(Vcurrent,r[:,None]*R))*ws
   axis=0 if sourcewise else None
   if oldg is None:direction=-grad;pr=0.
   else:
    pr=np.maximum(0.,np.sum(np.real(grad.conj()*(grad-oldg)),axis=axis)/np.maximum(np.sum(abs(oldg)**2,axis=axis),1e-300));direction=-grad+oldd*pr
   gd=np.sum(np.real(grad.conj()*direction),axis=axis)
   if sourcewise:
    bad=gd>=0;direction[:,bad]=-grad[:,bad]
   elif gd>=0:direction=-grad
   DD=apply(SU,direction);RR=bj*apply(U,direction)-r[:,None]*apply(Vcurrent,direction)
   num=-np.sum(np.real(grad.conj()*direction),axis=axis);den=wd*np.sum(abs(DD)**2,axis=0)+ws*np.sum(abs(RR)**2,axis=0)
   if not sourcewise:den=den.sum()
   alpha=np.divide(num,den,out=np.zeros_like(num),where=den>1e-300);A+=direction*alpha
   J=JD+apply(U,A);Q=E+GJD+apply(GU,A)
   if arm in ('a_material','a_plus'):
    fn=mobius_material_readout if arm=='a_plus' else onehop_material_readout
    m,inner=fn(blocks,J,Q,m,bm,Os,Ls,rho,(lo,hi),reg,config,cost)
   else:m,inner=_baseline_readout(J@Ls,Q@Ls,m,bm,np.ones(ns),rho,(lo,hi),reg,True,config,cost)
   oldg,oldd=grad,direction;cost['PR_CG_steps']+=1
   if it%25==0 or it==st['maxiter']-1:
    history.append(dict(stage=stageidx,iteration=it+1,aux_data=.5*float((wd*np.sum(abs(SJD-Y+apply(SU,A))**2,axis=0)).sum()),material_min=float(m.min()),material_max=float(m.max()),cpu_seconds=time.process_time()-start))
    if config.get('progress_path'):Path(config['progress_path']).write_text(json.dumps(history[-1]))
  validate(m,stageidx);stage_states.append(dict(m=m.copy(),stage=stageidx,objective=candidates[-1]['objective']))
  coeff={tuple(k):a.copy() for k,a in zip(freq,A)}
 best=min(candidates,key=lambda c:c['objective'])
 return dict(m=best['m'],last_joint_m=m,candidates=candidates,history=history,stage_states=stage_states,
  execution_receipt=dict(status='COMPLETED',cpu_seconds=time.process_time()-start,wall_seconds=time.perf_counter()-wall,
   operation_counts=dict(cost),receiver_ranks=ranks,selected_stage=best['stage'],
   method='source-specific mature PR CIE extension; established material readouts',
   qualification='No original independent NIE claim; no missing-data zero filling; per-source native DPCS/APCS'))
