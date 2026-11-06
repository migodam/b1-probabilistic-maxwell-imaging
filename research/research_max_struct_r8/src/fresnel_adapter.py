"""Measured Fresnel adapters; original measurements stay local.

Source-specific S_s exactly selects measured receiver channels. Equal-area
Richmond disk quadrature is an area-equivalent collocation rule, not disjoint
physical disks. Incident source multipoles are fitted ONLY to empty-target data.
"""
from pathlib import Path
from collections import Counter
import time
import numpy as np
from scipy.special import hankel1,jv
from scipy.linalg import lu_factor,lu_solve
from research.research_fbg_r3_experiments.src.physics import MaxwellModel
from research.research_max_struct_r8.src.paper_physics import MatrixFreeModel

def load_file(path,frequency=2.):
 a=np.loadtxt(path,comments='#')
 if a.ndim!=2 or a.shape[1]!=7 or not np.isfinite(a).all():raise ValueError('Original measurement file requires seven finite columns')
 a=a[np.isclose(a[:,2],frequency)]
 if not len(a):raise ValueError('Requested physical GHz frequency absent')
 ns=18 if 'Foam' in str(path) else 36;nr=360 if ns==18 else 72
 Y=np.full((nr,ns),np.nan+1j*np.nan);I=Y.copy()
 for row in a:
  j,k=int(row[1])-1,int(row[0])-1
  if np.isfinite(Y[j,k]):raise ValueError('Duplicate actual measurement')
  Y[j,k]=np.conj((row[3]+1j*row[4])-(row[5]+1j*row[6]));I[j,k]=np.conj(row[5]+1j*row[6])
 return a,Y,I


def operators(centers,sources,receivers,h,k,dense=True):
 a=h/np.sqrt(np.pi)
 def external(d):return .5j*np.pi*k*a*jv(1,k*a)*hankel1(0,k*d)
 out={}
 if dense:
  d=np.linalg.norm(centers[:,None]-centers[None,:],axis=2);np.fill_diagonal(d,1.)
  G=external(d);np.fill_diagonal(G,.5j*np.pi*k*a*hankel1(1,k*a)-1);out['G']=G
 return dict(**out,S=external(np.linalg.norm(receivers[:,None]-centers[None,:],axis=2)),radius=a,k=k)


def incident_fit(I,sources,receivers,centers,k,order=2):
 """Outgoing source multipole fit on near-axis empty-target receivers.

 Fit m=-2..2 outgoing cylindrical harmonics centered at each transmitting horn.
 Symmetric measured transmit/receive pattern is NOT separated. Hence restrict
 to +/-30deg about opposite receiver (15deg at transmitting horn), reserve
 alternating angles, and report empty-field validation. This is a source-free
 Helmholtz continuation in the object region; receive-horn contamination is
 retained as an explicit model uncertainty, not a calibrated antenna solution.
 """
 ns=len(sources);E=np.zeros((len(centers),ns),complex);rows=[];solutions=[]
 for s in range(ns):
  axis=-sources[s]/np.linalg.norm(sources[s]);normal=np.array([-axis[1],axis[0]])
  def basis(points):
   v=points-sources[s];dist=np.linalg.norm(v,axis=1);ang=np.arctan2(v@normal,v@axis)
   return np.array([.25j*hankel1(abs(m),k*dist)*np.exp(1j*m*ang) for m in range(-order,order+1)]).T
  angle=np.arccos(np.clip(receivers@axis/np.linalg.norm(receivers,axis=1),-1,1))
  ids=np.flatnonzero(np.isfinite(I[:,s])&(angle<=np.pi/6+1e-8));fit=ids[np.arange(len(ids))%3!=0];val=ids[np.arange(len(ids))%3==0]
  B=basis(receivers);scale=np.linalg.norm(B[fit],axis=0);Z=B/scale
  # Tikhonov chooses unique minimum-coefficient continuation from calibration,
  # fixed 1e-4 relative singular cutoff; no material/heldout scattered tuning.
  u,q,vh=np.linalg.svd(Z[fit],full_matrices=False);cut=q>1e-4*q[0]
  coef=vh[cut].conj().T@((u[:,cut].conj().T@I[fit,s])/q[cut])
  E[:,s]=(basis(centers)/scale)@coef
  solutions.append(coef/scale)
  pred=Z@coef
  # Retain cylindrical comparator for independently held-out empty channels.
  opposite=np.argmax(receivers@axis);base=I[opposite,s]/(.25j*hankel1(0,k*np.linalg.norm(receivers[opposite]-sources[s])))
  cyl=.25j*hankel1(0,k*np.linalg.norm(receivers-sources[s],axis=1))*base
  rows.append(dict(source=s+1,fit_channels=fit.tolist(),validation_channels=val.tolist(),rank=int(cut.sum()),
   fit_relative=float(np.linalg.norm((pred-I[:,s])[fit])/np.linalg.norm(I[fit,s])),
   validation_relative=float(np.linalg.norm((pred-I[:,s])[val])/np.linalg.norm(I[val,s])),
   cylindrical_validation_relative=float(np.linalg.norm((cyl-I[:,s])[val])/np.linalg.norm(I[val,s])),
   extrapolation_norm_ratio=float(np.linalg.norm(E[:,s])/max(np.linalg.norm(.25j*hankel1(0,k*np.linalg.norm(centers-sources[s],axis=1))*base),1e-300))))
 def evaluate(points,normals=None):
  result=[]
  from scipy.special import h1vp
  for s,coef in enumerate(solutions):
   axis=-sources[s]/np.linalg.norm(sources[s]);normal=np.array([-axis[1],axis[0]])
   v=np.asarray(points)-sources[s];dist=np.linalg.norm(v,axis=1);unit=v/dist[:,None];ang=np.arctan2(v@normal,v@axis)
   modes=np.arange(-order,order+1);H=hankel1(abs(modes)[None,:],k*dist[:,None]);phase=np.exp(1j*ang[:,None]*modes)
   B=.25j*H*phase
   if normals is not None:
    dr=np.sum(unit*normals,axis=1);tangent=np.c_[-unit[:,1],unit[:,0]];dt=np.sum(tangent*normals,axis=1)/dist
    B=.25j*phase*(k*h1vp(abs(modes)[None,:],k*dist[:,None])*dr[:,None]+1j*modes*H*dt[:,None])
   result.append(B@coef)
  return np.stack(result,axis=1)
 return E,rows,evaluate


class SourceReceiverModel(MaxwellModel):
 """Full original scalar-TM LU plus correctly gathered source receiver maps."""
 def __init__(self,public,source_indices,receiver_indices):
  super().__init__(dict(public,E0=public['E0'][:,source_indices]));self.source_indices=np.asarray(source_indices)
  self.receiver_indices=np.asarray(receiver_indices);self.Ss=np.stack([public['S'][idx] for idx in receiver_indices],axis=2)
  self.nr=self.Ss.shape[0];self._current=None
 def receive(self,J):return np.einsum('rns,ns->rs',self.Ss,J,optimize=True)
 def receive_adjoint(self,R):return np.einsum('rns,rs->ns',self.Ss.conj(),R,optimize=True)
 def state(self,m):
  state=super().state(m)
  if not state.get('source_receivers',False):state['Y']=self.receive(state['j']);state['source_receivers']=True
  return state
 def data_value_gradient(self,m,Y,beta=None,reference=None):
  if beta is not None:raise ValueError('Feedback objective is outside this measured adapter')
  state=self.state(m);R=state['Y']-Y;scale=float(np.vdot(Y,Y).real)
  with self._physics.operation('DATA_ADJOINT',rhs=self.ns,adjoint_RHS=self.ns):
   adj=lu_solve(state['lu'],self.receive_adjoint(R),trans=2,check_finite=False)
  reg,g,_=self.reg.evaluate(m)
  return float(.5*np.vdot(R,R).real/scale+reg),np.real(np.sum(adj.conj()*state['Q'],axis=1))/scale+g




class BoundedMatrixFreeModel(MatrixFreeModel):
 """Same VIE with bounded iterative effort then exact retained dense library LU.

 No scientific approximation or relaxed residual. A difficult first RHS
 triggers one LU for all source RHS and its matched analytic adjoint.
 """
 def solve(self,m,B,adjoint=False):
  from scipy.sparse.linalg import LinearOperator,gmres
  if self._dense_m is not None and np.array_equal(m,self._dense_m):return self.dense_fallback(m,B,adjoint)
  def mv(v):return v-self.green.adjoint((m*v)[:,None])[:,0] if adjoint else v-m*self.green.apply(v[:,None])[:,0]
  A=LinearOperator((self.n,self.n),matvec=mv,dtype=complex);diagonal=1-m*(np.conj(self.green.diagonal) if adjoint else self.green.diagonal)
  P=LinearOperator((self.n,self.n),matvec=lambda v:v/diagonal,dtype=complex);Z=np.empty_like(B)
  for j in range(B.shape[1]):
   old=self._solve_initial.get(adjoint);x0=None if old is None or old.shape!=B.shape else old[:,j];it=[]
   Z[:,j],info=gmres(A,B[:,j],M=P,x0=x0,rtol=1e-9,atol=0,restart=40,maxiter=5,callback=lambda v:it.append(v),callback_type='pr_norm')
   self.cost['GMRES_iterations']+=len(it);self.cost['actual_complex_RHS']+=1
   res=float(np.linalg.norm(mv(Z[:,j])-B[:,j])/max(np.linalg.norm(B[:,j]),1e-300))
   if info or res>2e-9:
    self.cost['failed_bounded_GMRES_attempts']+=1;self.cost['maximum_failed_GMRES_relative_residual']=max(self.cost.get('maximum_failed_GMRES_relative_residual',0),res)
    return self.dense_fallback(m,B,adjoint)
   self.cost['maximum_GMRES_relative_residual']=max(self.cost.get('maximum_GMRES_relative_residual',0),res)
  self._solve_initial[adjoint]=Z.copy();return Z


class SourceReceiverMatrixFreeModel(BoundedMatrixFreeModel):
 def __init__(self,public,source_indices,receiver_indices):
  super().__init__(dict(public,E0=public['E0'][:,source_indices]));self.source_indices=np.asarray(source_indices);self.receiver_indices=np.asarray(receiver_indices)
  self.Ss=np.stack([public['S'][idx] for idx in receiver_indices],axis=2);self._SsC=np.moveaxis(self.Ss,2,0).copy();self.nr=self.Ss.shape[0]
 def receive(self,J):return np.matmul(self._SsC,J.T[:,:,None])[:,:,0].T
 def receive_adjoint(self,R):return np.matmul(self._SsC.conj().transpose(0,2,1),R.T[:,:,None])[:,:,0].T
 def state(self,m):
  state=super().state(m)
  if not state.get('source_receivers',False):state['Y']=self.receive(state['j']);state['source_receivers']=True
  return state
 def data_value_gradient(self,m,Y):
  st=self.state(m);r=st['Y']-Y;scale=np.vdot(Y,Y).real;a=self.solve(m,self.receive_adjoint(r),adjoint=True)
  prior,gprior,_=self.reg.evaluate(m);self.cost['DATA_ADJOINT']+=1
  return .5*float(np.vdot(r,r).real/scale)+prior,np.real(np.sum(a.conj()*st['Q'],axis=1))/scale+gprior


def prepare(path,frequency=2.,side=24,width=.16,illumination='multipole',matrix_free=False):
 started=time.process_time();a,Y,I=load_file(path,frequency);k=2*np.pi*frequency*1e9/299792458.
 h=width/side;axis=(np.arange(side)+.5)*h-width/2;yy,xx=np.meshgrid(axis,axis,indexing='ij');centers=np.c_[xx.ravel(),yy.ravel()]
 nr,ns=Y.shape;foam=ns==18
 sa=2*np.pi*np.arange(ns)/ns if foam else np.arange(36)*np.pi/18
 ra=2*np.pi*np.arange(nr)/nr if foam else np.arange(72)*np.pi/36
 rs=1.67 if foam else .72;rr=1.67 if foam else .76
 sources=rs*np.c_[np.cos(sa),np.sin(sa)];receivers=rr*np.c_[np.cos(ra),np.sin(ra)]
 public=operators(centers,sources,receivers,h,k,dense=not matrix_free)
 if illumination=='multipole':E,cal,incident_evaluator=incident_fit(I,sources,receivers,centers,k)
 elif illumination=='cylindrical':
  opposite=(np.arange(ns)*(nr//ns)+nr//2)%nr;gain=I[opposite,np.arange(ns)]/(.25j*hankel1(0,k*(rs+rr)))
  E=.25j*hankel1(0,k*np.linalg.norm(centers[:,None]-sources[None,:],axis=2))*gain;cal=[]
 else:raise ValueError(illumination)
 public.update(E0=E,grid_shape=np.array([side,side]),centers=centers)
 held_sources=np.arange(0,ns,4);train_sources=np.setdiff1d(np.arange(ns),held_sources)
 train_mask=np.zeros(Y.shape,bool);receiver_indices=[]
 for s in train_sources:
  ids=np.flatnonzero(np.isfinite(Y[:,s]));# sorting relative angle retains rotating acquisition mask
  ids=ids[np.argsort((ra[ids]-sa[s])%(2*np.pi))];idx=ids[np.arange(len(ids))%5!=0]
  train_mask[idx,s]=True;receiver_indices.append(idx)
 model=(SourceReceiverMatrixFreeModel if matrix_free else SourceReceiverModel)(public,train_sources,receiver_indices)
 return dict(model=model,public=public,Y=Y,I=I,Y_train=np.stack([Y[idx,s] for idx,s in zip(receiver_indices,train_sources)],axis=1),
  train_sources=train_sources,train_mask=train_mask,withheld_mask=np.isfinite(Y)&~train_mask,held_sources=held_sources,centers=centers,
  calibration=cal,incident_evaluator=incident_evaluator if illumination=='multipole' else None,setup_cpu_seconds=time.process_time()-started,illumination=illumination)
