"""Paid, local target-material bias screen; not a nonlinear safety certificate."""
import numpy as np
from scipy.linalg import solve
from scipy.sparse import diags
from .fem_model import physical_lock

def _tv_hessian(a,smoothing):
 a=np.asarray(a,float);ny,nx=a.shape;n=a.size;H=np.zeros((n,n))
 for y in range(ny):
  for x in range(nx):
   i=y*nx+x;D=np.zeros((2,n))
   if x+1<nx:D[0,i]=-1;D[0,i+1]=1
   if y+1<ny:D[1,i]=-1;D[1,i+nx]=1
   v=D@a.ravel();m=np.sqrt(v@v+smoothing*smoothing)
   C=np.eye(2)/m-np.outer(v,v)/m**3;H+=D.T@C@D
 return H

def prior_hessian(model,z,o):
 nw,nt=model.nwall,model.ntarget;H=np.zeros((model.nparam,model.nparam))
 ny,nx=model.wall_shape;L=np.zeros((nw,nw))
 for y in range(ny):
  for x in range(nx):
   i=y*nx+x
   for j in ([i+1] if x+1<nx else [])+([i+nx] if y+1<ny else []):
    L[i,i]+=1;L[j,j]+=1;L[i,j]-=1;L[j,i]-=1
 H[:nw,:nw]=(np.eye(nw)+o.wall_spatial_weight*L)/o.wall_eps_sigma**2
 H[nw:2*nw,nw:2*nw]=(np.eye(nw)+o.wall_spatial_weight*L)/o.wall_loss_sigma**2
 for offset,scale in [(2*nw,o.target_chi_scale),(2*nw+nt,o.target_loss_scale)]:
  a=z[offset:offset+nt].reshape(model.target_shape)/scale
  H[offset:offset+nt,offset:offset+nt]=(o.tv_weight*_tv_hessian(a,o.tv_smoothing)+o.target_l2_weight*np.eye(nt))/scale**2
 return H

def full_material_jacobian(reference,z):
 """Receiver products reuse the already paid full material checkpoint LU."""
 model=reference.model;nw,nt=model.nwall,model.ntarget;pieces=[];data=[]
 for g in reference.greens:
  # The cache getter acquires its own lock only on a physical-state miss.
  state=g.physical_state_with_factor(z);p=g.p;P=g.P
  with physical_lock(model.book):
   adj=state['factor'].solve(g.ev.R.conj().T.toarray(),trans='H',label='target risk receiver adjoint16 cached full LU')
  field=state['field'];av=P@adj
  product=p.k**2*np.einsum('qr,qs,q->qsr',np.conj(av),field,g.weights,optimize=True)
  Jc=np.zeros((nw+nt,4,16),complex);np.add.at(Jc,g.cells,product)
  JJ=np.column_stack([Jc[:nw].reshape(nw,-1).T,(1j/p.freq*Jc[:nw]).reshape(nw,-1).T,Jc[nw:].reshape(nt,-1).T,(1j/p.freq*Jc[nw:]).reshape(nt,-1).T])
  pieces.append(JJ);data.append(state['data'])
 J=np.vstack(pieces);return J,np.asarray(data)

def target_bias_screen(reference,z,prior,candidates,sd_real,sd_loss,max_scaled_rms=.1):
 model=reference.model;nw,nt=model.nwall,model.ntarget;candidates=np.asarray(candidates,int)
 if not len(candidates):return candidates,dict(kind='TARGET_RISK_SCREEN',decision='KEEP_ACTIVE',reason='no preliminary candidates',additional_full_J_computed=False)
 J,data=full_material_jacobian(reference,z)
 J=np.vstack((J.real,J.imag))/model.sigma;Aw=J[:,:2*nw];B=J[:,2*nw:];Hprior=prior_hessian(model,z,prior)
 remaining=np.asarray([i for i in range(nw) if i not in candidates]);Ridx=np.r_[remaining,nw+remaining].astype(int)
 if len(Ridx):
  AR=Aw[:,Ridx];HR=AR.T@AR+Hprior[np.ix_(Ridx,Ridx)];BR=B.T@AR
  S=B.T@B+Hprior[2*nw:,2*nw:]-BR@solve(HR,BR.T,assume_a='pos')
 else:AR=None;S=B.T@B+Hprior[2*nw:,2*nw:]
 scale=np.r_[np.full(nt,prior.target_chi_scale),np.full(nt,prior.target_loss_scale)]
 rows=[];chosen=[]
 for cell in candidates:
  AA=Aw[:,[cell,nw+cell]];cross=B.T@AA
  if AR is not None:cross-=BR@solve(HR,AR.T@AA,assume_a='pos')
  bias=-solve(S,cross,assume_a='pos');sig=np.array([sd_real[cell],sd_loss[cell]])
  # Conditional diagonal surrogate only. No posterior/box/nonlinear guarantee.
  scaled=(bias*sig[None,:])/scale[:,None];rms=float(np.linalg.norm(scaled)/np.sqrt(len(scale)))
  ok=rms<=max_scaled_rms
  if ok:chosen.append(int(cell))
  rows.append(dict(cell=int(cell),local_target_material_scaled_RMS=rms,eligible=ok))
 ordered=sorted(rows,key=lambda r:(r['local_target_material_scaled_RMS'],r['cell']));chosen=[];joint_sq=0.
 for row in ordered:
  rr=row['local_target_material_scaled_RMS']
  if row['eligible'] and np.sqrt(joint_sq+rr*rr)<=max_scaled_rms:chosen.append(row['cell']);joint_sq+=rr*rr
 lo,hi,_=model.bounds_and_scale();active=int(np.sum((z<=lo+1e-10)|(z>=hi-1e-10)))
 event=dict(kind='TARGET_RISK_SCREEN',candidate_cells=candidates.tolist(),selected_cells=chosen,rows=rows,additional_full_J_computed=True,full_J_complex_shape=[J.shape[0]//2,J.shape[1]],local_target_schur_condition=float(np.linalg.cond(S)),active_box_coordinates=active,scope='LOCAL_UNCONSTRAINED_GAUSSIAN_TANGENT_WITH_CONDITIONAL_SD_PROXY; no nonlinear certificate or posterior claim',decision='ABSORB' if chosen else 'KEEP_ACTIVE',max_scaled_RMS=max_scaled_rms,joint_selected_scaled_RMS=float(np.sqrt(joint_sq)))
 model.book.event(model.phase,'paid local target bias screen',event_kind=event['kind'],**{k:v for k,v in event.items() if k!='kind'})
 return np.asarray(chosen,int),event
