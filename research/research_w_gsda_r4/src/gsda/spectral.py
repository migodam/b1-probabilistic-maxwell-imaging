"""Actual fixed-reference GD singular actions and projected weak-Schur audit.

The SVD is of the quadrature-L2 Green operator, not a material Jacobian.
Exact Schur elimination is exact only in the stated retained current space.
All omitted-space errors are separately compared with full FEM.
"""
import time
import numpy as np
from scipy.sparse.linalg import LinearOperator,svds
from scipy.linalg import solve

def weighted_operator(g):
 s=g.sqrt_weights
 def mv(v):return s*g.apply_GD(np.asarray(v).reshape(-1)/s)
 def rmv(v):return g.apply_GDH(s*np.asarray(v).reshape(-1))/s
 def mm(v):return s[:,None]*g.apply_GD(np.asarray(v)/s[:,None])
 def rmm(v):return g.apply_GDH(s[:,None]*np.asarray(v))/s[:,None]
 return LinearOperator((g.nq,g.nq),matvec=mv,rmatvec=rmv,matmat=mm,rmatmat=rmm,dtype=np.complex128)

def build_modes(g,rank=48):
 started=time.perf_counter();cpu=time.process_time();C=weighted_operator(g)
 rng=np.random.default_rng(44103001);v0=rng.normal(size=g.nq)+1j*rng.normal(size=g.nq)
 U,s,Vh=svds(C,k=rank,which='LM',tol=1e-6,maxiter=600,v0=v0)
 order=np.argsort(s)[::-1];U,s,V=U[:,order],s[order],Vh[order].conj().T
 # The right singular basis is a fixed current chart; no truth field enters.
 Q=V;CQ=C.matmat(Q)
 return dict(Q=Q,CQ=CQ,singular_values=s,left_modes=U,
             elapsed_seconds=time.perf_counter()-started,cpu_seconds=time.process_time()-cpu,
             basis='right singular GD modes in physical quadrature L2',rank=rank)

def compare_state(g,z,modes,sigma=.01,strong_ranks=(8,16,32)):
 chi=g.chi(z);Q=modes['Q'];CQ=modes['CQ'];n=Q.shape[1]
 M=np.eye(n,dtype=complex)-Q.conj().T@(chi[:,None]*CQ)
 b=Q.conj().T@(chi[:,None]*(g.sqrt_weights[:,None]*g.inc))
 p=solve(M,b);jt=Q@p;j=jt/g.sqrt_weights[:,None]
 exact=g.physical_state(z);d=g.ev.direct+g.dA+g.apply_GS(j)
 state_residual=j-chi[:,None]*(g.inc+g.apply_GD(j))
 fieldrel=np.linalg.norm(g.inc+g.apply_GD(j)-exact['field'])/max(np.linalg.norm(exact['field']),1e-30)
 # Full current/GD truth is used only for this independent operator audit.
 out=[]
 for h in strong_ranks:
  if h>=n:continue
  HH,HL,LH,LL=M[:h,:h],M[:h,h:],M[h:,:h],M[h:,h:]
  ilh=solve(LL,LH);ilb=solve(LL,b[h:]);Schur=HH-HL@ilh;rhs=b[:h]-HL@ilb
  ph=solve(Schur,rhs);pl=ilb-ilh@ph;ps=np.vstack((ph,pl))
  dropped_h=solve(HH,b[:h]);dropped_l=solve(LL,b[h:]-LH@dropped_h)
  jd=Q@np.vstack((dropped_h,dropped_l))/g.sqrt_weights[:,None]
  dd=g.ev.direct+g.dA+g.apply_GS(jd)
  out.append(dict(strong_rank=h,weak_rank=n-h,schur_relative_to_full_projected=float(np.linalg.norm(ps-p)/max(np.linalg.norm(p),1e-30)),
      weak_block_condition=float(np.linalg.cond(LL)),strong_schur_condition=float(np.linalg.cond(Schur)),
      cross_feedback_norm=float(np.linalg.norm(HL@ilh)),
      dropped_feedback_white_error_to_projected=float(np.linalg.norm(dd-d)/sigma),
      dropped_feedback_white_error_to_fullFEM=float(np.linalg.norm(dd-exact['data'])/sigma)))
 return dict(frequency=g.p.freq,nq=g.nq,rank=n,projected_white_error_to_fullFEM=float(np.linalg.norm(d-exact['data'])/sigma),
       projected_field_relative_error=float(fieldrel),current_equation_L2_residual=float(np.linalg.norm(g.sqrt_weights[:,None]*state_residual)),
       current_equation_relative_residual=float(np.linalg.norm(g.sqrt_weights[:,None]*state_residual)/max(np.linalg.norm(g.sqrt_weights[:,None]*j),1e-30)),
       schur_ablation=out,scope='actual GD SVD; exact weak Schur within retained current space; full-space approximation separately audited'),dict(projected_data=d,full_data=exact['data'],M=M,b=b,p=p,singular_values=modes['singular_values'],full_current=exact['current'],projected_current=j)
