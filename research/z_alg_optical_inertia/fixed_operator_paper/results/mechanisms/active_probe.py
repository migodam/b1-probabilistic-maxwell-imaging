from pathlib import Path
import numpy as np
from scipy.special import hankel1
from scipy.linalg import eigvalsh
np.set_printoptions(precision=6,suppress=True)
rng=np.random.default_rng(20261126)
nside=4; lam=1.; k=2*np.pi; side=.7; step=side/nside; area=step**2
axis=(np.arange(nside)-(nside-1)/2)*step
pos=np.array([(x,y) for y in axis for x in axis]); n=len(pos)
dist=np.linalg.norm(pos[:,None]-pos[None,:],axis=2); dist[np.diag_indices(n)]=1
G=1j*k*k*area/4*hankel1(0,k*dist)
aeq=np.sqrt(area/np.pi); np.fill_diagonal(G,1j*np.pi*k*aeq/2*hankel1(1,k*aeq)-1)
p=12; recang=np.linspace(-np.pi/2,np.pi/2,p); rec=2*np.column_stack((np.cos(recang),np.sin(recang)))
S=1j*k*k*area/4*hankel1(0,k*np.linalg.norm(rec[:,None]-pos[None,:],axis=2))
d=20; ang=np.arange(d)*2*np.pi/d
Bd=np.exp(1j*k*pos@np.array([np.cos(ang),np.sin(ang)]))/np.sqrt(n)
def herm(A):return (A+A.conj().T)/2
def form(lo,hi,mu):
 m=(lo+hi)/2; h=(hi-lo)/2; U=np.eye(n)-np.diag(m)@G
 K=herm(S.conj().T@S+mu*U.conj().T@U-mu*G.conj().T@np.diag(h*h)@G)
 ev=eigvalsh(K); nu=np.sum(ev<0)
 N=mu*(U.conj().T@np.diag(m)+G.conj().T@np.diag(h*h))@Bd
 Cb=mu*Bd.conj().T@np.diag(m*m-h*h)@Bd
 D0=np.zeros((d+p,d+p),complex);D0[:d,:d]=Cb;D0[d:,d:]=np.eye(p)
 L0=np.column_stack((N,S.conj().T));Q=herm(D0-L0.conj().T@np.linalg.solve(K,L0))
 return K,Q,int(nu)
def response(x):return S@np.linalg.solve(np.eye(n)-np.diag(x)@G,np.diag(x)@Bd)
found=[]
for u in [1,2,3,4,6,8,10,12]:
 lo=np.zeros(n); hi=np.full(n,u)
 for mu in [1e-4,1e-3,1e-2,.1,1,10]:
  K,Q,nu=form(lo,hi,mu)
  if nu==0 or nu>8:continue
  for trial in range(15):
   x=u*rng.uniform(.1,.9,n);F=response(x)
   s=min(nu+2,10); olds=rng.choice(d,s,replace=False); V=np.eye(d)[:,olds];Y=F[:,olds]
   Z=np.vstack((V,Y)); H=herm(Z.conj().T@Q@Z); t=1e-8
   vals=eigvalsh(H-t*np.eye(s)); idx=np.sum(vals>0)
   if idx!=nu or min(abs(vals))<1e-10:continue
   C=herm(Q-Q@Z@np.linalg.solve(H-t*np.eye(s),Z.conj().T@Q))
   Cyy=C[d:,d:];cvals=eigvalsh(Cyy)
   if cvals[0]<=1e-6:continue
   new=next(j for j in range(d) if j not in olds); aa=np.eye(d)[:,new]
   f=C[d:,:d]@aa; center=-np.linalg.solve(Cyy,f)
   beta=np.real(np.vdot(f,np.linalg.solve(Cyy,f))-np.vdot(aa,C[:d,:d]@aa)+t)
   residual=np.real(np.vdot(F[:,new]-center,Cyy@(F[:,new]-center)))
   if beta>0:
    found.append(dict(u=u,mu=mu,nu=nu,old_indices=olds.tolist(),new=new,minK=float(eigvalsh(K)[0]),minCyy=float(cvals[0]),beta=float(beta),true_quad=float(residual),ellipsoid_max_radius=float(np.sqrt(beta/cvals[0])),true_norm=float(np.linalg.norm(F[:,new])),margin=float(beta-residual),min_old_gap=float(min(abs(vals)))))
    if u==3 and mu==1: np.savez(Path(__file__).resolve().parent / 'active_tm_example.npz',G=G,S=S,Bd=Bd,K=K,Q=Q,x=x,F=F,olds=olds,new=new,t=t,C=C,center=center,beta=beta)
    break
print(found[:12]);print('count',len(found))
