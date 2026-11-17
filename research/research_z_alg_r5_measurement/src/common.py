"""R5 infrastructure. Historical R2 modules are imported read-only."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'): os.environ[key]='1'
import sys,json,time,hashlib,resource
from pathlib import Path
from contextlib import contextmanager
from collections import Counter
import numpy as np
from scipy.linalg import lu_factor,lu_solve,null_space
from scipy.special import j0
ROOT=Path(__file__).resolve().parents[1]
R2=ROOT.parent/'research_z_alg_r2'
sys.path.insert(0,str(R2))
from src.physics import Model,ReducedKernel
from src.general_quotient_solver import enumerate_general_quotient,filter_general_roots,GeneralBilinearDualSystem
from src.quotient_solver import enumerate_quotient
BOX=(.5,30.)
EPS=1e-4
MODEL_ALLOWANCE=2e-5
class BudgetExceeded(RuntimeError):pass
class Ledger:
 def __init__(self,label):
  self.label=label;self.events=[];self.counts=Counter();self.start=time.process_time();self.wall=time.monotonic()
 @contextmanager
 def op(self,kind,rhs=0,**meta):
  a=time.process_time();b=time.monotonic()
  self.counts[kind]+=1;self.counts['actual_complex_RHS']+=rhs
  for k,v in meta.items():
   if isinstance(v,(int,float)) and not isinstance(v,bool):self.counts[k]+=v
  try:yield None
  finally:self.events.append(dict(kind=kind,rhs=rhs,CPU_seconds=time.process_time()-a,wall_seconds=time.monotonic()-b,**meta))
 def receipt(self):return dict(label=self.label,CPU_seconds=time.process_time()-self.start,wall_seconds=time.monotonic()-self.wall,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,counts=dict(self.counts),events=self.events)
def serial(x):
 if isinstance(x,np.ndarray):return serial(x.tolist())
 if isinstance(x,(complex,np.complexfloating)):return dict(real=float(x.real),imag=float(x.imag))
 if isinstance(x,np.generic):return x.item()
 if isinstance(x,dict):return {str(k):serial(v) for k,v in x.items() if k!='system'}
 if isinstance(x,(list,tuple)):return [serial(v) for v in x]
 return x
def save(path,x):
 path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
 if path.exists():raise FileExistsError(path)
 path.write_text(json.dumps(serial(x),indent=2,allow_nan=True))
def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def base_geometry(n=6):
 g=json.loads((R2/'experiments/GENERAL6_CONTRACT.json').read_text())['geometry']
 if n==4:g['centers']=[g['centers'][i] for i in (0,1,3,4)];g['sources']=g['sources'][:2]
 return g
def geometry(case):
 n=4 if case in ('T1','T5') else 6;g=base_geometry(n)
 if case=='T3':
  p=np.array(g['sources'][0]);g['sources']=[(p+np.array([0.,j*.002])).tolist() for j in range(3)]
 if case=='T5':g['sources']=(30*np.asarray(g['sources'])/np.linalg.norm(g['sources'],axis=1)[:,None]).tolist()
 return g
def dictionary(case,g):
 s=len(g['sources']);theta=np.linspace(0,2*np.pi,12,endpoint=False)+.137
 if case=='T4':theta=np.linspace(-2.43,-2.41,12)
 if case=='T5':r=30.
 else:r=1.05
 positions=np.vstack([g['sources'],np.c_[r*np.cos(theta),r*np.sin(theta)]])
 Q=j0(g['k']*np.linalg.norm(positions[:,None]-positions[None,:],axis=2))
 return positions,Q
class Physics:
 def __init__(self,g,positions,book):
  self.book=book;self.model=Model(book,centers=g['centers'],sources=positions,receivers=g['receivers'],radius=g['radius'],k=g['k'])
  self.n=self.model.n;self.L=len(positions);self.s=len(g['sources']);self.H=self.model.G;self.B=self.model.E0;self.R=self.model.S
  # All cells active in Gate1, vacuum host. Independent nonzero-host tests cover host reduction.
  self.host=ReducedKernel(self.model,np.zeros(self.n))
 def response(self,m,A,jac=False,independent=False):
  A=np.asarray(A);A=A[:,None] if A.ndim==1 else A
  E=self.B@A;n=self.n
  with self.book.op('FULL_INDEPENDENT' if independent else 'REDUCED_FORWARD',rhs=A.shape[1] if independent else 0,**({'LU':1} if independent else {'reduced_LU':1,'reduced_complex_RHS':A.shape[1]+(1 if jac else 0)})):
   M=np.eye(n)-m[:,None]*self.H
   if independent:Jcur=np.linalg.solve(M,m[:,None]*E)
   else:
    lu=lu_factor(M);Jcur=lu_solve(lu,m[:,None]*E)
   y=(self.R@Jcur).reshape(-1)
   if not jac:return y
   transfer=lu_solve(lu,self.R.T,trans=1).reshape(-1)
   field=E+self.H@Jcur;J=(transfer[:,None]*field).T
   return y,np.stack([J.real,J.imag],axis=1).reshape(-1,n)
 def dual(self,m):
  with self.book.op('DUAL_ROW',reduced_LU=1,reduced_complex_RHS=1):
   return np.linalg.solve((np.eye(self.n)-m[:,None]*self.H).T,self.R.T).reshape(-1)*m
 def old_kernel(self):
  # Share identical G/E/S and host coefficients; avoid repeated host setup.
  import copy
  model=copy.copy(self.model);model.ns=self.s;model.sources=model.sources[:self.s];model.E0=model.E0[:,:self.s]
  kernel=copy.copy(self.host);kernel.model=model;kernel.q=kernel.q[:,:self.s];kernel.y0=kernel.y0[:,:self.s];kernel.j0=kernel.j0[:,:self.s]
  return kernel
 def normalize(self,a,Q):
  p=np.real(np.vdot(a,Q@a))
  if p<=1e-12:return None
  a=a/np.sqrt(p)
  if np.linalg.norm(a)>4+1e-10:return None
  return a
 def menu(self,Q,roots,gn):
  L=self.L;s=self.s;items=[]
  def add(label,a):
   a=self.normalize(a,Q)
   if a is not None:items.append(dict(label=label,a=a,energy=float(np.real(np.vdot(a,Q@a))),feed_norm=float(np.linalg.norm(a))))
  a=np.zeros(L,complex);a[:s]=1;add('old_sum',a)
  a=np.zeros(L,complex);a[:s]=np.exp(1j*np.arange(s));add('old_phase',a)
  for j in range(s,L):add('primitive_%02d'%(j-s),np.eye(L,dtype=complex)[:,j])
  rng=np.random.default_rng(18112026)
  for j in range(12):
   a=np.zeros(L,complex);a[s:]=rng.normal(size=L-s)+1j*rng.normal(size=L-s);add('random_%02d'%j,a)
  # Pair matched filters are conventional power-Gram discriminating directions.
  for i in range(len(roots)):
   for j in range(i):
    if len(items)>=44:break
    d=(self.dual(roots[i])-self.dual(roots[j]))@self.B
    a=np.linalg.pinv(Q,rcond=1e-10)@d.conj();add('pair_%d_%d'%(j,i),a)
  if gn is not None:
   # Local weak material direction from old-data J; realizable matched source.
   _,J=self.response(gn,np.eye(L)[:,:s],True);v=np.linalg.svd(J,full_matrices=False)[2][-1]
   _,Jall=self.response(gn,np.eye(L),True);Jc=Jall.reshape(L,2,self.n);d=(Jc[:,0]+1j*Jc[:,1])@v
   add('fisher_weak_mode',np.linalg.pinv(Q,rcond=1e-10)@d.conj())
  return items[:48]
def pack(y):return np.stack([np.real(y),np.imag(y)],axis=-1).reshape(-1)
def starts(n):
 from scipy.stats import qmc
 points=qmc.Sobol(n,scramble=False).random_base2(7)
 pts=[np.full(n,15.25)]
 for p in points[1:]:
  m=.5+29.5*p
  if not any(np.array_equal(m,x) for x in pts):pts.append(m)
  if len(pts)==64:break
 return pts
