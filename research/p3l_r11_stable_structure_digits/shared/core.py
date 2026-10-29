"""R11 shared exact physics, original source-major whitening, append-only receipts."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'): os.environ[k]='1'
import json,time,sys,hashlib
from pathlib import Path
import numpy as np
from scipy.linalg import lu_factor,lu_solve
P=Path(__file__).resolve().parents[1]; REPO=P.parents[1]; B1=REPO.parent
R10=B1/'b1-p3l-r10-bilateral-structure/research/p3l_r10_bilateral_structure'
R9=B1/'b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace'
R4=B1/'b1-p3l-r4-conditional-transport'
SIGMA=.00119211798514519; HOST=.1+.04j; DELTA_GT=1.3+.105j
sys.path.insert(0,str(REPO));sys.path.insert(0,str(R4/'src'))
from p1_material.prior import MaterialPrior,saturation,saturation_prime

def plain(x):
 if isinstance(x,np.ndarray):return x.tolist()
 if isinstance(x,np.generic):return x.item()
 if isinstance(x,complex):return [x.real,x.imag]
 if isinstance(x,Path):return str(x)
 if isinstance(x,dict):return {str(k):plain(v) for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [plain(v) for v in x]
 return x

def write(p,obj):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(plain(obj),indent=2,ensure_ascii=False,allow_nan=False)+'\n')

def load_setup():
 path=R9/'results/stable_roi/SETUP.npz'
 a=dict(np.load(path));a['E']=a['ER'];return a

def prior():
 cfg=json.loads((R4/'research/p3l_r4_conditional_transport/configs/TARGET.json').read_text())
 bg=np.load(R4/'research/p3l_r4_conditional_transport/inputs/background/weak.npz')
 return MaterialPrior(cfg,bg['points'])

def pack(F):return np.concatenate((F.real.T,F.imag.T),axis=1).ravel()
def white(F):return np.sqrt(2)/SIGMA*pack(F)
def white_jac(J):
 # source-major complex rows:8 blocks32. Real32 then imaginary32 eachsource.
 q=np.asarray(J).reshape(8,32,-1);return np.concatenate((q.real,q.imag),axis=1).reshape(512,-1)*np.sqrt(2)/SIGMA

def phi(F,y):return float(np.linalg.norm(F-y)**2/SIGMA**2)
def blake(p):return hashlib.blake2b(Path(p).read_bytes(),digest_size=32).hexdigest()

class Ledger:
 def __init__(self,path,context=None,cap=330000):
  self.path=Path(path);self.path.parent.mkdir(parents=True,exist_ok=True)
  if self.path.exists():raise RuntimeError('ledger already exists; preserve and use distinctjob')
  self.f=self.path.open('x');self.context=context or {};self.rhs=0;self.lu=0;self.calls=0;self.wall=0.;self.cap=cap
 def __call__(self,d):
  rhs=int(d.get('RHS',0));self.rhs+=rhs;self.lu+=int(d.get('LU',0));self.calls+=1;self.wall+=d.get('wall_seconds',0.)
  self.f.write(json.dumps(plain(dict(utc=time.time(),sequence=self.calls,context=self.context,**d)))+'\n');self.f.flush()
  if self.rhs>self.cap:raise RuntimeError('newjobRHS cap')
 def event(self,event,**d):self(dict(event=event,**d))
 def close(self):self.f.close()

class ExactSolver:
 name='FULL256'
 def __init__(self,C,O,E,Fh,sigma=SIGMA,ledger=None):
  self.C=C;self.O=O;self.E=E;self.Fh=Fh;self.sigma=sigma;self.ledger=ledger;self.calls=0
 def evaluate(self,a,delta,need_jac=False,y=None):
  ts=time.perf_counter();a=np.asarray(a,float);d=complex(delta)*a;M=np.eye(256)-d[:,None]*self.C
  fac=lu_factor(M,check_finite=False);h=lu_solve(fac,d[:,None]*self.E,check_finite=False);F=self.Fh+self.O@h;out=dict(F=F)
  if need_jac:
   e=self.E+self.C@h;T=lu_solve(fac,self.O.T,trans=1,check_finite=False).T
   Jo=np.einsum('ri,is->sri',complex(delta)*T,e,optimize=True).reshape(256,256)
   Jd=(T@(a[:,None]*e)).T.reshape(-1);out.update(J_occ=Jo,J_delta=Jd)
  self.calls+=1;diag=dict(method=self.name,call=self.calls,LU=1,RHS=40 if need_jac else 8,RHS_by_order={'256':40 if need_jac else 8},forward_RHS=8,adjoint_RHS=32 if need_jac else 0,support_size=int(np.count_nonzero(a>0)),factor_order=256,wall_seconds=time.perf_counter()-ts,fallback=False)
  out['diagnostics']=diag
  if y is not None:out['phi']=phi(F,y)
  if self.ledger is not None:self.ledger(diag)
  return out
