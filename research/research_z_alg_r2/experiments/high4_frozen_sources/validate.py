"""Observation-only endpoint validation and uniform accounting; no hidden imports."""
from contextlib import contextmanager
from pathlib import Path
import json,time
import numpy as np
P=Path(__file__).resolve().parents[1]

def write(path,value):
 def encode(x):
  if isinstance(x,np.ndarray):return x.tolist()
  if isinstance(x,np.generic):return x.item()
  raise TypeError(type(x).__name__)
 path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(value,indent=2,default=encode)+'\n')

class CountingBook:
 def __init__(self,path=None):self.path=Path(path) if path else None;self.events=[];self.t=time.perf_counter();self.c=time.process_time();self.rhs=0
 @contextmanager
 def op(self,kind,rhs=0,formula=0,**metadata):
  if self.rhs+rhs>15000 or time.perf_counter()-self.t>120:raise RuntimeError('frozen per-case budget exceeded')
  self.rhs+=rhs;row={'kind':kind,'source_or_adjoint_rhs':rhs,'formula_equivalent':formula,**metadata};t=time.perf_counter();c=time.process_time()
  try:yield len(self.events)+1
  except BaseException as exc:row['error']=repr(exc);raise
  finally:
   row.update(wall_seconds=time.perf_counter()-t,CPU_seconds=time.process_time()-c);self.events.append(row)
 def append(self,row):pass
 def save(self,*args,**kwargs):pass
 def summary(self,ns=None,nr=None):
  out={'wall_seconds':time.perf_counter()-self.t,'CPU_seconds':time.process_time()-self.c,'total_complex_rhs':self.rhs,'events':self.events}
  for name,fragments in [('forward',('FULL_FORWARD',)),('adjoint_J',('FULL_MATERIAL_J_ADJOINT',)),('validation',('INDEPENDENT_FORWARD',)),('setup',('SETUP','RECEIVER_SVD'))]:
   ev=[r for r in self.events if any(x in r['kind'] for x in fragments)]
   out[name]={'calls':len(ev),'complex_rhs':sum(r['source_or_adjoint_rhs'] for r in ev),'wall_seconds':sum(r['wall_seconds'] for r in ev),'CPU_seconds':sum(r['CPU_seconds'] for r in ev)}
  if ns is not None:
   out['per_source_forward_RHS']=[out['forward']['calls']]*ns
   out['per_source_validation_RHS']=[out['validation']['calls']]*ns
  if nr is not None:out['per_receiver_adjoint_RHS']=[out['adjoint_J']['calls']]*nr
  out['adjoint_policy']='receiver transpose solves shared by all sources; J contractions include every source'
  return out

def validate_endpoint(model,m,y):
 state=model.forward(m,independent=True);res=float(np.linalg.norm(state['y']-y)/max(np.linalg.norm(y),1e-300));m=np.asarray(m)
 return {'response_relative':res,'state_relative':float(state['state_relative_residual']),'resolvent_sigma_min':float(state['resolvent_sigma_min']),'legal':bool(np.all((m>=.5-1e-10)&(m<=4+1e-10))),'response_fit_pass':bool(res<=1e-8 and state['state_relative_residual']<=1e-10 and state['resolvent_sigma_min']>=1e-6 and np.all((m>=.5-1e-10)&(m<=4+1e-10)))}

def load_case(case):
 contract=json.loads((P/'experiments/contract.json').read_text());manifest=json.loads((P/'experiments/public/manifest.json').read_text());row=next(r for r in manifest['cases'] if r['case']==case)
 return contract,contract['presets'][row['preset']],np.load(P/'experiments/public'/row['observation_file'])['y']
