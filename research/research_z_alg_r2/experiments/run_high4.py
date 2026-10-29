"""Separate frozen higher-contrast blind root enumeration; no private access."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[name]='1'
import sys,json,time
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
P=Path(__file__).resolve().parents[1];sys.path.insert(0,str(P/'src'))
from physics import Model,ReducedKernel
from algebraic_solver import validate_candidate
from branch_enumerator import enumerate_bilinear_dual,filter_bilinear_dual_roots
from validate import CountingBook

def encode(x):
 if isinstance(x,np.ndarray):return x.tolist()
 if isinstance(x,np.generic):return x.item()
 if isinstance(x,complex):return {'real':x.real,'imag':x.imag}
 raise TypeError(type(x).__name__)
def write(path,row):path.write_text(json.dumps(row,indent=2,default=encode)+'\n')
def qualify(model,y,candidates,box):
 out=[];seen=[]
 for root in candidates:
  m=np.asarray(root['m'])
  if any(np.linalg.norm(m-z)/2<=1e-7 for z in seen):continue
  seen.append(m);v=validate_candidate(model,m,y,box=box,response_rtol=1e-8)
  out.append({'m':m,'validation':v,'source_path':root.get('path'),'status':'QUALIFIED' if v['qualified'] else 'REJECTED_INDEPENDENT_FORWARD'})
 return out

def run(case,c,y,method):
 outpath=P/f'results/high4_{case}_{method}.json'
 if outpath.exists():raise RuntimeError('Frozen prediction exists; no automatic rerun')
 book=CountingBook();model=Model(book,**c['geometry']);t=time.perf_counter();cpu=time.process_time();box=tuple(c['material_box']);out={'case':case,'method':method,'status':'RUNNING','box':box}
 try:
  if method=='six':
   kernel=ReducedKernel(model,np.full(4,box[0]),active=[0,1,2,3]);result=enumerate_bilinear_dual(kernel,y,**c['homotopy']);filtered=filter_bilinear_dual_roots(result,box=box,validate=False)
   out.update(homotopy={k:v for k,v in result.items() if k!='system'},physical_filter=filtered,independently_validated=qualify(model,y,[r for r in filtered['roots'] if r['status']=='CANDIDATE'],box))
  else:
   target=model.pack(y);scale=max(np.linalg.norm(target),1e-300);starts=[];candidates=[]
   for i,x0 in enumerate(c['initialization']['FullJ8_starts']):
    cache={}
    def state(x):
     if 'x' not in cache or not np.array_equal(x,cache['x']):cache.update(x=x.copy(),state=model.forward(x))
     return cache['state']
    def fun(x):return (model.pack(state(x)['y'])-target)/scale
    def jac(x):return model.jacobian(state(x))/scale
    fit=least_squares(fun,x0,jac=jac,bounds=box,method='trf',max_nfev=120,ftol=1e-12,xtol=1e-12,gtol=1e-12,x_scale=1.,loss='linear');relative=float(np.linalg.norm(fit.fun));starts.append({'start':i,'x0':x0,'m':fit.x,'scipy_success':bool(fit.success),'scipy_status':fit.status,'nfev':fit.nfev,'njev':fit.njev,'response_relative':relative,'message':fit.message})
    if relative<=1e-8:candidates.append({'m':fit.x,'path':i})
   out.update(starts=starts,independently_validated=qualify(model,y,candidates,box))
  out['status']='COMPLETE'
 except BaseException as exc:out.update(status='FAILED_PAID',error=repr(exc))
 out.update(method_wall_seconds=time.perf_counter()-t,method_CPU_seconds=time.process_time()-cpu,accounting=book.summary(model.ns,model.nr),model_counts=dict(model.counts),model_timings=dict(model.timings));write(outpath,out)
 print(case,method,out['status'],'qualified',sum(r['status']=='QUALIFIED' for r in out.get('independently_validated',[])))
 return out

def main():
 c=json.loads((P/'experiments/HIGH4_CONTRACT.json').read_text());manifest=json.loads((P/'experiments/high4_public/manifest.json').read_text())
 for row in manifest['cases']:
  y=np.load(P/'experiments/high4_public'/row['observation_file'])['y']
  for method in ('six','FullJ8'):run(row['case'],c,y,method)
if __name__=='__main__':main()
