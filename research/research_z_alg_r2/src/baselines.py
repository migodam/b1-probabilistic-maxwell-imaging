"""Frozen observation-only midpoint TRF baselines. Imports no private module."""
import argparse,json,time
import numpy as np
from scipy.optimize import least_squares
from physics import Model
from validate import P,CountingBook,validate_endpoint,load_case,write

def solve(case,frozen=False):
 contract,g,y=load_case(case);method='frozen_J_TRF' if frozen else 'FullJ_adjoint_TRF';output_path=P/'results'/f'blind_{case}_{method}.json'
 if output_path.exists():raise RuntimeError('Frozen result already exists; preserve the original run')
 book=CountingBook();model=Model(book,**g);target=model.pack(y);m0=np.full(model.n,2.25)
 cache={};scale=max(float(np.linalg.norm(target)),1e-12)
 def state(x):
  if 'm' not in cache or not np.array_equal(x,cache['m']):cache.update(m=x.copy(),s=model.forward(x))
  return cache['s']
 def fun(x):return (model.pack(state(x)['y'])-target)/scale
 def jac(x):
  s=state(x)
  return model.jacobian(s)/scale
 J0=jac(m0).copy() if frozen else None
 t=time.perf_counter();cpu=time.process_time()
 try:
  fit=least_squares(fun,m0,jac=(lambda x:J0) if frozen else jac,method='trf',bounds=(.5,4),max_nfev=120,ftol=1e-12,xtol=1e-12,gtol=1e-12,x_scale=1.,loss='linear')
  m=fit.x;output={'case':case,'method':method,'status':'COMPLETE','m':m,'scipy_success':bool(fit.success),'scipy_status':int(fit.status),'message':fit.message,'nfev':fit.nfev,'njev':fit.njev,'optimization_wall_seconds':time.perf_counter()-t,'optimization_CPU_seconds':time.process_time()-cpu,'validation':validate_endpoint(model,m,y)}
 except BaseException as e:output={'case':case,'method':method,'status':'FAILED','error':repr(e)}
 output['accounting']=book.summary(model.ns,model.nr);output['model_counts']=dict(model.counts);output['model_timings']=dict(model.timings);write(output_path,output);return output
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--case');a.add_argument('--frozen',action='store_true');args=a.parse_args()
 cases=[args.case] if args.case else [r['case'] for r in json.loads((P/'experiments/public/manifest.json').read_text())['cases']]
 for case in cases:
  out=solve(case,args.frozen);print(case,out['method'],out['status'],out.get('validation',{}).get('response_relative'))
