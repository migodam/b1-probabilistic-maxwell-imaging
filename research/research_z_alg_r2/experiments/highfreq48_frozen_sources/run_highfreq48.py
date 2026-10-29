"""18 frozen observation-only runs, shared geometry; no private imports."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[name]='1'
import sys,json,time
from pathlib import Path
from collections import Counter
import numpy as np
from scipy.optimize import least_squares
P=Path(__file__).resolve().parents[1];sys.path.insert(0,str(P/'src'))
from algebraic_solver import reconstruct_identifiable
from validate import CountingBook,validate_endpoint,write
from highfreq_model import FastModel

def accounting(book,model):
 a=book.summary(model.ns,model.nr);events=book.events;je=[e for e in events if e['kind'] in ['FULL_MATERIAL_J_TANGENT','FULL_MATERIAL_J_ADJOINT']]
 a['exact_J']={'calls':len(je),'complex_rhs':sum(e['source_or_adjoint_rhs'] for e in je),'wall_seconds':sum(e['wall_seconds'] for e in je),'CPU_seconds':sum(e['CPU_seconds'] for e in je)};a['per_tangent_coordinate_RHS']=[sum(e['kind']=='FULL_MATERIAL_J_TANGENT' for e in events)]*model.n
 return a

def main():
 c=json.loads((P/'experiments/HIGHFREQ48_CONTRACT.json').read_text());manifest=json.loads((P/'experiments/highfreq48_public/manifest.json').read_text());setup_path=P/'results/highfreq48_setup.json'
 if setup_path.exists() or list((P/'results').glob('highfreq48_HF48*.json')):raise RuntimeError('Frozen highfreq run exists; preserve paid outputs')
 sb=CountingBook();model=FastModel(sb,**c['geometry']);setup=accounting(sb,model);setup['model_counts']=dict(model.counts);setup['SVD_required_by']=['exactSVD','SVD_FullJ'];setup['baseline_setup_policy']='FullJ midpoint charged Green setup only; does not require receiverSVD';write(setup_path,setup)
 for row in manifest['cases']:
  y=np.load(P/'experiments/highfreq48_public'/row['observation_file'])['y'];target=model.pack(y);scale=max(np.linalg.norm(target),1e-300)
  for arm in ('FullJ_midpoint','exactSVD','SVD_FullJ'):
   book=CountingBook();model.book=book;model.counts=Counter();model.timings=Counter();model.state=model.factor=None;t=time.perf_counter();cpu=time.process_time();output={'case':row['case'],'input_id':row['input_id'],'family':row['family'],'noise_mode':row['noise_mode'],'method':arm,'status':'RUNNING'};latest=None
   try:
    if arm in ('exactSVD','SVD_FullJ'):
     answer=reconstruct_identifiable(model,y,**c['SVD_parameters'],box=tuple(c['material_box']));raw=answer['candidate'];output['SVD_status']=answer['status'];output['SVD_diagnostics']={k:v for k,v in answer.items() if k not in ('candidate','reconstructed_current','node_field_power')}
     if raw is None:raise RuntimeError('SVD returned no candidate')
     latest=np.asarray(raw,dtype=float)
    if arm!='exactSVD':
     x0=np.full(model.n,2.25) if arm=='FullJ_midpoint' else np.clip(latest,.5+1e-10,4-1e-10);output['initialization']=x0.copy();cache={}
     def state(x):
      nonlocal latest
      latest=x.copy()
      if 'm' not in cache or not np.array_equal(x,cache['m']):cache.update(m=x.copy(),state=model.forward(x))
      return cache['state']
     def fun(x):return (model.pack(state(x)['y'])-target)/scale
     def jac(x):return model.jacobian(state(x))/scale
     fit=least_squares(fun,x0,jac=jac,method='trf',bounds=(.5,4.),max_nfev=120,ftol=1e-12,xtol=1e-12,gtol=1e-12,x_scale=1.,loss='linear');latest=fit.x
     output.update(scipy_success=bool(fit.success),scipy_status=int(fit.status),message=fit.message,nfev=fit.nfev,njev=fit.njev,optimization_response_relative=float(np.linalg.norm(fit.fun)))
    output.update(m=latest.copy(),validation=validate_endpoint(model,latest,y),status='COMPLETE')
   except BaseException as exc:
    output.update(status='FAILED_PAID',error=repr(exc))
    if latest is not None:output['m']=latest.copy()
   output.update(incremental_wall_seconds=time.perf_counter()-t,incremental_CPU_seconds=time.process_time()-cpu,accounting=accounting(book,model),model_counts=dict(model.counts),model_timings=dict(model.timings));write(P/f"results/highfreq48_{row['input_id']}_{arm}.json",output)
   print(row['input_id'],arm,output['status'],'data',output.get('validation',{}).get('response_relative'),'nfev',output.get('nfev'),flush=True)
if __name__=='__main__':main()
