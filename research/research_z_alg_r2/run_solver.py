#!/usr/bin/env python3
"""Observation-only public CLI for algebraic Maxwell material inversion."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[name]='1'
import argparse,json,pathlib,sys,time
import numpy as np
from scipy.optimize import least_squares
P=pathlib.Path(__file__).resolve().parent
sys.path[:0]=[str(P/'src'),str(P/'experiments')]
from physics import ReducedKernel
from highfreq_model import FastModel
from algebraic_solver import reconstruct_identifiable,validate_candidate
from branch_enumerator import enumerate_bilinear_dual,filter_bilinear_dual_roots
from quotient_solver import enumerate_quotient
from validate import CountingBook

def encode(x):
 if isinstance(x,np.ndarray):return x.tolist()
 if isinstance(x,np.generic):return x.item()
 if isinstance(x,complex):return {'real':x.real,'imag':x.imag}
 if isinstance(x,pathlib.Path):return str(x)
 raise TypeError(type(x).__name__)

def main():
 ap=argparse.ArgumentParser(description=__doc__)
 ap.add_argument('--geometry',required=True,help='JSON with geometry={centers,sources,receivers,radius,k}, or that dictionary itself')
 ap.add_argument('--observations',required=True,help='NPZ containing complex y[receiver,source]')
 ap.add_argument('--output',required=True,help='New JSON path; existing files never overwritten')
 ap.add_argument('--method',choices=['auto','quotient','six','svd','fullj','hybrid'],default='auto')
 ap.add_argument('--box',nargs=2,type=float,default=[.5,4.])
 ap.add_argument('--active',nargs='+',type=int,help='Exactly4 independent active cells for quotient; other cells are fixed by --reference')
 ap.add_argument('--reference',help='JSON array containing known fixed/reference materials; default box midpoint')
 ap.add_argument('--relative-noise',type=float,default=0.,help='Only changes reported fit qualification, never accesses GT')
 args=ap.parse_args();out=pathlib.Path(args.output)
 if out.exists():raise FileExistsError('Preserving existing output: '+str(out))
 config=json.loads(pathlib.Path(args.geometry).read_text());geometry=config.get('geometry',config)
 y=np.load(args.observations)['y'];book=CountingBook();model=FastModel(book,**geometry);start=time.perf_counter();box=tuple(args.box)
 if y.shape!=(model.nr,model.ns):raise ValueError('y shape must be (receivers,sources)')
 reference=np.full(model.n,sum(box)/2) if not args.reference else np.asarray(json.loads(pathlib.Path(args.reference).read_text()))
 method=args.method
 if method=='auto':method='quotient' if model.n==4 and model.nr==1 and model.ns==2 else ('hybrid' if model.nr>=model.n else 'fullj')
 result={'method':method,'observations':args.observations,'geometry':args.geometry,'box':box,'scope':'finite discrete scalar-TM pulse VIE; no continuous-PDE or global material coverage claim','status':'RUNNING','candidates':[]}
 try:
  if method in ('quotient','six'):
   active=np.arange(model.n) if args.active is None else np.asarray(args.active)
   if len(active)!=4:raise ValueError('This quotient requires exactly4 active cells; other dimensions use an explicit declared slice')
   kernel=ReducedKernel(model,reference,active=active)
   if method=='quotient':
    try:
     r=enumerate_quotient(kernel,y)
     if not r['all_paths_numerically_resolved']:raise ArithmeticError('Unresolved quotient endpoints')
    except (ValueError,ArithmeticError,np.linalg.LinAlgError) as exc:
     result['fallback_reason']=repr(exc);r=enumerate_bilinear_dual(kernel,y)
   else:r=enumerate_bilinear_dual(kernel,y)
   f=filter_bilinear_dual_roots(r,box=box,validate=False)
   result['algebraic_result']={k:v for k,v in r.items() if k!='system'};result['physical_filter']=f
   for row in f['roots']:
    if row['status']=='CANDIDATE':result['candidates'].append({'m':row['m'],'validation':validate_candidate(model,row['m'],y,box=box,response_rtol=1e-8)})
  else:
   m0=reference.copy()
   if method in ('svd','hybrid'):
    direct=reconstruct_identifiable(model,y,box=box);result['direct_status']=direct['status']
    if direct['candidate'] is None:
     if method=='svd':raise ValueError(direct['status'])
    else:m0=np.asarray(direct['candidate'])
   if method!='svd':
    cache={};scale=max(np.linalg.norm(y),1e-300);target=model.pack(y)
    def state(x):
     if 'x' not in cache or not np.array_equal(cache['x'],x):cache.update(x=x.copy(),s=model.forward(x))
     return cache['s']
    fit=least_squares(lambda x:(model.pack(state(x)['y'])-target)/scale,np.clip(m0,*box),jac=lambda x:model.jacobian(state(x))/scale,bounds=box,method='trf',max_nfev=120,ftol=1e-12,xtol=1e-12,gtol=1e-12)
    m=fit.x;result['optimization']={'nfev':fit.nfev,'njev':fit.njev,'status':fit.status,'message':fit.message}
   else:m=m0
   result['candidates']=[{'m':m,'validation':validate_candidate(model,m,y,box=box,response_rtol=max(1e-8,2*args.relative_noise))}]
  result['status']='COMPLETED_NUMERICAL'
 except Exception as exc:result.update(status='FAILED_PAID',error=repr(exc))
 result['wall_seconds_after_geometry']=time.perf_counter()-start;result['accounting']=book.summary(model.ns,model.nr);result['model_counts']=dict(model.counts)
 out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(result,indent=2,default=encode)+'\n')
 print(json.dumps({'status':result['status'],'qualified':sum(c['validation']['qualified'] for c in result['candidates']),'RHS':result['accounting']['total_complex_rhs'],'output':str(out)}))
 if result['status']=='FAILED_PAID':sys.exit(1)
if __name__=='__main__':main()
