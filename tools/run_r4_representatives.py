#!/usr/bin/env python3
"""Prepare or EXPLICITLY run bounded Output3 material representatives.

--prepare alone writes the interface, no physics. --prepare --manifest freezes
online estimator raw/clipped arrays and qualified anchor metadata only. Actual
physics requires a separate --run --prepared-request ... --job ... invocation.
No GT file/key is imported. Offline GT evaluation is deliberately separate.
"""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[name]='1'
from pathlib import Path
import sys,json,argparse,time,traceback,gc,resource,re
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from p3r4.core import OUT,Book,Target,BudgetStop,write,save,plain,to_probit,from_probit
from p3r4.representative import PHI_THRESHOLD,GRADIENT_MODES,LocalBudgetStop,solve_representative,material_objective_gradient
PROTOCOL=OUT/'configs/REPRESENTATIVE_PROTOCOL.json'
ALLOWED_MU_KEYS=('prediction_raw','prediction_clipped')


def protocol():
 return dict(status='IMPLEMENTED_NOT_PHYSICALLY_RUN',output='Output3 candidate representative, not posterior mean/sample',
  objective='mean full-grid |exact original decoder(x)-frozen mu|^2',
  constraint='original Full nonlinear Phi(x)<=282.873792, strict final check',Phi_threshold=PHI_THRESHOLD,
  coordinates='full66 probit; exact native Dg times du/deta normal pdf; no prior clipping or new prior',
  initialization='existing Full-qualified online anchor only; exact finite saved probit context preferred',
  optimizer='SLSQP; material objective/gradient decoder-only; Full likelihood and explicit gradient query path',
  gradient_modes=GRADIENT_MODES,default_gradient='source_adjoint8',
  cached_anchor_gradient='reuse already-paid receiver T/J by matrix composition when available; otherwise explicitly charged gradient',
  per_parent_per_arm=dict(max_new_states=16,max_complex_RHS=256),
  solver_iterations_max=40,restoration_scales=[.999,.5],restoration_budget='reserve at most2 states inside SAME16/256 cap',
  final_selection='original anchor or accepted SLSQP feasible iterates; optional separately Full-verified accepted restoration; rejected/infeasible trials never final',
  fallback='retain original qualified anchor if no feasible improvement; expose fallback and optimizer status',
  KKT='numeric lambda>=0, stationarity, feasibility, complementarity, active flag; diagnostic only',
  ordinary_control='independent same-cap extra Phi minimization, no mu in proposals or selection; shared once per parent across raw/clipped mu',
  allowed_mu_keys=ALLOWED_MU_KEYS,mu='use frozen provided raw/clipped outputs exactly; never silently clip inside optimization',
  material_GT_import=False,new_physics_in_prepare=False,offline_evaluation='separate future module',
  required_manifest_fields=['request_id','estimator_id','predictions_npz','mu_layout','parents','outputs','frozen_estimator','no_GT'],
  mu_layouts=['complex_N576','realimag_N2_24_24'],
  actual_launch='separate explicit --run --prepared-request PATH --job UNIQUE; no default launch')


def append(path,row):
 with Path(path).open('a') as f:f.write(json.dumps(plain(row),allow_nan=False)+'\n');f.flush();os.fsync(f.fileno())


def read_request(path):
 req=json.loads(Path(path).read_text())
 required=set(protocol()['required_manifest_fields'])
 if not required.issubset(req) or req['frozen_estimator'] is not True or req['no_GT'] is not True:raise ValueError('frozen no-GT estimator request required')
 if not req['request_id'].replace('_','').replace('-','').isalnum():raise ValueError('simple unique request id required')
 if not req['outputs'] or any(k not in ALLOWED_MU_KEYS for k in req['outputs']) or len(set(req['outputs']))!=len(req['outputs']):raise ValueError('raw/clipped online prediction keys only')
 if len(set(req['parents']))!=len(req['parents']) or not req['parents']:raise ValueError('unique parent allowlist required')
 if any(not re.fullmatch(r'(R4DEV|R4L1|R4L2)_[0-9]{5}',p) for p in req['parents']):raise ValueError('declared R4 online parent ids required')
 if req.get('gradient_mode','source_adjoint8') not in GRADIENT_MODES:raise ValueError('explicit gradient mode')
 if req.get('state_cap',16)!=16 or req.get('RHS_cap',256)!=256:raise ValueError('frozen16state/256RHS caps')
 return req


def prepare(args):
 if PROTOCOL.exists():
  if json.loads(PROTOCOL.read_text())!=plain(protocol()):raise ValueError('existing representative protocol differs')
 else:write(PROTOCOL,protocol())
 job='representative_interface_prepare_v1' if args.manifest is None else 'representative_prepare_'+read_request(args.manifest)['request_id']
 b=Book('II',job,max_wall=900,analysis_only=True);status='COMPLETE';error=None
 audit=OUT/'AUDIT/representative_preparation';audit.mkdir(parents=True,exist_ok=True)
 try:
  if args.manifest is None:
   write(audit/'INTERFACE_READY.json',dict(status='INTERFACE_READY_NO_ESTIMATOR_REQUEST',protocol=str(PROTOCOL),
    module='src/p3r4/representative.py',runner='tools/run_r4_representatives.py',decoder_and_mock_tests=4,
    new_Maxwell_states=0,new_complex_RHS=0,new_physical_LU=0,actual_representative_solve='NOT_RUN'))
   print(str(PROTOCOL));return
  req=read_request(args.manifest);dest=OUT/'results/II/representative_preparation'/req['request_id']
  if dest.exists():raise FileExistsError('frozen representative request exists')
  dest.mkdir(parents=True);source=Path(req['predictions_npz']).expanduser().resolve()
  with np.load(source,allow_pickle=False) as z:
   sourceparents=z['parent'].astype(str).tolist()
   if len(set(sourceparents))!=len(sourceparents):raise ValueError('prediction parent duplicates')
   lookup={p:i for i,p in enumerate(sourceparents)}
   values={key:z[key].copy() for key in req['outputs']}
  for key,a in values.items():
   if req['mu_layout']=='complex_N576':
    if not np.iscomplexobj(a) or a.shape!=(len(sourceparents),576):raise ValueError('declared complex_N576 layout')
    values[key]=a.astype(complex)
   elif req['mu_layout']=='realimag_N2_24_24':
    if a.shape!=(len(sourceparents),2,24,24) or np.iscomplexobj(a):raise ValueError('declared real/imag layout')
    values[key]=(a[:,0]+1j*a[:,1]).reshape(len(a),576)
   else:raise ValueError('explicit mu layout required')
  cfg=json.loads((OUT/'configs/TARGET.json').read_text());qualified=[];skips=[]
  for parent in req['parents']:
   if parent not in lookup:skips.append(dict(parent=parent,reason='MISSING_FROZEN_ESTIMATOR_OUTPUT'));continue
   cache=OUT/'results/contexts/cache'/f'{parent}.npz';online=OUT/'inputs/online'/f'{parent}.npz'
   if not cache.exists() or not online.exists():skips.append(dict(parent=parent,reason='MISSING_EXISTING_ONLINE_ANCHOR_OR_OBSERVATION'));continue
   with np.load(cache,allow_pickle=False) as z:q0=z['q'].copy();prediction=z['prediction'].copy();backward=float(z['backward_residual'])
   with np.load(online,allow_pickle=False) as z:
    if set(z.files)!={'observation','sigma','target_id'} or str(z['target_id'])!=cfg['target_id'] or float(z['sigma'])!=cfg['noise']['sigma']:raise ValueError('online target contract')
    y=z['observation'].copy()
   difference=prediction-y;phi=float(np.vdot(difference,difference).real/cfg['noise']['sigma']**2)
   if not np.isfinite(phi) or phi>PHI_THRESHOLD or not np.isfinite(backward) or backward>1e-10:skips.append(dict(parent=parent,reason='ANCHOR_NOT_FULL_QUALIFIED',Phi=phi,backward_residual=backward));continue
   group='dev' if parent.startswith('R4DEV_') else 'locked' if parent.startswith('R4L1_') else 'confirmation' if parent.startswith('R4L2_') else None
   context=None if group is None else OUT/'results/contexts'/group/f'{parent}.npz'
   if context is not None and context.exists():
    with np.load(context,allow_pickle=False) as z:x0=z['x'].copy()
    x_source=str(context)
   else:x0=to_probit(q0);x_source='native interior ndtri; no inward clipping'
   if x0.shape!=(66,) or not np.isfinite(x0).all() or np.max(abs(from_probit(x0)-q0))>1e-12:raise ValueError('finite exact probit anchor context required')
   arrays=dict(q0=q0,x0=x0,anchor_Phi=np.array(phi))
   for key in req['outputs']:
    mu=values[key][lookup[parent]]
    if not np.isfinite(mu).all():raise ValueError('nonfinite frozen estimator material')
    arrays[key]=mu
   snapshot=dest/'frozen_mu'/f'{parent}.npz';save(snapshot,**arrays)
   qualified.append(dict(parent=parent,anchor_cache=str(cache),online_observation=str(online),x_source=x_source,
    anchor_Phi=phi,frozen_mu=str(snapshot),outputs=req['outputs']))
  prepared=dict(status='FROZEN_READY',request=req,source_predictions=dict(path=str(source),size=source.stat().st_size,mtime_ns=source.stat().st_mtime_ns),
    qualified=qualified,skips=skips,protocol=protocol(),new_physics_RHS=0,GT_read=False,
    ordinary_control_reuse='one independent original-anchor Phi optimization per parent; all mu distances evaluated afterward without proposal influence')
  write(dest/'FROZEN_REQUEST.json',prepared);print(str(dest/'FROZEN_REQUEST.json'))
 except BaseException as exc:status='ERROR';error=repr(exc);raise
 finally:b.close(status,error)


class ArmBudget:
 def __init__(self):self.cost=dict(states=0,RHS=0,LU=0)
 def require(self,*,states,RHS):
  if self.cost['states']+states>16 or self.cost['RHS']+RHS>256:raise LocalBudgetStop('representative per-parent/arm16states/256RHS hard cap')
 def used(self):return dict(self.cost)


def run(args):
 if not args.job or args.prepared_request is None:raise ValueError('explicit unique job and prepared frozen request required')
 prepared=json.loads(args.prepared_request.read_text())
 if prepared['status']!='FROZEN_READY' or prepared['protocol']!=plain(protocol()):raise ValueError('prepared representative protocol mismatch')
 if not prepared['qualified']:raise ValueError('no qualified existing online anchors; no paid job is started')
 for item in prepared['qualified']:
  with np.load(item['frozen_mu'],allow_pickle=False) as z:q0=z['q0'].copy()
  with np.load(item['anchor_cache'],allow_pickle=False) as z:
   if not np.array_equal(q0,z['q']):raise ValueError('prepared original anchor changed')
 dest=OUT/'results/II/representatives'/args.job
 if dest.exists():raise FileExistsError('representative job output exists')
 b=Book('II',args.job,max_wall=args.max_wall)
 active=[None];original=b._update;allcost=dict(states=0,RHS=0,LU=0)
 max_total_RHS=len(prepared['qualified'])*(1+len(prepared['request']['outputs']))*256
 def tracked(event,**kw):
  if allcost['RHS']+kw.get('rhs',0)>max_total_RHS:raise BudgetStop('representative prepared totalRHS cap')
  if active[0] is not None:active[0].require(states=kw.get('states',0),RHS=kw.get('rhs',0))
  usage=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)
  if usage>3*1024**3:raise BudgetStop('representative3GB RSS cap')
  result=original(event,**kw)
  for name,arg in [('states','states'),('RHS','rhs'),('LU','lu')]:
   amount=kw.get(arg,0);allcost[name]+=amount
   if active[0] is not None:active[0].cost[name]+=amount
  return result
 b._update=tracked;status='COMPLETE';error=None;results=[];started=time.perf_counter()
 try:
  t=Target(b,keep_fields=True);dest.mkdir(parents=True)
  write(dest/'FROZEN_REQUEST.json',prepared)
  for item in prepared['qualified']:
   parent=item['parent'];folder=dest/parent;folder.mkdir()
   with np.load(item['frozen_mu'],allow_pickle=False) as z:x0=z['x0'].copy();mus={key:z[key].copy() for key in item['outputs']}
   y=t.observation(parent)
   cases=[('ordinary_likelihood_shared','ordinary_likelihood',next(iter(mus)))]+[(key,'material_representative',key) for key in mus]
   for arm,kind,key in cases:
    case=folder/arm;case.mkdir();guard=ArmBudget();active[0]=guard
    anchor=t.load_cache(item['anchor_cache'],owner='existing qualified original online context; historical acquisition retained')
    event=lambda row:append(case/'TRIAL_EVENTS.jsonl',row)
    result,point=solve_representative(t,anchor,y,mus[key],x0=x0,
      gradient_mode=prepared['request'].get('gradient_mode','source_adjoint8'),require=guard.require,used=guard.used,event=event,kind=kind)
    mu_distances={mu_key:material_objective_gradient(t.prior,result['x'],mu)[0] for mu_key,mu in mus.items()}
    result.update(parent=parent,arm=arm,estimator_id=prepared['request']['estimator_id'],mu_distances=mu_distances,
      context_cost='old paid anchor/receiver acquisition retained, not cold-free',GT_read=False)
    save(case/'FINAL.npz',x=result['x'],q=point.q,chi=point.state.chi,prediction=point.state.data,
      current=point.state.current,field=point.state.field,Phi=np.array(result['Phi']),material_target=mus[key],
      fallback=np.array(result['fallback_to_original_anchor']))
    write(case/'RESULTS.json',result);results.append(result)
    write(dest/'PROGRESS.json',dict(completed_cases=len(results),actual_cost=allcost,last={k:v for k,v in result.items() if k not in ('x','q')}))
    active[0]=None;del anchor,point; t.backend()._transfers.clear();gc.collect()
  write(dest/'RESULTS.json',dict(status='COMPLETE',cases=results,actual_cost=allcost,elapsed_seconds=time.perf_counter()-started,
    GT_read=False,output='Output3 candidate representative/control only',posterior_mean=False,posterior_sample=False))
 except BaseException as exc:
  status='STOPPED' if isinstance(exc,BudgetStop) else 'ERROR';error=repr(exc)
  write(dest/'ERROR.json',dict(status=status,error=error,traceback=traceback.format_exc(),actual_cost=allcost));raise
 finally:b.close(status,error)


def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--prepare',action='store_true');parser.add_argument('--run',action='store_true');parser.add_argument('--manifest',type=Path);parser.add_argument('--prepared-request',type=Path);parser.add_argument('--job');parser.add_argument('--max-wall',type=float,default=900)
 args=parser.parse_args()
 if args.prepare==args.run:parser.error('choose --prepare or explicit --run')
 if not 0<args.max_wall<=3600:parser.error('controlled wall must be<=3600s')
 if args.prepare:prepare(args)
 else:run(args)

if __name__=='__main__':main()
