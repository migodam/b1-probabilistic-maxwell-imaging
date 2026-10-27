"""Read-only existing D anchor recovery audit; no truth/Target/solver imports."""
from pathlib import Path
import importlib.util,json,collections,time
import numpy as np
ROOT=Path.cwd();job='D_anchor_R6_DLOW_0000_full';folder=ROOT/'raw/jobs'/job;dest=ROOT/'track_d/contexts/DLOW/R6_DLOW_0000/native96_probit64'
closed=json.loads((folder/'CLOSED.json').read_text());start=json.loads((folder/'START.json').read_text());report=json.loads((dest/'ANCHOR_REPORT.json').read_text());snapshot=json.loads((ROOT/'RESOURCE_LEDGER_SNAPSHOT.json').read_text());sums=collections.Counter();jobs=collections.defaultdict(collections.Counter);activej=set();activea=set();errors=[]
with (ROOT/'RESOURCE_LEDGER.jsonl').open() as f:
 for line in f:
  r=json.loads(line)
  if r['stage']!='D':continue
  sums['RHS']+=r['RHS'];sums['states']+=r['states'];sums['LU']+=r['LU'];sums['forward_RHS' if r['kind']=='forward' else 'derivative_RHS']+=r['RHS'];jobs[r['job']]['RHS']+=r['RHS'];jobs[r['job']]['states']+=r['states'];jobs[r['job']]['forward_RHS' if r['kind']=='forward' else 'derivative_RHS']+=r['RHS']
  if r['event']=='JOB_BEGIN':activej.add(r['job'])
  if r['event']=='JOB_END':activej.remove(r['job'])
  if r['event']=='PHYSICS_BEGIN':activea.add(r['action'])
  if r['event']=='PHYSICS_END':activea.remove(r['action'])
# Only online observation/context plus already saved anchor endpoints/cache.
with np.load(dest/'CONTEXT.npz',allow_pickle=False) as z:context={k:z[k].copy() for k in z.files}
cachepath=dest/'original_optimizer/probit/RAW_MODE_CACHE.npz'
with np.load(cachepath,allow_pickle=False) as z:cache={k:z[k].copy() for k in z.files}
with np.load(ROOT/'raw/online/R6_DLOW_0000.npz',allow_pickle=False) as z:online={k:z[k].copy() for k in z.files}
for ca,co in [('q','q0'),('chi','chi0'),('J','J'),('prediction','prediction')]:
 if not np.array_equal(cache[ca],context[co]):errors.append('context_cache_'+ca)
if not np.array_equal(context['observation'],online['observation']):errors.append('online_observation')
if not np.array_equal(context['residual'],online['observation']-cache['prediction']):errors.append('residual')
# Static decoder copy only.
src=ROOT/'inputs/Research_P3LR5T/decoder_sources';spec=importlib.util.spec_from_file_location('static_prior_recovery',src/'prior.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);cfg=json.loads((src/'ONLINE_TARGET.json').read_text());points=json.loads((src/'geometry.json').read_text())['points'];prior=mod.MaterialPrior(cfg,points)
chi=prior.decode(cache['q'])[0]
if not prior.legal(cache['q']) or not np.array_equal(chi,cache['chi']):errors.append('decoder')
background=ROOT.parents[2]/'b1-p3l-r4-conditional-transport/research/p3l_r4_conditional_transport/inputs/background/weak.npz'
with np.load(background,allow_pickle=False) as z:G,S,incident=(z[k].copy() for k in ('G','S','incident'))
current=cache['current'];field=cache['field'];T=cache['receiver_T'];rhs=chi[:,None]*incident
relative=lambda a,b:float(np.linalg.norm(a-b)/max(1e-30,np.linalg.norm(b)))
algebra=dict(current_backward=relative(current-chi[:,None]*(G@current),rhs),field=relative(field,incident+G@current),prediction=relative(cache['prediction'],S@current),receiver_transfer=relative(T-(T*chi[None,:])@G,S))
C=np.einsum('qp,ps->qsp',T,field).reshape(256,576);a=C@prior.derivative(cache['q']);sigma=float(online['sigma']);rebuilt=np.column_stack([np.concatenate((a[:,k].reshape(32,8).real.T,a[:,k].reshape(32,8).imag.T),axis=1).ravel()*np.sqrt(2)/sigma for k in range(66)]);algebra['J_relative']=relative(rebuilt,cache['J'])
phi=float(np.vdot(cache['prediction']-online['observation'],cache['prediction']-online['observation']).real/sigma**2);algebra['phi_absolute_difference']=abs(phi-float(context['phi']))
endpoint_res=[]
for p in (folder/'endpoints').glob('*.npz'):
 with np.load(p,allow_pickle=False) as z:endpoint_res.append(float(z['backward_residual']))
native=report['original_native_report'];probit=report['original_probit_report'];trace=(ROOT/'track_d/EXECUTION.log').read_text();expectedjob=jobs[job]
for key in ['RHS','states','LU']:
 if sums[key]!=snapshot['stages']['D'][key]:errors.append('D_ledger_snapshot_'+key)
if activej or activea or snapshot['active_jobs'] or snapshot['active_actions']:errors.append('live_action_or_job')
if expectedjob['RHS']!=5120 or sums['RHS']!=9216 or len(endpoint_res)!=156:errors.append('unexpected_physics_count')
if max(algebra[k] for k in ('current_backward','field','prediction','receiver_transfer','J_relative'))>1e-10:errors.append('saved_cache_equation')
if report.get('path_complete') is not True or report.get('GT_read') is not False:errors.append('report_completion')
if native['fresh_forward_state_cap']!=96 or probit['fresh_state_attempts']>64:errors.append('optimizer_caps')
result=dict(status='RECOVERABLE_WITH_ZERO_NEW_PHYSICS' if not errors else 'AUDIT_ERROR',errors=errors,serialization_failure_confirmed='TypeError: Object of type ndarray is not JSON serializable' in trace,original_CLOSED=closed,original_START=start,
 original_CLOSED_must_remain_FAILED=True,original_ANCHOR_REPORT=dict(path_complete=report['path_complete'],GT_read=report['GT_read'],attempted_forward_states=report['attempted_forward_states'],path_resource_cost=report['path_resource_cost'],native_states=native['fresh_forward_states'],native_cap=native['fresh_forward_state_cap'],native_status=native['status'],probit_states=probit['fresh_state_attempts'],probit_optimizer=probit['optimizer']),
 D_stage_counters=dict(sums),D_jobs={k:dict(v) for k,v in jobs.items()},active_D_jobs=list(activej),active_D_actions=list(activea),endpoint_count=len(endpoint_res),max_saved_endpoint_backward_residual=max(endpoint_res),saved_cache_algebra=algebra,context_keys_shapes={k:list(v.shape) for k,v in context.items()},
 recovery=dict(new_RHS_required=0,may_reconstruct_outer_REPORT_from='existing ANCHOR_REPORT + immutable CONTEXT/online observation + retained original CLOSED failure receipt',must_not_replay_optimizer=True,must_not_change_GT_or_statistics=True,
 timing_qualification='Exact pre-serialization cleanup/end_to_end expression was not persisted. Existing ANCHOR_REPORT wall excludes subsequent cleanup; CLOSED wall is conservative including serialization failure and Book overhead. A recovered timing must explicitly bound/qualify this parent or retain conservative CLOSED wall; never invent exact lost cleanup time.'),GT_arrays_opened=False,new_Maxwell_RHS=0,new_training=0)
(ROOT/'audit/D_SERIALIZATION_RECOVERY_QA.json').write_text(json.dumps(result,indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else None));print(json.dumps(dict(status=result['status'],errors=errors,D=dict(sums),job=dict(expectedjob),algebra=algebra,max_backward=max(endpoint_res)),indent=2))
