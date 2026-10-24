#!/usr/bin/env python3
"""Fixed primarySpline_s1/R2 raw+clipped representatives; one sharedFull control.

prepare is observation-only; run requires separate explicit physical release;
score is independent after actual Book closure. No command runs automatically.
"""
from pathlib import Path
import argparse,json,os,sys,time,resource,traceback,gc,csv
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from scipy.special import ndtri,ndtr
from p3r4.evaluation import checkpoint_identity,write_new_json,material_losses,paired_bootstrap
S=ROOT/'research/p3l_r4_conditional_transport';KEYS=('spline_raw','spline_clipped','r2_raw','r2_clipped')
def read(p):return json.loads(Path(p).read_text())
def check(v):
 if checkpoint_identity(v['path'])!=v:raise ValueError('frozen source changed')
def online_check(manifest):
 m=read(manifest);check(m['online_gate_identity']);g=read(m['online_gate'])
 if m['parents']!=[f'R4L1_{i:05d}' for i in range(64)] or m['no_GT'] is not True or g['truth_read'] is not False or g['status']!='ONLINE_COMPLETE' or g.get('canonical_batch_size')!=64:raise ValueError('fixed first64 canonical noGT online parents required')
 if m['component_identity']!=g['freeze'].get('component_identity') or m['top_freeze_identity']!=g['freeze'].get('top_freeze_identity'):raise ValueError('representative manifest is not bound to actual parent risk freeze')
 sealed=g['sealed_artifact_identities']
 for item in m['estimators']:
  check(item['identity'])
  if item['identity'] not in sealed:raise ValueError('source prediction not sealed by actual inference')
 if [(e['label'],e['id']) for e in m['estimators']]!=[('spline','fusion_spline_s2026102501'),('r2','frozen_r2')]:raise ValueError('fixed estimator seeds required')
 if m['online_gate_identity']['mtime_ns']>=Path(manifest).stat().st_mtime_ns:raise ValueError('request must follow online completion')
 return m,g

def prepare(a):
 m,g=online_check(a.manifest)
 if a.stage!='IV' or m.get('stage')!='IV':raise ValueError('fresh L1 representative request must explicitly charge StageIV')
 cfg=read(S/'configs/TARGET.json');threshold=282.873792
 if a.out.exists():raise FileExistsError('fresh prepared request directory required')
 a.out.mkdir(parents=True);mu={}
 for e in m['estimators']:
  with np.load(e['predictions_npz'],allow_pickle=False) as z:
   parents=z['parent'].astype(str).tolist();lookup={p:i for i,p in enumerate(parents)}
   for action in ('raw','clipped'):mu[e['label']+'_'+action]={p:z['prediction_'+action][lookup[p]].copy() for p in m['parents']}
 qualified=[];skips=[];sources=[]
 for parent in m['parents']:
  try:
   cache=S/'results/contexts/cache'/(parent+'.npz');obs=S/'inputs/online'/(parent+'.npz');ctx=S/'results/contexts/locked'/(parent+'.npz')
   for p in [cache,obs,ctx]:sources.append(checkpoint_identity(p))
   with np.load(cache,allow_pickle=False) as z:q=z['q'];pred=z['prediction'];backward=float(z['backward_residual'])
   with np.load(obs,allow_pickle=False) as z:y=z['observation']
   with np.load(ctx,allow_pickle=False) as z:x=z['x']
   phi=float(np.linalg.norm(pred-y)**2/cfg['noise']['sigma']**2)
   if q.shape!=(66,) or not np.isfinite(q).all() or not np.isfinite(backward) or not np.isfinite(phi) or phi>threshold or backward>1e-10:raise ValueError('ANCHOR_NOT_FULL_QUALIFIED')
   if x.shape!=(66,) or not np.isfinite(x).all() or np.max(abs(np.r_[x[:64],ndtr(x[64:])]-q))>1e-12:raise ValueError('exact native/probit anchor mismatch')
   arrays={k:v[parent] for k,v in mu.items()}
   if any(v.shape!=(576,) or not np.isfinite(v).all() for v in arrays.values()):raise ValueError('INVALID_FROZEN_MATERIAL_ACTION')
   dest=a.out/'frozen_mu'/(parent+'.npz');dest.parent.mkdir(exist_ok=True)
   with dest.open('xb') as stream:np.savez(stream,q0=q,x0=x,anchor_Phi=np.array(phi),**arrays)
   qualified.append(dict(parent=parent,anchor_cache=str(cache.resolve()),online_observation=str(obs.resolve()),frozen_mu=str(dest.resolve()),mu_identity=checkpoint_identity(dest),anchor_Phi=phi))
  except Exception as exc:skips.append(dict(parent=parent,reason=repr(exc),status='NOT_EVALUATED_UNQUALIFIED_OR_MISSING'))
 for v in sources:check(v)
 req=dict(status='FROZEN_READY',stage=a.stage,manifest=m,manifest_identity=checkpoint_identity(a.manifest),qualified=qualified,skips=skips,source_identities=sources,GT_read=False,
  outputs=KEYS,ordinary='one ordinary_likelihood_shared per parent; no mu in its selection/proposals',state_cap_per_arm=16,RHS_cap_per_arm=256,maximum_total_RHS=64*5*256,
  arms_per_qualified_parent=5,shared_cost='ordinaryFull control charged once; old anchor/receiver cache acquisition retained and disclosed; no free setup claim',
  scope='fixed Spline2026102501 andR2seed33 actions, not ensemble image or posterior object',created_time_ns=time.time_ns())
 write_new_json(a.out/'FROZEN_REQUEST.json',req)

def run(a):
 req=read(a.prepared_request);release=read(a.release)
 request_stage=req.get('stage','II')
 if a.stage!=request_stage or release.get('stage','II')!=request_stage:raise ValueError('explicit request/release/CLI stage mismatch')
 if any(p.startswith('R4L1_') for p in req.get('manifest',{}).get('parents',[])) and request_stage!='IV':raise ValueError('fresh L1 physical validation requires StageIV')
 if release.get('released') is not True or release.get('experiment')!='fresh_locked64_representatives_v1' or release.get('prepared_request_identity')!=checkpoint_identity(a.prepared_request):raise ValueError('explicit parent-owned physical release required')
 if not a.job or a.max_wall is None or not 0<a.max_wall<=3600 or a.out.exists():raise ValueError('bounded unique physical run required')
 online_check(req['manifest_identity']['path']);check(req['manifest_identity'])
 for v in req['source_identities']:check(v)
 for item in req['qualified']:check(item['mu_identity'])
 # The only physical imports/calls are inside this separately released branch.
 from p3r4.core import Book,Target,BudgetStop,save,write
 from p3r4.representative import solve_representative
 from run_r4_representatives import ArmBudget,append
 b=Book(request_stage,a.job,max_wall=a.max_wall);cost=dict(RHS=0,states=0,LU=0);original=b._update;active=[None];status='COMPLETE';error=None;results=[]
 def tracked(event,**kw):
  if cost['RHS']+kw.get('rhs',0)>req['maximum_total_RHS']:raise BudgetStop('frozen representative totalRHS cap')
  if active[0] is not None:active[0].require(states=kw.get('states',0),RHS=kw.get('rhs',0))
  if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)>3*1024**3:raise BudgetStop('3GB RSS cap')
  value=original(event,**kw)
  for key,arg in [('RHS','rhs'),('states','states'),('LU','lu')]:
   amount=kw.get(arg,0);cost[key]+=amount
   if active[0] is not None:active[0].cost[key]+=amount
  return value
 b._update=tracked
 try:
  a.out.mkdir(parents=True);write(a.out/'FROZEN_REQUEST.json',req);t=Target(b,keep_fields=True)
  for item in req['qualified']:
   parent=item['parent'];folder=a.out/parent;folder.mkdir();y=t.observation(parent)
   with np.load(item['frozen_mu'],allow_pickle=False) as z:x0=z['x0'];mus={k:z[k] for k in KEYS}
   for arm in ['ordinary_likelihood_shared',*KEYS]:
    case=folder/arm;case.mkdir();guard=ArmBudget();active[0]=guard
    anchor=t.load_cache(item['anchor_cache'],owner='original paid online native/probit context; setup not free')
    target=mus[KEYS[0] if arm=='ordinary_likelihood_shared' else arm]
    result,point=solve_representative(t,anchor,y,target,x0=x0,gradient_mode='source_adjoint8',require=guard.require,used=guard.used,event=lambda row:append(case/'TRIAL_EVENTS.jsonl',row),kind='ordinary_likelihood' if arm=='ordinary_likelihood_shared' else 'material_representative')
    result.update(parent=parent,arm=arm,GT_read=False,actual_cost=guard.used(),posterior_mean=False,posterior_sample=False)
    save(case/'FINAL.npz',x=result['x'],q=point.q,chi=point.state.chi,current=point.state.current,field=point.state.field,prediction=point.state.data,Phi=np.array(result['Phi']),fallback=np.array(result['fallback_to_original_anchor']))
    result['endpoint_identity']=checkpoint_identity(case/'FINAL.npz');write(case/'RESULTS.json',result);results.append(result);active[0]=None;del anchor,point;t.backend()._transfers.clear();gc.collect()
    write(a.out/'PROGRESS.json',dict(cases=results,actual_cost=cost))
  write(a.out/'RESULTS.json',dict(status='COMPLETE',cases=results,actual_cost=cost,GT_read=False,shared_ordinary_cost_charged_once=True,skips=req['skips'],posterior_mean=False,posterior_sample=False,job=a.job,stage=request_stage))
 except BaseException as exc:
  status='ERROR';error=repr(exc);write(a.out/'ERROR.json',dict(error=error,traceback=traceback.format_exc(),actual_cost=cost));raise
 finally:b.close(status,error)

def score(a):
 req=read(a.prepared_request);physical=read(a.physical/'RESULTS.json');closed=read(a.receipt);start=read(a.receipt.parent/'START.json');m,g=online_check(req['manifest_identity']['path'])
 request_stage=req.get('stage','II')
 if a.stage!=request_stage or physical.get('stage','II')!=request_stage:raise ValueError('score must use explicit request stage')
 if a.out.exists() or closed['status']!='COMPLETE' or physical['status']!='COMPLETE' or physical['GT_read'] is not False or start['kind']!='physics' or start['stage']!=request_stage or start['job']!=physical['job'] or a.receipt.parent.name!=physical['job'] or start['pid']==os.getpid() or g['pid']==os.getpid():raise ValueError('independent actual closed physical/online gates required')
 if read(a.physical/'FROZEN_REQUEST.json')!=req or not g['completed_time_ns']<=req['created_time_ns']<=int(start['start_unix']*1e9):raise ValueError('online/request/physics chronology failed')
 if not (a.physical/'RESULTS.json').stat().st_mtime_ns<=a.receipt.stat().st_mtime_ns or 'offline_only' not in a.truth.resolve().parts:raise ValueError('closed result and isolated truth required')
 check(req['manifest_identity'])
 for v in req['source_identities']:check(v)
 for item in req['qualified']:check(item['mu_identity'])
 for case in physical['cases']:
  check(case['endpoint_identity'])
  if read(a.physical/case['parent']/case['arm']/'RESULTS.json')!=case:raise ValueError('perarm result differs from closed aggregate')
 # All checks precede the first offline-only truth read.
 from p3r4.probability_eval import decoder_from_geometry
 cfg=read(S/'configs/TARGET.json');points=read(ROOT/'research/p1_material_target_20261018/inputs/geometry.json')['points'];geom=decoder_from_geometry(cfg,points);roi=geom['roi'];A=geom['A']
 parents=m['parents'];actions=[*KEYS,'ordinary_extra_Phi',*['representative_'+k for k in KEYS],'original_anchor'];losses={k:[] for k in actions};records=[];qualified={v['parent']:v for v in req['qualified']};a.out.mkdir(parents=True)
 source_mu={}
 for e in m['estimators']:
  with np.load(e['predictions_npz'],allow_pickle=False) as z:
   ids=z['parent'].astype(str).tolist()
   for action in ['raw','clipped']:source_mu[e['label']+'_'+action]={p:z['prediction_'+action][ids.index(p)].copy() for p in parents}
 for parent in parents:
  truth_error=None
  try:
   with np.load(a.truth/(parent+'.npz'),allow_pickle=False) as z:truth=z['chi'];tq=z['q']
   if truth.shape!=(576,) or tq.shape!=(66,) or not np.iscomplexobj(truth) or not np.isfinite(truth).all() or not np.isfinite(tq).all():raise ValueError('invalid truth arrays')
  except Exception as exc:truth=None;tq=None;truth_error=repr(exc)
  truephase=None if tq is None else (cfg['field_mean']+A[roi]@tq[:64])/cfg['phase_transition_halfwidth'];arrays={};metadata={}
  try:
   with np.load(S/'results/contexts/cache'/(parent+'.npz'),allow_pickle=False) as z:anchor=z['chi'];aq=z['q'];pred=z['prediction']
   with np.load(S/'inputs/online'/(parent+'.npz'),allow_pickle=False) as z:y=z['observation']
   anchorphase=(cfg['field_mean']+A[roi]@aq[:64])/cfg['phase_transition_halfwidth'];arrays['original_anchor']=anchor;metadata['original_anchor']=dict(Phi=float(np.linalg.norm(pred-y)**2/cfg['noise']['sigma']**2),prior=float(.5*np.r_[aq[:64],ndtri(aq[64:])]@np.r_[aq[:64],ndtri(aq[64:])]),source_kind='original_paid_cache')
  except Exception:anchor=anchorphase=None
  arrays.update({k:source_mu[k][parent] for k in KEYS})
  for arm in ['ordinary_likelihood_shared',*KEYS]:
   action='ordinary_extra_Phi' if arm=='ordinary_likelihood_shared' else 'representative_'+arm
   try:
    if parent not in qualified:raise ValueError('UNQUALIFIED_SKIP')
    case=a.physical/parent/arm;r=read(case/'RESULTS.json')
    if r['GT_read'] is not False:raise ValueError('physicalGT contract')
    with np.load(case/'FINAL.npz',allow_pickle=False) as z:chi=z['chi'];x=z['x'];phi=float(z['Phi']);fallback=bool(z['fallback'])
    arrays[action]=chi;metadata[action]=dict(Phi=phi,prior=float(.5*x@x),fallback=fallback,actual_cost=r['actual_cost'],source_kind='new_shared_ordinary' if arm=='ordinary_likelihood_shared' else 'new_representative')
   except Exception as exc:metadata[action]=dict(status='NOT_EVALUATED',reason=repr(exc))
  rows={};predictions=[]
  for action in actions:
   value=arrays.get(action,np.full(576,np.nan+1j*np.nan));predictions.append(value)
   if np.isfinite(value).all() and truth is not None:
    v=material_losses(value[None],truth[None],roi,None if anchorphase is None else anchorphase[None],None if truephase is None else truephase[None]);rows[action]={k:(float(x[0]) if np.isfinite(x[0]) else None) for k,x in v.items()};losses[action].append(float(v['full_mse'][0]));metadata.setdefault(action,{}).update(status='EVALUATED',RMS_from_anchor=None if anchor is None else float(np.sqrt(np.mean(abs(value-anchor)**2))))
   else:
    rows[action]=dict(full_mse=None);losses[action].append(np.nan);metadata.setdefault(action,{}).update(status='NOT_EVALUATED',reason=truth_error or metadata.get(action,{}).get('reason','MISSING_ACTION'))
  write_new_json(a.out/(parent+'.json'),dict(parent=parent,losses=rows,metadata=metadata));np.savez(a.out/(parent+'.npz'),actions=np.array(actions),predictions=np.array(predictions));records.append(dict(parent=parent,losses=rows,metadata=metadata))
 comparisons=[]
 for key in KEYS:
  for ref in [key,'ordinary_extra_Phi','original_anchor']:
   ca=np.array(losses['representative_'+key]);re=np.array(losses[ref]);valid=np.isfinite(ca)&np.isfinite(re)
   row=dict(candidate='representative_'+key,reference=ref,MSE_change_candidate_minus_reference=[float(v) if np.isfinite(v) else None for v in ca-re],declared_parents=64,evaluable_pairs=int(valid.sum()),status='EVALUATED' if valid.all() else 'NOT_EVALUATED_INCOMPLETE_NO_DROPPING')
   if valid.all():row.update(paired_bootstrap(ca,re))
   comparisons.append(row)
 csvrows=[dict(parent=r['parent'],action=k,**{z:r['losses'][k].get(z) for z in ['full_mse','full_mae','roi_mse','roi_mae']},**{z:r['metadata'].get(k,{}).get(z) for z in ['status','Phi','prior','RMS_from_anchor','fallback','source_kind','reason']}) for r in records for k in actions]
 with (a.out/'PER_PARENT.csv').open('x') as stream:
  writer=csv.DictWriter(stream,fieldnames=list(csvrows[0]));writer.writeheader();writer.writerows(csvrows)
 np.savez(a.out/'PER_PARENT.npz',parent=np.array(parents),**{k:np.array(v) for k,v in losses.items()})
 write_new_json(a.out/'RESULTS.json',dict(status='OFFLINE_REPRESENTATIVE_SCORE_COMPLETE',parents=records,comparisons=comparisons,physical_cost=physical['actual_cost'],shared_ordinary_cost_charged_once=True,GT_access='separate process after closed gates',physics_calls=0,posterior_claim=False))

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('command',choices=['prepare','run','score'])
 for k in ['manifest','prepared-request','release','physical','receipt','truth','out']:p.add_argument('--'+k,type=Path)
 p.add_argument('--stage',choices=['II','IV'],default='II');p.add_argument('--job');p.add_argument('--max-wall',type=float);a=p.parse_args();{'prepare':prepare,'run':run,'score':score}[a.command](a)
if __name__=='__main__':main()
