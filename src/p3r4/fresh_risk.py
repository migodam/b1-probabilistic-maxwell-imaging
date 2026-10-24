"""Canonical batch64 fresh-cohort risk inference. No physics or truth imports."""
from pathlib import Path
import json,os,time,resource,sys
import numpy as np
from .evaluation import checkpoint_identity,write_new_json,read_online_context,read_anchor_fallback


def canonical_infer(frozen,contexts,native,features,out,check_frozen,device='cpu'):
 from .risk import load_r2_baseline
 from .risk_anchorbin import load_risk_compatible
 from .fusion_provider_compat import load_fusion_compatible
 from .risk_fusion import read_moments
 from p3r2.risk_models import observation_features
 from p3r2.fusion import passive_box
 import torch
 wall,cpu=time.monotonic(),time.process_time();torch.set_num_threads(1)
 if frozen.get('canonical_batch_size')!=64 or frozen.get('array_codec')!='stored':raise ValueError('canonical stored batch64 required')
 if device!='cpu':raise ValueError('frozen CPU inference; do not alter arithmetic backend')
 out=Path(out);features=Path(features);contexts=Path(contexts);native=Path(native)
 if out.exists():raise FileExistsError('new canonical online directory required')
 featuregate=json.loads((features/'FEATURES_COMPLETE.json').read_text())
 if featuregate['freeze']!=frozen:raise ValueError('feature freeze mismatch')
 sourceids=featuregate['gates'].copy()+featuregate['observation_sources']
 for provider in frozen['providers']:sourceids.extend(json.loads((features/provider/'MOMENTS_COMPLETE.json').read_text())['artifacts'])
 for v in sourceids:
  if checkpoint_identity(v['path'])!=v:raise ValueError('feature source changed')
 models={};load_errors={};geometry=None
 for m in frozen['methods']:
  if m['loader'] in ('anchor','native_saved'):continue
  try:
   model=load_fusion_compatible(m['checkpoint'],device) if m['loader']=='fusion' else load_r2_baseline(m['checkpoint'],device) if m['loader']=='r2' else load_risk_compatible(m['checkpoint'],device)
   models[m['id']]=model;g=model.p['geometry'] if m['loader']=='r2' else model.geometry
   if geometry is None:geometry=g
   elif not np.array_equal(geometry['roi'],g['roi']):raise ValueError('ROI mismatch')
  except Exception as exc:load_errors[m['id']]=repr(exc)
 if geometry is None:raise ValueError('no frozen geometry loadable')
 out.mkdir(parents=True);parents=frozen['parents'];preds={m['id']:[] for m in frozen['methods']};clips={m['id']:[] for m in frozen['methods']};flags={m['id']:[] for m in frozen['methods']};failures={m['id']:[] for m in frozen['methods']};allanchors=[];allqs=[];artifacts=[];batches=[]
 for start in range(0,len(parents),64):
  ids=parents[start:start+64]
  if len(ids)!=64:raise ValueError('cohort must consist of complete canonical64 blocks')
  online=[];anchors=[];qs=[];context_errors=[]
  for parent in ids:
   path=contexts/(parent+'.npz');sourceids.append(checkpoint_identity(path));err=None
   try:on=read_online_context(path,parent);anchor,q=on['chi0'],on['q0']
   except Exception as exc:
    err=repr(exc);on=None;anchor,q=read_anchor_fallback(path,parent)
   online.append(on);anchors.append(anchor);qs.append(q);context_errors.append(err)
  anchor=np.concatenate(anchors);q=np.concatenate(qs);allanchors.append(anchor);allqs.append(q)
  # Failed contexts still occupy their original batch positions. Fill only
  # disallowed-feature rows with zero placeholders; final action is anchor.
  valid=next((v for v in online if v is not None),None)
  if valid is None:context=None
  else:context={k:np.concatenate([v[k] if v is not None else np.zeros_like(valid[k]) for v in online]) for k in valid}
  for m in frozen['methods']:
   errs=[None]*64;fallback=np.array([e is not None for e in context_errors],bool)
   try:
    if m['id'] in load_errors:raise ValueError(load_errors[m['id']])
    if m['loader']=='anchor':prediction=anchor.copy();fallback[:]=False
    elif m['loader']=='native_saved':
     prediction=anchor.copy();fallback[:]=False
     for i,parent in enumerate(ids):
      try:
       path=native/parent/'NATIVE.npz';sourceids.append(checkpoint_identity(path))
       with np.load(path,allow_pickle=False) as z:prediction[i]=z['chi']
       if not np.isfinite(prediction[i]).all():raise ValueError('nonfinite native cached material')
      except Exception as exc:fallback[i]=True;errs[i]=repr(exc)
    else:
     if context is None:raise ValueError('all64 contexts invalid; anchors retained')
     model=models[m['id']]
     if m['loader']=='fusion':
      extra=read_moments(features/m['moment_provider'],ids,model.p['moment_checkpoint_identity'],m['moment_provider'],m['moment_K'])
      prediction=model.predict(observation_features(context['observation']),context,moments=extra,parents=ids,batch_size=64,moment_provider=m['moment_provider'],moment_K=m['moment_K'])
     else:prediction=model.predict_full(context,batch_size=64) if m['loader']=='r2' else model.predict(observation_features(context['observation']),context,batch_size=64)
    if prediction.shape!=(64,576) or not np.iscomplexobj(prediction):raise ValueError('canonical64 full576 complex required')
    bad=~np.isfinite(prediction).all(1);fallback|=bad
    for i in np.flatnonzero(bad):errs[i]='nonfinite prediction'
   except Exception as exc:
    prediction=anchor.copy();fallback[:]=True;errs=[repr(exc)]*64
   prediction[fallback]=anchor[fallback]
   for i in np.flatnonzero(fallback):failures[m['id']].append(dict(parent=ids[i],error=errs[i] or context_errors[i],fallback='anchor',canonical_batch_start=start))
   clipped=passive_box(prediction,geometry['config'],geometry['roi']);dest=out/'batches'/m['id']/f'batch_{start//64:03d}.npz';dest.parent.mkdir(parents=True,exist_ok=True)
   with dest.open('xb') as stream:np.savez(stream,parent=np.array(ids),prediction_raw=prediction,prediction_clipped=clipped,fallback_anchor=fallback)
   artifacts.append(checkpoint_identity(dest));preds[m['id']].append(prediction.copy());clips[m['id']].append(clipped.copy());flags[m['id']].append(fallback)
  batches.append(dict(start=start,count=64,status='PERSISTED_CANONICAL_BATCH',arithmetic='CPU float32 neural forward; exact arrays, no tolerance substitution'))
  print(json.dumps(dict(canonical_completed_parents=start+64)),flush=True)
 for m in frozen['methods']:
  dest=out/(m['id']+'.npz')
  with dest.open('xb') as stream:np.savez(stream,parent=np.array(parents),prediction_raw=np.concatenate(preds[m['id']]),prediction_clipped=np.concatenate(clips[m['id']]),fallback_anchor=np.concatenate(flags[m['id']]))
  artifacts.append(checkpoint_identity(dest))
 dest=out/'online_anchor.npz'
 with dest.open('xb') as stream:np.savez(stream,parent=np.array(parents),chi0=np.concatenate(allanchors),q0=np.concatenate(allqs))
 artifacts.append(checkpoint_identity(dest));check_frozen(frozen)
 for v in sourceids:
  if checkpoint_identity(v['path'])!=v:raise ValueError('source changed during canonical inference')
 write_new_json(out/'ONLINE_COMPLETE.json',dict(status='ONLINE_COMPLETE',pid=os.getpid(),completed_time_ns=time.time_ns(),freeze=frozen,truth_read=False,physics_calls=0,training_updates=0,parents=parents,
  geometry=dict(config=geometry['config'],roi=np.asarray(geometry['roi']).tolist(),A=np.asarray(geometry['A']).tolist()),anchor_context_used=True,failures=failures,
  sealed_artifact_identities=artifacts,source_identities=sourceids,canonical_batches=batches,array_codec='stored',canonical_batch_size=64,
  score_definition='losses of exactly persisted batch64 predictions; never recompute via tolerance replay; average seed losses, not images',
  wall_seconds=time.monotonic()-wall,cpu_seconds=time.process_time()-cpu,peak_rss_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)*(1 if sys.platform=='darwin' else 1024)))
