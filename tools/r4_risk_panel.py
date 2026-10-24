#!/usr/bin/env python3
"""Reusable prefeature freeze / observation-only features / online risk panel.

No truth argument in freeze/features/infer. Offline score delegates in a new
process only after all model/provider/context/prediction identities are checked.
LOCKED operation requires a separately authorized caller; this tool never
chooses parents, checkpoints, actions, or residuals from data.
"""
from pathlib import Path
import argparse,json,sys,os,time,resource,subprocess
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from p3r4.evaluation import checkpoint_identity,write_new_json,read_online_context,read_anchor_fallback,CI_LEVEL,BOOTSTRAP_SEED

def read(p):return json.loads(Path(p).read_text())
def identity_check(v):
 if checkpoint_identity(v['path'])!=v:raise ValueError('frozen source changed: '+v['path'])
def safe(p):
 if 'offline_only' in Path(p).resolve().parts:raise ValueError('offline_only prohibited online')
def frozen_check(f):
 if f['status']!='FROZEN_BEFORE_LOCKED_TRUTH':raise ValueError('freeze status')
 for m in f['methods']:
  if 'checkpoint_identity' in m:identity_check(m['checkpoint_identity'])
  if 'compatibility_metadata_identity' in m:identity_check(m['compatibility_metadata_identity'])
 for p in f['providers'].values():identity_check(p['checkpoint_identity'])

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('command',choices=('freeze','features','infer','score'))
 for k in ('spec','freeze','observations','contexts','native','online','truth','out'):ap.add_argument('--'+k,type=Path)
 ap.add_argument('--device',default='cpu',choices=('cpu','mps'));a=ap.parse_args();wall,cpu=time.monotonic(),time.process_time()
 if a.command=='freeze':
  import torch
  from p3r4.fusion_provider_compat import load_fusion_compatible
  torch.set_num_threads(1);f=read(a.spec);ids=[m['id'] for m in f['methods']]
  if len(set(ids))!=len(ids) or any(not i.replace('_','').isalnum() for i in ids):raise ValueError('safe unique methods')
  if len(set(f['parents']))!=len(f['parents']) or any('/' in p or '..' in p for p in f['parents']):raise ValueError('safe unique parents')
  providers={}
  for m in f['methods']:
   if m['loader'] not in ('r4','r2','fusion','anchor','native_saved'):raise ValueError('panel loader')
   if m['loader'] in ('anchor','native_saved'):continue
   m['checkpoint_identity']=checkpoint_identity(m['checkpoint']);m['checkpoint']=m['checkpoint_identity']['path']
   if m['loader']=='fusion':
    model=load_fusion_compatible(m['checkpoint'],'cpu');provider=model.moment_provider
    if hasattr(model,'metadata_recovery'):
     from p3r4.fusion_provider_compat import AUDIT
     m['compatibility_metadata_identity']=checkpoint_identity(AUDIT)
    p=dict(provider=provider,K=model.moment_K,checkpoint_identity=model.p['moment_checkpoint_identity'],sample_seed=20261025121)
    identity_check(p['checkpoint_identity'])
    if provider in providers and providers[provider]!=p:raise ValueError('provider name cannot conceal different sources')
    providers[provider]=p;m.update(moment_provider=provider,moment_K=model.moment_K)
  f.update(status='FROZEN_BEFORE_LOCKED_TRUTH',panel_contract='prefeature_provider_sealed_v1',providers=providers,
   freeze_time_ns=time.time_ns(),freeze_pid=os.getpid(),truth_read=False,physics_calls=0,
   statistics=dict(draws=4000,seed=BOOTSTRAP_SEED,primary_CI=CI_LEVEL,unit='independent parent; average seed losses before bootstrap'),
   selection='scientific selection deferred; no residual/sample/GT selection',fallback='anchor; missing original anchor is fatal, never drop parent')
  write_new_json(a.out,f);print(json.dumps(dict(status='FROZEN',methods=len(ids),providers=list(providers))));return
 f=read(a.freeze);frozen_check(f)
 if a.command=='features':
  if a.truth is not None:raise ValueError('truth prohibited')
  safe(a.observations);a.out.mkdir(parents=True,exist_ok=True)
  observation_sources=[checkpoint_identity(a.observations/(p+'.npz')) for p in f['parents']] if a.observations.is_dir() else [checkpoint_identity(a.observations)]
  for provider,p in f['providers'].items():
   dest=a.out/provider
   subprocess.run([sys.executable,str(ROOT/'tools/r4_fusion_moments.py'),'--checkpoint',p['checkpoint_identity']['path'],'--observations',str(a.observations),'--parents',str(a.freeze),'--out',str(dest),'--provider',provider,'--device',a.device,'--sample-seed',str(p['sample_seed'])],check=True)
  frozen_check(f)
  for v in observation_sources:identity_check(v)
  write_new_json(a.out/'FEATURES_COMPLETE.json',dict(status='FEATURES_COMPLETE',freeze=f,truth_read=False,physics_calls=0,completed_time_ns=time.time_ns(),observation_sources=observation_sources,gates=[checkpoint_identity(a.out/p/'MOMENTS_COMPLETE.json') for p in f['providers']]))
  return
 if a.command=='score':
  gate=read(a.online/'ONLINE_COMPLETE.json')
  if gate['freeze']!=f or gate['truth_read'] is not False or gate['status']!='ONLINE_COMPLETE' or gate['pid']==os.getpid() or gate['completed_time_ns']<=f['freeze_time_ns']:raise ValueError('independent closed online chronology required')
  for v in gate['sealed_artifact_identities']+gate['source_identities']:identity_check(v)
  subprocess.run([sys.executable,str(ROOT/'tools/evaluate_r4_risk.py'),'score','--freeze',str(a.freeze),'--online',str(a.online),'--truth',str(a.truth),'--out',str(a.out)],check=True);return
 # Inference only. No imports or access to truth/target solver.
 if a.truth is not None:raise ValueError('truth prohibited')
 safe(a.contexts);safe(a.native) if a.native else None
 from p3r4.risk import load_r2_baseline
 from p3r4.risk_anchorbin import load_risk_compatible
 from p3r4.fusion_provider_compat import load_fusion_compatible
 from p3r4.risk_fusion import read_moments
 from p3r2.risk_models import observation_features
 from p3r2.fusion import passive_box
 import torch
 torch.set_num_threads(1);a.online.mkdir(parents=True,exist_ok=True)
 featuregate=read(a.out/'FEATURES_COMPLETE.json')
 if featuregate['freeze']!=f:raise ValueError('feature freeze mismatch')
 sourceids=featuregate['gates'].copy()+featuregate['observation_sources']
 for provider in f['providers']:sourceids.extend(read(a.out/provider/'MOMENTS_COMPLETE.json')['artifacts'])
 for v in sourceids:identity_check(v)
 models={};errors={};geometry=None
 for m in f['methods']:
  if m['loader'] in ('anchor','native_saved'):continue
  try:
   model=load_fusion_compatible(m['checkpoint'],a.device) if m['loader']=='fusion' else load_r2_baseline(m['checkpoint'],a.device) if m['loader']=='r2' else load_risk_compatible(m['checkpoint'],a.device)
   models[m['id']]=model;g=model.p['geometry'] if m['loader']=='r2' else model.geometry
   if geometry is None:geometry=g
   elif not np.array_equal(g['roi'],geometry['roi']):raise ValueError('geometry mismatch')
  except Exception as e:errors[m['id']]=repr(e)
 if geometry is None:raise ValueError('no valid geometry')
 parents=f['parents'];preds={m['id']:[] for m in f['methods']};flags={m['id']:[] for m in f['methods']};failures={m['id']:[] for m in f['methods']};anchors=[];qs=[]
 for parent in parents:
  path=a.contexts/f'{parent}.npz';sourceids.append(checkpoint_identity(path));context_error=None
  try:on=read_online_context(path,parent);anchor,q=on['chi0'],on['q0']
  except Exception as e:context_error=repr(e);on=None;anchor,q=read_anchor_fallback(path,parent)
  anchors.append(anchor[0]);qs.append(q[0])
  for m in f['methods']:
   fallback=False
   try:
    if m['id'] in errors:raise ValueError(errors[m['id']])
    if m['loader']=='anchor':prediction=anchor.copy()
    elif m['loader']=='native_saved':
     path=a.native/parent/'NATIVE.npz';sourceids.append(checkpoint_identity(path))
     with np.load(path,allow_pickle=False) as z:prediction=z['chi'][None]
    else:
     if context_error:raise ValueError(context_error)
     model=models[m['id']]
     if m['loader']=='fusion':
      provider=m['moment_provider'];extra=read_moments(a.out/provider,[parent],model.p['moment_checkpoint_identity'],provider,m['moment_K'])
      prediction=model.predict(observation_features(on['observation']),on,moments=extra,parents=[parent],moment_provider=provider,moment_K=m['moment_K'])
     else:prediction=model.predict_full(on) if m['loader']=='r2' else model.predict(observation_features(on['observation']),on)
    if prediction.shape!=(1,576) or not np.iscomplexobj(prediction) or not np.isfinite(prediction).all():raise ValueError('finite full576 required')
   except Exception as e:
    prediction=anchor.copy();fallback=True;failures[m['id']].append(dict(parent=parent,error=repr(e),fallback='anchor'))
   preds[m['id']].append(prediction[0]);flags[m['id']].append(fallback)
  if (len(anchors)%16)==0:print(json.dumps(dict(inferred_parents=len(anchors))),flush=True)
 artifacts=[]
 for m in f['methods']:
  raw=np.array(preds[m['id']]);dest=a.online/(m['id']+'.npz')
  with dest.open('xb') as stream:np.savez_compressed(stream,parent=np.array(parents),prediction_raw=raw,prediction_clipped=passive_box(raw,geometry['config'],geometry['roi']),fallback_anchor=np.array(flags[m['id']],bool))
  artifacts.append(checkpoint_identity(dest))
 dest=a.online/'online_anchor.npz'
 with dest.open('xb') as stream:np.savez_compressed(stream,parent=np.array(parents),chi0=np.array(anchors),q0=np.array(qs))
 artifacts.append(checkpoint_identity(dest));frozen_check(f)
 for v in sourceids:identity_check(v)
 write_new_json(a.online/'ONLINE_COMPLETE.json',dict(status='ONLINE_COMPLETE',pid=os.getpid(),completed_time_ns=time.time_ns(),freeze=f,truth_read=False,physics_calls=0,training_updates=0,parents=parents,
  geometry=dict(config=geometry['config'],roi=np.asarray(geometry['roi']).tolist(),A=np.asarray(geometry['A']).tolist()),anchor_context_used=True,failures=failures,
  source_identities=sourceids,sealed_artifact_identities=artifacts,wall_seconds=time.monotonic()-wall,cpu_seconds=time.process_time()-cpu,
  peak_rss_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)*(1 if sys.platform=='darwin' else 1024),
  contract='original online observation-only provider generation; anchor/J/residual model allowlist; all parents/fallbacks retained; native_saved already-paid physics output only'))
 print(json.dumps(dict(status='ONLINE_COMPLETE',parents=len(parents),failures={k:len(v) for k,v in failures.items()})))
if __name__=='__main__':main()
