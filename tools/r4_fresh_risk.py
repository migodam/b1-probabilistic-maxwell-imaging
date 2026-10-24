#!/usr/bin/env python3
"""Prepared fresh risk pipeline: bind/features/infer/score/export, never generate.

Every operational command requires a parent-owned top-level freeze that binds
the sealed risk component. No observation, truth or Maxwell generation here.
"""
from pathlib import Path
import argparse,json,sys,os,time,subprocess
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from p3r4.evaluation import checkpoint_identity,write_new_json
from r4_risk_panel import frozen_check,identity_check
S=ROOT/'research/p3l_r4_conditional_transport'

def read(p):return json.loads(Path(p).read_text())
def component_check(c):
 if c['status']!='PREPARED_RISK_COMPONENT_FROZEN' or len(c['methods'])!=45:raise ValueError('frozen45-action component required')
 for v in c['source_identities']+c['feature_state_identities']:identity_check(v)
 for m in c['methods']:
  if 'checkpoint_identity' in m:identity_check(m['checkpoint_identity'])
  if 'compatibility_metadata_identity' in m:identity_check(m['compatibility_metadata_identity'])

def bound(a):
 c=read(a.component);component_check(c);top=read(a.top_freeze)
 if top.get('released') is not True:raise ValueError('top-level template is not a released parent freeze')
 if top.get('risk_component_identity')!=checkpoint_identity(a.component) or top.get('risk_rule')!=c['risk_rule'] or top.get('cohort')!=a.cohort:raise ValueError('parent-owned matching top-level freeze required; no implicit release')
 if top.get('methods_frozen') is not True or top.get('adapt_after_L1') is not False:raise ValueError('same risk rule without postL1 adaptation required')
 cohort=c['cohorts'][a.cohort]
 if top.get('parents')!=cohort['parents'] or top.get('generation_seed')!=cohort['generation_seed']:raise ValueError('top-level parent/seed binding mismatch')
 if a.cohort=='confirmation':
  if top.get('rule_frozen_before_L1_generation') is not True:raise ValueError('confirmation rule must be committed beforeL1 generation')
 return c,top

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('command',choices=('bind','features','infer','score','export'))
 for k in ('component','top-freeze','freeze','observations','contexts','native','features','online','truth','out'):ap.add_argument('--'+k,type=Path)
 ap.add_argument('--cohort',choices=('locked','confirmation'),required=True);a=ap.parse_args();c,top=bound(a)
 if a.command=='bind':
  # Binding is metadata-only; the parent performs its own generation scheduling.
  f={k:c[k] for k in ['methods','providers','comparisons','statistics']};f.update(status='FROZEN_BEFORE_LOCKED_TRUTH',parents=c['cohorts'][a.cohort]['parents'],cohort=a.cohort,risk_rule=c['risk_rule'],component_identity=checkpoint_identity(a.component),top_freeze_identity=checkpoint_identity(a.top_freeze),freeze_time_ns=time.time_ns(),freeze_pid=os.getpid(),truth_read=False,physics_calls=0,canonical_batch_size=64,array_codec='stored',adapt_after_L1=False)
  # Parent freeze must already precede any genuinely fresh observation.
  online=S/'inputs/online'
  existing=[p for p in f['parents'] if (online/(p+'.npz')).exists()]
  if existing:raise ValueError('bind must precede fresh observations; refuse retrospective freeze')
  write_new_json(a.out,f);return
 f=read(a.freeze);frozen_check(f)
 if f.get('component_identity')!=checkpoint_identity(a.component) or f.get('top_freeze_identity')!=checkpoint_identity(a.top_freeze) or f['parents']!=c['cohorts'][a.cohort]['parents']:raise ValueError('bound risk freeze mismatch')
 if a.command=='features':
  if a.truth:raise ValueError('truth prohibited')
  subprocess.run([sys.executable,str(ROOT/'tools/r4_risk_panel.py'),'features','--freeze',str(a.freeze),'--observations',str(a.observations),'--out',str(a.out)],check=True);return
 if a.command=='infer':
  if a.truth:raise ValueError('truth prohibited')
  from p3r4.fresh_risk import canonical_infer
  canonical_infer(f,a.contexts,a.native,a.features,a.out,frozen_check);component_check(c);return
 gate=read(a.online/'ONLINE_COMPLETE.json')
 if gate['freeze']!=f or gate['truth_read'] is not False or gate['canonical_batch_size']!=64 or gate['array_codec']!='stored' or gate['pid']==os.getpid():raise ValueError('closed independent canonical batch64 inference required')
 for v in gate['sealed_artifact_identities']+gate['source_identities']:identity_check(v)
 if a.command=='score':
  if 'offline_only' not in a.truth.resolve().parts:raise ValueError('isolated offline_only truth required')
  subprocess.run([sys.executable,str(ROOT/'tools/r4_risk_panel.py'),'score','--freeze',str(a.freeze),'--online',str(a.online),'--truth',str(a.truth),'--out',str(a.out)],check=True);return
 if a.truth:raise ValueError('truth prohibited in online exports')
 if a.out.exists():raise FileExistsError('new export directory required')
 a.out.mkdir(parents=True)
 # Arbitrary material diagnostic all parents: no anchor-qualification selection.
 selected=[('spline','fusion_spline_s2026102501'),('r2','frozen_r2'),('strong','strong_spectral_s2026102501')];allpred={};sources=[]
 for tag,mid in selected:
  path=a.online/(mid+'.npz');sources.append(checkpoint_identity(path))
  with np.load(path,allow_pickle=False) as z:
   if z['parent'].astype(str).tolist()!=f['parents']:raise ValueError('prediction parent identity')
   for act in ['raw','clipped']:allpred[tag+'_'+act]=z['prediction_'+act].copy()
 rows=[]
 for i,parent in enumerate(f['parents']):
  dest=a.out/'raw_mu'/(parent+'.npz');dest.parent.mkdir(exist_ok=True)
  with dest.open('xb') as stream:np.savez(stream,**{k:v[i] for k,v in allpred.items()})
  obs=a.observations/(parent+'.npz');identity=checkpoint_identity(dest)
  for key in allpred:rows.append(dict(parent=parent,action=key,frozen_mu=str(dest.resolve()),mu_identity=identity,observation=str(obs.resolve()),observation_identity=checkpoint_identity(obs)))
 write_new_json(a.out/'RAW_FORWARD_FROZEN.json',dict(status='FROZEN_ARBITRARY_MATERIAL_FORWARD',stage='IV',experiment='raw_material_forward_diagnostic_v1',rows=rows,planned_RHS=8*len(rows),GT_read=False,declared_parents=f['parents'],predeclared_unqualified_skips=[],scope='all declared parents; no prior feasibility selection',decoder_realizability='NOT_PROVEN',P1_legal_world='NOT_CLAIMED',arbitrary_materials_preserved_no_clipping=True,component_identity=checkpoint_identity(a.component),online_gate_identity=checkpoint_identity(a.online/'ONLINE_COMPLETE.json'),source_predictions=sources))
 if a.cohort=='locked':
  write_new_json(a.out/'REPRESENTATIVE_MANIFEST.json',dict(status='ONLINE_FROZEN_REPRESENTATIVE_MANIFEST',stage='IV',request_id='LOCKED64_PRIMARY_SPLINE_S1_AND_FROZEN_R2',parents=f['parents'][:64],estimators=[dict(id=mid,label=tag,predictions_npz=str((a.online/(mid+'.npz')).resolve()),identity=checkpoint_identity(a.online/(mid+'.npz'))) for tag,mid in selected[:2]],outputs=['prediction_raw','prediction_clipped'],seed_rule='fixed Spline2026102501 and frozenR2seed33; no seed-ensemble image',online_gate=str((a.online/'ONLINE_COMPLETE.json').resolve()),online_gate_identity=checkpoint_identity(a.online/'ONLINE_COMPLETE.json'),component_identity=checkpoint_identity(a.component),top_freeze_identity=checkpoint_identity(a.top_freeze),no_GT=True,frozen_estimator=True,per_arm=dict(max_states=16,max_RHS=256),ordinary_control='one shared same-budget Full-Phi arm per parent across all4material targets; charged once, never free',maximum_total_RHS=64*5*256))
 print(json.dumps(dict(status='ONLINE_EXPORTS_FROZEN',parents=len(f['parents']),raw_forward_RHS=8*len(rows),representative_eligible=a.cohort=='locked')))
if __name__=='__main__':main()
