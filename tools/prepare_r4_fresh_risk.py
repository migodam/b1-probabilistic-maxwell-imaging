#!/usr/bin/env python3
"""Seal all45 existing actions and feature provenance; no observations/fits/physics."""
from pathlib import Path
import sys,json,time,copy,importlib.util
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np,torch
from p3r4.evaluation import checkpoint_identity,write_new_json
from p3r4.fusion_provider_compat import load_fusion_compatible,AUDIT
S=ROOT/'research/p3l_r4_conditional_transport';D=S/'results/development_risk_strong_panel_v1';A=S/'AUDIT/risk/fresh_freeze_v3'
def main():
 torch.set_num_threads(1)
 if A.exists():raise FileExistsError('do not replace an existing freeze')
 dev=json.loads((D/'FROZEN.json').read_text());protocol=json.loads((S/'configs/PROTOCOL.json').read_text());methods=copy.deepcopy(dev['methods']);providers=copy.deepcopy(dev['providers']);A.mkdir(parents=True)
 featureids=[];validation=[];sources=[]
 def extract(v,prefix,arrays):
  if isinstance(v,dict):
   for k,x in v.items():extract(x,prefix+'__'+str(k),arrays)
  elif isinstance(v,np.ndarray):arrays[prefix]=v
  elif isinstance(v,torch.Tensor):arrays[prefix]=v.detach().cpu().numpy()
 for m in methods:
  if m['loader'] in ['anchor','native_saved']:continue
  if checkpoint_identity(m['checkpoint'])!=m['checkpoint_identity']:raise ValueError('DEV-frozen checkpoint changed')
  payload=torch.load(m['checkpoint'],map_location='cpu',weights_only=False);arrays={}
  for key in ['normalization','feature_normalization','geometry','moment_mean','moment_scale']:
   if key in payload:extract(payload[key],key,arrays)
  path=A/'feature_state'/(m['id']+'.npz');path.parent.mkdir(exist_ok=True)
  with path.open('xb') as stream:np.savez(stream,**arrays)
  featureids.append(checkpoint_identity(path));m['input_feature_state_identity']=checkpoint_identity(path)
  result=Path(m['checkpoint']).parent/'RESULT.json';vp=result.parent/'validation.npz'
  if result.exists():
   r=json.loads(result.read_text());validation.append(dict(id=m['id'],family=m['family'],result_identity=checkpoint_identity(result),validation_identity=checkpoint_identity(vp) if vp.exists() else None,validation_parents=r.get('validation_size'),mean_clipped_full_complex_mse=r.get('validation_clipped_full_complex_mse')))
   sources.append(checkpoint_identity(result))
   if vp.exists():sources.append(checkpoint_identity(vp))
  if m['loader']=='fusion':
   model=load_fusion_compatible(m['checkpoint']);m['moment_provider']=model.moment_provider;m['moment_K']=model.moment_K
   if hasattr(model,'metadata_recovery'):m['compatibility_metadata_identity']=checkpoint_identity(AUDIT)
 # Source identities refer to feature mathematics/parameters, never mutable VI
 # corrected-transport/core scheduler code. New observation bytes do not yet exist.
 files=['tools/r4_fresh_risk.py','tools/r4_fresh_representatives.py','tools/r4_risk_panel.py','tools/r4_fusion_moments.py','tools/evaluate_r4_risk.py','tools/prepare_r4_fresh_risk.py','tools/queue_r4_fresh_risk.py','src/p3r4/fresh_risk.py','src/p3r4/risk.py','src/p3r4/risk_bin.py','src/p3r4/risk_anchorbin.py','src/p3r4/risk_fusion.py','src/p3r4/fusion_provider_compat.py','src/p3r4/evaluation.py','src/p3r4/probability_eval.py','src/p3r4/density.py','src/p3r4/field_moments.py','src/p3r4/representative.py','tools/run_r4_representatives.py','tools/r4_raw_material_forward.py']
 sources.extend(checkpoint_identity(ROOT/p) for p in files)
 for module in ['p3r2.residual_risk','p3r2.risk_models','p3r2.fusion','p1_material.prior']:
  spec=importlib.util.find_spec(module)
  if spec is None or spec.origin is None:raise ValueError('source module missing '+module)
  sources.append(checkpoint_identity(spec.origin))
 for path in [S/'configs/TARGET.json',S/'configs/PROTOCOL.json',S/'configs/REPRESENTATIVE_PROTOCOL.json',ROOT/'research/p1_material_target_20261018/inputs/geometry.json',S/'inputs/strong_context_1024_256/anchor_train.npz',S/'inputs/strong_context_1024_256/anchor_validation.npz',AUDIT]:sources.append(checkpoint_identity(path))
 for provider,directory in [('spline','fusion_moments'),('diag_gaussian','fusion_diag_gaussian'),('field_moments','fusion_field_moments')]:
  for fold in ['train','validation']:
   gate=S/'inputs'/directory/fold/'MOMENTS_COMPLETE.json';sources.append(checkpoint_identity(gate));g=json.loads(gate.read_text())
   if g['checkpoint_identity']!=providers[provider]['checkpoint_identity'] or g.get('provider','spline')!=provider or g['K']!=providers[provider]['K']:raise ValueError('frozen provider source contract mismatch')
  sources.append(providers[provider]['checkpoint_identity'])
 cohorts={}
 for group,key in [('locked','new_locked'),('confirmation','independent_confirmation')]:
  p=protocol[key];parents=[f"{p['prefix']}_{i:05d}" for i in range(p['n'])]
  if any((S/'inputs/online'/(p+'.npz')).exists() for p in parents):raise ValueError('risk component must precede fresh observations')
  cohorts[group]=dict(parents=parents,n=p['n'],generation_seed=p['seed'],parent_rng='SeedSequence([seed,parent_index,0])',noise_rng='SeedSequence([seed,parent_index,1])',online_observation_shape=[32,8],context='native96 + probit64 original full66 online anchor; no newJ unless parent physically releases context',context_directory=str((S/'results/contexts'/group).resolve()))
 def meanval(family):
  rows=[v for v in validation if v['family']==family]
  if len(rows)!=3 or any(v['validation_parents']!=256 for v in rows):raise ValueError('three-seed R2VAL256 selection source required')
  return float(np.mean([v['mean_clipped_full_complex_mse'] for v in rows],dtype=np.float64))
 spline,diag=meanval('fusion_spline'),meanval('fusion_diag')
 if not spline<diag:raise ValueError('parent primary selection conflicts with fixed validation sources')
 comparisons=copy.deepcopy(dev['comparisons'])
 for action in ['raw','clipped']:comparisons.append(dict(id='primary_spline_vs_frozen_r2_'+action,candidate_family='fusion_spline',reference_family='frozen_r2',action=action))
 rule=dict(primary_family='fusion_spline',primary_action='clipped',training_seeds=[2026102501,2026102502,2026102503],primary_selection='smallest mean clipped R2VAL256 risk among frozen stronger risk families; selected by validation, not DEV; near-tie acknowledged',validation_spline=spline,validation_diag=diag,primary_references=['frozen_r2','strong_spectral','fusion_diag','fusion_permuted'],seed_aggregation='average per-parent losses across3 independent training seeds; not an ensemble image',representative_fixed_models=['fusion_spline_s2026102501','frozen_r2'],adapt_after_L1=False,identical_rule_both_cohorts=True,phase_thresholds='unchanged original decoder saturation/bin/phase source contracts',prior_and_Target='unchanged; no physics VI dependency')
 comp=dict(status='PREPARED_RISK_COMPONENT_FROZEN',created_time_ns=time.time_ns(),before_new_observations=True,methods=methods,providers=providers,comparisons=comparisons,statistics=dev['statistics'],risk_rule=rule,cohorts=cohorts,validation_sources=validation,source_identities=sources,feature_state_identities=featureids,
  downstream_physics_budget=dict(stage='IV',representative_per_arm_states=16,representative_per_arm_RHS=256,representative_max_total_RHS=81920,shared_ordinary_count_per_parent=1,raw_L1_RHS=24576,raw_confirmation_RHS=12288,all_release_and_receipt_stages_must_match=True,legacy_default='II for nonfresh existing callers only'),canonical_inference=dict(batch_size=64,device='cpu',array_codec='stored',definition='exact persisted batch64 predictions define actual risk scores; no tolerance replay replacement'),
  feature_provenance='all checkpoint/provider/TRAIN normalization/frame/sourcecode identities sealed now; future raw observation bytes cannot exist before generation and will be sealed after acquisition before feature use',
  historical_replay=dict(original_2e6_gate='FAIL_RETAINED',original_receipt=checkpoint_identity(S/'I/jobs/common_dev128_risk_summary_v1/CLOSED.json'),v2_1e5='technical heuristic replay only; not originalgatePASS, scientific qualification, or scientific threshold relaxation',v2_receipt=checkpoint_identity(S/'I/jobs/common_dev128_risk_summary_v2/CLOSED.json')),
  fallback='all45 methods/allparents retained; prediction failures explicit originalanchor fallback; missing originalanchor is fatal before offline truth',new_RHS=0,new_fits=0,new_observations=0,GT_read=False,supersedes='unbound fresh_freeze_v2 retained; only downstream explicit stage-aware request/release/receipt guards updated to StageIV before observations; same45models/providers/validation/riskrule, unchanged P1Target/thresholds/hardRHS caps; no data-dependent revision',parent_top_level_freeze='NOT_BOUND; parent binds both FROZEN_LOCKED andFROZEN_CONFIRMATION afterVIselection, before anyL1generation; thiscomponent grants no physicalrelease')
 write_new_json(A/'RISK_COMPONENT_FROZEN.json',comp)
 # Release templates are deliberately not active generation gate filenames.
 for group,v in cohorts.items():
  top=dict(status='TEMPLATE_NOT_RELEASED',cohort=group,risk_component_identity=checkpoint_identity(A/'RISK_COMPONENT_FROZEN.json'),risk_rule=rule,parents=v['parents'],generation_seed=v['generation_seed'],methods_frozen=True,adapt_after_L1=False,rule_frozen_before_L1_generation=True,parent_VI_component_binding='PENDING_PARENT_SELECTION',released=False)
  write_new_json(A/('TOP_LEVEL_'+group.upper()+'_TEMPLATE.json'),top)
 print(json.dumps(dict(status=comp['status'],path=str(A/'RISK_COMPONENT_FROZEN.json'),primary=rule['primary_family'],validation_spline=spline,validation_diag=diag)))
if __name__=='__main__':main()
