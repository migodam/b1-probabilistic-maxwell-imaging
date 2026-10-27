"""Read-only source/epoch and saved-C-context KKT-bound audit, zero physics."""
from pathlib import Path
import importlib.util,json,hashlib,datetime,collections
import numpy as np
ROOT=Path.cwd();OUT=ROOT/'audit/r6_evidence';freeze=json.loads((ROOT/'track_d/D_FREEZE.json').read_text());errors=[]
def ident(path):
 p=Path(path);h=hashlib.blake2b()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return dict(path=str(p.resolve()),bytes=p.stat().st_size,blake2b=h.hexdigest())
checks={}
for key in ('manifest_identities','frozen_checkpoint_identities','frozen_source_identities','development_evidence'):
 values=list(freeze[key].values()) if isinstance(freeze[key],dict) else freeze[key];checks[key]=[dict(path=x['path'],matches_initial_freeze=ident(x['path'])==x,current=ident(x['path'])) for x in values]
created=freeze['created_utc'];cfg=json.loads((ROOT/'track_c/kkt_closure/CONFIG_FROZEN.json').read_text());weights=json.loads((ROOT/'track_c/kkt_closure/SELECTED_WEIGHTS_IDENTITY.json').read_text());stamp=lambda s:datetime.datetime.fromisoformat(s).timestamp()
starts={}
for job in ('generate_DLOW','generate_DHIGH','D_anchor_R6_DLOW_0000_full'):
 p=ROOT/'raw/jobs'/job/'START.json';starts[job]=json.loads(p.read_text())
chronology=dict(KKT_config_utc=stamp(cfg['freeze_UTC']),selected_weights_record_utc=stamp(weights['created_UTC']),D_root_freeze_utc=created,D_job_start_times={k:v['utc'] for k,v in starts.items()},KKT_metadata_before_D_freeze=stamp(weights['created_UTC'])<created,all_D_job_starts_after_freeze=all(v['utc']>created for v in starts.values()),initial_truth_absence_evidence='Root freeze explicitly records D_truth_exists_before_freeze=false; current audit verifies job chronology, not past truth bytes or independent past filesystem absence')
if not chronology['KKT_metadata_before_D_freeze'] or not chronology['all_D_job_starts_after_freeze']:errors.append('epoch_order')
expectedstats=dict(bootstrap=20000,rng_seed=202610270699,familywise_alpha_global=.05,conditions=2,familywise_alpha_per_condition=.025,comparisons_per_condition=4,comparisons_global=8,per_comparison_alpha=.00625,NI_fraction=.01,wall_fraction=.8)
if freeze['statistics']!=expectedstats:errors.append('stats')
expectedarms={'full_deterministic':dict(original='deterministic_matched',context='full'),'selected_deterministic':dict(original='deterministic_matched',context='low',feature_mode='kkt_q'),'low_actual_deterministic':dict(original='deterministic_matched',context='low'),'full_kkt_deterministic':dict(original='deterministic_matched',context='full',feature_mode='kkt_q'),'full_distill':dict(original='distilled_w128_n32768',context='full'),'full_strong':dict(original='strongCNN',context='full')}
if freeze['estimator_arms']!=expectedarms:errors.append('arm_labels')
# Saved already-exposed CDEV ONLINE contexts, never D truth/context or model inference.
src=ROOT/'inputs/Research_P3LR5T/decoder_sources';spec=importlib.util.spec_from_file_location('prior_kkt_audit',src/'prior.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);targetcfg=json.loads((src/'ONLINE_TARGET.json').read_text());points=json.loads((src/'geometry.json').read_text())['points'];prior=module.MaterialPrior(targetcfg,points);sigma=targetcfg['noise']['sigma']
rows=[]
for config in ('native96_probit64','native24_probit64'):
 for i in range(8):
  parent=f'R6_CDEV_{i:04d}';p=ROOT/'track_c/contexts/CDEV'/parent/config/'CONTEXT.npz'
  with np.load(p,allow_pickle=False) as z:q=z['q0'].copy();J=z['J'].copy();res=z['residual'].copy()
  packed=np.concatenate((res.real.T,res.imag.T),axis=1).ravel()*np.sqrt(2)/sigma;g=J[:,:64].T@packed;zfield=q[:64];chi,info=prior.decode(q);df=module.saturation_prime(info['field']/prior.temperature)/prior.temperature*prior.mask;L=info['inclusion_increment']*df[:,None]*prior.A;op=float(abs(info['inclusion_increment'])*np.linalg.norm(df[:,None]*prior.A,2))/24
  project=lambda v:v/max(1.,np.linalg.norm(v))
  for r in (.5,1.,2.):
   # Algebraic rescale of this fixed C point only: not a new-noise physical test.
   gf=r*g;zf=r*zfield;e=float(np.linalg.norm(zf-gf));gap=float(np.linalg.norm(L@(project(gf)-project(zf))))/24;bound=op*min(e,2.);slack=bound-gap
   if slack<-1e-12:errors.append('raw_feature_bound')
   rows.append(dict(parent=parent,config=config,synthetic_feature_noise_ratio=r,stationarity_e_in_feature_units=e,raw_material_RMS_gap=gap,operator_norm_RMS=op,bound=bound,slack=slack))
result=dict(status='READ_ONLY_EPOCH_AND_RAW_BOUND_AUDIT' if not errors else 'AUDIT_ERROR',errors=errors,chronology=chronology,freeze_statistics=freeze['statistics'],freeze_arm_labels=freeze['estimator_arms'],identities=checks,
 source_differences=[x for x in checks['frozen_source_identities'] if not x['matches_initial_freeze']],source_difference_policy='Root-reported genuine ndarray serializer fix is being versioned; any initial identity mismatch must remain explicitly recorded, not silently erased.',
 KKT_bound_saved_CDEV_checks=rows,min_bound_slack=min(x['slack'] for x in rows),max_raw_feature_RMS_gap=max(x['raw_material_RMS_gap'] for x in rows),
 proof_qualifications=['g=+z uses residual y-F and Gaussian-field first-order stationarity, not constrained property-gradient stationarity.','For shifted sigma, e_feature=r*z-g_feature; theorem raw bound uses the same feature units on both sides.','Closed-unit-ball projection is globally1-Lipschitz and diameter2, so raw feature bound holds for any e; it does not require actual stationarity.','Raw material RMS bound does not directly bound normalized neural features or output risk; image_scale/head sensitivity and finite-precision errors would need separate propagation.','Same weights are unchanged but qclosure is a new feature algorithm; cost NI cannot establish mechanism superiority without same-point low-actual comparison.','Saved CDEV points and algebraic r-rescaling are reused development diagnostics; no new noise or D physical evidence.'],D_truth_opened=False,GT_scoring_run=False,new_Maxwell_RHS=0,new_model_inference=0)
(ROOT/'audit/KKT_EPOCH_PROVENANCE_QA.json').write_text(json.dumps(result,indent=2));print(json.dumps(dict(errors=errors,chronology=chronology,min_bound_slack=result['min_bound_slack'],source_differences=[x['path'] for x in result['source_differences']]),indent=2))
