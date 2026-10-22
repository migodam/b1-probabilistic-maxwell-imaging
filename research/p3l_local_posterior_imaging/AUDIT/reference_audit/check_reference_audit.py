"""Supplied synthetic chains only; no measured observations or physics."""
from pathlib import Path
import json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT/'src'))
from p1_material.prior import MaterialPrior
from p3l.reference_audit import audit_reference,decode_material_batches
OUT=Path(__file__).resolve().parent
cfg=json.loads((ROOT/'research/p1_material_target_20261018/configs/ONLINE_TARGET.json').read_text())
points=json.loads((ROOT/'research/p1_material_target_20261018/inputs/geometry.json').read_text())['points']
prior=MaterialPrior(cfg,points)
rng=np.random.default_rng(2026102250)
checks={}
def check(name,ok,**detail):
    checks[name]=dict(passed=bool(ok),**detail)
    if not ok:raise AssertionError(name)
q=np.c_[rng.normal(size=(17,64)),rng.uniform(size=(17,2))]
b=decode_material_batches(q,prior)
expected=np.stack([prior.decode(row)[0] for row in q])
check('vectorized_original_decoder_agreement',np.allclose(b,expected,rtol=0,atol=1e-14),maximum_absolute_error=float(np.max(abs(b-expected))))
paths=[];n=384;warmup=8
for c in range(4):
    d=OUT/'supplied_chains'/f'chain_{c}';d.mkdir(parents=True,exist_ok=True)
    q=np.c_[rng.normal(size=(n+warmup,64)),rng.uniform(size=(n+warmup,2))]
    q[:warmup,:64]=10. # Distinct warmup to expose removal errors.
    path=d/'CHAIN.npz';paths.append(path)
    np.savez_compressed(path,q=q,Phi=rng.gamma(256,size=n+warmup),warmup=warmup,target_family=np.array('FULL'))
s=audit_reference(paths,prior,warmup,OUT/'independent_draw_control')
with np.load(OUT/'independent_draw_control/REFERENCE_ARRAYS.npz') as a:
    q=a['production_q'];full=decode_material_batches(q,prior)
    err=float(np.max(abs(full.mean(axis=(0,1))-a['material_mean'])))
    check('mean_variance_quantiles_from_full_draws',err<1e-13 and np.allclose(np.mean(abs(full-full.mean(axis=(0,1)))**2,axis=(0,1)),a['material_variance']) and np.allclose(np.quantile(full.real,[.025,.5,.975],axis=(0,1)),a['quantile_real']),mean_error=err)
    check('fixed_material_function_manifest',len(a['material_function_names'])==52 and a['selected_roi_indices'].shape==(16,) and a['projection_weights'].shape==(8,256) and int(a['projection_seed'])==2026102250)
check('warmup_removed',s['production_draws_per_chain']==n and s['warmup_removed']==warmup)
check('cross_mechanism_always_separate',s['reference_qualified'] is False and s['cross_mechanism_status']=='NOT_EVALUATED_SEPARATE_INDEPENDENT_MECHANISM_REQUIRED')
check('host_constants_explicit',s['fixed_host_pixel_count']==320 and s['outside_raw_constant_component_count']==640)
check('strict_json_has_no_nan',not any(t in (OUT/'independent_draw_control/REFERENCE_STATS.json').read_text() for t in [': NaN',': Infinity',': -Infinity']))
# Identical chains create unassessed latent/functions and ROI constants. They
# cannot pass by collecting more fixed host zeros.
constant_paths=[]
for c in range(4):
    folder=OUT/'supplied_constant_chains'/f'chain_{c}';folder.mkdir(parents=True,exist_ok=True)
    path=folder/'CHAIN.npz';constant_paths.append(path)
    q=np.tile(np.r_[np.zeros(64),.5,.5],(n+warmup,1))
    np.savez_compressed(path,q=q,Phi=np.ones(n+warmup),warmup=warmup,target_family=np.array('FULL'))
s=audit_reference(constant_paths,prior,warmup,OUT/'constant_chain_control')
check('constant_roi_cannot_hide_failure',s['status']=='UNQUALIFIED_WITHIN_MECHANISM_DIAGNOSTICS' and s['roi_constant_component_count']==512 and not s['gates']['all_material_function_bulk_and_tail_ess_at_least_minimum'] and not s['gates']['all_roi_components_assessed'])
check('nonfinite_diagnostics_are_null',s['named_diagnostics']['latent']['rhat'][0] is None and s['named_diagnostics']['latent']['status'][0]=='CONSTANT_UNASSESSED')
try:audit_reference(paths[:3],prior,warmup,OUT/'should_not_exist')
except ValueError:check('four_distinct_chains_required',True)
else:check('four_distinct_chains_required',False)
try:audit_reference(paths,prior,warmup+1,OUT/'should_not_exist')
except ValueError:check('warmup_metadata_mismatch_rejected',True)
else:check('warmup_metadata_mismatch_rejected',False)
report=dict(status='PASS_SYNTHETIC_SAVED_CHAIN_AUDIT_ONLY',checks=checks,new_rhs=0,new_physics=0,truth_access=False)
(OUT/'results.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
