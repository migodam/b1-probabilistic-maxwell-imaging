"""Read historical VAL256 only; never new confirmation truth or physics."""
import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
import sys,json,hashlib,platform,time
from pathlib import Path
sys.dont_write_bytecode=True
OUT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(OUT/'src'))
import numpy as np,torch,scipy
from model_actions import Action,HistoricalAction,HIST,S,R4,SEEDS,ONLINE_KEYS
from p3r2.residual_risk import read_anchor
from p3r4.risk import full_material
from p3r2.fusion import passive_box
torch.set_num_threads(1);torch.set_num_interop_threads(1)
a=read_anchor(S/'inputs/strong_context_1024_256/anchor_validation.npz')
with np.load(HIST/'prepared_validation.npz',allow_pickle=False) as z:d={k:z[k] for k in z.files}
assert np.array_equal(a['parent'],d['parent'])
sealed=json.loads((HIST/'RESULTS.json').read_text());summary={r['family']:r for r in sealed['summary']}
records=[];aggregate=[];paths=set();momcheck=[];start=time.monotonic()
for method in ('physics','strongCNN','original_K32_frozen','deterministic_matched','distilled_w128_n32768','det_frozen'):
 preds=[]
 for seed in ([SEEDS[0]] if method=='physics' else SEEDS):
  rng=torch.get_rng_state().clone();act=Action(method,seed);ctor_rng_equal=torch.equal(rng,torch.get_rng_state());rng=torch.get_rng_state().clone()
  pred=act.predict(d['y'],a,d['parent'],64);inference_rng_equal=torch.equal(rng,torch.get_rng_state());act2=Action(method,seed);repeat=act2.predict(d['y'],a,d['parent'],64)
  oldseed=seed-100;reference=None
  if method=='strongCNN':reference=full_material(d[f'strong_s{oldseed}'],act.geometry);paths.add(S/f'models/risk/strongrisk_spectral_s{oldseed}/BEST.pt')
  elif method in ('original_K32_frozen','det_frozen'):
   reference=full_material(d[f'{"K32" if method=="original_K32_frozen" else "det_frozen"}_s{oldseed}'],act.geometry);paths.add(S/f'models/risk/riskfusion_diag_gaussian_s{oldseed}/BEST.pt')
  elif method in ('deterministic_matched','distilled_w128_n32768'):
   folder=HIST/('experiments' if method=='deterministic_matched' else 'learned_context')/f'{method}_s{seed}'
   with np.load(folder/'validation.npz') as z:reference=z['prediction']
   paths.add(folder/('BEST.pt' if method=='deterministic_matched' else 'HEAD.pt'))
   if method.startswith('distilled'):paths.add(folder/'PROVIDER.pt')
  if method=='original_K32_frozen':
   chunks=[act.moments(d['y'][i:i+64],d['parent'][i:i+64]) for i in range(0,256,64)];mom=np.concatenate(chunks);momcheck.append(dict(seed=seed,max_abs_vs_stored_moments=float(np.max(abs(mom-d['K32']))),tolerance=2e-6))
   # Historical sampling path versus vectorized one, all VAL rows.
   hist=HistoricalAction(method,seed);orig=hist.predict(d['y'],a,d['parent'],64)
   momcheck[-1].update(max_abs_prediction_vs_unoptimized=float(np.max(abs(pred-orig))))
  risk=np.mean(abs(pred-d['truth'])**2,axis=1);preds.append(pred)
  rec=dict(method=method,seed=seed,parent_count=256,mean_full_complex_MSE=float(risk.mean()),exact_repeated_CPU_predictions=bool(np.array_equal(pred,repeat)),max_abs_repeat=float(np.max(abs(pred-repeat))),constructor_preserves_global_torch_rng=ctor_rng_equal,inference_preserves_global_torch_rng=inference_rng_equal)
  if reference is not None:rec.update(max_abs_vs_historical=float(np.max(abs(pred-reference))),historical_MSE=float(np.mean(abs(reference-d['truth'])**2)))
  records.append(rec);print(json.dumps(rec),flush=True)
  np.savez_compressed(OUT/'audit/old_val'/f'{method}_s{seed}.npz',parent=d['parent'],prediction=pred,risk=risk)
 family={'strongCNN':'strong_CNN','det_frozen':'deterministic_frozen'}.get(method,method)
 risks=np.stack([np.mean(abs(p-d['truth'])**2,axis=1) for p in preds]);agg=dict(method=method,seed_count=len(preds),seed_mean_parent_MSE=float(risks.mean()),ensemble_image_MSE=float(np.mean(abs(np.mean(preds,axis=0)-d['truth'])**2)))
 if family in summary:agg.update(sealed_seed_mean_parent_MSE=summary[family]['seed_mean_parent_MSE'],absolute_MSE_difference=abs(agg['seed_mean_parent_MSE']-summary[family]['seed_mean_parent_MSE']),sealed_MSE_tolerance=1e-8)
 aggregate.append(agg)
paths.add(S/'models/density/extended_diag_gaussian_s2026102501/BEST.pt')
manifest=[]
for path in sorted(paths):
 h=hashlib.blake2b()
 with path.open('rb') as f:
  for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
 st=path.stat();manifest.append(dict(path=str(path),bytes=st.st_size,mtime_ns=st.st_mtime_ns,blake2b_512=h.hexdigest()))
(OUT/'audit/model_replay/CHECKPOINT_IDENTITIES_BLAKE2B.json').write_text(json.dumps(dict(algorithm='BLAKE2b-512',checkpoints=manifest),indent=2)+'\n')
original_cross=sealed['inference_entry_replay']['crossdevice']
report=dict(status='PASS' if all(r['exact_repeated_CPU_predictions'] and r['constructor_preserves_global_torch_rng'] and r['inference_preserves_global_torch_rng'] for r in records) and all(x.get('absolute_MSE_difference',0)<=1e-8 for x in aggregate) and all(x['max_abs_vs_stored_moments']<=2e-6 for x in momcheck) else 'FAIL',scope='existing exposed VAL256 checkpoint recovery only; not independent confirmation',canonical_backend='CPU1',records=records,aggregate=aggregate,K32_equivalence=momcheck,original_CPU_MPS_gate=original_cross,original_CPU_MPS_gate_retained='FAIL; original tolerance 3e-6 unchanged',physics_RHS=0,training_updates=0,new_test_truth_read=False,environment=dict(python=sys.version,torch=torch.__version__,numpy=np.__version__,scipy=scipy.__version__,platform=platform.platform(),threads=torch.get_num_threads(),interop_threads=torch.get_num_interop_threads()),normalization_audit=dict(observation='512 source-major Re32 then Im32; y supplied as frozen historical raw measurement',feature='frozen TRAIN feature_normalization through historical feature_image',moments='frozen moment_mean/moment_scale per head',distilled='frozen provider y_mean/y_scale and target_mean/target_scale; nonnegative variance clamp',K32='normalized Gaussian then x=(sx*x_scale+x_mean)@frame_inverse; four chunk8 draws; seed parent_seed(parent,20261025121); ddof=0',per_parent_rng='local CPU torch.Generator; no parent/global RNG changes',allowed_online=list(ONLINE_KEYS)),historical_dependencies=[str(HIST/'inference.py'),str(HIST/'material_moments.py'),str(HIST/'distill_compact_w128.py'),str(R4/'src/p3r4/fusion_provider_compat.py')],wall_seconds=time.monotonic()-start)
(OUT/'CHECKPOINT_REPRODUCTION.json').write_text(json.dumps(report,indent=2)+'\n')
lines=['# Frozen checkpoint reproduction','',f"Status: **{report['status']}**. Existing VAL256 only; CPU1 canonical. No new truth, model training, or Maxwell RHS.",'','| Action | CPU seed-mean MSE | Sealed MSE | Absolute difference |','|---|---:|---:|---:|']
for row in aggregate:lines.append(f"| {row['method']} | {row['seed_mean_parent_MSE']:.12g} | {row.get('sealed_seed_mean_parent_MSE','N/A')} | {row.get('absolute_MSE_difference','N/A')} |")
lines+=['','All action/seed pairs use two fresh loads across all 256 VAL parents; exact equality and global torch RNG preservation are individually recorded. MSE reproduction gate is 1e-8; optimized K32 stored-moment equivalence gate remains 2e-6.','', 'The original CPU-versus-MPS replay remains **FAIL**, at the original 3e-6 tolerance. Exact CPU repeatability is a separate claim and does not repair that failure.','', 'Normalization, versions, allowed observation keys, K32 sampling contract, and source dependency paths are in CHECKPOINT_REPRODUCTION.json. BLAKE2b-512 checkpoint identities are in audit/model_replay/CHECKPOINT_IDENTITIES_BLAKE2B.json. Compact predictions and per-parent losses are in audit/old_val/.']
(OUT/'CHECKPOINT_REPRODUCTION.md').write_text('\n'.join(lines)+'\n')
print(report['status'],json.dumps(aggregate),flush=True)
