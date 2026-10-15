"""Zero-physics prior reconstruction, TRAIN replay, and fresh material draws."""
import csv,json,time,sys,subprocess
from pathlib import Path
import numpy as np
from scipy import ndimage
ROOT=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT/'src'))
from prior_constrained import prior as p
OUT=Path(__file__).resolve().parent
start=time.perf_counter(); t0=time.process_time()

def write(name,x):
 (OUT/name).write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n')

def comparison(a,b):
 sa,sb=a['smoothed'],b['smoothed'];x=sa-sa.mean(axis=(-2,-1),keepdims=True);y=sb-sb.mean(axis=(-2,-1),keepdims=True)
 disagreement=np.mean(a['masks']!=b['masks'],axis=(-2,-1))
 return np.array([np.sqrt(np.mean((x-y)**2)),np.mean(disagreement),np.mean(a['boundaries']!=b['boundaries']),np.sum(np.mean(disagreement,axis=0)>=.05)],float)

pts=p.geometry_points();c=p.config();rows=[];errors=[];train_rows=[]
# Explicit TRAIN-only metadata selection: no validation/locked geometry is read.
with np.load(p.INPUTS/'train.npz',allow_pickle=False) as z:
 ids=z['parent_id']; bg=z['background_id']; amp=z['amplitude']; chi=z['chi']
 assert set(ids)==set(range(250001,250065)) and set(z['split'])=={'train'}
 for fi in range(4):
  for local in range(16):
   pid=250001+16*fi+local;row=p.draw_train_row(251008,fi,local,parent_id=pid);train_rows.append(row)
   indices=np.where(ids==pid)[0]
   for k in indices:
    replay=p.decode(row,background_id=str(bg[k]),amplitude=float(amp[k]),points=pts)
    error=float(np.max(abs(replay-chi[k])));errors.append(error)
    rows.append(dict(parent_id=pid,family=row['family'],background_id=str(bg[k]),amplitude=float(amp[k]),max_error=error,prior_status='P0_EXACT' if error<=1e-12 else 'REPLAY_FAILURE'))
 train_material=np.stack([p.decode(r,points=pts)-p.background('weak',pts) for r in train_rows])
 train_rho=np.real(train_material*np.conj(p.V))
with (OUT/'TRAIN_REPLAY.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

# Fresh independent master_seed; exact original parameter support/draw order.
fresh=[]
for fi in range(4):
 for local in range(64):
  fresh.append(p.draw_train_row(20261015,fi,local,parent_id=f'pc20261015_{fi}_{local:03d}'))
assert len({tuple(r['parent_seed']) for r in fresh})==256
assert not {tuple(r['parent_seed']) for r in fresh}&{tuple(r['parent_seed']) for r in train_rows}
chi_fresh={b:np.stack([p.decode(r,background_id=b,points=pts) for r in fresh]) for b in ['weak','ring']}
rho=np.real((chi_fresh['weak']-p.background('weak',pts))*np.conj(p.V))
np.savez_compressed(OUT/'MATERIAL_ONLY_256.npz',rho=rho,weak=chi_fresh['weak'],ring=chi_fresh['ring'],points=pts,amplitude=np.asarray(1.2),parent_id=np.asarray([r['parent_id'] for r in fresh]))
write('MATERIAL_ONLY_PROVENANCE.json',{'status':'P0_EXACT','master_seed':20261015,'count':256,'family_counts':{f:64 for f in p.FAMILIES},'rows':fresh,'geometry_law':'Original A25 TRAIN branch; balanced four-family mixture; no topology/hole-count conditioning','paired_Maxwell_labels':0,'backgrounds':['weak','ring'],'amplitude':1.2,'generation_parent_source':'TRAIN law only; no development or locked geometry'})

# Freeze structural scales from TRAIN only, independent of measurements/MAE.
desc=[p.structure_descriptor(x) for x in train_rho]
metrics=[]
for i in range(64):
 for j in range(i+1,64):
  metrics.append((i,j,comparison(desc[i],desc[j])))
nearest=[]
for i in range(64):
 candidates=[v for a,b,v in metrics if a==i or b==i]
 # joint mean rank is avoided: nearest by centered material variation only.
 nearest.append(min(candidates,key=lambda v:v[0]))
nearest=np.asarray(nearest)
# Grid robustness transforms are not P0 candidates and never enter reservoir.
jitter=[]
for x,d in zip(train_rho,desc):
 for dx,dy in [(.5,0),(-.5,0),(0,.5),(0,-.5)]:
  moved=ndimage.shift(x.reshape(24,24),(dx,dy),order=1,mode='nearest',prefilter=False)
  jitter.append(comparison(d,p.structure_descriptor(moved.ravel())))
jitter=np.asarray(jitter)
low=np.quantile(nearest[:,:3],.1,axis=0);robust=np.quantile(jitter[:,:3],.95,axis=0)
freeze=np.maximum(low,robust)
cal={'status':'FROZEN_BEFORE_RESPONSE_SCREEN','source':'64 distinct TRAIN parents; no y, clean response, dev truth, MAE, or test metadata used',
 'thresholds_normalized_rho':list(p.THRESHOLDS),'gaussian_scales_cells':list(p.SCALES),'foreground_connectivity':4,'background_connectivity':8,'small_component_cutoff_cells':4,
 'metric_order':['shape_multiscale_rms','occupancy_multithreshold_mean','boundary_multithreshold_mean'],
 'TRAIN_nearest_neighbor_p10':low.tolist(),'half_cell_interpolation_jitter_p95':robust.tolist(),
 'recommended_minima_max_of_scales':dict(zip(['shape_multiscale_rms','occupancy_multithreshold_mean','boundary_multithreshold_mean'],freeze.tolist())),
 'persistent_structure_rule':'At least three of five normalized excursion thresholds have mean-over-scales occupancy disagreement >=0.05. Report components/holes separately; no forced topology-count difference.',
 'dc_policy':'DC is reported independently and removed from continuous shape metric; DC alone cannot qualify novelty.',
 'jitter_status':'256 existing TRAIN image diagnostic transforms, non-generator interpolation, not P0_EXACT and excluded from proposal reservoir',
 'threshold_role':'Parent must freeze exact joint acceptance semantics before E1; these are data-independent calibration scales, not an existence certificate'}
write('STRUCTURE_CALIBRATION.json',cal)
np.savez_compressed(OUT/'STRUCTURE_CALIBRATION_ARRAYS.npz',train_rho=train_rho,nearest=nearest,jitter=jitter)
write('TRAIN_ROWS_REPLAYED.json',{'rows':train_rows,'read_scope':'TRAIN-only regenerated geometry; no generator_manifest.json and no locked arrays read'})
write('SUMMARY.json',{'TRAIN_independent_parents':64,'TRAIN_condition_rows_checked':len(errors),'TRAIN_generator_replay_max_error':max(errors),'TRAIN_exact_replay_rows':sum(e<=1e-12 for e in errors),'fresh_material_only_parents':256,'fresh_family_counts':{f:64 for f in p.FAMILIES},'amplitude_index_1':1.2,'increment_unit_phase':[p.V.real,p.V.imag],'weak_background':[.1,.04],'ring_background':[2.,.1],'ring_radii':[.2,.35],'phase_acceptance_on_fresh':{b:sum(p.phase_audit(x,p.background(b,pts))['phase_valid'] for x in arr) for b,arr in chi_fresh.items()},'max_new_material_sample_budget':2048,'new_generator_samples_spent':256,'diagnostic_interpolated_transforms':256,'conservative_total_material_arrays_if_transforms_counted':512,'new_Maxwell_states':0,'new_Maxwell_RHS':0,'GPU_updates':0,'wall_seconds':time.perf_counter()-start,'process_cpu_seconds':time.process_time()-t0})
write('ZERO_PHYSICS_RECEIPT.json',{'entrypoint':'make_assets.py','new_Maxwell_states':0,'new_Maxwell_RHS':0,'forward_entrypoints_imported':False,'fresh_material_draws':256,'diagnostic_image_transforms':256,'original_training_parents_replayed':64,'original_training_conditions_replayed':512,'paired_labels_created':0,'wall_seconds':time.perf_counter()-start,'process_cpu_seconds':time.process_time()-t0})
print(json.dumps(json.loads((OUT/'SUMMARY.json').read_text()),ensure_ascii=False,indent=2))
