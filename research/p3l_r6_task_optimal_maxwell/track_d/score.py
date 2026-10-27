"""Separate root-unlocked zero-shot-noise evaluator; no physics/model/training imports."""
import os
for _key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[_key]='1'
from pathlib import Path
import csv,json,sys,time,collections
import numpy as np
TRACK=Path(__file__).resolve().parent;OUT=TRACK.parent;sys.path.insert(0,str(OUT))
from track_d.common import freeze,parents,identity,write,ROOT,GROUPS,SEEDS
PHASES=('outsideROI','host','transition','pure')

def validate_global(d):
 globalpath=ROOT/'ALL_512_PREDICTIONS_SEALED.json';seal=json.loads(globalpath.read_text())
 if seal['status']!='ALL_512_PREDICTIONS_SEALED_BEFORE_D_GT' or seal['GT_read'] is not False:raise ValueError('global512 seal required')
 if identity(ROOT/'D_FREEZE.json')!=seal['freeze_identity']:raise ValueError('D freeze changed after seal')
 unlock=json.loads((ROOT/'GT_UNLOCK.json').read_text())
 if unlock.get('status')!='ROOT_UNLOCKED_AFTER_SEAL' or unlock.get('global_seal_identity')!=identity(globalpath):raise RuntimeError('root must unlock GT against exact global seal')
 for item in seal['group_seals']:
  if identity(item['path'])!=item:raise ValueError('group seal changed')
 result={}
 # Validate ALL predictions AND immutable context/report inputs of both groups
 # before the first independent-test truth is opened.
 for group in GROUPS:
  p=ROOT/'predictions'/group/'PREDICTIONS_SEALED.json';g=json.loads(p.read_text());ids=parents(group)
  if g['parents']!=ids or g['freeze_identity']!=seal['freeze_identity']:raise ValueError('group parent order/freeze identity')
  expected={(name,seed) for name in d['estimator_arms'] for seed in SEEDS};actual=[(r['method'],r['seed']) for r in g['records']]
  if len(actual)!=len(expected) or set(actual)!=expected:raise ValueError('all frozen methods/seeds required before GT')
  for rec in g['records']:
   if identity(rec['identity']['path'])!=rec['identity']:raise ValueError('prediction changed after seal')
  if identity(g['phases_identity']['path'])!=g['phases_identity']:raise ValueError('full-anchor phases changed')
  expectedctx={(parent,kind) for parent in ids for kind in ('full','low')};actualctx=[(c['parent'],c['context']) for c in g['context_artifacts']]
  if len(actualctx)!=len(expectedctx) or set(actualctx)!=expectedctx:raise ValueError('full/low context input family incomplete')
  for ctx in g['context_artifacts']:
   for key in ('context_identity','report_identity'):
    if identity(ctx[key]['path'])!=ctx[key]:raise ValueError('sealed context/report changed')
   report=json.loads(Path(ctx['report_identity']['path']).read_text())
   if report.get('context_identity')!=ctx['context_identity'] or report.get('freeze_identity')!=seal['freeze_identity'] or report.get('GT_read') is not False or report.get('path_complete') is not True:raise ValueError('context report identity/isolation contract')
  result[group]=g
 return result

def bootstrap(contrasts,seed,statistics,destination):
 rng=np.random.default_rng(seed);n=len(contrasts);draws=[];indices=[]
 for start in range(0,statistics['bootstrap'],500):
  ix=rng.integers(n,size=(min(500,statistics['bootstrap']-start),n));indices.append(ix.astype(np.uint16));draws.append(contrasts[ix].mean(1))
 draws=np.concatenate(draws);indices=np.concatenate(indices);alpha=statistics['per_comparison_alpha']
 np.savez_compressed(destination,parent_indices=indices,contrast_mean_draws=draws,bootstrap_seed=np.array(seed),original_parent_contrasts=contrasts)
 return [dict(estimate=float(contrasts[:,j].mean()),one_sided_upper=float(np.quantile(draws[:,j],1-alpha,method='linear')),
    one_sided_confidence=1-alpha,descriptive_two_sided_confidence=1-alpha,
    descriptive_two_sided_interval=[float(x) for x in np.quantile(draws[:,j],[alpha/2,1-alpha/2],method='linear')]) for j in range(contrasts.shape[1])]

def passive_clip(pred,cfg,roi):
 # Exact fixed p3r2.fusion.passive_box rule, copied to avoid model imports.
 a=np.asarray(pred).copy();host=complex(*cfg['host_chi'])
 a[...,roi]=np.clip(a[...,roi].real,host.real,host.real+cfg['inclusion_increment_real_interval'][1])+1j*np.clip(a[...,roi].imag,host.imag,host.imag+cfg['inclusion_increment_imag_interval'][1])
 a[...,~roi]=host;return a

def metrics(pred,truth,roi,phase):
 absolute=np.abs(pred-truth);sq=absolute**2
 return dict(mse_full=sq.mean(-1),mae_full=absolute.mean(-1),mse_roi=sq[...,roi].mean(-1),mae_roi=absolute[...,roi].mean(-1),
   phase_squared_error_sums=np.stack([(sq*(phase==k)).sum(-1) for k in (-1,0,1,2)],-1))

def summarize(metric,counts):
 answer={k:float(v.mean()) for k,v in metric.items() if k!='phase_squared_error_sums'};sums=metric['phase_squared_error_sums'];answer['anchor_phase']={}
 for pi,name in enumerate(PHASES):
  count=int(counts[:,pi].sum());answer['anchor_phase'][name]=dict(cell_count=count,cell_weighted_mse=float(sums[:,pi].sum()/count) if count else None,full_grid_mse_contribution=float(sums[:,pi].mean()/576))
 return answer

def main():
 d=freeze();seals=validate_global(d)
 if (ROOT/'scores').exists():raise FileExistsError('scores exist; no implicit repeat/overwrite')
 (ROOT/'scores').mkdir();write(ROOT/'scores/GT_OPEN_STARTED.json',dict(utc=time.time(),global_seal=identity(ROOT/'ALL_512_PREDICTIONS_SEALED.json'),all512_predictions_and_contexts_verified_before_truth=True))
 # Static original target parameters only; does not construct a solver/prior.
 cfg=json.loads((OUT/'inputs/Research_P3LR5T/decoder_sources/ONLINE_TARGET.json').read_text())
 for gi,(group,g) in enumerate(seals.items()):
  ids=parents(group);truth=[]
  for parent in ids:
   with np.load(OUT/'raw/offline_truth'/f'{parent}.npz',allow_pickle=False) as z:
    if str(z['parent'])!=parent:raise ValueError('truth parent')
    truth.append(z['chi'].copy())
  truth=np.array(truth);losses={};timing={};rawmetrics={};clipmetrics={};predictions={};clipped={};seedrows=[]
  with np.load(g['phases_identity']['path'],allow_pickle=False) as z:
   if z['parent'].astype(str).tolist()!=ids:raise ValueError('phase parent')
   phasecode=z['phase_code'].copy()
  roi=phasecode[0]>=0
  if roi.sum()!=256 or not np.all((phasecode>=0)==roi[None,:]):raise ValueError('fixed256 ROI phase mask required')
  counts=np.stack([(phasecode==k).sum(1) for k in (-1,0,1,2)],1)
  for rec in g['records']:
   with np.load(rec['identity']['path'],allow_pickle=False) as z:
    if z['parent'].astype(str).tolist()!=ids:raise ValueError('prediction parent')
    pred=z['prediction'].copy()
   if pred.shape!=(256,576) or not np.isfinite(pred).all():raise ValueError('prediction shape/finite')
   clip=passive_clip(pred,cfg,roi);name,seed=rec['method'],rec['seed'];raw=metrics(pred,truth,roi,phasecode);clp=metrics(clip,truth,roi,phasecode)
   losses.setdefault(name,{})[seed]=raw['mse_full'];rawmetrics.setdefault(name,{})[seed]=raw;clipmetrics.setdefault(name,{})[seed]=clp
   predictions.setdefault(name,{})[seed]=pred;clipped.setdefault(name,{})[seed]=clip
   timing.setdefault(name,{})[seed]=np.array(rec['per_parent_cold_anchor_seconds'])+rec['per_parent_warm_learning_seconds']
   for i,parent in enumerate(ids):
    row=dict(group=group,parent=parent,method=name,seed=seed)
    for label,met in [('raw',raw),('fixed_passive_clip',clp)]:
     for key in ('mse_full','mae_full','mse_roi','mae_roi'):row[label+'_'+key]=float(met[key][i])
    for pi,phase in enumerate(PHASES):row[phase+'_cells']=int(counts[i,pi])
    seedrows.append(row)
  for name in losses:
   if set(losses[name])!=set(SEEDS):raise ValueError('three original seed losses required')
  risk={name:np.stack([a[s] for s in SEEDS],1).mean(1) for name,a in losses.items()};seconds={name:np.stack([a[s] for s in SEEDS],1).mean(1) for name,a in timing.items()}
  full,selected,distill=d['primary_full'],d['primary_selected'],d['distill_reference'];actual=d['actualgradient_reference']
  contrast=np.column_stack([risk[selected]-1.01*risk[full],seconds[selected]-.8*seconds[full],risk[selected]-risk[actual],risk[selected]-risk[distill]])
  names=['risk_NI_1pct_same_frozen_weights','whole_system_wall_20pct','KKT_vs_same_low_actual_gradient_risk','zero_shot_risk_vs_full_Distill_SECONDARY'];bound=bootstrap(contrast,202610270699+gi,d['statistics'],ROOT/'scores'/f'{group}_BOOTSTRAP_DRAWS.npz')
  comparisons={name:dict(b,gate=b['one_sided_upper']<0) for name,b in zip(names,bound)}
  summaries={};ensemble={}
  for name in losses:
   summaries[name]={}
   for label,store in [('raw',rawmetrics),('fixed_passive_clip',clipmetrics)]:
    averaged={key:np.stack([store[name][s][key] for s in SEEDS],1).mean(1) for key in store[name][SEEDS[0]]};summaries[name][label]=summarize(averaged,counts)
   ensemble[name]=dict(raw=summarize(metrics(np.stack([predictions[name][s] for s in SEEDS],1).mean(1),truth,roi,phasecode),counts),fixed_passive_clip=summarize(metrics(np.stack([clipped[name][s] for s in SEEDS],1).mean(1),truth,roi,phasecode),counts),qualification='prediction ensemble separate descriptive action; never primary seed-mean loss')
  physical={};optimizer={}
  for kind in ('full','low'):
   rows=[c for c in g['context_artifacts'] if c['context']==kind];rows.sort(key=lambda c:ids.index(c['parent']));material=[];phis=[];native_status=collections.Counter();probit_status=collections.Counter();native_success=[];probit_success=[]
   for ctx in rows:
    with np.load(ctx['context_identity']['path'],allow_pickle=False) as z:material.append(z['chi0'].copy());phis.append(float(z['phi']))
    r=json.loads(Path(ctx['report_identity']['path']).read_text());nr=r['original_native_report'];pr=r['original_probit_report'];native_status[nr['status']]+=1;probit_status[pr['optimizer']['status']]+=1
    native_success.append(all(p.get('success',False) for p in nr.get('phases',[])));probit_success.append(bool(pr['optimizer'].get('success',False)))
   physical[kind]=summarize(metrics(np.array(material),truth,roi,phasecode),counts)
   optimizer[kind]=dict(native_status_counts=dict(native_status),probit_status_counts=dict(probit_status),native_all_phases_scipy_success_count=sum(native_success),probit_scipy_success_count=sum(probit_success),parents=256,
     phi_summary=dict(min=float(np.min(phis)),mean=float(np.mean(phis)),median=float(np.median(phis)),max=float(np.max(phis))),pipeline_fallback_count=0,fallback_policy='any missing/nonfinite/error stops; no pipeline fallback',qualification='finite best endpoints retained as original optimizer returns; success is not posterior/global-optimum certification')
  with (ROOT/'scores'/f'{group}_PER_PARENT_SEED_METRICS.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(seedrows[0]));w.writeheader();w.writerows(seedrows)
  np.savez_compressed(ROOT/'scores'/f'{group}_PARENT_LOSSES.npz',parent=np.array(ids),phase_counts=counts,**{name:np.stack([a[s] for s in SEEDS],1) for name,a in losses.items()})
  record=dict(status='COMPLETE_SEPARATE_NOISE_CONDITION',group=group,n=256,sigma_scale=GROUPS[group],experiment_scope='zero-shot noise-scale robustness only; no new acquisition/source experiment',
    primary_estimator='parent mean of three frozen seed raw full-grid MSE losses',mean_raw_full_grid_complex_MSE={k:float(v.mean()) for k,v in risk.items()},mean_anchor_plus_learning_seconds={k:float(v.mean()) for k,v in seconds.items()},comparisons=comparisons,
    secondary_seed_mean_metrics=summaries,secondary_ensemble_metrics=ensemble,physical_point_baselines=physical,optimizer_and_fallback=optimizer,
    phase_order=list(PHASES),phase_full_grid_counts=counts.sum(0).tolist(),phase_ROI_counts_host_transition_pure=counts[:,1:].sum(0).tolist(),phase_definition='frozen full96/64 anchor field; never truth',fixed_passive_clip='original p3r2 passive_box host-to-host+increment_max rectangle; outsideROI fixed host; no fitting',
    statistics=d['statistics'],bootstrap_seed=202610270699+gi,statistical_qualification='finite-sample paired percentile bootstrap; global8 Bonferroni over2 noise conditions; conditions never pooled',
    probability_claim=False,new_training=0,new_physics_RHS=0,overall_pass=all(comparisons[name]['gate'] for name in names[:2]),cost_GO=all(comparisons[name]['gate'] for name in names[:2]),mechanism_GO=all(comparisons[name]['gate'] for name in names[:3]),actualgradient_reference=actual,primary_full=full,primary_selected=selected,secondary_reference=distill,secondary_not_used_for_champion_selection=True,
    wall_qualification='full96/64 and low24/64 separately paid physical paths plus measured original learning chain',timing_repeat_qualification='one batch64 timing per frozen seed; parent bootstrap does not estimate repeated hardware timing variability')
  write(ROOT/'scores'/f'{group}_RESULTS.json',record)
 write(ROOT/'scores/SCORING_COMPLETE.json',dict(status='COMPLETE',conditions=list(GROUPS),global_family_size=8,no_noise_pooling=True,scope='zero-shot noise robustness only',new_physics_RHS=0,probability_claim=False))
if __name__=='__main__':main()
