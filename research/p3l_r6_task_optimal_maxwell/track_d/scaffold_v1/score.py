"""Separate root-unlocked offline D scoring; no physics, model or training imports."""
import os
for _key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[_key]='1'
from pathlib import Path
import json,sys,time
import numpy as np
TRACK=Path(__file__).resolve().parent;OUT=TRACK.parent;sys.path.insert(0,str(OUT))
from track_d.common import freeze,parents,identity,write,ROOT,GROUPS,SEEDS

def validate_global(d):
 globalpath=ROOT/'ALL_512_PREDICTIONS_SEALED.json';seal=json.loads(globalpath.read_text())
 if seal['status']!='ALL_512_PREDICTIONS_SEALED_BEFORE_D_GT' or seal['GT_read'] is not False:raise ValueError('global512 seal required')
 if identity(ROOT/'D_FREEZE.json')!=seal['freeze_identity']:raise ValueError('D freeze changed after seal')
 unlock=json.loads((ROOT/'GT_UNLOCK.json').read_text())
 if unlock.get('status')!='ROOT_UNLOCKED_AFTER_SEAL' or unlock.get('global_seal_identity')!=identity(globalpath):raise RuntimeError('root must unlock GT against exact global seal')
 for item in seal['group_seals']:
  if identity(item['path'])!=item:raise ValueError('group seal changed')
 result={}
 # Validate every prediction of both conditions before first truth is opened.
 for group in GROUPS:
  p=ROOT/'predictions'/group/'PREDICTIONS_SEALED.json';g=json.loads(p.read_text())
  if g['parents']!=parents(group) or g['freeze_identity']!=seal['freeze_identity']:raise ValueError('group parent order/freeze identity')
  expected={(name,seed) for name in ('strongCNN','distilled_w128_n32768','deterministic_matched',*d['matched_arms']) for seed in SEEDS}
  actual=[(r['method'],r['seed']) for r in g['records']]
  if len(actual)!=len(expected) or set(actual)!=expected:raise ValueError('all frozen methods/seeds required before GT')
  for rec in g['records']:
   if identity(rec['identity']['path'])!=rec['identity']:raise ValueError('prediction changed after seal')
  if identity(g['phases_identity']['path'])!=g['phases_identity']:raise ValueError('full-anchor phases changed')
  result[group]=g
 return result

def bootstrap(contrasts,seed):
 # Rows are independent parents; same draws for all comparisons. No seed resampling.
 rng=np.random.default_rng(seed);n=contrasts.shape[0];draws=[]
 for start in range(0,20000,500):
  ix=rng.integers(n,size=(min(500,20000-start),n));draws.append(contrasts[ix].mean(1))
 draws=np.concatenate(draws);alpha=.05/3
 return [dict(estimate=float(contrasts[:,j].mean()),one_sided_upper_98_333333=float(np.quantile(draws[:,j],1-alpha,method='linear')),
    descriptive_two_sided_98_333333=[float(x) for x in np.quantile(draws[:,j],[alpha/2,1-alpha/2],method='linear')]) for j in range(contrasts.shape[1])]

def main():
 d=freeze();seals=validate_global(d)
 if (ROOT/'scores').exists():raise FileExistsError('scores exist; no silent repeat or overwrite')
 (ROOT/'scores').mkdir();write(ROOT/'scores/GT_OPEN_STARTED.json',dict(utc=time.time(),global_seal=identity(ROOT/'ALL_512_PREDICTIONS_SEALED.json'),all512_verified_before_truth=True))
 for gi,(group,g) in enumerate(seals.items()):
  ids=parents(group);truth=[]
  for parent in ids:
   with np.load(OUT/'raw/offline_truth'/f'{parent}.npz',allow_pickle=False) as z:
    if str(z['parent'])!=parent:raise ValueError('truth parent')
    truth.append(z['chi'].copy())
  truth=np.array(truth);losses={};timing={};phase_risk={}
  with np.load(g['phases_identity']['path'],allow_pickle=False) as z:
   if z['parent'].astype(str).tolist()!=ids:raise ValueError('phase parent')
   phasecode=z['phase_code'].copy()
  for rec in g['records']:
   with np.load(rec['identity']['path'],allow_pickle=False) as z:
    if z['parent'].astype(str).tolist()!=ids:raise ValueError('prediction parent')
    pred=z['prediction'].copy()
   sq=abs(pred-truth)**2;name=rec['method'];losses.setdefault(name,{})[rec['seed']]=sq.mean(1)
   timing.setdefault(name,{})[rec['seed']]=np.array(rec['per_parent_cold_anchor_seconds'])+rec['per_parent_warm_learning_seconds']
   phase_risk.setdefault(name,{})[rec['seed']]=np.stack([(sq*(phasecode==k)).sum(1)/576 for k in (-1,0,1,2)],axis=1)
  for name in losses:
   if set(losses[name])!=set(SEEDS):raise ValueError('three frozen seed losses required')
  risk={name:np.stack([a[s] for s in SEEDS],1).mean(1) for name,a in losses.items()}
  seconds={name:np.stack([a[s] for s in SEEDS],1).mean(1) for name,a in timing.items()}
  # Negative supports each predeclared endpoint. Ratio is never the NI test.
  contrast=np.column_stack([risk['matched_selected']-1.01*risk['matched_full'],
      risk['matched_selected']-risk['distilled_w128_n32768'],seconds['matched_selected']-.8*seconds['matched_full']])
  bounds=bootstrap(contrast,202610270699+gi);names=['risk_NI_1pct_same_estimator','risk_vs_historical_distill','whole_system_wall_20pct']
  comparisons={name:dict(bound,gate=bound['one_sided_upper_98_333333']<0) for name,bound in zip(names,bounds)}
  # Timing gate uses a paid Native-prefix counterfactual if Native-selected.
  native=any(r['method']=='matched_selected' and r['context']=='native' for r in g['records'])
  qualification='Native cold cost estimated from paid prefix; not separately measured standalone process' if native else 'paid full and low anchor paths plus amortized measured learning-chain inference'
  record=dict(status='COMPLETE_SEPARATE_NOISE_CONDITION',group=group,n=256,
    sigma_scale=GROUPS[group],primary_estimator='parent mean of three frozen seed losses; no prediction ensemble',
    mean_raw_full_grid_complex_MSE={k:float(v.mean()) for k,v in risk.items()},
    mean_anchor_plus_learning_seconds={k:float(v.mean()) for k,v in seconds.items()},comparisons=comparisons,
    full_anchor_phase_MSE_contributions={name:np.stack([a[s] for s in SEEDS],1).mean((0,1)).tolist() for name,a in phase_risk.items()},
    phase_order=['outsideROI','host','transition','pure'],phase_definition='frozen full96/64 anchor field; never truth',
    statistics=d['statistics'],bootstrap_seed=202610270699+gi,statistical_qualification='finite-sample paired percentile bootstrap approximation; families separate per noise condition; no pooled-noise or global FWER claim',
    wall_qualification=qualification,probability_claim=False,new_training=0,new_physics_RHS=0,
    overall_pass=all(x['gate'] for x in comparisons.values()),standalone_speedup_confirmed=False if native else comparisons['whole_system_wall_20pct']['gate'])
  np.savez_compressed(ROOT/'scores'/f'{group}_PARENT_LOSSES.npz',parent=np.array(ids),**{name:np.stack([a[s] for s in SEEDS],1) for name,a in losses.items()})
  write(ROOT/'scores'/f'{group}_RESULTS.json',record)
 write(ROOT/'scores/SCORING_COMPLETE.json',dict(status='COMPLETE',conditions=list(GROUPS),no_noise_pooling=True,new_physics_RHS=0,probability_claim=False))
if __name__=='__main__':main()
