#!/usr/bin/env python3
"""DEV16-only saved-query additive selector; no L1 command, no GT/physics."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import sys,json,argparse,collections
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from p3r4.core import OUT,Book,write,save,from_probit
from p3r4.best_query_candidate import ordinary_best_query,material_best_query
from p1_material.prior import MaterialPrior

PROTOCOL=OUT/'configs/BEST_QUERY_DEV16_ADDENDUM_PROTOCOL.json';DEST=OUT/'results/II/best_query_DEV16_addendum_v1'
SOURCE=OUT/'results/II/representatives/representatives_DEV_REOPT_CNN';JOB=OUT/'II/jobs/representatives_DEV_REOPT_CNN'
PARENTS=tuple(f'R4DEV_{i:05d}' for i in range(16))
def protocol():return dict(status='PREREGISTERED_ADDITIVE_QUERY_SELECTOR_DEV_ONLY',parents=PARENTS,source=str(SOURCE),ordinary='minimumvalidatedPhi,UIDtie',
 material='minimummean|chi-mu|² amongoriginalpriorlegalFullPhi<=282.873792,UIDtie',include_anchor=True,include_events=['EVALUATED','GRADIENT_FAILED'],
 validation='originaldecoderchi/q,finiteFullpred/current/field,backward<=1e-10;gradientfailuredoesnotinvalidatevalidFullworld',
 callback_flags_used_for_selection=False,role='BEST_VALIDATED_QUERY_CANDIDATE',weights=None,posterior_mean=False,posterior_sample=False,
 old_outputs_unchanged=True,alloriginalquerycostownersretained=True,GT_read=False,new_RHS=0,new_L1_selector_authorized=False,retroactive_primary_claim=False)

def run():
 if DEST.exists():raise FileExistsError('preserveaddendum')
 if json.loads(PROTOCOL.read_text())!=json.loads(json.dumps(protocol())) or json.loads((JOB/'CLOSED.json').read_text())['status']!='COMPLETE':raise ValueError('frozenaddendumandclosedsource required')
 b=Book('II','best_query_DEV16_addendum_v1',max_wall=600,analysis_only=True);status='COMPLETE';error=None
 try:
  DEST.mkdir();cfg=json.loads((OUT/'configs/TARGET.json').read_text());sigma=cfg['noise']['sigma'];prior=MaterialPrior(cfg,json.loads((ROOT/'research/p1_material_target_20261018/inputs/geometry.json').read_text())['points'])
  req=json.loads((SOURCE/'FROZEN_REQUEST.json').read_text());items={r['parent']:r for r in req['qualified']};records=[];inventory=[]
  for parent in PARENTS:
   if parent not in items:raise ValueError('fixedDEVrosterincomplete')
   item=items[parent]
   with np.load(item['online_observation'],allow_pickle=False) as z:y=z['observation'].copy()
   with np.load(item['frozen_mu'],allow_pickle=False) as z:mus={key:z[key].copy() for key in item['outputs']}
   for arm in ('ordinary_likelihood_shared',*item['outputs']):
    folder=SOURCE/parent/arm;original=json.loads((folder/'RESULTS.json').read_text());events=[json.loads(s) for s in (folder/'TRIAL_EVENTS.jsonl').read_text().splitlines()]
    callback_ids={r['endpoint'] for r in events if r['event']=='ACCEPTED_ITERATE'};queries=[];seen=set()
    for e in events:
     if e['event'] not in ('EVALUATED','GRADIENT_FAILED') or e['endpoint'] in seen:continue
     seen.add(e['endpoint']);path=Path(item['anchor_cache']) if e['endpoint']==parent else JOB/'endpoints'/f"{e['endpoint']}.npz"
     r=dict(UID=e['endpoint'],event=e['event'],source=str(path),prior_legal=False,forward_valid=False,Phi=None,was_callback=e['endpoint'] in callback_ids)
     try:
      with np.load(path,allow_pickle=False) as z:a={k:z[k].copy() for k in ('q','chi','prediction','current','field')};back=float(z['backward_residual'])
      legal=prior.legal(a['q']);decoded=prior.decode(a['q'])[0] if legal else None;valid=legal and np.array_equal(decoded,a['chi']) and all(np.isfinite(v).all() for v in a.values()) and np.isfinite(back) and back<=1e-10
      phi=float(np.vdot(a['prediction']-y,a['prediction']-y).real/sigma**2);r.update(prior_legal=legal,forward_valid=valid,Phi=phi,chi=a['chi'],arrays=a)
     except Exception as exc:r['error']=repr(exc)
     inventory.append(dict(parent=parent,arm=arm,**{k:v for k,v in r.items() if k not in ('chi','arrays')}));queries.append(r)
    chosen=ordinary_best_query(queries) if arm=='ordinary_likelihood_shared' else material_best_query(queries,mus[arm])
    dest=DEST/parent/arm;dest.mkdir(parents=True);save(dest/'BEST_QUERY_CANDIDATE.npz',**chosen['arrays'],Phi=np.array(chosen['Phi']),UID=np.array(chosen['UID']))
    mu=mus[next(iter(mus)) if arm=='ordinary_likelihood_shared' else arm]
    diff=chosen['chi']-mu;value=float(np.vdot(diff,diff).real/576)
    result=dict(parent=parent,arm=arm,role='BEST_VALIDATED_QUERY_CANDIDATE',UID=chosen['UID'],source=chosen['source'],was_original_callback=chosen['was_callback'],
     original_returned_endpoint=original['endpoint'],original_Phi=original['Phi'],best_query_Phi=chosen['Phi'],original_material_objective=original['material_objective'],
     selected_material_objective=None if arm=='ordinary_likelihood_shared' else value,old_query_cost=original['actual_cost'],source_cost_receipt=str(JOB/'CLOSED.json'),
     new_RHS=0,GT_read=False,weights=None,candidate_only=True,old_primary_modified=False,new_L1_prediction_published=False)
    write(dest/'RESULTS.json',result);records.append(result)
  write(DEST/'QUERY_VALIDATION_INVENTORY.json',dict(rows=inventory,GT_read=False))
  result=dict(status='DEV16_SAVED_QUERY_ADDENDUM_COMPLETE',cases=records,total_cases=len(records),original_cases_unchanged=True,alloriginalqueryownersretained=True,
   ordinary_parents_improved=sum(r['best_query_Phi']<r['original_Phi']-1e-12 for r in records if r['arm']=='ordinary_likelihood_shared'),
   mu_cases_improved=sum(r['selected_material_objective']<r['original_material_objective']-1e-14 for r in records if r['selected_material_objective'] is not None),
   new_RHS=0,new_states=0,new_LU=0,new_training=0,GT_read=False,L1_selector_executed=False,not_retroactive_primary_claim=True)
  write(DEST/'RESULTS.json',result);print(json.dumps({k:v for k,v in result.items() if k!='cases'}),flush=True)
 except BaseException as exc:status='ERROR';error=repr(exc);raise
 finally:b.close(status,error)

def main():
 p=argparse.ArgumentParser();p.add_argument('action',choices=('prepare','run-dev'));a=p.parse_args()
 if a.action=='prepare':
  if PROTOCOL.exists():raise FileExistsError('preservepre-registration')
  write(PROTOCOL,protocol())
 else:run()
if __name__=='__main__':main()
