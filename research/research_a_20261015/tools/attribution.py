import json,numpy as np
from b1_joint.accounting import dump
from research_a.budget import Budget
from research_a.target import OUT,Target,TAU
from research_a.search import correct,endpoint,EDIT_BOXES
from research_a.material_metrics import distance,structural_difference
ledger=Budget(OUT)
with ledger.job('protected-independent-restart-attribution'):
 target=Target(ledger);dev=json.loads((OUT/'results/DEV1_CAMPAIGN.json').read_text());rng=np.random.default_rng(261015114);rows=[]
 spec=dict(scope='post-confirmation development-only attribution; same target andsameprotectedgeometry, independentcompensationinitialization; no primarypaneltuning',cases=['solid','annular_cavity'],steps=[0,1,3,5],new_full_attempts=8,seed=261015114)
 dump(OUT/'configs/ATTRIBUTION_FROZEN.json',spec)
 for o in [dev[0],dev[2]]:
  original=next(m for m in o['methods'] if m['arm']=='finite_static')
  for step in spec['steps']:
   c=original['attempts'][step];p=np.array(c['history']['proposal_p']);p[2:]=[rng.uniform(2.2,2.8),rng.uniform(.04,.16)]
   p,r=correct(target,p,np.array(o['observation']),np.array(o['sigma']),EDIT_BOXES[step])
   new=endpoint(target,p,np.array(o['observation']),np.array(o['sigma']),'PROTECTED_INDEPENDENT_RESTART_ATTRIBUTION',{})
   rows.append(dict(parent=o['id'],step=step,protected_restart=new,finite_edit_endpoint=c['p'],material_l2_difference=distance(p,c['p']),full_residual_difference=new['residual']-c['residual']))
 dump(OUT/'results/STRATIFICATION_ATTRIBUTION.json',dict(protocol=spec,records=rows,all_fullvalidated=True,maximum_material_endpoint_difference=max(r['material_l2_difference'] for r in rows)))
 # Exact observed-prediction pair distances, not resoftmax probabilities.
 campaign=json.loads((OUT/'results/CONFIRM_CAMPAIGN.json').read_text());pairs=[]
 for o in campaign:
  method=next(m for m in o['methods'] if m['arm']=='finite_static');a=method['archive'];sd=np.array(o['sigma'])
  for i in range(len(a)):
   for j in range(i+1,len(a)):
    if structural_difference(a[i]['p'],a[j]['p']):
     dy=float(np.linalg.norm((np.array(a[i]['fitted_prediction'])-np.array(a[j]['fitted_prediction']))/sd))
     pairs.append(dict(parent=o['id'],candidate_i=i,candidate_j=j,p_i=a[i]['p'],p_j=a[j]['p'],material_distance=distance(a[i]['p'],a[j]['p']),whitened_pair_distance=dy,each_observation_residual=[a[i]['residual'],a[j]['residual']],nuisance_gain=[a[i]['gain'],a[j]['gain']],measurement_semantics='two fitted mean responses with respectiveboundedrealgain; notposteriorodds'))
 qualifying=[p for p in pairs if p['whitened_pair_distance']<=1]
 dump(OUT/'results/PAIRWISE_AMBIGUITY.json',dict(structural_pairs=pairs,pairs_at_noise_distance_le1=len(qualifying),first_fixed_order_witness=qualifying[0] if qualifying else None,total_pairs=len(pairs),all_unexplored_structures='unknown'))
print(ledger.summary())
