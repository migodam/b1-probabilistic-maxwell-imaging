import json,time
from pathlib import Path
import numpy as np
from b1_joint.accounting import dump
from b1_joint.material import geometry
from research_a.budget import Budget
from research_a.target import Target,OUT,profile,TAU
ledger=Budget(OUT)
with ledger.job('E1-paid-controls'):
 target=Target(ledger);campaign=json.loads((OUT/'results/CONFIRM_CAMPAIGN.json').read_text());truth=json.loads((OUT/'results/offline_truth/CONFIRM_TRUTH.json').read_text());rng=np.random.default_rng(261015113);rows=[]
 # Fixed first8, base and maximum material-distance displayed static candidate.
 from research_a.material_metrics import distance
 for i in range(8):
  o=campaign[i];tr=truth[i];ftruth=target.forward(tr['p'],2.4,kind='FIXED_FIRST8_TRUTH_SECOND_FREQUENCY')[0]
  y2=tr['gain']*ftruth+target.sigma[32:]*rng.normal(size=32);yy=np.r_[o['observation'],y2]
  arm=next(a for a in o['methods'] if a['arm']=='finite_static');candidates=[o['base'],max(arm['archive'],key=lambda c:distance(c['p'],o['base']['p']))]
  for role,c in zip(['base','farthest_material'],candidates):
   f2=target.forward(c['p'],2.4,kind='FIXED_FIRST8_SELECTED_SECOND_FREQUENCY')[0];f=np.r_[c['prediction'],f2];g,r=profile(f,yy,target.sigma)
   rows.append(dict(parent=o['id'],role=role,material=c['p'],single_residual=c['residual'],two_frequency_residual=float(np.linalg.norm(r)),shared_two_frequency_gain=float(g),tau64=float(np.sqrt(__import__('scipy.stats',fromlist=['chi2']).chi2.ppf(.95,64))),new_attempts=1,noise_tuning=False,truth_only_used_to_generate_independent_new_measurement=True))
 dump(OUT/'results/SECOND_FREQUENCY_CONTROL.json',dict(cases=rows,new_material_frequency_attempts=24,selection='fixedfirst8,baseandfarthestmaterial,nottruthMAE',new_noise_seed=261015113))
 # Two exact raw proposals per developmentobject, one withcentral andonewithring;
 # same material proposal as savedfirsttrial, not a changed intervention.
 dev=json.loads((OUT/'results/DEV1_CAMPAIGN.json').read_text());no=[]
 for o in dev:
  arm=next(a for a in o['methods'] if a['arm']=='finite_static')
  for step in [1,4]:
   c=arm['attempts'][step];p=c['history']['proposal_p'];f=target.forward(p,kind='SAME_PROPOSAL_NO_COMPENSATION')[0];g,r=profile(f,np.asarray(o['observation']),np.asarray(o['sigma']))
   no.append(dict(parent=o['id'],step=step,proposal=p,raw_full_residual=float(np.linalg.norm(r)),corrected_full_residual=c['residual'],raw_compatible=bool(np.linalg.norm(r)<=TAU),corrected_compatible=c['valid'],nuisance_profile_before_and_after='sameboundedrealgain',compensated_epsilon_change=(np.array(c['p'])[2:]-np.array(p)[2:]).tolist(),full_nonlinear_secant=(np.asarray(c['prediction'])-f).tolist(),raw_gain=float(g),corrected_gain=c['gain']))
 dump(OUT/'results/NO_COMPENSATION_ABLATION.json',dict(records=no,new_attempts=8,probabilities=None))
print(ledger.summary())
