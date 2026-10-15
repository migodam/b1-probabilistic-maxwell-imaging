import json,sys,time,numpy as np
from b1_joint.accounting import dump
from research_a.budget import Budget
from research_a.target import OUT,ROOT,Target
from research_a.search import search
ledger=Budget(OUT)
with ledger.job('warm-cache-search-replay'):
 data=json.loads((OUT/'results/CONFIRM_CAMPAIGN.json').read_text());target=Target(ledger);cache={}
 for o in data:
  for c in [o['base']]+[c for a in o['methods'] for c in a['attempts']]:cache[np.asarray(c['p'],dtype='float64').tobytes()]=c
 hits=[0]
 def cached(p,*args,**kwargs):
  key=np.asarray(p,dtype='float64').tobytes()
  if key not in cache:raise RuntimeError('EXACT_MATERIAL_CACHE_MISS: warm audit forbids anynewsolve')
  c=cache[key];hits[0]+=1
  return np.array(c['prediction']),None,None,dict(Qabs=c['Qabs'],boundary_backward_residual=c['boundary_residual'],new_physics_attempts=0,dense_boundary_solves=0)
 target.forward=cached
 sys.path.insert(0,str(ROOT.parents[1]/'C1/c1-joint-p2/p2'))
 from b1_joint.proposal_correction import load_conditioned_proposal,measurement_condition
 model=ROOT.parents[1]/'C1/c1-joint-p2/p2/results/models/full_distill-26101321.pt'
 arms=['finite_static','multistart','population','mdn','adaptive'];rows=[]
 for i,o in enumerate(data[:4]):
  y=np.array(o['observation']);sigma=np.array(o['sigma']);p,receipt,density=load_conditioned_proposal(model,measurement_condition(y,sigma))
  for rep in range(5):
   for order in range(5):
    j=(order+rep)%5;arm=arms[j];rng=np.random.default_rng(261015101+1000*i+j+100000)
    t=time.perf_counter();draw=density.sample_parameters(p,12,rng)['latent'][0] if arm=='mdn' else None
    result=search(target,y,sigma,o['base'],arm,12,rng,revision=3,mdn_draws=draw)
    wall=time.perf_counter()-t;cold=next(m for m in o['methods'] if m['arm']==arm)
    assert len(result['archive'])==len(cold['archive'])
    rows.append(dict(parent=o['id'],arm=arm,repetition=rep,rotated_order=order,warm_wall_seconds=wall,cold_initial_search_seconds=cold['wall_seconds'],new_physics=0,endpoint_cache_requests=12,matching_material_keys='exactfloat64physicalcoordinates',scope='rerunproxyproposal/correction/selection withfrozenbase andexactvalidfullendpointcache; notnewmeasurements orOS-cacheflush; MDNconditionload separate'))
 dump(OUT/'results/COLD_WARM_SEARCH_TIMING.json',dict(records=rows,cache_hits=hits[0],new_physics=0,first_fixed_parents=4,rotated_repetitions=5,baseline_setup_separate=True,CPU_only=True))
print(ledger.summary())
