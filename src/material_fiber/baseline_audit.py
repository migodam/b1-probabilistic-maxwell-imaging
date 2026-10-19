"""Exposed, falsification-only baseline completion; no truth inputs or retuning."""
import json
import numpy as np
from .core import OUT,write,BlockStop
from .explore import geometry,direction,feasible_direction,constrained_fit
from p1_material.proposals import frozen_latent_step

def run(t):
    parents=json.loads((OUT/'configs/PARENTS.json').read_text());anchors=json.loads((OUT/'A/ANCHORS.json').read_text());results=[]
    selected=[r for r in parents if r['split']!='train' or r['index']<4]
    for r in selected:
        parent=r['parent'];p=t.load(anchors[parent]);t.parent=parent;y=t.observation(parent)
        if not t.qualified(p,y):results.append(dict(parent=parent,status='INVALID_BASE'));continue
        for seed in ([2026101910,2026101911] if r['split']=='heldout' and r['index']<4 else [2026101910]):
            for method in ['fullJ_sphere_only','hard_prior_sphere_only']:
                before=t.book.used['B'].copy();rr=[];rng=np.random.default_rng(seed)
                with t.book.scope(4,128):
                    for s in [.02,.04]:
                        g=geometry(t,p,s)
                        if method.startswith('hard'):v,w=direction(t,p,g,'hard_prior',rng)
                        else:v=rng.normal(size=66);v[64:]=0;v,w=feasible_direction(t,p,v)
                        if v is None:continue
                        v[64:]=0;v,w=feasible_direction(t,p,v);q,ok=frozen_latent_step(t.prior,p.q,v,s)
                        if not ok:continue
                        q,fit=constrained_fit(t,p,p,w,s,q,progress_required=False)
                        if not fit['radius_matched']:continue
                        a=t.evaluate(q,f'{parent}_{method}_{seed}_{s}_predictor');rr.append(t.record(a,parent,method,p,target_RMS=s,seed=seed,role='sphere_only_predictor',post_scoring_baseline_audit=True))
                        try:
                            t.jacobian(a);qq,fit=constrained_fit(t,a,p,w,s,a.q,progress_required=False)
                            if fit['radius_matched'] and np.linalg.norm(qq-a.q)>1e-10:
                                b=t.evaluate(qq,f'{parent}_{method}_{seed}_{s}_corrector');rr.append(t.record(b,parent,method,p,target_RMS=s,seed=seed,role='sphere_only_corrector',post_scoring_baseline_audit=True))
                        except BlockStop:break
                results.append(dict(parent=parent,method=method,seed=seed,rows=rr,states=t.book.used['B']['STATE']-before['STATE'],RHS=t.book.used['B']['RHS']-before['RHS']))
    write(OUT/'B/BASELINE_COMPLETION.json',dict(results=results,reason='Remove random angular-progress restriction from the ordinary Full-J comparator; radius alone prevents return-to-anchor',new_target=False,heldout_status='EXPOSED audit after primary scoring; no new blind claim',original_results_preserved=True,probabilities=None))
