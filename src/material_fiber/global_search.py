"""Independent/global candidate search and a distinct full-likelihood MH kernel."""
import json
import numpy as np
from .core import OUT,write,save,BlockStop
from .optimize import reconstruct
from .explore import block

def independent(t,parent,anchor,method,seed,states=32,version='v1'):
    if t.book.recovery_id:version=version+'_recovery'+t.book.recovery_id
    rng=np.random.default_rng(seed);before=t.book.used['C'].copy();records=[];points=[]
    with t.book.scope(states,40*states+64):
        for j in range(2 if states>=16 else 1):
            cap=states//(2 if states>=16 else 1)
            q=t.prior.draw(rng)[0] if method=='multistart' else t.prior.local_prior_reversible(anchor.q,rng,beta=.8,property_sd=.2)
            start=t.evaluate(q,f'{parent}_{method}_{version}_independent_start_{j}')
            p,r=reconstruct(t,parent,'coarse_TR',cap-1,start,run=f'{parent}_{method}_{version}_{seed}_{j}')
            r=t.record(p,parent,method,anchor,seed=seed,version=version,role='global_corrected_endpoint');records.append(r);points.append(p)
    return points,dict(parent=parent,method=method,version=version,records=records,states=t.book.used['C']['STATE']-before['STATE'],RHS=t.book.used['C']['RHS']-before['RHS'])

def tempered(t,parent,anchor,seed):
    rng=np.random.default_rng(seed);y=t.observation(parent);temps=np.array([1.,.25,.0625,0.]);betas=[.02,.05,.1,.3]
    chain=[anchor];accept=np.zeros(4,int);trials=np.zeros(4,int);swaps=0;swaptries=0;trace=[];before=t.book.used['C'].copy()
    with t.book.scope(64,1024):
        for i in range(1,4):chain.append(t.evaluate(t.prior.draw(rng)[0],f'{parent}_tempered_init_{i}'))
        for sweep in range(16):
            if sweep:
                for i in range(4):
                    q=t.prior.local_prior_reversible(chain[i].q,rng,beta=betas[i],property_sd=.01 if i==0 else .04)
                    proposal=t.evaluate(q,f'{parent}_tempered_t{i}_sweep{sweep}')
                    logratio=-temps[i]*(t.phi(proposal,y)-t.phi(chain[i],y));trials[i]+=1
                    yes=bool(np.log(rng.uniform())<min(0.,logratio))
                    if yes:chain[i]=proposal;accept[i]+=1
                    t.record(proposal,parent,'tempered_pCN',anchor,temperature=float(temps[i]),sweep=sweep,MH_accepted=yes,log_MH_ratio=float(logratio),prior_reversible_proposal=True,sampler_role='proposal, not an independent posterior draw')
            for i in range(sweep%2,3,2):
                a,b=chain[i],chain[i+1];lr=(temps[i]-temps[i+1])*(t.phi(a,y)-t.phi(b,y));swaptries+=1
                if np.log(rng.uniform())<min(0.,lr):chain[i],chain[i+1]=b,a;swaps+=1
            trace.append(dict(sweep=sweep,keys=[p.key for p in chain],Phi=[t.phi(p,y) for p in chain],logprior=[t.prior.logprior(p.q) for p in chain],q=np.array([p.q for p in chain])))
        for p in chain:t.persist_cache(p)
    save(OUT/f'C/chains/{parent}.npz',q=np.array([r['q'] for r in trace]),Phi=np.array([r['Phi'] for r in trace]),temperature=temps,accept=accept,trials=trials)
    result=dict(parent=parent,method='tempered_pCN',acceptance=(accept/np.maximum(trials,1)).tolist(),accepted=accept.tolist(),trials=trials.tolist(),swap_acceptance=swaps/max(1,swaptries),target='p0(q) exp(-temperature Phi(q))',proposal='Gaussian pCN x symmetric reflected-uniform-property proposal',samples_are_correlated=True,posterior_convergence='UNASSESSED_SHORT_CHAIN',credible_intervals=None,trace=[{k:v for k,v in r.items() if k!='q'} for r in trace],states=t.book.used['C']['STATE']-before['STATE'],RHS=t.book.used['C']['RHS']-before['RHS'])
    write(OUT/f'C/chains/{parent}.json',result)
    return chain[0],result

def stage_C(t):
    parents=json.loads((OUT/'configs/PARENTS.json').read_text());anchors=json.loads((OUT/'A/ANCHORS.json').read_text())
    panel=[r for r in parents if (r['split']=='train' and r['index']<4) or (r['split']=='heldout' and r['index']<4)]
    reports=[];best={};winner='multistart'
    for ix,r in enumerate(panel):
        parent=r['parent'];anchor=t.load(anchors[parent]);y=t.observation(parent);allpoints=[anchor]
        if r['split']=='heldout' and not (OUT/'C/GLOBAL_METHOD_FROZEN.json').exists():
            rr=[v for v in reports if v['method'] in ['multistart','prior_jump']]
            def score(m):
                a=[c for z in rr if z['method']==m for c in z['records']]
                return (-sum(x['data_feasible'] for x in a),float(np.mean([x['objective'] for x in a])))
            winner=min(['multistart','prior_jump'],key=score)
            write(OUT/'C/GLOBAL_METHOD_FROZEN.json',dict(winner=winner,selection='first four train only; no GT',heldout_GT_read=False))
        for m in ['multistart','prior_jump']:
            pp,result=independent(t,parent,anchor,m,2026101920+ix);allpoints.extend(pp);reports.append(result)
        cp,result=tempered(t,parent,anchor,2026101930+ix);allpoints.append(cp);reports.append(result)
        best[parent]=min(allpoints,key=lambda p:(not t.qualified(p,y),t.objective(p,y)))
        # Second, fixed paired test: transport from best discovered point vs one extra global restart.
        p=best[parent]
        if t.qualified(p,y):
            for seed in [2026101940,2026101941]:reports.append(block(t,parent,p,'soft_adaptive',seed,version='C_transport_v2'))
        else:reports.append(dict(parent=parent,method='soft_adaptive',status='INVALID_BASE_AFTER_GLOBAL'))
        selected=winner
        if r['split']=='train':
            a=[z for z in reports if z.get('parent')==parent and z['method'] in ['multistart','prior_jump']]
            selected=min(a,key=lambda z:min(c['objective'] for c in z['records']))['method']
        pp,result=independent(t,parent,p,selected,2026101950+ix,states=8,version='v2');reports.append(result)
        for p2 in pp:
            if (not t.qualified(p2,y),t.objective(p2,y))<(not t.qualified(best[parent],y),t.objective(best[parent],y)):best[parent]=p2
        t.persist_cache(best[parent]);print(json.dumps(dict(stage='C',parent=parent,consumed=dict(t.book.used['C']))),flush=True)
    write(OUT/'C/BEST_POINTS.json',{k:v.key for k,v in best.items()});write(OUT/'C/SUMMARY.json',dict(panel=panel,runs=reports,posterior_reference=False,GT_online=False))
