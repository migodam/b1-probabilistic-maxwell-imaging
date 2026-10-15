"""Zero-solve exact TRAIN-cache screen; no inverse algorithm sees truth."""
import csv,json,time
import numpy as np
from . import prior
from .common import OUT,dump,descriptor,separation,plain

def run():
    start=time.perf_counter(); dest=OUT/'cached_pairs'; dest.mkdir(exist_ok=True)
    with np.load(prior.INPUTS/'train.npz',allow_pickle=False) as z:
        a={k:z[k] for k in ['chi','clean','sigma','case_id','parent_id','background_id','amplitude','noise_db']}
    rows=[]; best=[]
    for bg in ['weak','ring']:
        b=prior.background(bg)
        for amp in [.3,1.2]:
            ix=np.flatnonzero((a['background_id']==bg)&(a['amplitude']==amp)&(a['noise_db']==20))
            ds=[descriptor(a['chi'][i],b,amp) for i in ix]
            candidates=[]
            for j in range(len(ix)):
                for k in range(j+1,len(ix)):
                    i,l=ix[j],ix[k]; sigma=float(a['sigma'][i])
                    s=separation(ds[j],ds[k])
                    r=dict(background_id=bg,amplitude=amp,parent_a=int(a['parent_id'][i]),parent_b=int(a['parent_id'][l]),
                        case_a=str(a['case_id'][i]),case_b=str(a['case_id'][l]),sigma_fixed=sigma,
                        d_pred=float(np.linalg.norm(a['clean'][i]-a['clean'][l])*np.sqrt(2)/sigma),
                        source='TRAIN64 cached full clean; exact generator replay separately verified',prior_status='P0_EXACT',**s)
                    rows.append(r)
                    if s['structure_qualified']: candidates.append(r)
            best.extend(sorted(candidates,key=lambda r:r['d_pred'])[:8])
    with (OUT/'CACHED_PAIR_SCREEN.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    summary=dict(comparisons=len(rows),independent_TRAIN_parents=64,amplitude_index_1=1.2,
        sigma_policy='For each pair the lower-ID first endpoint defines one fixed original sigma, used for both endpoints; no pairwise rescaling',
        structural_qualified=sum(r['structure_qualified'] for r in rows),
        qualified_d_le_2=sum(r['structure_qualified'] and r['d_pred']<=2 for r in rows),
        minimum_by_background_amplitude=[min((r for r in rows if r['background_id']==b and r['amplitude']==a),key=lambda r:r['d_pred']) for b in ['weak','ring'] for a in [.3,1.2]],
        closest_structurally_qualified=sorted(best,key=lambda r:r['d_pred']),
        new_states=0,new_RHS=0,wall_seconds=time.perf_counter()-start)
    dump(dest/'SUMMARY.json',summary); print(json.dumps(plain({k:v for k,v in summary.items() if k!='closest_structurally_qualified'})))
if __name__=='__main__': run()
