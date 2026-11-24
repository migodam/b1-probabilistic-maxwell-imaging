"""Combine already outward-enclosed bounds on each ordered residual eigenvalue."""
import json
from pathlib import Path
from time import perf_counter
from flint import arb
import numpy as np
HERE=Path(__file__).resolve().parent
ARRAY=HERE.parents[2]/'physical_revision'/'array'

def main():
    start=perf_counter()
    cohorts={}
    for p in (20,32):
        for prefix in ['native_p','high_grid_p','high_grid_pREPLACE_recovery','gap_refinement_p']:
            name=prefix.replace('REPLACE',str(p)) if 'REPLACE' in prefix else prefix+str(p)
            f=ARRAY/(name+'.json')
            for arm in json.loads(f.read_text())['profiles']:
                for row in arm.get('observations',[]):
                    if 'source_eigenvalues_lower' not in row: continue
                    vals=row['source_eigenvalues_lower']
                    k=row['negative_state_count']
                    tail=np.r_[np.maximum(vals[k:],0),np.zeros(k)]
                    key=(arm['cap'],row['id'])
                    group=cohorts.setdefault(key,dict(cap=arm['cap'],id=row['id'],n=len(vals),r=row['radius'],bounds=[],names=[]))
                    group['bounds'].append(tail);group['names'].append(arm['name'])
    result=[]
    for group in cohorts.values():
        a=np.stack(group.pop('bounds'));names=group.pop('names');r=arb(group['r'])
        choice=np.argmax(a,axis=0);els=a[choice,np.arange(a.shape[1])]
        robust=sum((arb(float(x)) for x in els),arb(0)).sqrt()/r
        old=max(float(sum((arb(float(x)) for x in row),arb(0)).sqrt()/r) for row in a)
        result.append(dict(**group,best_single=old,lower=float(np.nextafter(float(robust.lower()),-np.inf)) if robust>0 else 0.,winner_by_index=[names[i] for i in choice],source_eigenvalue_lowers=els.tolist()))
    (HERE/'spectral_envelope.json').write_text(json.dumps(dict(controls=result,wall_seconds=perf_counter()-start),indent=2))
    print([(x['id'],round(x['best_single'],5),round(x['lower'],5)) for x in result if x['cap']==1 and '32x32__pair2pct' in x['id']])
    print('strict new excluded',[(x['id'],x['cap']) for x in result if x['best_single']<=1<x['lower']])
if __name__=='__main__':main()
