"""What fails when no central old root is returned? Replay only the 20 algebraic paths."""
import json
import numpy as np
from spectral_online import ROOT, find_roots

def main():
    layouts={x['label']:x for x in json.loads((ROOT/'results/spectral_new_layouts.json').read_text())}
    trials=json.loads((ROOT/'results/spectral_online.json').read_text())
    costs=[]
    for trial in trials:
        y=np.asarray(trial['old_ri'])@np.array([1,1j])
        s=find_roots(layouts[trial['label']]['geometry'],y)
        old=np.asarray(trial['design']['central_roots'])
        new=np.asarray(s['validated_roots'])
        dist=None
        if len(old) and len(new):
            d=np.linalg.norm(old[:,None,:]-new[None,:,:],axis=-1)
            dist=max(np.max(np.min(d,axis=0)),np.max(np.min(d,axis=1)))
        costs.append(dict(label=trial['label'],side=trial['truth_side'],rep=trial['rep'],original_count=len(old),
                          replay_count=len(new),maximum_nearest_distance=dist,
                          **{k:s[k] for k in ('paths_attempted','solve_status','path_classifications',
                            'path_jacobian_evaluations','tracking_cpu_seconds','total_wall_seconds')}))
    (ROOT/'results/online_root_diagnostics.json').write_text(json.dumps(costs,indent=2)+'\n')

if __name__=='__main__':
    main()
