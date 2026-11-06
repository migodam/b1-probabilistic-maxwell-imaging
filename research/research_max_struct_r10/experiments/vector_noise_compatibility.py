"""Postseal vector noise/model separation, using only already saved fields."""
from pathlib import Path
import json,time
import numpy as np
from scipy.stats import chi2
ROOT=Path(__file__).resolve().parents[1]

def main():
    wall,cpu=time.perf_counter(),time.process_time()
    completed=json.loads((ROOT/'delegated/vector_score/scored/RECEIPT.json').read_text())
    if completed['status']!='COMPLETE':raise ValueError('Primary root scoring must already be complete')
    rows=[]
    for folder in sorted((ROOT/'data_public/vector').iterdir()):
        if not (folder/'observations.npz').exists():continue
        case=folder.name.split('_')[0]
        with np.load(folder/'observations.npz') as z:Y=z['Y']
        sigma=json.loads((folder/'PUBLIC.json').read_text())['sigma']
        limit=float(chi2.ppf(.99,2*Y.size)/2)
        with np.load(ROOT/'sealed_truth/vector'/case/'grid12.npz') as z:
            coarse=z['Y_clean']
        generating=coarse
        if 'cross' in folder.name:
            with np.load(ROOT/'sealed_truth/vector'/case/'grid18.npz') as z:generating=z['Y_clean']
        norm=lambda a:float(np.vdot(a,a).real/sigma**2)
        row=dict(public_id=folder.name,CN99_squared_radius=limit,
            generating_GT_whitened_squared=norm(generating-Y),
            coarse_reference_whitened_squared=norm(coarse-Y),
            generating_GT_data_compatible=norm(generating-Y)<=limit,
            coarse_reference_data_compatible=norm(coarse-Y)<=limit,arms={})
        for arm in ('sourceCIE','GlobalDual','AplusGlobal'):
            task=case+'_'+arm if 'cross' not in folder.name else folder.name+'_'+arm
            with np.load(ROOT/'RESULTS/vector'/task/'prediction.npz') as z:
                value=norm(z['prediction']-Y);m=z['material']
                row['arms'][arm]=dict(whitened_squared=value,data_compatible=value<=limit,
                    box_legal=bool((m>=-1e-12).all() and (m<=8+1e-12).all()))
        rows.append(row)
    result=dict(status='POSTSEAL_DIAGNOSTIC_COMPLETE',rows=rows,
        scope='Per-case CN99 on declared iid scalar sigma; not a joint99 event or the source-normalized optimization objective',
        finite_noise_pair_scope='Matched-grid: generating GT and estimate share one forward class. Cross-grid: fine generating field and coarse reference are different models; do not call them two solutions of the same inverse map.',
        new_forward_calls=0,new_adjoint_calls=0,CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall)
    (ROOT/'delegated/vector_score/scored/NOISE_MODEL_COMPATIBILITY.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))

if __name__=='__main__':main()
