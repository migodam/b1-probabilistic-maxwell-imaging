"""Score the fixed v04 public m80 continuation after all three outputs seal."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import argparse,json,time,sys
import numpy as np
from scipy.stats import chi2
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT.parents[1]))
from research.research_max_struct_r10.experiments.score_vector_comparison import material_metrics,table,dump

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--authorize-postseal',required=True);args=ap.parse_args()
    if args.authorize_postseal!='ROOT_EXPLICIT_POSTSEAL_AUTHORIZATION':raise ValueError('Root authorization required')
    wall,cpu=time.perf_counter(),time.process_time();base=ROOT/'RESULTS/vector_continuation_v04'
    receipt=json.loads((base/'PROCESS_COST.json').read_text())
    if receipt['status']!='COMPLETE':raise ValueError('All three diagnostics must finish')
    arms=('sourceCIE','GlobalDual','AplusGlobal')
    for arm in arms:
        seal=json.loads((base/arm/'SEALED.json').read_text())
        if seal['status']!='SEALED_BEFORE_ROOT_OFFLINE_DIAGNOSTIC_SCORE':raise ValueError('Unsealed arm')
    out=ROOT/'delegated/vector_continuation_score';out.mkdir(parents=True,exist_ok=False)
    dump(out/'ALL3_PUBLIC_SEAL.json',dict(scope='EXPOSED_POSTHOC_BUDGET_DIAGNOSTIC',arms=arms,
        prior_truth_access='Primary18 truth was already postseal scored; producer never reads GT',new_prediction_scoring='All3 verified before scoring'))
    with np.load(ROOT/'sealed_truth/vector/v04/grid12.npz') as z:truth=z['material']
    with np.load(ROOT/'data_public/vector/v04_matched12/observations.npz') as z:Y=z['Y']
    sigma=json.loads((ROOT/'data_public/vector/v04_matched12/PUBLIC.json').read_text())['sigma']
    limit=float(chi2.ppf(.99,2*Y.size)/2);rows=[];volumes={'GT':truth}
    for arm in arms:
        extra=json.loads((base/arm/'RECEIPT.json').read_text())
        for phase,folder,name in [('at80',ROOT/'RESULTS/vector'/('v04_'+arm),'physical_accepted80'),('continued',base/arm,'prediction')]:
            with np.load(folder/(name+'.npz')) as z:m=z['material'];prediction=z['prediction']
            meta=json.loads((folder/(name+'.json')).read_text());residual=prediction-Y
            white=float(np.vdot(residual,residual).real/sigma**2)
            row=dict(arm=arm,phase=phase,**material_metrics(m,truth),Phi=meta['Phi'],
                relative_data_residual=float(np.linalg.norm(residual)/np.linalg.norm(Y)),
                whitened_squared=white,CN99_squared_radius=limit,CN99_data_compatible=white<=limit,
                extra_accepted=extra['physical_corrector']['accepted'] if phase=='continued' else 0,
                extra_requests=extra['physical_corrector']['requests'] if phase=='continued' else 0,
                extra_CPU_nested=extra['CPU_seconds'] if phase=='continued' else 0,
                extra_wall_nested=extra['wall_seconds'] if phase=='continued' else 0,
                projected_gradient_inf=meta.get('PG_inf'),
                extra_stop=extra['physical_corrector']['stop'] if phase=='continued' else None)
            rows.append(row);volumes[(arm,phase)]=m
    table(out/'SCORES.csv',rows)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    from matplotlib.cm import ScalarMappable
    fig,axes=plt.subplots(2,4,figsize=(12,6.5),constrained_layout=True)
    norm=Normalize(0,8);cmap='viridis'
    for i,phase in enumerate(('at80','continued')):
        for j,arm in enumerate(('GT',)+arms):
            m=volumes['GT'] if arm=='GT' else volumes[(arm,phase)]
            axes[i,j].imshow(m.reshape(12,12,12)[:,:,6].T,origin='lower',extent=(-.3,.3,-.3,.3),norm=norm,cmap=cmap)
            label=arm
            if arm!='GT':label+=' | chi error %.3f'%next(x['chi_relative_L2'] for x in rows if x['arm']==arm and x['phase']==phase)
            axes[i,j].set_title(label)
            if j==0:axes[i,j].set_ylabel(phase+' | z=+.025 lambda')
    fig.colorbar(ScalarMappable(norm=norm,cmap=cmap),ax=axes,shrink=.75,label='contrast chi')
    fig.suptitle('v04 public m80 + fixed extra80 physical steps; all three reset history')
    fig.savefig(out/'MATERIALS.png',dpi=180);fig.savefig(out/'MATERIALS.pdf');plt.close(fig)
    result=dict(status='COMPLETE',scope='EXPOSED_POSTHOC_BUDGET_DIAGNOSTIC',rows=6,
        new_forward_calls=0,new_adjoint_calls=0,CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,
        primary18_modified=False,extra_inverse_cost_source='RESULTS/vector_continuation_v04/PROCESS_COST.json',
        limitation='One previously scored object, same reset-history and extra budget; not three new materials or proof of convergence. Nested perarm extra times must not be added to the whole process again.')
    dump(out/'RECEIPT.json',result);print(json.dumps(result))

if __name__=='__main__':main()
