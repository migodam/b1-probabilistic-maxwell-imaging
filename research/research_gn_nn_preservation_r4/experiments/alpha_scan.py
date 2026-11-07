"""VAL-only posthoc final-current soft release, not a trained-soft Cascade."""
import argparse,json,sys,time
from pathlib import Path
S=Path(__file__).resolve().parents[1];sys.path.insert(0,str(S.parents[1]))
import numpy as np
from research.research_gn_nn_preservation_r4.experiments.run_evaluation import read_npz
p=argparse.ArgumentParser();p.add_argument('--source-run',required=True);p.add_argument('--alphas',nargs='+',type=float,default=[0,.1,.3,.6,1]);a=p.parse_args()
records=json.loads((S/'data/manifests/online.json').read_text())['records']
ops=read_npz(S/'data/operators.npz');basis=read_npz(S/'cache/projectors.npz')['W'];cost=[]
for r in records:
    if r['split']!='val':continue
    obs=r['obs_id'];source=S/'predictions'/a.source_run/f'{obs}.npz'
    if not source.exists():continue
    pred=read_npz(source);info=json.loads(source.with_suffix('.json').read_text());anchor=read_npz(S/'anchors'/f'{obs}.npz')
    for alpha in a.alphas:
        run=f'{a.source_run}_posthoc_a{alpha:g}';dest=S/'predictions'/run;dest.mkdir(parents=True,exist_ok=True)
        if (dest/f'{obs}.json').exists():continue
        tick=time.perf_counter();z=np.asarray(pred['delta_current'],complex)
        bright=basis@(basis.conj().T@z);delta=z-(1-alpha)*bright
        j=np.asarray(anchor['j0'])+delta;q=ops['E0']+ops['G']@j
        m=np.clip(np.real(np.sum(np.conj(q)*j,axis=1))/np.sum(abs(q)**2,axis=1),.02,20.)
        extra=time.perf_counter()-tick;prep=dict(info['preparation_cost']);prep['wall_seconds']+=extra
        prep.update(posthoc_projection_seconds=extra,posthoc_alpha=alpha,selection_only=True)
        np.savez_compressed(dest/f'{obs}.npz',m_proposal=m,j_aux=j,q_aux=q,delta_current=delta,
                            preparation_cost_json=json.dumps(prep))
        out=dict(info,arm='soft_dark_posthoc',preparation_cost=prep,
                 scope='VAL alpha selection from fixed J-full final delta; not trained-soft stage trajectory')
        (dest/f'{obs}.json').write_text(json.dumps(out,indent=2)+'\n')
        cost.append({'obs_id':obs,'run':run,'projection_seconds':extra})
(S/'audit/alpha_scan_public_preparation.json').write_text(json.dumps(cost,indent=2)+'\n')
print(json.dumps({'status':'VAL_ONLY_PREPARED','proposals':len(cost)}))
