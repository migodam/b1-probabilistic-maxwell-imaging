"""Post-seal public-only finite-noise material-branch comparisons.

No truth, optimization, new field solve, rank or uniqueness inference. Two
actually audited legal materials in the same observation ball demonstrate
ambiguity of this noisy data/prior class, not exact PDE nonuniqueness.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
import argparse,json,itertools,time
from pathlib import Path
from seal_results import read,validate_seal


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--seal',type=Path,action='append',required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    cpu,wall=time.process_time(),time.perf_counter();root=a.root.resolve()
    seals=[read(s) for s in a.seal]
    for seal in seals:
        validate_seal(seal,True)
        if Path(seal['root']).resolve()!=root:raise ValueError('Wrong study root')
    import numpy as np
    from score_results import table
    out=a.out.resolve();out.mkdir(parents=True,exist_ok=False)
    bycase={};attempts=[]
    for seal in seals:
        campaign=seal['campaign']
        for task in seal['manifest_snapshot']['tasks']:
            cid=task['case_id'];folder=root/'RESULTS'/campaign/task['task_id']
            with np.load(root/task['input_file'].replace('\\','/'),allow_pickle=False) as z:
                y=z['Y'].copy();lower=float(z['lower']) if 'lower' in z else .02;upper=float(z['upper']) if 'upper' in z else 20.
            entry=bycase.setdefault(cid,dict(Y=y,materials=[],lower=lower,upper=upper))
            if not np.array_equal(y,entry['Y']):raise ValueError('Same case ID with different observations: '+cid)
            for view in ['raw','corrected','primary','secondary']:
                meta=folder/(view+'.json');arrays=folder/(view+'.npz')
                if not meta.is_file() or not arrays.is_file():continue
                h=read(meta)
                if not(h.get('eligible') and h.get('legal') and h.get('noise_feasible')):continue
                with np.load(arrays,allow_pickle=False) as z:m=z['m'].copy().ravel();pred=z['Y_prediction'].copy()
                ref=dict(campaign=campaign,task_id=task['task_id'],method=task['method'],view=view,candidate=h['candidate'])
                old=next((v for v in entry['materials'] if np.array_equal(v['m'],m)),None)
                if old is not None:
                    old['references'].append(ref);continue
                rr=float(np.linalg.norm(pred-y));ratio=float(h['residual_over_delta']);delta=rr/ratio if ratio>0 else None
                entry['materials'].append(dict(m=m,pred=pred,residual=rr,delta=delta,Phi=h['Phi'],D=h['D'],R=h['R'],references=[ref]))
    pairs=[];summary=[]
    for cid,entry in bycase.items():
        materials=entry['materials'];this=[]
        for i,j in itertools.combinations(range(len(materials)),2):
            x,y=materials[i],materials[j];distance=float(np.linalg.norm(x['m']-y['m'])/np.sqrt(x['m'].size));dy=float(np.linalg.norm(x['pred']-y['pred']));deltas=[v['delta'] for v in (x,y) if v['delta'] is not None]
            row=dict(case_id=cid,pair_i=i,pair_j=j,material_RMS_distance=distance,public_box_normalized_RMS_distance=distance/(entry['upper']-entry['lower']),data_distance=dy,data_distance_over_delta=dy/max(deltas) if deltas else None,data_triangle_bound=sum(deltas) if len(deltas)==2 else None,first_R=x['R'],second_R=y['R'],first_Phi=x['Phi'],second_Phi=y['Phi'],first_references=x['references'],second_references=y['references'],interpretation='Two audited legal materials explain the actual same observation under the declared raw discrepancy balls. Not exact-data or continuous-PDE nonuniqueness.')
            this.append(row);pairs.append(row)
        largest=max(this,key=lambda v:v['material_RMS_distance']) if this else None
        summary.append(dict(case_id=cid,distinct_feasible_materials=len(materials),qualified_pairs=len(this),largest_material_RMS_distance=largest['material_RMS_distance'] if largest else None,largest_box_normalized_distance=largest['public_box_normalized_RMS_distance'] if largest else None,largest_pair_i=largest['pair_i'] if largest else None,largest_pair_j=largest['pair_j'] if largest else None,qualification='Numerical finite-discrepancy evidence only. Fresnel/classical delta is a model diagnostic; crossgrid mismatch is not included in CN radius.'))
    table(out,'FEASIBLE_MATERIAL_PAIRS',pairs);table(out,'PER_OBJECT_BRANCH_SUMMARY',summary)
    receipt=dict(status='PUBLIC_POST_SEAL_ANALYSIS',truth_access='NONE',new_Maxwell_LU=0,new_adjoint_RHS=0,new_optimization=0,objects=len(summary),pair_rows=len(pairs),CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,source_seals=[str(v.resolve()) for v in a.seal],warning='Pair count is not an independent-material count; deterministic minPhi selection need not choose the closest-to-GT candidate. No new inference of global uniqueness or nonuniqueness.')
    (out/'RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt))
if __name__=='__main__':main()
