#!/usr/bin/env python3
"""Zero-physics replay from already-disclosed, seedless truth and frozen predictions.

This is a retrospective reproducibility check, never a new blind experiment.
It does not import the solver, run Maxwell, or open private_truth.
"""
import argparse,json
from pathlib import Path
import numpy as np
from experiments.score_predictions import metrics

ROOT=Path(__file__).resolve().parent
METHODS=['fullj_native','fullj_reciprocal','globaldual']
FIELDS=['RMSE','MAE','foreground_RMS','background_RMS','IoU','SSIM','oracle_best_at_K_RMSE']

def replay(split):
    truthroot=ROOT/'truth_after_freeze'/split
    mapping=json.loads((truthroot/'mapping_seedless.json').read_text())['records']
    frozen=ROOT/'experiments/frozen_predictions'/split
    reference=json.loads((ROOT/'results'/f'{split}_scores.json').read_text())
    expected={(r['observation_id'],r['method']):r for r in reference['rows']}
    rows=[];differences=[]
    for case in mapping:
        ident=case['observation_id']
        with np.load(truthroot/f'{ident}.npz',allow_pickle=False) as z:t,mask=z['m'],z['mask'].astype(bool)
        for method in (['data_adjoint'] if split=='data_adjoint' else METHODS):
            r=json.loads((frozen/method/f'{ident}.json').read_text());m=r.get('m')
            if m is None:
                rows.append({'observation_id':ident,'method':method,'status':'NO_FINITE_ENDPOINT'});continue
            m=np.asarray(m,float);row={'observation_id':ident,'method':method,**metrics(m,t,mask)}
            candidates=r.get('candidates') or [{'m':m.tolist()}]
            row['oracle_best_at_K_RMSE']=min(float(np.sqrt(np.mean((np.asarray(c['m'])-t)**2))) for c in candidates)
            for key in FIELDS:
                differences.append(abs(row[key]-expected[(ident,method)][key]))
            # The online selected material must be a physical-score minimum, never GT-selected.
            scores=[c.get('objective',c.get('score')) for c in candidates]
            if scores and all(x is not None for x in scores):
                winner=int(np.argmin(scores));row['online_top1_matches_frozen_score_minimum']=bool(np.array_equal(m,np.asarray(candidates[winner]['m'])))
            rows.append(row)
    maxdiff=max(differences,default=0.)
    return {'scope':'Retrospective replay of disclosed truth; not a new blind trial',
            'split':split,'new_Maxwell_LU':0,'new_complex_RHS':0,'rows':rows,
            'maximum_metric_difference':maxdiff,'metric_replay_pass':bool(maxdiff<=1e-12),
            'truth_source':'truth_after_freeze only; generator entropy not used'}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--split',choices=['dev','dev_v2','holdout','budget_filled','data_adjoint','noise_extension'],default='holdout');ap.add_argument('--output',required=True);a=ap.parse_args()
    target=Path(a.output)
    if target.exists():raise FileExistsError('Preserve existing result '+str(target))
    result=replay(a.split);target.parent.mkdir(parents=True,exist_ok=True);target.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='rows'}))

if __name__=='__main__':main()
