"""Separate exposed-DEV scoring after fixed reference projection/Full checks.

Reads truth only after the exact32 online actions and evaluator configuration
are frozen. No solver/Target/trainer import, selection or statistical target
change. No locked parent is accessed.
"""
from pathlib import Path
import json
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'research/p3l_local_posterior_imaging'
PARENTS=[f'DEV_{i:03d}' for i in range(32)]


def array(path,key):
    with np.load(path,allow_pickle=False) as data:value=data[key].copy()
    if value.shape!=(576,) or not np.isfinite(value).all():raise ValueError('finite full576 material required')
    return value


def main():
    folder=OUT/'B/legal_reference_projection'
    result_path=folder/'OFFLINE_RISK.json'
    if result_path.exists():raise FileExistsError('no scoring overwrite')
    closed=json.loads((OUT/'B/jobs/validate_legal_reference_projection_all32/CLOSED.json').read_text())
    validation=json.loads((folder/'FULL_VALIDATION.json').read_text())
    if closed['status']!='COMPLETE' or [r['parent'] for r in validation['rows']]!=PARENTS:
        raise RuntimeError('all32 online actions must finish first')
    config=dict(parents=PARENTS,exposed_development=True,locked_methods_changed=False,
        methods=['optimized_anchor','supplied_Gaussian_mean','legal_reference_projection'],
        primary_loss='full_grid complex MSE',secondary_loss='full_grid complex MAE',
        bootstrap_seed=2026102356,replicates=2000,unit='paired32parents',
        failure='retain parent and missing raw aggregate; no outcome-based fallback or omission',
        probability='UNQUALIFIED_REFERENCE_ONLY_NOT_FULL_BAYES_ACTION')
    freeze=folder/'SCORING_FROZEN.json'
    if freeze.exists():raise FileExistsError('scoring config already exists')
    freeze.write_text(json.dumps(config,indent=2)+'\n')
    estimates={p:{'optimized_anchor':array(OUT/f'A/anchors/{p}.npz','chi'),
        'supplied_Gaussian_mean':array(folder/p/'INPUTS.npz','supplied_mean'),
        'legal_reference_projection':array(folder/p/'RESULT.npz','material')} for p in PARENTS}
    # All online materials/config are fixed before any GT material is read.
    rows=[]
    for parent in PARENTS:
        row=dict(parent=parent,validation=next(v for v in validation['rows'] if v['parent']==parent))
        try:
            truth=array(OUT/f'inputs/offline_only/{parent}.npz','chi')
            row['risk']={name:dict(MSE=float(np.mean(np.abs(value-truth)**2)),
                MAE=float(np.mean(np.abs(value-truth)))) for name,value in estimates[parent].items()}
        except Exception as exc:row['error']=repr(exc)
        rows.append(row)
    aggregate={};paired={}
    if all('risk' in r for r in rows):
        rng=np.random.default_rng(config['bootstrap_seed']);indices=rng.integers(0,32,(2000,32))
        for method in config['methods']:
            aggregate[method]={metric:float(np.mean([r['risk'][method][metric] for r in rows])) for metric in ['MSE','MAE']}
        for baseline in ['optimized_anchor','supplied_Gaussian_mean']:
            a=np.array([r['risk']['legal_reference_projection']['MSE'] for r in rows])
            b=np.array([r['risk'][baseline]['MSE'] for r in rows]);delta=a-b
            paired['projection_minus_'+baseline]=dict(mean_difference=float(delta.mean()),
                bootstrap95_difference=np.quantile(delta[indices].mean(1),[.025,.975]).tolist(),
                relative_MSE_reduction=float(1-a.mean()/b.mean()),
                interpretation='EXPOSED_DESCRIPTIVE_COMPARISON_NOT_HELDOUT_CONFIRMATION')
    report=dict(config=config,status='EXPOSED_DEV32_OFFLINE_EVALUATED_NOT_POSTERIOR_QUALIFIED',
        rows=rows,mean_risks=aggregate,paired=paired,all32_retained=True,
        missing_truth=[r['parent'] for r in rows if 'risk' not in r],
        noise_ball_qualified=sum(r['validation'].get('noise_ball',False) for r in rows),
        source_reference_qualified=False,new_physics_RHS=0)
    result_path.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(status=report['status'],mean_risks=aggregate,paired=paired,
        noise_ball_qualified=report['noise_ball_qualified'])))


if __name__=='__main__':main()
