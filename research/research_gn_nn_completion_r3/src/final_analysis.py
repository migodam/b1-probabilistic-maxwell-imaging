"""Post-seal report diagnostics; no training, candidate selection, or new physics."""
from pathlib import Path
from collections import defaultdict
import json, sys
import numpy as np

S=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(S.parents[1]))
from research.research_gn_nn_completion_r3.src.score import material_metrics

def avg(xs):
    a=[x for x in xs if x is not None]
    return float(np.mean(a)) if a else None

def summarize(rows):
    top=[r for r in rows if r.get('top1')]
    return dict(n=len(rows),parents=len({r['parent_id'] for r in rows}),
                selected=len(top),joint=sum(r['joint_success'] for r in rows),
                strict_budget_joint=sum(r['joint_success'] and r['accounted_within_budget'] for r in rows),
                within_budget=sum(r['accounted_within_budget'] for r in rows),
                nrmse=avg(r['top1']['NRMSE'] for r in top),
                iou=avg(r['top1']['IoU_Otsu'] for r in top),
                raw_nrmse=avg(next((c['NRMSE'] for c in r['candidates'] if c['label']=='raw_1'),None) for r in rows),
                wall_mean=avg(r['wall_seconds'] for r in rows),
                LU_mean=avg(r['costcounts'].get('LU',0) for r in rows),
                RHS_mean=avg(r['costcounts'].get('actual_complex_RHS',0) for r in rows),
                fallback=sum(r['top1']['kind']=='anchor' for r in top))

def main():
    scores=json.loads((S/'evaluation/CONFIRMATION_SCORES.json').read_text())
    rows=scores['rows'];runs=list(dict.fromkeys(r['run'] for r in rows))
    out=dict(scope='POSTHOC_REPORT_DIAGNOSTIC_AFTER_TEST_SEAL; no new physics or online choice',
             by_cohort=[],by_family=[],oracle=[],training=[],GT_segmentation_calibration=[])
    for run in runs:
        for level in ['low','stress']:
            for noise in [0.,.001]:
                rr=[r for r in rows if r['run']==run and r['level']==level and r['noise']==noise]
                out['by_cohort'].append(dict(run=run,level=level,noise=noise,**summarize(rr)))
                if level=='low':
                    for family in ['digits','primitive_composites','islands','smooth']:
                        ff=[r for r in rr if r['family']==family]
                        out['by_family'].append(dict(run=run,noise=noise,family=family,**summarize(ff)))
    for level in ['low','stress']:
        for noise in [0.,.001]:
            rr=[r for r in scores['oracle_rows'] if r['level']==level and r['noise']==noise]
            out['oracle'].append(dict(level=level,noise=noise,n=len(rr),
                nrmse=avg(r['NRMSE'] for r in rr),physical=sum(r['physical_data_feasible'] for r in rr),
                joint=sum(r['joint_success'] for r in rr),
                median_constitutive_defect=float(np.median([r['constitutive_relative_defect'] for r in rr])),
                median_physical_residual=float(np.median([r['physical_data_relative_residual'] for r in rr]))))
    for p in sorted((S/'models/confirm1024').glob('*/summary.json')):
        d=json.loads(p.read_text());h=json.loads((p.parent/'history.json').read_text());b=h[d['best_epoch']-1]
        out['training'].append(dict(model=p.parent.name,epochs_run=d['epochs_run'],best_epoch=d['best_epoch'],
            best_loss=d['best_loss'],val_material_nrmse=b['val']['material_amplitude_nrmse'],
            training_seconds=d['training_seconds'],parameter_cost=d['config']['arithmetic_estimate']))
    # Calibration only: Otsu of a continuous GT can disagree with the frozen alpha>=.5 mask.
    manifest=json.loads((S/'data/manifests/evaluation.json').read_text())['records']
    for r in manifest:
        if r['split']!='test' or r['noise']!=0 or r['level']!='low':continue
        with np.load(S/'data'/r['truth_path']) as t:
            d=material_metrics(t['m'],t['m'],t['amplitude'],t['alpha'])
        out['GT_segmentation_calibration'].append(dict(parent_id=r['parent_id'],family=r['family'],
            GT_self_Otsu_IoU=d['IoU_Otsu'],GT_self_passes_frozen_shape_gate=d['IoU_Otsu']>=.7))
    cal=out['GT_segmentation_calibration']
    out['GT_self_shape_gate_count']=sum(r['GT_self_passes_frozen_shape_gate'] for r in cal)
    out['GT_segmentation_warning']='Posthoc calibration; frozen scores/gates remain unchanged. Otsu(alpha) need not equal alpha>=.5 for smooth/nonbinary materials. This is not a hard upper bound on reconstruction IoU.'
    (S/'results/FINAL_ANALYSIS.json').write_text(json.dumps(out,indent=2))
    print(json.dumps(dict(output='results/FINAL_ANALYSIS.json',GT_self_shape_gate_count=out['GT_self_shape_gate_count'],n=len(cal))))

if __name__=='__main__':main()
