"""Offline mechanism tables from sealed scores; no new physics or selection."""
from pathlib import Path
from collections import defaultdict
import argparse,json,csv
import numpy as np

def stat(values,median=False):
    a=[v for v in values if v is not None and np.isfinite(v)]
    return float(np.median(a)if median else np.mean(a))if a else None

def main():
    p=argparse.ArgumentParser();p.add_argument('--scores',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    scores=json.loads(a.scores.read_text());groups=defaultdict(list)
    for row in scores['rows']:
        groups[(row['run'],row.get('level','low'),row['noise'])].append(row)
    result=[]
    for (run,level,noise),rows in groups.items():
        raw=[next((c for c in r.get('candidates',[])if c['label']=='raw_1'),{})for r in rows]
        ret=[x for r in rows for x in r.get('correction_retention',[])if x['alpha']==1.]
        cur=[r['current_metrics']for r in rows if r.get('current_metrics')]
        selected=[r for r in rows if r.get('top1')]
        result.append(dict(run=run,level=level,noise=noise,n=len(rows),
            independent_parents=len({r['parent_id']for r in rows}),
            anchor_NRMSE=stat([r.get('anchor',{}).get('NRMSE')for r in rows]),
            raw_NRMSE=stat([c.get('NRMSE')for c in raw]),
            top1_NRMSE=stat([r['top1']['NRMSE']for r in selected]),
            qualified_bestK_NRMSE=stat([(r.get('qualified_bestK_OFFLINE')or{}).get('NRMSE')for r in rows]),
            all_bestK_NRMSE=stat([(r.get('all_bestK_OFFLINE')or{}).get('NRMSE')for r in rows]),
            raw_physical_count=sum(bool(c.get('physical_data_feasible'))for c in raw),
            top1_physical_count=len(selected),
            top1_joint_count=sum(r.get('joint_success',False)for r in rows),
            raw_full_finite_step_data_residual=stat([c.get('relative_residual')for c in raw]),
            raw_full_constitutive_defect=stat([(r.get('auxiliary')or{}).get('constitutive_relative_defect')for r in rows],True),
            full_correction_count=len(ret),
            correction_displacement_projection_retention_median=stat([r.get('displacement_projection_retention')for r in ret],True),
            correction_displacement_norm_retention_median=stat([r.get('displacement_norm_retention')for r in ret],True),
            positive_raw_material_gain_count=sum(r.get('raw_had_positive_material_gain',False)for r in ret),
            positive_material_gain_retention_median=stat([r.get('material_gain_retention')for r in ret],True),
            selected_RMS_displacement_mean=stat([r.get('selected_RMS_displacement_from_anchor')for r in selected]),
            current_projected_target_relative_error=stat([r.get('projected_target_current_error_relative')for r in cur]),
            current_auxiliary_relative_error=stat([r.get('auxiliary_current_relative_L2_error')for r in cur]),
            current_constitutive_defect_median=stat([r.get('constitutive_relative_defect')for r in cur],True),
            current_increment_visible_fraction=stat([r.get('current_increment_visible_fraction')for r in cur]),
            current_GT_error_null_fraction=stat([r.get('anchor_current_error_null_energy',0)/max(r.get('anchor_current_error_total_energy',0),1e-300)for r in cur]),
            finite_distinct_candidates_mean=stat([r.get('candidate_distinct_material_count')for r in rows])))
    a.output.mkdir(parents=True,exist_ok=True)
    out=dict(source_scores=str(a.scores),scope='Post-seal diagnostic, no GT selection and no new experiment',
             groups=result,interpretation='Current L2, auxiliary constitutive defect, material readout and full-VIE qualification are different quantities. Retention is conditional on a completed full-scale correction; positive gain retention uses only proposals with positive raw gain.')
    (a.output/'MECHANISM_SUMMARY.json').write_text(json.dumps(out,indent=2))
    with (a.output/'MECHANISM_SUMMARY.csv').open('w',newline='')as f:
        w=csv.DictWriter(f,fieldnames=list(result[0]));w.writeheader();w.writerows(result)
    print(json.dumps(dict(groups=len(result),output=str(a.output))))

if __name__=='__main__':main()
