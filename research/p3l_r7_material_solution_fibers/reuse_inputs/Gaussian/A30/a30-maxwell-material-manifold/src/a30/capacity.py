"""Complete exposed-cohort capacity run. Target-bearing fitting is OFFLINE only."""
from pathlib import Path
from time import perf_counter
import json
import numpy as np
from .data import online_e0,train_deltas,offline_targets
from .representations import prepare_train,freeze_case,fit_capacity,METHODS,BUDGETS
from .frozen_a26_metrics import image_metrics
from .hull_audit import csv_write,topology_metrics
from .statistics import paired,paired_difference,COMPLEX


def serial(x):
    if isinstance(x,np.ndarray):return serial(x.tolist())
    if isinstance(x,np.generic):return serial(x.item())
    if isinstance(x,complex):return [x.real,x.imag]
    if isinstance(x,dict):return {str(k):serial(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [serial(v) for v in x]
    return x


def run(root):
    root=Path(root).resolve();out=root/'results/a30/capacity';out.mkdir(parents=True,exist_ok=True)
    rows,n3,bg=online_e0(root);train,train_meta=train_deltas(root)
    model=prepare_train(train);config=json.loads((root/'a30_config.json').read_text())['E0']
    settings=dict(amplitude_scale=1.2,max_starts=config['starts_max'],max_iterations=config['iterations_per_start_max'])
    # Every discrete template/support decision is completed without a target argument.
    decisions=[]
    for i,meta in enumerate(rows):
        for method in METHODS:
            for budget in BUDGETS:
                f=freeze_case(model,method,budget,(n3[i]-bg[i])/1.2,bg[i],settings)
                starts=[{k:s[k] for k in ('id','eligible','reasons') if k in s} for s in f['starts']]
                decisions.append(dict(**meta,method=method,budget=budget,starts=starts,settings=f['settings']))
    (out/'TRAIN_AND_CHOICES.json').write_text(json.dumps(serial(dict(train=train_meta,model=model.metadata,
        decisions=decisions,target_read_for_choices=False)),indent=2)+'\n')
    truth,target_bg,target_n3=offline_targets(root,[x['case_id'] for x in rows])
    if not np.array_equal(target_bg,bg) or not np.array_equal(target_n3,n3):raise ValueError('BASELINE_BACKGROUND_MISMATCH')
    arrays=dict(case_id=np.asarray([x['case_id'] for x in rows]),truth=truth,background=bg,N3=n3)
    for method in METHODS:
        for budget in BUDGETS:arrays[f'{method}_{budget}']=np.full_like(truth,np.nan+1j*np.nan)
    metrics=[];fit_dir=out/'fits';fit_dir.mkdir(exist_ok=True)
    for i,meta in enumerate(rows):
        initial=image_metrics(n3[i],truth[i],bg[i]);top0=topology_metrics(n3[i],truth[i],bg[i])
        for method in METHODS:
            for budget in BUDGETS:
                start=perf_counter();fit=fit_capacity(model,method,budget,(n3[i]-bg[i])/1.2,bg[i],truth[i],settings)
                image=fit.get('image');valid=image is not None and fit['box']['feasible']
                rec=dict(**meta,method=method,budget=budget,status=fit['status'],valid=valid,
                    N3_nrmse=initial['delta_nrmse'],N3_edge=initial['edge_f1'],N3_iou=initial['support_iou'],
                    N3_component=top0['component_abs_error_25'],N3_hole=top0['hole_abs_error_25'],
                    fitwall_seconds=perf_counter()-start,optimizer_label=fit['global_vs_local'],
                    certificate=bool(fit.get('certificate',{}).get('certified',False)),
                    optimizer_iterations=fit['cost']['material_optimizer_iterations'],starts=fit['cost']['starts'])
                if valid:
                    arrays[f'{method}_{budget}'][i]=image
                    rec.update(image_metrics(image,truth[i],bg[i]),**topology_metrics(image,truth[i],bg[i]))
                else:
                    rec.update(delta_nrmse=None,edge_f1=None,support_iou=None)
                    with (root/'results/a30/FAILURE_LEDGER.jsonl').open('a') as f:
                        f.write(json.dumps(dict(stage='E0_capacity',case_id=meta['case_id'],method=method,budget=budget,status=fit['status']))+'\n')
                # Save actual trajectories without redundant full material arrays per start.
                detail={k:fit[k] for k in ('status','scope','method','budget','settings','best_start_id','params','losses',
                    'box','certificate','dimensions','global_vs_local','cost','offline_best_s_tail',
                    'fixed_support_box_vs_best_s_mse_gap') if k in fit}
                detail['fits']=[{k:v for k,v in s.items() if k not in ('image','delta')} for s in fit['fits']]
                (fit_dir/f"{meta['case_id']}__{method}_{budget}.json").write_text(json.dumps(serial(detail))+'\n')
                metrics.append(rec);csv_write(out/'CAPACITY_CASES.csv',metrics)
        np.savez_compressed(out/'CAPACITY_IMAGES.npz',**arrays)
        print(json.dumps(dict(completed_conditions=i+1,total_conditions=len(rows),case=meta['case_id'])),flush=True)
    summary={}
    for method in METHODS:
        for budget in BUDGETS:
            rr=[r for r in metrics if r['method']==method and r['budget']==budget]
            key=f'{method}_{budget}'
            if not all(r['valid'] for r in rr):
                summary[key]=dict(status='INCOMPLETE',valid=sum(r['valid'] for r in rr),total=len(rr));continue
            summary[key]=dict(status='COMPLETE',complex=paired(rr,'delta_nrmse'),
                family={f:paired(rr,'delta_nrmse',families={f}) for f in sorted(COMPLEX|{'smooth'})},
                edge=paired_difference(rr,'edge_f1','N3_edge'),iou=paired_difference(rr,'support_iou','N3_iou'),
                component=paired_difference(rr,'component_abs_error_25','N3_component'),
                holes=paired_difference(rr,'hole_abs_error_25','N3_hole'))
    (out/'CAPACITY_SUMMARY.json').write_text(json.dumps(serial(summary),indent=2)+'\n')
    # The original protocol explicitly admits full hull as a legal representation.
    hull=json.loads((root/'results/a30/hull/HULL_AUDIT_SUMMARY.json').read_text())['full_hull']
    passed=(hull['complex']['point']>=.1 and hull['complex']['ci95'][0]>0 and
        all(hull['family'][f]['point']>0 for f in COMPLEX) and hull['edge']['point']>0 and hull['iou']['point']>0)
    gate=dict(status='PASS' if passed else 'NO_GO',basis='ORIGINAL_A30_E0_FULL_HULL_CAPACITY_GATE',
        evidence='OFFLINE_ONLY_NOT_DEPLOYABLE',continuous_capacity_complete=len(metrics)==384,
        full_hull_complex_gain=hull['complex'],edge=hull['edge'],iou=hull['iou'],
        reason='Useful legal convex combinations exist; actual saved L weights do not recover them.' if passed else 'Original capacity/structure requirements were not met.',
        E1_authorized=bool(passed),E2_authorized=False)
    (root/'results/a30/E0_GATE.json').write_text(json.dumps(gate,indent=2)+'\n')
    print(json.dumps(gate),flush=True)


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path('.'));run(p.parse_args().root)
