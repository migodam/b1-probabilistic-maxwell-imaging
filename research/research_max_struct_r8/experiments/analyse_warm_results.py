"""Read-only warm-start analysis. Default: public validation development only.

No solver, sealing, scoring, model loading, or test-truth access in default mode.
Confirmation requires existing complete 750-task seal and scored outputs. Private
image access additionally requires --authorized-private-material-images supplied
only after explicit owner authorization. Never select a task or iterate by GT.
"""
from pathlib import Path
import argparse
import csv
import json
import time
PROCESS_BEGIN_CPU=time.process_time()
PROCESS_BEGIN_WALL=time.perf_counter()
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[1]
OUT = ROOT / 'delegated/warm_analysis'
FIG = ROOT / 'FIGURES'
COSTS = ('solver_cpu_seconds', 'warm_deployment_cpu_seconds',
         'actual_workflow_cpu_seconds', 'cold_deployment_cpu_seconds', 'LU', 'complex_RHS')
METRICS = ('material_nrmse', 'otsu_iou', 'original_Phi', 'relative_residual')


def load(path):
    return json.loads(Path(path).read_text())


def save(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def table(name, rows):
    path = OUT / name
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as f:
        fields = list(dict.fromkeys(k for r in rows for k in r))
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def finite(x):
    return x is not None and np.isfinite(x)


def summary(rows, include_costs=True):
    out = dict(objects=len({r['case_id'] for r in rows}), tasks=len(rows),
               failed=sum(r.get('status') == 'FAILED' for r in rows),
               material_successes=sum(bool(r.get('material_image_success')) for r in rows))
    for k in METRICS + (COSTS if include_costs else ()):
        vals = [r[k] for r in rows if finite(r.get(k))]
        out[k + '_observed'] = len(vals)
        out[k + '_mean'] = float(np.mean(vals)) if vals else None
        out[k + '_median'] = float(np.median(vals)) if vals else None
    return out


def bootstrap_pairs(rows, reference, resamples=10000):
    """One row per object/method. Resample objects, never iterations or sources."""
    by = {(r['case_id'], r['method']): r for r in rows}
    groups = defaultdict(list)
    for r in rows:
        groups[(r['cohort'], r['method'])].append(r)
    out = []
    rng = np.random.default_rng(20261106)
    for (cohort, method), rr in sorted(groups.items()):
        if method == reference:
            continue
        for metric in METRICS:
            pairs = [(r, by.get((r['case_id'], reference))) for r in rr]
            d = np.array([a[metric] - b[metric] for a, b in pairs
                          if b and finite(a.get(metric)) and finite(b.get(metric))])
            boots = np.mean(d[rng.integers(0, len(d), (resamples, len(d)))], axis=1) if len(d) else []
            ci = np.quantile(boots, [.025, .975]).tolist() if len(d) else [None, None]
            out.append(dict(cohort=cohort, method=method, reference=reference,
                            metric=metric, available_objects=len(rr), paired_objects=len(d),
                            missing_or_failed_pairs=len(rr)-len(d),
                            mean_paired_difference=float(np.mean(d)) if len(d) else None,
                            percentile95_low=ci[0], percentile95_high=ci[1],
                            resamples=resamples, unit='independent object',
                            interpretation='descriptive paired effect; no p value or significance declaration'))
    return out


def grouped_tables(rows, label):
    groups = defaultdict(list)
    for r in rows:
        groups[(r['cohort'], r['method'])].append(r)
    table(label + '_cohort_summary.csv', [dict(cohort=c, method=m, **summary(rr, include_costs=False)) for (c,m),rr in sorted(groups.items())])
    groups = defaultdict(list)
    for r in rows:
        groups[(r['cohort'], r['hardware'], r['method'])].append(r)
    table(label + '_absolute_cost_by_hardware.csv', [dict(cohort=c, hardware=h, method=m, **summary(rr)) for (c,h,m),rr in sorted(groups.items())])
    table(label + '_object_rows.csv', rows)
    if any(r['method']=='physical_lbfgs / uniform' for r in rows):
        table(label + '_paired_bootstrap.csv', bootstrap_pairs(rows, 'physical_lbfgs / uniform'))


def training_cost_accounting():
    """Use outer process totals once: nested forward/features are diagnostic only."""
    parts=[]
    for split, expected in (('train',1536),('validation',192)):
        path=ROOT/f'audit/generation_{split}.json'; rec=load(path)
        if rec['objects']!=expected:
            raise ValueError('Frozen original physics training split size changed')
        parts.append(dict(source=str(path.relative_to(ROOT)), scope=split+' generation + features + setup + writing',
                          cpu_seconds=rec['cpu_seconds'], wall_seconds=rec['wall_seconds'],
                          objects=rec['objects'], nested_costs_not_added_again=True))
    path=ROOT/'RESULTS/training/physics/RECEIPT.json'; rec=load(path)
    if (rec['epochs'],rec['training_objects'],rec['validation_objects'])!=(120,1536,192):
        raise ValueError('Original physics-network training contract changed')
    parts.append(dict(source=str(path.relative_to(ROOT)),scope='physics network full 120 epochs including validation selection',
                      cpu_seconds=rec['cpu_seconds'],wall_seconds=rec['wall_seconds']))
    original=dict(hardware='Mac',cpu_seconds=sum(p['cpu_seconds'] for p in parts),
                  wall_seconds=sum(p['wall_seconds'] for p in parts),components=parts,
                  excludes='Born-network ablation, domain repair, test generation, test inference and inverse solves',
                  deployment_semantics='Original physics NN numerator only; no hypothetical retraining per object')
    # Broader scientific setup ledger: each independent outer receipt exactly once.
    allparts=[]
    paths=list((ROOT/'audit').glob('generation_*.json'))+list((ROOT/'RESULTS/training').glob('*/RECEIPT.json'))
    for path in sorted(paths):
        rec=load(path)
        if finite(rec.get('cpu_seconds')):
            allparts.append(dict(source=str(path.relative_to(ROOT)),cpu_seconds=rec['cpu_seconds'],
                                 wall_seconds=rec.get('wall_seconds'),status=rec.get('status')))
    save(OUT/'TRAINING_COST_ACCOUNTING.json',dict(original_physics_network=original,
         research_all_generation_and_training=dict(cpu_seconds=sum(p['cpu_seconds'] for p in allparts),
             components=allparts,scope='All available generation and all network training receipts including Born/domain; setup subtotal, excludes inverse solves and unisolated orchestration',
             hardware='Mac',not_a_deployment_amortization_numerator=True),
         field_semantics=dict(feature_cpu_seconds='Original actual Mac preprocessing time',
             inference_cpu_seconds='Original actual Mac inference time',
             deployment_prefix_cpu_seconds='Independent measured preprocessing + inference on the solver hardware',
             actual_workflow_cpu_seconds='Original Mac preprocessing + actual solver CPU; mixed-machine work when solver is XINAN')))
    return original


def amortization(rows):
    """Expected cost accounting from paired object MEAN CPU savings, not median."""
    original=training_cost_accounting(); total=original['cpu_seconds']
    methods=('physical_lbfgs','cie_single','cie_dual','a_plus_sourcewise','a_plus')
    by={(r['case_id'],r['arm'],r['initialization']):r for r in rows}
    out=[];rng=np.random.default_rng(2026110607)
    for arm in methods:
        groups=defaultdict(list)
        for r in rows:
            if r['arm']==arm and r['initialization']=='nn':
                groups[(r['cohort'],r['hardware'])].append(r)
        for (cohort,hardware),rr in sorted(groups.items()):
            pairs=[(r,by.get((r['case_id'],arm,'uniform'))) for r in rr]
            good=[(a,b) for a,b in pairs if b and b['hardware']==hardware
                  and finite(a.get('warm_deployment_cpu_seconds')) and finite(b.get('warm_deployment_cpu_seconds'))]
            savings=np.array([b['warm_deployment_cpu_seconds']-a['warm_deployment_cpu_seconds'] for a,b in good])
            mean=float(np.mean(savings)) if len(savings) else None
            ci=[None,None]
            if len(savings):
                boot=np.mean(savings[rng.integers(0,len(savings),(10000,len(savings)))],axis=1)
                ci=np.quantile(boot,[.025,.975]).tolist()
            r=dict(cohort=cohort,hardware=hardware,arm=arm,comparison='paired uniform minus original physics NN',
                available_objects=len(rr),paired_cost_objects=len(good),missing_cost_pairs=len(rr)-len(good),
                pairs_with_failed_branch=sum(a['status']=='FAILED' or b['status']=='FAILED' for a,b in good),
                mean_paired_saved_deployment_CPU_seconds=mean,
                mean_savings_percentile95_low=ci[0],mean_savings_percentile95_high=ci[1],
                original_Mac_training_and_generation_CPU_seconds=total,
                statistical_unit='independent material object; arithmetic mean paired savings',
                cost_scope='Same-machine independently measured warm prefix plus full actual solver CPU; training separate',
                CPU_savings_do_not_establish_material_quality_advantage=True)
            if mean is None:
                r.update(status='UNAVAILABLE_NO_PAIRED_COSTS',K=None)
            elif mean<=0:
                r.update(status='NO_FINITE_AMORTIZATION_ADVANTAGE',K=None,
                         reason='Mean NN deployment CPU saving is nonpositive')
            elif hardware=='Mac':
                r.update(status='SAME_MACHINE_CPU_BOOKKEEPING_ONLY',K=total/mean,
                         whole_objects_to_repay_CPU=int(np.ceil(total/mean)),
                         interpretation='Mac CPU-accounting amortization; no monetary prices specified')
            else:
                r.update(status='PRICE_RATIO_REQUIRED',K=None,
                         K_coefficient_per_cMac_over_cXINAN=total/mean,
                         K_formula='K = (cMac * Ttrain_Mac_CPU) / (cXINAN * mean_paired_saved_XINAN_CPU)',
                         interpretation='Hardware CPU seconds are not economically interchangeable; cMac/cXINAN unspecified, so no numeric economic break-even')
            out.append(r)
    table('PRIMARY_training_amortization_mean_paired_CPU.csv',out)
    save(OUT/'PRIMARY_training_amortization.json',dict(rows=out,original_physics_training=original,
         assumptions='Expected savings use object-paired arithmetic means separately by cohort/hardware. No median-based total-cost claim; no monetary break-even without prices. Failed paid work retained when finite.',
         field_semantics='feature_cpu_seconds/inference_cpu_seconds are original Mac measurements; deployment_prefix_cpu_seconds is independently measured on each solver hardware.'))


def training_development():
    training_cost_accounting()
    FIG.mkdir(exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(11,4))
    receipts = {}
    for name in ('physics', 'born'):
        d = ROOT / 'RESULTS/training' / name
        if not (d/'history.json').exists():
            continue
        h, rec = load(d/'history.json'), load(d/'RECEIPT.json')
        receipts[name] = rec
        for ax, key in zip(axes, ('train_normalized_MSE','validation_normalized_MSE')):
            ax.plot([v['epoch'] for v in h], [v[key] for v in h], label=name)
            ax.set(xlabel='Epoch', ylabel=key.replace('_',' '), yscale='log')
            ax.axvline(rec['best_epoch'], alpha=.25, linestyle=':')
            ax.legend()
    fig.suptitle('Training / public validation only; epoch selected by validation MSE')
    fig.tight_layout(); fig.savefig(FIG/'warm_dev_training.png', dpi=170); plt.close(fig)
    rows = load(ROOT/'RESULTS/warm_development_analysis/SCORES.json')['rows']
    cases = sorted({r['case_id'] for r in rows})
    methods = sorted({r['arm']+' / '+r['initialization'] for r in rows})
    fig, axes = plt.subplots(2,1,figsize=(12,9))
    for method in methods:
        rr = {r['case_id']: r for r in rows if r['arm']+' / '+r['initialization']==method}
        for ax,key in zip(axes,('material_nrmse','cpu_seconds')):
            ax.plot(range(len(cases)), [rr.get(c,{}).get(key,np.nan) for c in cases], 'o-', label=method)
            ax.set(xticks=range(len(cases)), xticklabels=[c[-4:] for c in cases], ylabel=key)
    axes[0].set_yscale('log');axes[0].legend(ncol=2,fontsize=8)
    axes[1].set_xlabel('Predefined development validation object')
    fig.suptitle('8 public development objects: unequal final budgets; exploratory only')
    fig.tight_layout();fig.savefig(FIG/'warm_dev_objects.png',dpi=170);plt.close(fig)
    table('DEV_object_scores.csv', rows)
    save(OUT/'DEV_RECEIPT.json',dict(status='PUBLIC_DEVELOPMENT_PLOTTED',objects=len(cases),
         rows=len(rows),training=receipts,truth_access='NONE: reused existing validation metric rows',
         new_physical_calls=0,test_gt_read=False))
    (OUT/'DEV_REPORT.md').write_text(
        '# Warm-start development figures\n\n'
        'Only the eight predefined public validation objects and existing training histories were read. '
        'These are exploratory development observations, not confirmation statistics. No test truth, '
        'new inverse solve, sealing, scoring, or new model training was performed.\n\n'
        '- `FIGURES/warm_dev_training.png`: full recorded epochs, training and validation MSE; dotted lines indicate validation-selected epochs.\n'
        '- `FIGURES/warm_dev_objects.png`: all eight objects, all six existing method/initialization combinations, material NRMSE and actual solver CPU. '
        'Budgets differ (direct physical permits 6000 steps, auxiliary corrector 3000), so the final-point panel is not an equal-budget claim.\n\n'
        'The NN direct physical branch retains useful development trajectories on objects 0002/0003/0005/0006 '
        'that the auxiliary pre-corrector selection can discard. Object 0004 remains a failure; '
        '0007 is already well recovered from uniform. The future public-Phi-selected two-branch comparison '
        'must charge both complete solvers and one shared NN prefix. These development observations do not '
        'establish confirmation-set performance or novelty.\n\n'
        'Cost semantics: feature_cpu_seconds and inference_cpu_seconds retain the original actual Mac timing. '
        'deployment_prefix_cpu_seconds is the separate independently measured preprocessing/inference time '
        'on the solver machine. These quantities cannot be interchanged. TRAINING_COST_ACCOUNTING.json '
        'uses the outer 1536-train/192-validation generation receipts once (including their features/setup), '
        'plus all 120 original physics-network epochs. Born and domain-repair costs appear only in the '
        'broader scientific setup ledger. Confirmation amortization will use independent-object mean '
        'paired uniform-minus-NN deployment CPU savings, separately by cohort/hardware. Mac CPU bookkeeping '
        'can give K; XINAN requires K=(cMac*Ttrain)/(cXINAN*mean_saved_CPU), with no numeric economic '
        'break-even until the price ratio is supplied. A nonpositive mean saving gives no finite advantage.\n')


def portfolio_rows(rows):
    by = {r['task_id']:r for r in rows}
    out = []
    for p in load(ROOT/'audit/WARM_PORTFOLIO_PUBLIC_SELECTION.json')['records']:
        rr = [by[x] for x in p['candidate_ids']]
        good = [r for r in rr if r['status']=='SEALED_PUBLIC_PREDICTION']
        expected = min(good, key=lambda r:(r['original_Phi'],r['task_id']))['task_id'] if good else None
        if expected != p['selected_id']:
            raise ValueError('Public portfolio selection mismatch')
        selected = by[p['selected_id']] if p['selected_id'] else rr[0]
        r = dict(selected, method='portfolio / '+p['method'],task_id='portfolio / '+p['case_id']+' / '+p['method'],
                 candidate_ids=p['candidate_ids'], selected_id=p['selected_id'],
                 status='SEALED_PUBLIC_PREDICTION' if good else 'FAILED')
        hosts = sorted({s['hardware'] for s in rr});r['hardware']='+'.join(hosts)
        # Sum actual work including failed branches; never hide an unavailable cost.
        for key in COSTS:
            r[key] = sum(s[key] for s in rr) if all(finite(s.get(key)) for s in rr) else None
        # Same NN prefix is shared exactly once. Never subtract a Mac prefix
        # from a same-machine XINAN deployment estimate.
        shared_nn = all(s['initialization'].startswith('nn') for s in rr)
        for key,prefixkey in [('warm_deployment_cpu_seconds','deployment_prefix_cpu_seconds'),
                              ('actual_workflow_cpu_seconds','actual_workflow_prefix_cpu_seconds'),
                              ('cold_deployment_cpu_seconds','cold_prefix_cpu_seconds')]:
            if shared_nn and r[key] is not None:
                vals=[s.get(prefixkey) for s in rr]
                if key=='cold_deployment_cpu_seconds':
                    vals=[s['cold_deployment_cpu_seconds']-s['solver_cpu_seconds']
                          if finite(s.get('cold_deployment_cpu_seconds')) and finite(s.get('solver_cpu_seconds')) else None for s in rr]
                if key=='actual_workflow_cpu_seconds':
                    vals=[s.get('feature_cpu_seconds',0)+s.get('inference_cpu_seconds',0) for s in rr]
                if len(hosts)==1 and all(finite(v) for v in vals) and np.allclose(vals, vals[0]):
                    r[key]-=vals[0]
                else:
                    r[key]=None
                    r[key+'_unavailable_reason']='Shared prefix timing absent/inconsistent or mixed machines; no fabricated deduction'
        r['cost_rule']='Both complete branches charged; NN prefix shared once only with explicit measured matching prefix'
        out.append(r)
    return out


def budget_curves(rows, budgets):
    """Fixed common budgets; latest paid original-Phi incumbent. No GT choice."""
    meta={r['task_id']:r for r in rows}
    groups=defaultdict(lambda:defaultdict(list))
    for b in budgets:
        if b['task_id'] in meta:
            m=meta[b['task_id']]
            groups[(m['cohort'],m['hardware'],m['method'])][b['task_id']].append(b)
    records=[]; omissions=[]
    for (cohort,host,method),tasks in sorted(groups.items()):
        for axis in ('CPU_seconds','solver_CPU_seconds','LU','complex_RHS'):
            def cost(b):
                if axis not in ('CPU_seconds','solver_CPU_seconds'): return b.get(axis)
                m=meta[b['task_id']]
                if m.get('deployment_preprocessing_hardware')==host and finite(m.get('deployment_prefix_cpu_seconds')):
                    return b['CPU_seconds']-(m['deployment_prefix_cpu_seconds'] if axis=='solver_CPU_seconds' else 0)
                # Accept only explicitly documented same-hardware deployment CPU.
                if finite(b.get('same_hardware_CPU_seconds')): return b['same_hardware_CPU_seconds']
                if b.get('cpu_hardware')==host and finite(b.get('CPU_seconds')): return b['CPU_seconds']
                return None
            tracks={t:sorted(bs,key=lambda b:b['index']) for t,bs in tasks.items()}
            if any(not all(finite(cost(b)) for b in bs) for bs in tracks.values()):
                omissions.append(dict(cohort=cohort,hardware=host,method=method,axis=axis,
                                      reason='Missing explicitly same-hardware paid CPU coordinates'));continue
            # A common grid across all methods in this cohort/hardware/axis.
            peers=[bs for (cc,hh,_),tt in groups.items() if (cc,hh)==(cohort,host) for bs in tt.values()]
            peers=[bs for bs in peers if all(finite(cost(b)) for b in bs)]
            lo=max(min(cost(b) for b in bs) for bs in peers)
            hi=min(max(cost(b) for b in bs) for bs in peers)
            if hi<lo:
                omissions.append(dict(cohort=cohort,hardware=host,method=method,axis=axis,reason='No common paid-budget interval across methods'))
                continue
            for budget in np.linspace(lo,hi,50):
                selected=[]
                for bs in tracks.values():
                    paid=[b for b in bs if cost(b)<=budget]
                    # Original Phi determines incumbent; latest index breaks an exact tie.
                    selected.append(min(paid,key=lambda b:(b['original_Phi'],-b['index'])))
                for metric in ('original_Phi','material_nrmse','otsu_iou'):
                    vals=[b[metric] for b in selected if finite(b.get(metric))]
                    records.append(dict(cohort=cohort,hardware=host,method=method,axis=axis,budget=float(budget),
                        metric=metric,objects_with_checkpoints=len(tracks),
                        all_objects=len({r['case_id'] for r in rows if r['cohort']==cohort and r['hardware']==host and r['method']==method}),
                        mean=float(np.mean(vals)) if vals else None,
                        median=float(np.median(vals)) if vals else None,
                        selection='Minimum paid original Phi; latest index tie-break; no interpolation or GT selection'))
    table('PRIMARY_equal_paid_budget_curves.csv',records);save(OUT/'BUDGET_OMISSIONS.json',omissions)
    for metric in ('original_Phi','material_nrmse'):
        keys=sorted({(r['cohort'],r['hardware'],r['axis']) for r in records})
        with PdfPages(FIG/f'warm_primary_budget_{metric}.pdf') as pdf:
            for cohort,host,axis in keys:
                fig,ax=plt.subplots(figsize=(12,6))
                for method in sorted({r['method'] for r in records}):
                    rr=[r for r in records if (r['cohort'],r['hardware'],r['axis'],r['metric'],r['method'])==(cohort,host,axis,metric,method)]
                    if rr:ax.plot([r['budget'] for r in rr],[r['median'] for r in rr],label=method)
                ax.set(xlabel='Paid '+axis,ylabel='Object median '+metric,title=f'{cohort} / {host}; sparse paid checkpoints')
                ax.legend(fontsize=7,ncol=2);fig.tight_layout();pdf.savefig(fig);plt.close(fig)


def paired_budget_endpoints(rows, budgets):
    """Per-object common allowance min(final A,final B); sparse paid incumbents.

    The allowance is identical. Actual checkpoint spend can differ, so unused
    allowance and the next checkpoint gap are explicit. A stopped branch never
    receives invented iterations. Comparisons are stratified by actual hardware.
    """
    meta={r['task_id']:r for r in rows}
    by={(r['case_id'],r['method']):r for r in rows}
    tracks=defaultdict(list)
    for b in budgets:
        if b['task_id'] in meta:tracks[b['task_id']].append(b)
    for bs in tracks.values():bs.sort(key=lambda b:b['index'])
    contrasts=[]
    for arm in sorted({r['arm'] for r in rows if r['initialization']=='nn'}):
        contrasts.append((arm+' / nn',arm+' / uniform'))
    contrasts += [('a_plus_sourcewise / nn','cie_dual / nn'),
                  ('a_plus / nn','global_dual / nn'),
                  ('global_dual / nn','cie_dual / nn'),
                  ('a_plus / nn','a_plus_sourcewise / nn'),
                  ('a_plus / nn','physical_lbfgs / nn'),
                  ('a_plus_sourcewise / nn','physical_lbfgs / nn')]
    axes={'CPU_seconds':'warm_deployment_cpu_seconds',
          'solver_CPU_seconds':'solver_cpu_seconds','LU':'LU','complex_RHS':'complex_RHS'}
    def coordinate(b,axis):
        if axis not in ('CPU_seconds','solver_CPU_seconds'):return b.get(axis)
        r=meta[b['task_id']]
        if r.get('deployment_preprocessing_hardware')!=r['hardware'] or not finite(r.get('deployment_prefix_cpu_seconds')):
            return None
        return b['CPU_seconds']-(r['deployment_prefix_cpu_seconds'] if axis=='solver_CPU_seconds' else 0)
    def choose(rr,axis,allowance):
        bs=tracks[rr['task_id']]
        if not bs or any(not finite(coordinate(b,axis)) for b in bs):return None
        paid=[b for b in bs if coordinate(b,axis)<=allowance]
        if not paid:return None
        last=max(paid,key=lambda b:b['index'])
        incumbent=min(paid,key=lambda b:(b['original_Phi'],-b['index']))
        nextcosts=[coordinate(b,axis) for b in bs if coordinate(b,axis)>allowance]
        return dict(incumbent=incumbent,last_paid=coordinate(last,axis),
                    unused=max(0.,allowance-coordinate(last,axis)),
                    next_checkpoint_cost=min(nextcosts) if nextcosts else None,
                    gap=(min(nextcosts)-coordinate(last,axis)) if nextcosts else None,
                    age=coordinate(last,axis)-coordinate(incumbent,axis),
                    at_terminal_checkpoint=last['index']==max(b['index'] for b in bs))
    result=[]
    for method,reference in contrasts:
        for case in sorted({r['case_id'] for r in rows}):
            a,b=by.get((case,method)),by.get((case,reference))
            if not a or not b:continue
            for axis,finalkey in axes.items():
                row=dict(case_id=case,cohort=a['cohort'],hardware=a['hardware'],method=method,
                         reference=reference,axis=axis,method_status=a['status'],reference_status=b['status'],
                         selection='Original-Phi incumbent among checkpoints <= common allowance; latest index tie-break; never GT',
                         interpretation='Same at-most allowance, not necessarily equal actual expenditure')
                if a['hardware']!=b['hardware']:
                    row.update(status='MIXED_HARDWARE_PAIR_OMITTED',reference_hardware=b['hardware']);result.append(row);continue
                if not finite(a.get(finalkey)) or not finite(b.get(finalkey)):
                    row.update(status='MISSING_FINAL_COST');result.append(row);continue
                allowance=min(a[finalkey],b[finalkey]);row['common_allowance']=allowance
                aa,bb=choose(a,axis,allowance),choose(b,axis,allowance)
                if aa is None or bb is None:
                    row.update(status='MISSING_PAID_CHECKPOINT_AT_ALLOWANCE');result.append(row);continue
                row['status']='PAIRED_PAID_INCUMBENTS'
                for prefix,state in [('method',aa),('reference',bb)]:
                    row[prefix+'_last_paid_checkpoint_cost']=state['last_paid']
                    row[prefix+'_unused_allowance']=state['unused']
                    row[prefix+'_next_checkpoint_cost']=state['next_checkpoint_cost']
                    row[prefix+'_checkpoint_bracketing_gap']=state['gap']
                    row[prefix+'_incumbent_age_in_paid_cost']=state['age']
                    row[prefix+'_terminal_checkpoint_reached']=state['at_terminal_checkpoint']
                    row[prefix+'_incumbent_index']=state['incumbent']['index']
                    for metric in ('original_Phi','material_nrmse','otsu_iou'):
                        row[prefix+'_'+metric]=state['incumbent'].get(metric)
                for metric in ('original_Phi','material_nrmse','otsu_iou'):
                    va,vb=row.get('method_'+metric),row.get('reference_'+metric)
                    row['difference_'+metric]=va-vb if finite(va) and finite(vb) else None
                result.append(row)
    table('PRIMARY_object_matched_at_most_budget_endpoints.csv',result)
    groups=defaultdict(list)
    for r in result:groups[(r['cohort'],r['hardware'],r['method'],r['reference'],r['axis'])].append(r)
    effects=[];rng=np.random.default_rng(2026110608)
    for (cohort,hardware,method,reference,axis),rr in sorted(groups.items()):
        for metric in ('original_Phi','material_nrmse','otsu_iou'):
            valid=[r for r in rr if finite(r.get('difference_'+metric))]
            delta=np.array([r['difference_'+metric] for r in valid])
            ci=[None,None]
            if len(delta):
                boot=np.mean(delta[rng.integers(0,len(delta),(10000,len(delta)))],axis=1)
                ci=np.quantile(boot,[.025,.975]).tolist()
            effects.append(dict(cohort=cohort,hardware=hardware,method=method,reference=reference,axis=axis,
                metric=metric,objects=len(rr),paired_objects=len(valid),missing_pairs=len(rr)-len(valid),
                mean_paired_difference=float(np.mean(delta)) if len(delta) else None,
                percentile95_low=ci[0],percentile95_high=ci[1],
                mean_method_unused_allowance=float(np.mean([r['method_unused_allowance'] for r in valid])) if valid else None,
                mean_reference_unused_allowance=float(np.mean([r['reference_unused_allowance'] for r in valid])) if valid else None,
                sampling_unit='independent object; 10000 descriptive paired resamples',
                budget_rule='Per-object min(final actual cost of the two complete methods); same allowance, sparse checkpoints, no interpolation'))
    table('PRIMARY_object_matched_budget_paired_bootstrap.csv',effects)


def material_gallery(rows, tasks):
    """Authorized offline gallery, common full dynamic range; no score-based subset."""
    methods=sorted({r['method'] for r in rows})
    fixed_png={f'r8_test_{i:04d}' for i in (0,3,24,27,36,39)}
    with PdfPages(FIG/'warm_material_all_42.pdf') as allpdf:
        for cohort in sorted({r['cohort'] for r in rows}):
            cases=sorted({r['case_id'] for r in rows if r['cohort']==cohort})
            with PdfPages(FIG/f'warm_material_all_{cohort}.pdf') as pdf:
                for case in cases:
                    with np.load(ROOT/f'sealed_truth/warm_test/{case}.npz',allow_pickle=False) as z:
                        truth=z['m'].reshape(24,24).copy()
                    material={}; case_rows={r['method']:r for r in rows if r['case_id']==case}
                    for method,r in case_rows.items():
                        if r['status']=='SEALED_PUBLIC_PREDICTION':
                            p=Path(tasks[r['task_id']]['output_dir']);p=p if p.is_absolute() else REPO/p
                            with np.load(p/'prediction.npz',allow_pickle=False) as z:
                                material[method]=z['m'].reshape(24,24).copy()
                    images=[truth]+list(material.values())
                    vmin=min(.02,*(float(m.min()) for m in images))
                    vmax=max(.02,*(float(m.max()) for m in images))
                    fig,axes=plt.subplots(4,4,figsize=(14,14));axes=axes.ravel()
                    for ax in axes:ax.axis('off')
                    im=axes[0].imshow(truth,vmin=vmin,vmax=vmax,cmap='viridis')
                    axes[0].set_title(f'GT (offline authorized)\nmin={truth.min():.3g}; max={truth.max():.3g}',fontsize=9)
                    for ax,method in zip(axes[1:],methods):
                        r=case_rows[method]; m=material.get(method)
                        score=f"NRMSE={r['material_nrmse']:.3g}" if finite(r.get('material_nrmse')) else 'FAILED / no score'
                        limits=f'min={m.min():.3g}; max={m.max():.3g}' if m is not None else 'FAILED / no image'
                        ax.set_title(method+'\n'+score+'; '+limits,fontsize=8)
                        if m is not None:ax.imshow(m,vmin=vmin,vmax=vmax,cmap='viridis')
                    fig.suptitle(f'{cohort} / {case}; all 15 primary arms; failures retained\nFull shared object range [{vmin:g}, {vmax:g}], including all prediction peaks')
                    fig.subplots_adjust(top=.91,bottom=.025,left=.025,right=.93,hspace=.31,wspace=.09)
                    cax=fig.add_axes([.95,.12,.012,.7]);fig.colorbar(im,cax=cax,label='Material contrast, shared full range')
                    pdf.savefig(fig);allpdf.savefig(fig)
                    if case in fixed_png:fig.savefig(FIG/f'warm_material_fixed_{case}.png',dpi=160)
                    plt.close(fig)


def confirmation(images=False):
    seal=load(ROOT/'audit/WARM_ALL_PUBLIC_SEAL.json')
    score=load(ROOT/'RESULTS/warm_scoring/SCORES.json')
    receipt=load(ROOT/'RESULTS/warm_scoring/RECEIPT.json')
    if seal['status']!='ALL_FROZEN_WARM_TASKS_SEALED' or seal['tasks']!=750 or receipt['status']!='SCORED_ALL_FROZEN_PREDICTIONS':
        raise ValueError('Complete 750-task public seal and independent scoring required')
    rows=score['rows'];tasks={t['task_id']:t for t in load(ROOT/'configs/warm_all_manifest.json')['tasks']}
    primaryids={t['task_id'] for t in load(ROOT/'configs/warm_primary_manifest.json')['tasks']}
    if len(rows)!=750 or len({r['task_id'] for r in rows})!=750 or set(tasks)!={r['task_id'] for r in rows}:
        raise ValueError('Expected exact frozen 750 identities')
    primary=[r for r in rows if r['task_id'] in primaryids];extra=[r for r in rows if r['task_id'] not in primaryids]
    if (len(primary)!=630 or len(extra)!=120 or len({r['case_id'] for r in extra})!=8
        or len({r['case_id'] for r in primary})!=42
        or any(sum(r['case_id']==c for r in primary)!=15 for c in {r['case_id'] for r in primary})):
        raise ValueError('Expected 630=42x15 main rows, 120=8x15 extra rows')
    grouped_tables(primary,'PRIMARY_630')
    amortization(primary)
    paired=[]
    for arm in sorted({r['arm'] for r in primary}):
        paired+=bootstrap_pairs([r for r in primary if r['arm']==arm and r['initialization'] in ('uniform','nn')],arm+' / uniform')
    table('PRIMARY_NN_minus_uniform_object_paired.csv',paired)
    # Additional arms and their matched primary controls only within the fixed 8 objects.
    subset={r['case_id'] for r in extra}
    grouped_tables([r for r in rows if r['case_id'] in subset],'FIXED8_120_EXTRA_WITH_MATCHED_CONTROLS')
    portfolios=portfolio_rows(primary);grouped_tables(portfolios+[r for r in primary if r['arm']=='physical_lbfgs'],'PAID_PORTFOLIOS')
    # Pair neural/uniform portfolios directly, retaining object clusters and failures.
    effects=[]
    for method in load(ROOT/'configs/warm_portfolio_protocol.json')['methods']:
        effects+=bootstrap_pairs([r for r in portfolios if r['method'] in ('portfolio / '+method+' / uniform / retained_physics','portfolio / '+method+' / nn / retained_physics')], 'portfolio / '+method+' / uniform / retained_physics')
    table('PORTFOLIO_paired_bootstrap_vs_equal_two_candidate_uniform.csv',effects)
    factorial=[]
    cells={'cie_dual':('sourcewise','ordinary'),'global_dual':('global','ordinary'),
           'a_plus_sourcewise':('sourcewise','A+'),'a_plus':('global','A+')}
    for r in primary:
        if r['initialization']=='nn' and r['arm'] in cells:
            factorial.append(dict(r,current_scope=cells[r['arm']][0],readout=cells[r['arm']][1]))
    table('NN_2x2_object_rows.csv',factorial);grouped_tables(factorial,'NN_2x2')
    budget_rows=load(ROOT/'RESULTS/warm_scoring/BUDGET_SCORES.json')['rows']
    budget_curves(primary,budget_rows)
    paired_budget_endpoints(primary,budget_rows)
    if images:material_gallery(primary,tasks)
    (OUT/'CONFIRMATION_REPORT.md').write_text(
        '# Warm confirmation analysis\n\n'
        'Primary comparisons use 630 = 42 x 15 tasks, with ID/OOD/cross-grid cohorts separate. '
        'The additional 120 tasks and their matched controls are reported only on the fixed eight-object subset. '
        'Failures remain in denominators; available continuous pairs are explicitly counted.\n\n'
        'Cost fields: feature_cpu_seconds and inference_cpu_seconds are the original actual Mac timings. '
        'deployment_prefix_cpu_seconds is the independent preprocessing/inference measurement on each '
        'solver hardware. Solver CPU, same-machine warm deployment CPU, and actual mixed-machine workflow '
        'CPU have separate columns. Absolute CPU comparisons are stratified by hardware.\n\n'
        'Training amortization uses the original 1536 training / 192 validation generation outer receipts '
        '(including feature/setup work once) plus the full physics 120-epoch training receipt. Born/domain '
        'costs enter the broader scientific setup ledger, not this deployment numerator. K uses arithmetic '
        'mean independent-object paired uniform-minus-NN warm deployment CPU savings, never a median. '
        'Mac results are CPU bookkeeping. XINAN results require K=(cMac*Ttrain)/(cXINAN*mean_saved_CPU); '
        'the unspecified price ratio prevents a numeric economic break-even. Nonpositive mean savings '
        'give no finite amortization advantage. Costs of failed branches remain reported; saving time '
        'does not establish equivalent material quality.\n\n'
        'Two-branch portfolios use the already frozen public original-Phi selection and charge both full '
        'solvers with a shared NN prefix once. No unfrozen partial-branch budget schedule is invented. '
        'The strict all-method common-overlap curves are retained as a diagnostic and can end early. '
        'Additional object-paired endpoints use min(final_cost_A, final_cost_B) as the common allowance '
        'for each object and hardware, with last-paid original-Phi incumbents. Unused allowance and '
        'checkpoint bracketing gaps are explicit: equal at-most allowance does not imply equal actual spend. '
        'Primary budget curves use sparse paid original-Phi incumbents, not GT-picked iterations. '
        'Object-paired bootstrap intervals are descriptive and carry no significance declaration.\n')
    save(OUT/'CONFIRMATION_RECEIPT.json',dict(status='POST_SCORING_ANALYSIS',primary_tasks=630,extra_tasks=120,
        fixed_extra_objects=8,cohort_denominators='Separate ID 24 / OOD 12 / cross-grid 6; never pooled',
        bootstrap='10000 object-paired resamples; percentile 95% interval; no significance test',
        CPU='Compare absolute CPU only within hardware; solver, same-machine warm deployment and mixed-machine actual workflow are distinct',
        portfolio_budget_scope='Only final total full-branch costs; no unfrozen partial-branch execution schedule was invented',
        failures='Retained in task denominators; continuous metrics report observed and missing pairs',
        feature_CPU_semantics='feature_cpu_seconds/inference_cpu_seconds are original Mac actual timings; deployment_prefix_cpu_seconds is independent same-solver-machine measurement',
        truth_access='EXPLICITLY_AUTHORIZED_OFFLINE_MATERIAL_GALLERY' if images else 'NONE: existing score rows only',
        new_physical_calls=0,private_material_images=images))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=('dev','confirmation'),nargs='?',default='dev')
    parser.add_argument('--authorized-private-material-images',action='store_true',
                        help='Only after explicit owner authorization and existing complete seal+scoring')
    args=parser.parse_args();OUT.mkdir(parents=True,exist_ok=True);FIG.mkdir(exist_ok=True)
    if args.mode=='dev':
        if args.authorized_private_material_images:parser.error('Private images are not permitted in development mode')
        training_development()
    else:
        error=None
        try:confirmation(args.authorized_private_material_images)
        except Exception as exc:
            error=type(exc).__name__+': '+str(exc)
            raise
        finally:
            save(OUT/'CONFIRMATION_EXECUTION.json',dict(status='FAILED' if error else 'COMPLETE',error=error,
                cpu_seconds=time.process_time()-PROCESS_BEGIN_CPU,wall_seconds=time.perf_counter()-PROCESS_BEGIN_WALL,
                LU=0,complex_RHS=0,cost_measurement='Actual outer analysis including imports; solver and training receipts are reused, not charged again'))
