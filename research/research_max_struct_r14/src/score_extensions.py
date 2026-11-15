"""Separately sealed offline extension cohorts; no forward solve, NN or optimizer.

--make-seal always writes authorized_scoring=false. Parent approval is required
before private truth/family labels can be imported or read.
"""
from __future__ import annotations
import argparse
from collections import Counter
from copy import deepcopy
import html
import importlib.util
import json
from pathlib import Path
import sys
import time
from datetime import datetime,timezone
import score_dualmetric as base

CONTROL_IDS=(2,5,7,16,17,24,29,37,38,42,44,45,8,20,32,46)
SCHEMAS={
    'independent24_v1':(24,('pair00','pair01','pair10','pair11','pair01_to11','cie','dual_cie','g19','dr','global_dual','direct')),
    'self_only16_v1':(16,('pair01','pair11','global_dual')),
    'stage_control16_v1':(16,('pair01_to11','pair11_to01')),
}


def configuration(root,path):
    path=base.checked_path(root,path);cfg=base.read(path);name=cfg['campaign']
    if name not in SCHEMAS:raise ValueError('No3D or implicit arbitrary-cohort scoring is supported')
    count,methods=SCHEMAS[name]
    cases=tuple(f'r14_new_{i:04d}' for i in range(24)) if name=='independent24_v1' else tuple(f'r13_fresh_{i:04d}' for i in CONTROL_IDS)
    expected={(c,m) for c in cases for m in methods};actual={(t['case_id'],t['method']) for t in cfg['tasks']}
    if actual!=expected or len(cfg['tasks'])!=count*len(methods) or len({t['task_id'] for t in cfg['tasks']})!=len(expected):
        raise ValueError('Frozen complete expected cohort/method matrix required')
    for t in cfg['tasks']:
        if t['campaign']!=name or not t['task_id'].startswith(t['case_id']+'__') or any(c in t['task_id'] for c in '/\\'):
            raise ValueError('Task identity/campaign mismatch')
        base.checked_path(root,t['input_file'])
        if name=='self_only16_v1' and t['method'] in ('pair01','pair11') and t.get('onehop_self_only') is not True:
            raise ValueError('Self-only task lost frozen Goff-zero control')
        if name=='stage_control16_v1':
            wanted=['01','01','01','11'] if t['method']=='pair01_to11' else ['11','11','11','01']
            if t.get('metric_schedule')!=wanted:raise ValueError('Stage control schedule changed')
    return cfg,cases,methods


def task_inventory(root,task):
    folder=root/'RESULTS'/task['campaign']/task['task_id'];rp=folder/'RECEIPT.json';receipt=None;error=None
    if rp.is_file():
        try:receipt=base.read(rp)
        except Exception as exc:error=repr(exc)
    pool_new=folder/'CANDIDATE_POOL.npz';pool_old=folder/'CANDIDATES.npz'
    pool=pool_new if pool_new.is_file() else pool_old if pool_old.is_file() else pool_new
    return dict(task_id=task['task_id'],case_id=task['case_id'],method=task['method'],
        receipt_path=str(rp),receipt_file=base.stamp(rp),status=receipt.get('status') if receipt else None,
        receipt_error=error,predictions={v:dict(**base.stamp(folder/(v+'.npz')),
            metadata=base.stamp(folder/(v+'.json'))) for v in base.VIEWS},
        history_path=str(folder/'HISTORY.json'),history_file=base.stamp(folder/'HISTORY.json'),
        candidate_pool_path=str(pool),candidate_pool_file=base.stamp(pool),
        candidate_pool_alternatives=[base.stamp(pool_new),base.stamp(pool_old)],
        native_history_file=base.stamp(folder/'NATIVE_HISTORY.json'))


def make_seal(root,manifest,seal_path,source_note):
    cfg,cases,methods=configuration(root,manifest)
    queue=root/'queues'/cfg['campaign']/'RECEIPT.json';q=base.read(queue) if queue.is_file() else None
    selector=root/'configs/public_selection_independent24_v1.json'
    seal=dict(schema_version='r14_extension_prediction_seal_v1',root=str(root),campaign=cfg['campaign'],
        expected_cases=len(cases),expected_tasks=len(cfg['tasks']),manifest_path=str(Path(manifest).resolve()),
        manifest_snapshot=cfg,manifest_file=base.stamp(manifest),fixed_predictions=list(base.VIEWS),
        tasks=[task_inventory(root,t) for t in cfg['tasks']],authorized_scoring=False,
        queue_receipt=base.stamp(queue),queue_status=q.get('status') if q else None,
        created_at_utc=datetime.now(timezone.utc).isoformat(),source_note=source_note,
        public_selector_file=base.stamp(selector) if cfg['campaign']=='independent24_v1' else None,
        public_selector_snapshot=base.read(selector) if cfg['campaign']=='independent24_v1' else None,
        authorization_note='Parent alone approves after complete fixed-prediction accounting; failed/missing rows retained',
        truth_access='NONE',provenance='Existence/size only,no SHA or mtime checks; same-size replacement undetectable')
    if seal_path.is_relative_to(root/'RESULTS') or seal_path.is_relative_to(root/'queues'):raise ValueError('Seal outside predictions/queues required')
    seal_path.parent.mkdir(parents=True,exist_ok=True);base.write_json(seal_path,seal)
    return seal


def authorize(root,seal_path):
    seal=base.read(seal_path)
    if seal.get('schema_version')!='r14_extension_prediction_seal_v1' or seal.get('authorized_scoring') is not True:
        raise PermissionError('Separate extension seal must be parent-approved before truth access')
    if Path(seal['root']).resolve()!=root:raise ValueError('Seal root mismatch')
    cfg,cases,methods=configuration(root,seal['manifest_path'])
    if cfg!=seal['manifest_snapshot'] or cfg['campaign']!=seal['campaign'] or len(cfg['tasks'])!=seal['expected_tasks'] or len(cases)!=seal['expected_cases']:
        raise ValueError('Frozen extension manifest changed')
    if base.validate_stamp(root,seal['manifest_file'])!=Path(seal['manifest_path']).resolve():raise ValueError('Manifest path mismatch')
    if set(seal['fixed_predictions'])!=set(base.VIEWS):raise ValueError('All five views must be explicitly sealed')
    qpath=base.validate_stamp(root,seal['queue_receipt'])
    if qpath!=root/'queues'/cfg['campaign']/'RECEIPT.json':raise ValueError('Queue identity mismatch')
    if qpath.is_file():
        q=base.read(qpath)
        if q.get('status')!=seal['queue_status'] or q.get('status')=='RUNNING' or q.get('running') or q.get('pending_tasks',q.get('pending',0)):
            raise ValueError('Queue remains active or changed after sealing')
    if cfg['campaign']=='independent24_v1':
        path=base.validate_stamp(root,seal['public_selector_file'])
        if path!=root/'configs/public_selection_independent24_v1.json' or base.read(path)!=seal['public_selector_snapshot']:
            raise ValueError('Public selector changed after seal')
        validate_selector(seal['public_selector_snapshot'])
    account={}
    for row in seal['tasks']:
        ident=row['task_id'];folder=root/'RESULTS'/cfg['campaign']/ident
        if ident in account:raise ValueError('Duplicate sealed task')
        for field,name in [('receipt_file','RECEIPT.json'),('history_file','HISTORY.json'),('native_history_file','NATIVE_HISTORY.json')]:
            if base.validate_stamp(root,row[field])!=folder/name:raise ValueError('Task context mismatch')
        rp=Path(row['receipt_path']).resolve()
        if rp!=folder/'RECEIPT.json' or Path(row['history_path']).resolve()!=folder/'HISTORY.json':raise ValueError('Context alias mismatch')
        if rp.is_file() and row.get('receipt_error') is None:
            receipt=base.read(rp)
            if receipt.get('status')!=row['status'] or receipt.get('task',{}).get('task_id')!=ident or receipt.get('truth_access')!='NONE':
                raise ValueError('Receipt status/task/truth declaration changed')
        pool=base.validate_stamp(root,row['candidate_pool_file'])
        if pool.parent!=folder or pool.name not in ('CANDIDATE_POOL.npz','CANDIDATES.npz') or Path(row['candidate_pool_path']).resolve()!=pool:
            raise ValueError('Pool identity mismatch')
        for item in row['candidate_pool_alternatives']:
            if base.validate_stamp(root,item).parent!=folder:raise ValueError('Pool alternatives outside task')
        if set(row['predictions'])!=set(base.VIEWS):raise ValueError('Missing fixed view accounting')
        for view,item in row['predictions'].items():
            if base.validate_stamp(root,item)!=folder/(view+'.npz') or base.validate_stamp(root,item['metadata'])!=folder/(view+'.json'):
                raise ValueError('Fixed view identity mismatch')
        account[ident]=row
    if set(account)!={t['task_id'] for t in cfg['tasks']}:raise ValueError('Every expected task must be sealed, including failures')
    tasks=deepcopy(cfg['tasks'])
    for t in tasks:
        a=account[t['task_id']]
        if (a['case_id'],a['method'])!=(t['case_id'],t['method']):raise ValueError('Case/method identity mismatch')
        t['solver_label']=t['method']+('_self_only' if t.get('onehop_self_only') else '')
        t['pair']=t['solver_label'] # unique plot key,including stage/baseline controls
    return seal,cfg,cases,methods,tasks,account


def validate_selector(cfg):
    if cfg.get('order')!=['pair01','g19','dr'] or cfg.get('r14_BASE')!='pair01':raise ValueError('Unexpected public selector order/base')
    if cfg.get('stop_exact')!='if bool(incumbent.noise_feasible): break' or cfg.get('replace_exact')!='if float(trial.Phi) < float(incumbent.Phi): replace incumbent':
        raise ValueError('Unexpected frozen public selector rules')
    if cfg.get('child_initialization')!='Independent uniform-start stored results, never terminal warmrestart':raise ValueError('Selector initialization changed')


def costs_and_histories(root,cfg,tasks,caches,physical_counts,issues):
    queue=root/'queues'/cfg['campaign']/'RECEIPT.json'
    jobs={j['task_id']:j for j in base.read(queue).get('jobs',[])} if queue.is_file() else {}
    costs=[];history=[];native=[]
    for task in tasks:
        ident=task['task_id'];c=caches[ident];r=c['receipt'] or {};j=jobs.get(ident,{})
        ext=(j.get('external') or {}).get('child_process_CPU_seconds',j.get('child_process_CPU_seconds'))
        costs.append(dict(task_id=ident,case_id=task['case_id'],solver=task['solver_label'],run_status=c['item']['status'],
            external_child_CPU_seconds=ext,internal_final_process_CPU_seconds=r.get('CPU_seconds',r.get('whole_process_CPU_seconds')),
            external_wall_seconds=j.get('wall_seconds'),internal_wall_seconds=r.get('wall_seconds'),
            physical_cost=physical_counts(r.get('physical_cost',{})),auxiliary_cost=r.get('auxiliary_cost',r.get('auxiliary_receipt')),
            actual_corrector_policy=cfg['protocol'].get('corrector_policy'),core=j.get('cpu'),core_type=j.get('core_type'),
            all_spent_CPU_charged=True,no_pool_cutoff_cost_reduction_claim=True))
        for row in c['history']:history.append(dict(task_id=ident,solver=task['solver_label'],case_id=task['case_id'],**row))
        path=c['item']['native_history_file']
        if path['exists']:
            try:
                for row in base.read(path['path']):native.append(dict(task_id=ident,solver=task['solver_label'],case_id=task['case_id'],**row))
            except Exception as exc:issues.append(dict(task_id=ident,issue='NATIVE_HISTORY_ERROR',error=repr(exc)))
    return costs,history,native


def selector_public_plans(plans,costs,cases,cfg):
    """Exactly the frozen saved-uniform replay; no material metrics are inputs."""
    validate_selector(cfg)
    primary={(p['case_id'],p['solver']):p for p in plans if p['view']=='primary'}
    paid={(r['case_id'],r['solver']):r.get('external_child_CPU_seconds') for r in costs}
    output=[]
    for cid in cases:
        selected=None;invoked=[];trace=[];status='SELECTED';amounts=[]
        for solver in cfg['order']:
            if selected is not None and bool((selected.get('metadata') or {}).get('noise_feasible')):break
            trial=primary[(cid,solver)];meta=trial.get('metadata') or {};invoked.append(solver);amounts.append(paid.get((cid,solver)))
            if not trial['source_exists'] or not base.finite(meta.get('Phi')):
                status='UNAVAILABLE_INVOKED_CHILD';trace.append(dict(solver=solver,available=False));break
            replace=selected is None or float(meta['Phi'])<float(selected['metadata']['Phi'])
            if replace:selected=trial
            trace.append(dict(solver=solver,Phi=meta.get('Phi'),noise_feasible=meta.get('noise_feasible'),replace=replace,
                invoked_external_CPU_seconds=amounts[-1]))
        output.append(dict(case_id=cid,selector='frozen_public_saved_uniform',status=status,
            selected_solver=selected['solver'] if selected is not None else None,
            selected_task_id=selected['task_id'] if selected is not None else None,
            retained_incumbent_after_missing_invocation=selected is not None and status!='SELECTED',
            replay_complete=status=='SELECTED',
            invoked_solvers=invoked,invoked_full_child_CPU_seconds=sum(amounts) if all(base.finite(x) for x in amounts) else None,
            measured_invoked_CPU_subtotal=sum(x for x in amounts if base.finite(x)),
            missing_cost_solvers=[s for s,x in zip(invoked,amounts) if not base.finite(x)],trace=trace,
            scope='Frozen saved-result replay; independent uniform children; not measured conditional deployment,not a120s total-budget method'))
    return output


def summarize(rows,cases):
    result=base.summaries(rows)
    for row in result:
        row['expected_cases']=len(cases);row['recovery_denominator']=len(cases)
    return result


def pairwise(rows,cases,methods,np):
    lookup={(r['case_id'],r['solver']):r for r in rows if r['view']=='primary'}
    differences=[];memberships=[];boot=[]
    rng=np.random.default_rng(141515);draws=2000
    for method in methods:
        for ref in methods:
            if method==ref:continue
            gain=[];loss=[];both=[];neither=[];missing=[];metric_diffs={k:[] for k in ('NRMSE','IoU','joint_recovery')}
            for cid in cases:
                a,b=lookup[(cid,method)],lookup[(cid,ref)];available=base.accuracy_available(a) and base.accuracy_available(b)
                row=dict(case_id=cid,solver=method,reference=ref,paired_available=available)
                for key in base.METRICS:
                    av,bv=a.get(key),b.get(key);row['candidate_'+key]=av;row['reference_'+key]=bv
                    row['delta_'+key]=av-bv if available and base.finite(av) and base.finite(bv) else None
                differences.append(row)
                if not available:missing.append(cid);continue
                ar,br=base.recovered(a),base.recovered(b)
                (gain if ar and not br else loss if br and not ar else both if ar else neither).append(cid)
                metric_diffs['NRMSE'].append(a['NRMSE']-b['NRMSE']);metric_diffs['IoU'].append(a['IoU']-b['IoU'])
                metric_diffs['joint_recovery'].append(float(ar)-float(br))
            memberships.append(dict(solver=method,reference=ref,denominator=len(cases),gain_cases=gain,loss_cases=loss,
                both_cases=both,neither_cases=neither,unavailable_cases=missing,
                candidate_pass_count=sum(base.recovered(lookup[(c,method)]) for c in cases),
                reference_pass_count=sum(base.recovered(lookup[(c,ref)]) for c in cases),scope='Posthoc descriptive threshold,not selection'))
            for metric,values in metric_diffs.items():
                arr=np.asarray(values,float);low=high=mean=None
                if len(arr):
                    sampled=arr[rng.integers(0,len(arr),(draws,len(arr)))].mean(axis=1)
                    low,high=map(float,np.quantile(sampled,[.025,.975]));mean=float(arr.mean())
                boot.append(dict(solver=method,reference=ref,metric=metric,paired_available=len(arr),cohort_denominator=len(cases),
                    mean_difference=mean,descriptive_bootstrap95_low=low,descriptive_bootstrap95_high=high,
                    bootstrap_draws=draws,bootstrap_seed=141515,
                    scope='Paired object resampling of this fixed cohort; descriptive,unadjusted multiple comparisons,no population/novelty certificate'))
    return differences,memberships,boot


def corrector_effects(rows,np,issues):
    lookup={(r['case_id'],r['solver'],r['view']):r for r in rows if r['view'] in ('raw','corrected')}
    result=[]
    for cid,solver,view in list(lookup):
        if view!='raw' or solver=='direct':continue
        a,b=lookup[(cid,solver,'raw')],lookup[(cid,solver,'corrected')]
        row=dict(case_id=cid,solver=solver,raw_initial_fallback=a.get('raw_initial_fallback'),
            raw_source=a.get('source_file'),original_raw_source=a.get('raw_original_source'),
            corrected_available=b.get('available'),corrected_phase=b.get('phase'),
            actual_corrector_candidate=b.get('phase')=='corrector',material_RMS_change=None)
        for key in base.METRICS:
            av,bv=a.get(key),b.get(key);row['raw_'+key]=av;row['corrected_'+key]=bv
            row['delta_'+key]=bv-av if base.finite(av) and base.finite(bv) else None
        if a.get('available') and b.get('available'):
            try:
                with np.load(a['source_file'],allow_pickle=False) as za,np.load(b['source_file'],allow_pickle=False) as zb:
                    row['material_RMS_change']=float(np.sqrt(np.mean((zb['m']-za['m'])**2)))
            except Exception as exc:issues.append(dict(case_id=cid,solver=solver,issue='CORRECTOR_MATERIAL_MOVEMENT_ERROR',error=repr(exc)))
        result.append(row)
    return result


def figures(out,rows,cases,methods,pictures,truths,np):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    dest=out/'FIGURES';dest.mkdir();artifacts=[]
    panel_width=6 if len(methods)==11 else 5
    panels=[methods[i:i+panel_width] for i in range(0,len(methods),panel_width)]
    lookup={(r['case_id'],r['solver']):r for r in rows if r['view']=='primary'}
    with PdfPages(dest/'ALL_MATERIALS.pdf') as pdf:
        for offset in range(0,len(cases),4):
            subset=cases[offset:offset+4]
            scales={cid:max([float(np.max(m)) for m in ([truths[cid]['m']] if truths.get(cid) else [])+
                [pictures[(cid,s)] for s in methods if (cid,s) in pictures] if m.size]+[.1]) for cid in subset}
            for pi,panel in enumerate(panels):
                fig,axes=plt.subplots(len(subset),len(panel)+1,figsize=(16,2.6*len(subset)),squeeze=False)
                for ri,cid in enumerate(subset):
                    arrays=[truths[cid]['m'] if truths.get(cid) else None]+[pictures.get((cid,s)) for s in panel];im=None
                    for ci,m in enumerate(arrays):
                        ax=axes[ri,ci];ax.set_xticks([]);ax.set_yticks([])
                        if m is not None and m.size==576:im=ax.imshow(m.reshape(24,24),cmap='viridis',vmin=.02,vmax=scales[cid])
                        else:ax.text(.5,.5,'MISSING',ha='center',va='center',transform=ax.transAxes)
                        if ri==0:ax.set_title('GT' if ci==0 else panel[ci-1],fontsize=9)
                        if ci==0:ax.set_ylabel(cid,fontsize=8)
                        elif base.accuracy_available(lookup[(cid,panel[ci-1])]):
                            r=lookup[(cid,panel[ci-1])];ax.set_xlabel(f"NRMSE {r['NRMSE']:.3f} | IoU {r['IoU']:.3f}",fontsize=8)
                    if im is not None:fig.colorbar(im,ax=axes[ri,:].tolist(),fraction=.012,pad=.01)
                fig.suptitle('Frozen primary — same case IDs and same per-case scale across both method panels')
                fig.subplots_adjust(top=.93,bottom=.06,left=.08,right=.92,hspace=.28,wspace=.08)
                path=dest/f'materials_{offset//4+1:02d}_panel{pi+1}.png';fig.savefig(path,dpi=160);pdf.savefig(fig);plt.close(fig);artifacts.append(str(path.relative_to(out)))
    artifacts.append('FIGURES/ALL_MATERIALS.pdf')
    membership=np.array([[1 if base.recovered(lookup[(c,s)]) else 0 if base.accuracy_available(lookup[(c,s)]) else -1 for s in methods] for c in cases])
    fig,ax=plt.subplots(figsize=(13,max(6,len(cases)*.3)));ax.imshow(membership,cmap='RdYlGn',vmin=-1,vmax=1,aspect='auto')
    ax.set_xticks(range(len(methods)),methods,rotation=45,ha='right');ax.set_yticks(range(len(cases)),cases,fontsize=7)
    ax.set_title('Posthoc joint threshold membership: unavailable−1,fail0,pass1; not online selection');fig.tight_layout()
    path=dest/'JOINT_MEMBERSHIP.png';fig.savefig(path,dpi=160);plt.close(fig);artifacts.append(str(path.relative_to(out)))
    summary=summarize(rows,cases)
    for joint in (False,True):
        fig,axes=plt.subplots(1,1 if joint else 2,figsize=(13,5),squeeze=False)
        for ai,metric in enumerate(['recovered_posthoc_count'] if joint else ['NRMSE_mean','IoU_mean']):
            ax=axes[0,ai]
            for method in methods:
                points=sorted([r for r in summary if r['view']=='paid_prefix' and r['solver']==method],key=lambda r:r['CPU_prefix'])
                ax.plot([r['CPU_prefix'] for r in points],[r.get(metric) if r.get(metric) is not None else np.nan for r in points],marker='o',label=method)
            ax.set_xlabel('Native whole-process CPU cutoff (seconds)');ax.set_ylabel(metric);ax.grid(alpha=.2)
        axes[0,0].legend(fontsize=7);fig.suptitle('Posthoc joint counts (separate descriptive threshold)' if joint else 'Continuous material errors (no threshold in selection)');fig.tight_layout()
        path=dest/('POSTHOC_JOINT_PREFIXES.png' if joint else 'CONTINUOUS_PREFIXES.png');fig.savefig(path,dpi=160);plt.close(fig);artifacts.append(str(path.relative_to(out)))
    return artifacts


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--history-root',type=Path);parser.add_argument('--manifest',type=Path)
    parser.add_argument('--campaign',choices=list(SCHEMAS));parser.add_argument('--seal',type=Path,required=True)
    parser.add_argument('--out',type=Path);parser.add_argument('--make-seal',action='store_true');parser.add_argument('--source-note',default='')
    args=parser.parse_args();root=args.root.resolve();seal_path=args.seal.resolve()
    if args.make_seal:
        manifest=(args.manifest or (root/'configs'/(args.campaign+'.json') if args.campaign else None))
        if manifest is None:parser.error('--manifest or --campaign required to make a seal')
        seal=make_seal(root,manifest.resolve(),seal_path,args.source_note)
        print(json.dumps(dict(status='SEALED_NOT_AUTHORIZED',campaign=seal['campaign'],tasks=seal['expected_tasks'],seal=str(seal_path))));return 0
    if args.out is None or args.history_root is None:parser.error('--history-root and --out required for scoring')
    start_cpu,start_wall=time.process_time(),time.perf_counter();history=args.history_root.resolve();out=args.out.resolve()
    seal,cfg,cases,declared_methods,tasks,account=authorize(root,seal_path)
    if out.exists() or out==root or out.is_relative_to(history) or any(out.is_relative_to(root/n) for n in ('RESULTS','queues','sealed_truth','data_public','src','configs','audit')):
        raise ValueError('New output directory outside historical reports,predictions,truth and source required')
    out.mkdir(parents=True,exist_ok=False);issues=[]
    # Task order may vary; derive method order from the actual frozen method set.
    methods=tuple(m+('_self_only' if cfg['campaign']=='self_only16_v1' and m in ('pair01','pair11') else '') for m in declared_methods)
    if set(cfg.get('methods',declared_methods))!=set(declared_methods):issues.append(dict(issue='REGISTERED_METHODS_METADATA_OMITS_ACTUAL_FROZEN_TASK_METHOD',declared_metadata=cfg.get('methods'),actual_methods=list(declared_methods),tasks_unchanged=True))
    plans,caches=base.public_plans(root,tasks,account,issues,cases=cases,expected_methods=len(methods))
    base.write_json(out/'PUBLIC_SELECTIONS_BEFORE_TRUTH.json',plans)
    raw_public=[r for r in plans if r['view']=='raw']
    base.write_json(out/'PUBLIC_RAW_INITIAL_COUNTS_BEFORE_TRUTH.json',dict(expected_tasks=len(tasks),audited_rows=len(raw_public),
        comparable=sum(r.get('raw_audit_comparable') is True for r in raw_public),
        initial_fallbacks=sum(r.get('raw_initial_fallback') is True for r in raw_public),
        unavailable_or_uncomparable=sum(r.get('raw_audit_comparable') is not True for r in raw_public),truth_access='NONE'))
    # Costs and selector decisions below read only public receipts/history. No GT import yet.
    costs,histories,native=costs_and_histories(root,cfg,tasks,caches,lambda c:c,issues)
    selector=selector_public_plans(plans,costs,cases,seal['public_selector_snapshot']) if cfg['campaign']=='independent24_v1' else []
    base.write_json(out/'PUBLIC_SELECTOR_REPLAY_BEFORE_TRUTH.json',selector)
    import numpy as np
    helper_path=history/'src/score_results.py';sys.path.insert(0,str(history/'src'))
    module=importlib.util.spec_from_file_location('r14_extension_offline_metrics',helper_path)
    helper=importlib.util.module_from_spec(module);module.loader.exec_module(helper)
    truths={};labels={};truth_sources=[]
    if cfg['campaign']=='independent24_v1':
        # This private read is intentionally after approved seal and frozen public selectors.
        choices=base.read(root/'sealed_truth/independent24/GENERATOR_MANIFEST.json')
        labels={r['case_id']:dict(family=r['family'],regime=r['contrast_regime_name']) for r in choices['cases']}
    for cid in cases:
        spec=dict(path=str(root/'sealed_truth/independent24'/(cid+'.npz')),key='material_coarse',support_rule='fixed_chi_gt_1',amplitude_rule='max_truth_minus_host',host=.1,object_id=cid) if cfg['campaign']=='independent24_v1' else helper.truth_spec(history,dict(case_id=cid))
        truth_sources.append(dict(case_id=cid,**spec))
        try:truths[cid]=helper.load_truth(spec)
        except Exception as exc:truths[cid]=None;issues.append(dict(case_id=cid,issue='AUTHORIZED_TRUTH_UNAVAILABLE',error=repr(exc)))
    rows,pictures=base.score_plans(plans,caches,truths,helper.material_metrics,helper.physical_counts,np,issues)
    for row in rows:
        row['provenance']='NEW_R14_'+cfg['campaign'];row.update(labels.get(row['case_id'],{}))
    primary={(r['case_id'],r['solver']):r for r in rows if r['view']=='primary'}
    selector_rows=[]
    for decision in selector:
        chosen=primary.get((decision['case_id'],decision['selected_solver']))
        selector_rows.append(dict(chosen or dict(case_id=decision['case_id'],available=False,availability_reason='PUBLIC_SELECTOR_UNAVAILABLE'),
            solver='frozen_public_saved_uniform',**decision))
    differences,membership,bootstrap=pairwise(rows,cases,methods,np)
    effects=corrector_effects(rows,np,issues);summary=summarize(rows,cases)
    direct_rows=[dict(r,direct_scope='Whole120CPU physical material optimization;raw is initial,corrected is best direct search;not auxiliarycorrector effect') for r in rows if r['solver']=='direct']
    family=[]
    for key in ('family','regime'):
        for value in sorted({r.get(key) for r in rows if r.get(key) is not None}):
            subset=[r for r in rows if r.get(key)==value and r['view']=='primary'];subcases=sorted({r['case_id'] for r in subset})
            for r in summarize(subset,subcases):family.append(dict(group_by=key,group=value,**r))
    raw_audit=[r for r in rows if r['view']=='raw']
    for name,data in [('ALL_VIEWS',[r for r in rows if r['view'] in base.VIEWS]),('PAID_PREFIXES',[r for r in rows if r['view']=='paid_prefix']),
        ('COMMON_CPU_PREFIXES',[r for r in rows if r['view']=='common_CPU_prefix']),
        ('STRICT48_PAID_POOL_BRACKETS',[r for r in rows if r['view'] in ('strict48_lower','strict48_upper')]),
        ('COMMON_PHYSICS_PAID_POOL',[r for r in rows if r['view']=='common_physics_pool']),
        ('RAW_INITIAL_AUDIT',raw_audit),('PUBLIC_SELECTOR_REPLAY',selector_rows),('PAIRED_METRIC_DIFFERENCES',differences),
        ('POSTHOC_GAIN_LOSS',membership),('DESCRIPTIVE_PAIRED_BOOTSTRAP',bootstrap),('CORRECTOR_EFFECTS',effects),
        ('DIRECT_PHYSICAL_SEARCH_VIEWS',direct_rows),('COST_LEDGER',costs),('PAID_CANDIDATE_HISTORY',histories),('NATIVE_STAGE_HISTORY',native),
        ('SUMMARY',summary),('FAMILY_REGIME_AFTER_SEAL',family),('TRUTH_SOURCES',truth_sources),('ISSUES',issues)]:base.table(out,name,data)
    if selector_rows:
        base.write_json(out/'PUBLIC_SELECTOR_COUNTS.json',dict(expected_cases=len(cases),rows=len(selector_rows),
            available=sum(r.get('available') is True for r in selector_rows),
            replay_complete=sum(r.get('replay_complete') is True for r in selector_rows),
            missing_fallback_retained=sum(r.get('retained_incumbent_after_missing_invocation') is True for r in selector_rows),
            unavailable=sum(r.get('available') is not True for r in selector_rows),
            recovered_posthoc=sum(base.recovered(r) for r in selector_rows),denominator=len(cases),
            rule='Missing invocation prevents complete policy replay;retain any prior validincumbent,never GT-select an alternative;all cases retained'))
    artifacts=figures(out,rows,cases,methods,pictures,truths,np)
    measured=[r['external_child_CPU_seconds'] for r in costs]
    receipt=dict(status='AUTHORIZED_EXTENSION_OFFLINE_SCORING_COMPLETE',campaign=cfg['campaign'],expected_cases=len(cases),
        expected_tasks=len(tasks),method_labels=list(methods),task_status_counts=dict(Counter(r['run_status'] or 'MISSING' for r in costs)),
        expected_view_rows=len(tasks)*5,actual_view_rows=sum(r['view'] in base.VIEWS for r in rows),
        expected_prefix_rows=len(tasks)*5,actual_prefix_rows=sum(r['view']=='paid_prefix' for r in rows),
        raw_initial_fallback_count=sum(r.get('raw_initial_fallback') is True for r in raw_audit),
        raw_comparable_count=sum(r.get('raw_audit_comparable') is True for r in raw_audit),
        raw_unavailable_or_uncomparable=sum(r.get('raw_audit_comparable') is not True for r in raw_audit),
        source_seal=str(seal_path),protocol_snapshot=cfg['protocol'],issues=len(issues),
        direct_scope='Direct physical L-BFGS gets whole120CPU without proposal/corrector hierarchy;strict48bracket NA;not included in auxiliarycorrector effects',
        budget='All independent24 methods actual absolute120 total;unused proposal reassigned;strict48 outputs are paid-pool cutoff brackets,not actual stopped runs;full spentCPU charged',
        equal_physics='Componentwise same-case minima of all actual final paid original Maxwell LU/RHS;select min eligiblePhi within both cumulative limits;auxiliaryFFT/block/SVD/CPU separate and charged;unknown counters null,not0',
        clocks='Common native final CPU_seconds/whole_process_CPU_seconds across all methods per case;external lifetime CPU separate,not exact per-candidate external clock',
        external_child_CPU_total=sum(measured) if all(base.finite(x) for x in measured) else None,
        measured_external_CPU_subtotal=sum(x for x in measured if base.finite(x)),
        missing_external_CPU_tasks=[r['task_id'] for r in costs if not base.finite(r['external_child_CPU_seconds'])],
        selector='Exact frozen pair01→g19→dr noise-feasible stop/Phi replacement;sum all invoked fullchildCPU,not120deployment',
        threshold='NRMSE<=.15 andIoU>=.70 descriptive after seal only;continuous errors separate',
        noise='Raw CN99 noise only;crossgrid not model-error qualified',
        source_scope=cfg['scope'],family_access='Independent24 private labels read only after authorization/public-selector freezing',
        new_forward=0,new_RHS=0,new_optimizer=0,new_NN=0,scientific_review='PARENT_REQUIRED',
        offline_CPU_seconds=time.process_time()-start_cpu,offline_wall_seconds=time.perf_counter()-start_wall,
        timer_scope='Through metric/figure generation;receipt/report write andexitexcluded',artifacts=artifacts)
    base.write_json(out/'RECEIPT.json',receipt)
    body=['<!doctype html><meta charset="utf-8"><title>R14 extension offline evidence</title>',
        '<style>body{font:16px system-ui;max-width:1500px;margin:auto;padding:30px}img{max-width:100%}td,th{padding:6px;border:1px solid #ccc}table{border-collapse:collapse}</style>',
        '<h1>'+html.escape(cfg['campaign'])+'</h1>',
        '<p>Separate parent-approved cohort; every missing/failed task retained. Continuous errors and posthoc joint threshold are separate. No forward/optimizer/NN is run.</p>',
        '<p>Actual2Dcorrector budget is total absolute120 with unused proposal reallocated. Strict48 pool cutoffs are lower/upper robustness brackets, not measured stopped deployments. All actual costs remain charged.</p>',
        '<p>Public selector replay sums every invoked independent full child cost; it is neither warm restart nor a120second deployment. Cross-grid discrepancy is not a noise certificate.</p>',
        '<table><tr><th>Method</th><th>Available</th><th>Posthoc pass</th><th>Denominator</th><th>Mean NRMSE</th><th>Mean IoU</th></tr>']
    for r in summary:
        if r['view']=='primary':body.append('<tr>'+''.join('<td>'+html.escape(str(v))+'</td>' for v in [r['solver'],r['available_cases'],r['recovered_posthoc_count'],len(cases),r.get('NRMSE_mean'),r.get('IoU_mean')])+'</tr>')
    body.append('</table>')
    for path in sorted(out.glob('*.csv')):body.append('<p><a href="'+path.name+'">'+path.name+'</a></p>')
    for path in artifacts:
        body.append('<p><a href="'+path+'">'+path+'</a></p>')
        if path.endswith('.png'):body.append('<img src="'+path+'">')
    with (out/'report.html').open('x',encoding='utf-8') as f:f.write('\n'.join(body))
    print(json.dumps(dict(status=receipt['status'],campaign=cfg['campaign'],out=str(out),issues=len(issues))));return 0


if __name__=='__main__':raise SystemExit(main())
