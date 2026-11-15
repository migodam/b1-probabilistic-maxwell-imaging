"""Offline observable19 evidence. --make-seal never authorizes scoring or reads GT.

--root is supplement_observable19; the parent must separately set authorized=true
in a copy of its terminal prediction seal before invoking --seal. No forward solve.
"""
from __future__ import annotations
import argparse,csv,html,importlib.util,json,sys,time
from datetime import datetime,timezone
from pathlib import Path
sys.dont_write_bytecode=True
R14=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(R14/'src'))
import score_dualmetric as base
import score_extensions as extension

CAMPAIGN='observable19_exposed24_v1'
METHODS=('g19_internal_galerkin','g19_observation_galerkin','g19_anchor_observation_galerkin',
         'g19_internal_right','g19_observation_right','g19_anchor_observation_right')
LABELS=('G / PGP','S / PGP','O0 / PGP','G / GP','S / GP','O0 / GP')
CASES=tuple(f'r14_new_{i:04d}' for i in range(24))
SCHEMA='observable19_prediction_seal_v1'


def configuration(root,path):
    path=base.checked_path(root,path);cfg=base.read(path)
    expected={(c,m) for c in CASES for m in METHODS};tasks=cfg['tasks']
    if cfg['campaign']!=CAMPAIGN or len(tasks)!=144 or {(t['case_id'],t['method']) for t in tasks}!=expected or len({t['task_id'] for t in tasks})!=144:
        raise ValueError('Exactly frozen24 x six arms required, including missing/failed tasks')
    if set(cfg['methods'])!=set(METHODS) or cfg.get('basis_rank')!=19 or cfg.get('exposure_status')!='EXPOSED24_FOLLOWUP':raise ValueError('Frozen method/rank/exposure contract changed')
    if float(cfg['protocol']['CPU_cap'])!=120 or float(cfg['protocol']['proposal_fraction'])!=.6:raise ValueError('Frozen120/72 contract required')
    for t in tasks:
        if t['campaign']!=CAMPAIGN or t['task_id']!=t['case_id']+'__'+t['method'] or t.get('algorithm')!='g19':raise ValueError('Unexpected task identity/dispatch')
        p=(root/t['input_file']).resolve()
        if p!=root.parent/'data_public/independent24'/(t['case_id']+'.npz'):raise ValueError('Unexpected public observation path')
    return cfg


def queue_paths(root,explicit):
    paths=[p.resolve() for p in explicit] if explicit else sorted((root/'queues'/CAMPAIGN).rglob('RECEIPT.json'))
    if not paths:raise ValueError('Terminal queue receipt required; supply --queue for an explicit attempt')
    for p in paths:
        base.checked_path(root,p)
        if not p.is_relative_to(root/'queues'/CAMPAIGN):raise ValueError('Queue outside campaign')
    return paths


def terminal_queues(root,paths,expected):
    snapshots=[];jobs={}
    for p in paths:
        q=base.read(p)
        if q.get('status') in (None,'RUNNING') or q.get('running') or q.get('pending_tasks',q.get('pending',0)) or q.get('pending_cases',0):raise ValueError('Queue missing terminal closure: '+str(p))
        for job in q.get('jobs',[]):
            if job.get('task_id') not in expected:raise ValueError('Unexpected task in queue')
            jobs[job['task_id']]=dict(job,queue_receipt_path=str(p))
        snapshots.append(dict(file=base.stamp(p),snapshot=q))
    if set(jobs)!=expected:raise ValueError('Terminal queue must account for all144 task IDs, failed launches included')
    return snapshots,jobs


def inventory(root,t):
    item=extension.task_inventory(root,t);folder=root/'RESULTS'/CAMPAIGN/t['task_id']
    # This cohort uses original C2 CANDIDATES, never the new runner pool schema.
    item['candidate_pool_path']=str(folder/'CANDIDATES.npz');item['candidate_pool_file']=base.stamp(folder/'CANDIDATES.npz')
    expected=[folder/(v+s) for v in base.VIEWS for s in ('.json','.npz')]
    expected += [folder/n for n in ('RECEIPT.json','ACTUAL_CONFIG.json','HISTORY.json','CANDIDATES.npz','DENOMINATOR_SETUP.json','DENOMINATOR_SETUP_FAILED.json','WRAPPER_PROVENANCE.json')]
    files=set(expected)|{p for p in folder.rglob('*') if p.is_file()}
    item['all_result_files']=[base.stamp(p) for p in sorted(files)]
    return item


def make_seal(root,manifest,path,queues,note):
    cfg=configuration(root,manifest);snapshots,jobs=terminal_queues(root,queue_paths(root,queues),{t['task_id'] for t in cfg['tasks']})
    seal=dict(schema_version=SCHEMA,root=str(root),campaign=CAMPAIGN,authorized=False,authorized_scoring=False,
        authorization_note='Parent must separately save authorized=true after terminal accounting; this CLI never authorizes.',
        manifest_file=base.stamp(manifest),manifest_snapshot=cfg,expected_tasks=144,expected_cases=24,
        fixed_predictions=list(base.VIEWS),queue_receipts=snapshots,tasks=[inventory(root,t) for t in cfg['tasks']],
        created_at_utc=datetime.now(timezone.utc).isoformat(),source_note=note,truth_access='NONE',
        provenance='Local paths/existence/bytes and public JSON only; no hash/mtime; same-size replacement is undetectable')
    path=base.checked_path(root,path)
    if any(path.is_relative_to(root/n) for n in ('RESULTS','queues','src','configs')):raise ValueError('Seal must be a new audit/predictions path')
    path.parent.mkdir(parents=True,exist_ok=True);base.write_json(path,seal)
    return seal


def authorize(root,path):
    seal=base.read(path)
    if seal.get('schema_version')!=SCHEMA or seal.get('authorized') is not True:raise PermissionError('Parent-authorized observable19 seal (authorized=true) required before any truth access')
    if Path(seal['root']).resolve()!=root or seal.get('expected_tasks')!=144 or seal.get('expected_cases')!=24 or set(seal.get('fixed_predictions',[]))!=set(base.VIEWS):raise ValueError('Seal identity/complete fixed-view accounting mismatch')
    manifest=base.validate_stamp(root,seal['manifest_file']);cfg=configuration(root,manifest)
    if cfg!=seal['manifest_snapshot']:raise ValueError('Manifest changed')
    paths=[]
    for q in seal['queue_receipts']:
        p=base.validate_stamp(root,q['file'])
        if base.read(p)!=q['snapshot']:raise ValueError('Terminal queue content changed')
        paths.append(p)
    _,jobs=terminal_queues(root,paths,{t['task_id'] for t in cfg['tasks']})
    tasks=[];account={};original={t['task_id']:t for t in cfg['tasks']}
    for item in seal['tasks']:
        ident=item['task_id']
        if ident not in original or ident in account:raise ValueError('Duplicate/unexpected sealed task')
        task=original[ident];folder=root/'RESULTS'/CAMPAIGN/ident
        if (item['case_id'],item['method'])!=(task['case_id'],task['method']):raise ValueError('Task case/method changed')
        current=inventory(root,task)
        if current!=item:raise ValueError('Sealed result inventory/status/bytes changed: '+ident)
        if item['receipt_file']['exists'] and item.get('receipt_error') is None:
            receipt=base.read(item['receipt_path'])
            if receipt.get('truth_access')!='NONE' or receipt.get('task',{}).get('task_id')!=ident:raise ValueError('Unexpected online truth declaration/task identity')
        if Path(item['candidate_pool_path'])!=folder/'CANDIDATES.npz':raise ValueError('Wrong C2 candidate-pool schema')
        account[ident]=item;tasks.append(dict(task,solver_label=task['method'],pair=task['method']))
    if set(account)!=set(original):raise ValueError('Every task including failures/missing must remain sealed')
    return seal,cfg,tasks,account,jobs


def audit_views(plans,caches,issues):
    # The shared planner implements paid CPU/physics selectors; corrected fallback
    # is original C2 phase_best.get('corrector',best), including zero-corrector cases.
    audits=[]
    for p in plans:
        if p['selection_kind']!='saved_view':continue
        c=caches[p['task_id']];hist=c['history'];valid=[r for r in hist if base.eligible(r)];view=p['view'];meta=p.get('metadata') or {}
        if view=='corrected':
            expected=base.pick([r for r in valid if r.get('phase')=='corrector'],120.)
            fallback=expected is None
            if fallback:expected=base.pick([r for r in valid if r.get('phase') in ('initial','proposal')],120.)
        elif view=='raw':expected=base.pick([r for r in valid if r.get('phase') in ('initial','proposal')],120.);fallback=False
        elif p.get('source_view')=='primary':expected=base.pick(hist,120.);fallback=False
        elif view=='initial':expected=next((r for r in hist if r.get('phase')=='initial'),None);fallback=False
        else:
            feasible=[r for r in valid if r.get('noise_feasible') is True and base.finite(r.get('R'))]
            expected=min(feasible,key=lambda r:(r['R'],r['Phi'],r['candidate'])) if feasible else None;fallback=False
        idx=expected.get('candidate') if expected else None
        consistent=bool(p['source_exists'] and expected is not None and meta.get('candidate')==idx)
        discrepancies=[k for k in ('Phi','D','R','CPU_seconds','phase','eligible','legal') if expected is not None and meta.get(k)!=expected.get(k)]
        consistent=consistent and not discrepancies
        p['expected_public_candidate']=idx;p['archived_selector_consistent']=consistent
        row=dict(task_id=p['task_id'],case_id=p['case_id'],solver=p['solver'],view=view,source_view=p.get('source_view'),
            archived_candidate=meta.get('candidate'),expected_candidate=idx,consistent=consistent,metadata_discrepancies=discrepancies,
            archived_phase=meta.get('phase'),corrected_uses_pre_corrector_fallback=view=='corrected' and fallback,
            raw_proposal_source_includes_initial=view=='raw',raw_original_source=p.get('raw_original_source'),
            unavailable_reason=None if p['source_exists'] else 'MISSING_SAVED_VIEW_OR_DEPLOYMENT_FALLBACK',
            array_consistency='checked separately against CANDIDATES.m by shared score_plans after seal')
        audits.append(row)
        if p['source_exists'] and not consistent:issues.append(dict(row,issue='ARCHIVED_SELECTOR_MISMATCH'))
    return audits


def cost_rows(root,cfg,tasks,caches,jobs):
    rows=[]
    for t in tasks:
        ident=t['task_id'];folder=root/'RESULTS'/CAMPAIGN/ident;r=caches[ident]['receipt'] or {};job=jobs.get(ident,{})
        external=job.get('external') or {};cpu=external.get('child_process_CPU_seconds',job.get('child_process_CPU_seconds'))
        wrapper=base.read(folder/'WRAPPER_PROVENANCE.json') if (folder/'WRAPPER_PROVENANCE.json').is_file() else {}
        setup_path=next((folder/n for n in ('DENOMINATOR_SETUP_FAILED.json','DENOMINATOR_SETUP.json') if (folder/n).is_file()),None)
        setup=base.read(setup_path) if setup_path else wrapper.get('actual_setup',{})
        aux=r.get('auxiliary_receipt') or {};actual=setup.get('setup_operation_counts');nominal=aux.get('denominator_setup_cost',aux.get('setup_cost'))
        if nominal is None:nominal={k:v for k,v in aux.items() if 'setup' in k}
        rows.append(dict(task_id=ident,case_id=t['case_id'],solver=t['method'],run_status=caches[ident]['item']['status'],queue_status=job.get('status'),
            returncode=job.get('returncode'),queue_receipt=job.get('queue_receipt_path'),external_child_CPU_seconds=cpu,
            external_CPU_null_reason=None if base.finite(cpu) else 'NO_EXTERNAL_WHOLE_CHILD_CPU_MEASUREMENT',external_CPU_scope=external.get('CPU_scope'),
            internal_final_process_CPU_seconds=r.get('whole_process_CPU_seconds',r.get('CPU_seconds')),
            internal_CPU_null_reason=None if base.finite(r.get('whole_process_CPU_seconds',r.get('CPU_seconds'))) else 'MISSING_CHILD_FINAL_CLOCK',
            external_wall_seconds=job.get('wall_seconds'),internal_wall_seconds=r.get('wall_seconds'),cpu=job.get('cpu'),core_type=job.get('core_type'),
            physical_cost_raw=r.get('physical_cost'),physical_cost=base.cumulative_physics(r.get('physical_cost')),
            original_R13_nominal_setup=nominal or None,original_R13_nominal_setup_null_reason=None if nominal else 'NO_NOMINAL_SETUP_FIELD',
            actual_factory_setup=setup or None,actual_factory_setup_counters=actual,
            actual_factory_setup_null_reason=None if actual is not None else 'NO_ACTUAL_FACTORY_COUNTER_RECEIPT',actual_factory_setup_source=str(setup_path) if setup_path else str(folder/'WRAPPER_PROVENANCE.json'),
            wrapper_actual_setup=wrapper.get('actual_setup'),auxiliary_receipt=aux or None,projector_cost=r.get('projector_cost'),
            O0_physical_delta_already_in_physical_aggregate=setup.get('anchor_physical_cost_delta'),
            setup_RHS_accounting='O0 receiver adjoint RHS already in original physical aggregate; never added twice. Dense Green/basis RHS are auxiliary actions, not physical solves.',
            all_spent_CPU_charged=True,failed_attempts_retained=True,wrapper_CPU_not_added_to_child=True))
    return rows


def historical(root):
    p=root.parent/'analysis/independent24_v1/ALL_VIEWS.csv';data=[]
    original=list(csv.DictReader(p.open())) if p.is_file() else []
    for cid in CASES:
        for method in ('cie','dual_cie','pair01','pair11','g19','dr'):
            found=next((x for x in original if x['case_id']==cid and x['solver']==method and x['view']=='primary'),None)
            row=dict(found or {},case_id=cid,solver=method,view='primary',source_table=str(p),provenance='FROZEN_HISTORICAL_R14_EXPOSED_REFERENCE',causal_equal_cost_comparison=False,reference_available=found is not None)
            if found is None:row['availability_reason']='MISSING_HISTORICAL_REFERENCE_ROW'
            data.append(row)
    return data


def figures(out,rows,pictures,truths,np,costs,effects):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    paths=[];folder=out/'figures';folder.mkdir()
    # Four cases/page: seven material columns + a dedicated colorbar column.
    for start in range(0,24,4):
        group=CASES[start:start+4];fig=plt.figure(figsize=(20,3.1*len(group)));gs=fig.add_gridspec(len(group),8,width_ratios=[1]*7+[.055],hspace=.35,wspace=.12)
        for ri,cid in enumerate(group):
            arrays=[truths[cid]['m'] if truths.get(cid) else None]+[pictures.get((cid,m)) for m in METHODS]
            valid=[a for a in arrays if a is not None];vmax=max([20e-3]+[float(np.max(a)) for a in valid]);im=None
            for col,(label,array) in enumerate(zip(('GT',)+LABELS,arrays)):
                ax=fig.add_subplot(gs[ri,col]);ax.set_xticks([]);ax.set_yticks([])
                if ri==0:ax.set_title(label,fontsize=10)
                if col==0:ax.set_ylabel(cid,fontsize=9)
                if array is None:ax.text(.5,.5,'MISSING',ha='center',va='center',transform=ax.transAxes);continue
                im=ax.imshow(array.reshape(24,24),origin='lower',vmin=.02,vmax=vmax,cmap='viridis')
            cax=fig.add_subplot(gs[ri,7])
            if im is not None:fig.colorbar(im,cax=cax)
            else:cax.axis('off')
        fig.suptitle('Exposed24: archived primary; shared scale within each case',y=.995)
        path=folder/f'materials_{start:02d}_{start+len(group)-1:02d}.png';fig.savefig(path,dpi=150,bbox_inches='tight');plt.close(fig);paths.append(str(path.relative_to(out)))
    lookup={(r['case_id'],r['solver']):r for r in rows if r['view']=='primary'}
    membership=np.full((24,6),np.nan)
    for i,cid in enumerate(CASES):
        for j,m in enumerate(METHODS):
            r=lookup[cid,m]
            if base.accuracy_available(r):membership[i,j]=float(base.recovered(r))
    fig,ax=plt.subplots(figsize=(10,10));im=ax.imshow(np.ma.masked_invalid(membership),vmin=0,vmax=1,cmap='Greys',aspect='auto');ax.set_xticks(range(6),LABELS,rotation=30,ha='right');ax.set_yticks(range(24),CASES);ax.set_title('Descriptive joint recovery; blank = unavailable, denominator24');fig.colorbar(im,ax=ax,ticks=[0,1]);fig.tight_layout();p=folder/'recovery_membership.png';fig.savefig(p,dpi=140);plt.close(fig);paths.append(str(p.relative_to(out)))
    fig,axes=plt.subplots(2,3,figsize=(17,9))
    for ax,metric in zip(axes[0],('NRMSE','IoU','Phi')):
        for method,label in zip(METHODS,LABELS):
            vals=[lookup[c,method].get(metric) for c in CASES];ax.plot(range(24),[v if base.finite(v) else np.nan for v in vals],'.-',label=label)
        ax.set_title('Primary '+metric);ax.set_xlabel('Fixed case index')
    for ax,metric in zip(axes[1],('NRMSE','IoU','Phi')):
        for method,label in zip(METHODS,LABELS):
            means=[]
            for cut in base.PREFIXES:
                vs=[r.get(metric) for r in rows if r['solver']==method and r['view']=='paid_prefix' and r['CPU_prefix']==cut and base.finite(r.get(metric))];means.append(sum(vs)/len(vs) if vs else np.nan)
            ax.plot(base.PREFIXES,means,'.-',label=label)
        ax.set_title('Saved paid-prefix '+metric+' (available means)');ax.set_xlabel('Native process CPU s; full costs retained')
    axes[0,0].legend(fontsize=8);fig.tight_layout();p=folder/'continuous_and_cost_curves.png';fig.savefig(p,dpi=150);plt.close(fig);paths.append(str(p.relative_to(out)))
    fig,axes=plt.subplots(1,2,figsize=(14,10))
    for ax,metric in zip(axes,('NRMSE','IoU')):
        mat=np.full((24,6),np.nan)
        for i,cid in enumerate(CASES):
            ref=lookup[cid,METHODS[0]].get(metric)
            for j,m in enumerate(METHODS):
                v=lookup[cid,m].get(metric)
                if base.finite(v) and base.finite(ref):mat[i,j]=v-ref
        bound=max(1e-12,float(np.nanmax(abs(mat)))) if np.isfinite(mat).any() else 1
        im=ax.imshow(np.ma.masked_invalid(mat),cmap='coolwarm',vmin=-bound,vmax=bound,aspect='auto');ax.set_xticks(range(6),LABELS,rotation=40,ha='right');ax.set_yticks(range(24),CASES);ax.set_title(metric+' minus G/PGP');fig.colorbar(im,ax=ax)
    fig.tight_layout();p=folder/'paired_metric_differences.png';fig.savefig(p,dpi=140);plt.close(fig);paths.append(str(p.relative_to(out)))
    fig,axes=plt.subplots(2,3,figsize=(17,9))
    for ax,metric in zip(axes.ravel()[:5],('NRMSE','IoU','Phi','D','R')):
        for j,(method,label) in enumerate(zip(METHODS,LABELS)):
            values=[e.get('delta_'+metric) for e in effects if e['solver']==method]
            values=[v for v in values if base.finite(v)]
            ax.scatter([j]*len(values),values,s=13,alpha=.7)
        ax.axhline(0,color='grey',linewidth=.6);ax.set_title('Corrected - raw '+metric);ax.set_xticks(range(6),LABELS,rotation=40,ha='right')
    ax=axes.ravel()[5]
    for j,(method,label) in enumerate(zip(METHODS,LABELS)):
        values=[c['external_child_CPU_seconds'] for c in costs if c['solver']==method and base.finite(c.get('external_child_CPU_seconds'))]
        ax.scatter([j]*len(values),values,s=13,alpha=.7)
    ax.set_title('Full external child CPU; missing excluded explicitly');ax.set_ylabel('Whole lifetime OS CPU seconds');ax.set_xticks(range(6),LABELS,rotation=40,ha='right')
    fig.tight_layout();p=folder/'corrector_effects_and_actual_costs.png';fig.savefig(p,dpi=140);plt.close(fig);paths.append(str(p.relative_to(out)))
    return paths


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]);p.add_argument('--manifest',type=Path);p.add_argument('--make-seal',type=Path);p.add_argument('--seal',type=Path);p.add_argument('--queue',type=Path,action='append',default=[]);p.add_argument('--out',type=Path);p.add_argument('--source-note',default='');a=p.parse_args()
    start=time.perf_counter();cpu=time.process_time();root=a.root.resolve();manifest=(a.manifest or root/'configs'/f'{CAMPAIGN}.json').resolve()
    if bool(a.make_seal)==bool(a.seal):p.error('Choose exactly one of --make-seal PATH or --seal APPROVED_PATH')
    if a.make_seal:
        make_seal(root,manifest,a.make_seal.resolve(),a.queue,a.source_note)
        receipt=dict(status='UNAUTHORIZED_SEAL_CREATED',truth_access='NONE',CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-start,seal=str(a.make_seal.resolve()))
        base.write_json(a.make_seal.resolve().with_suffix('.receipt.json'),receipt);print(json.dumps(receipt));return
    if not a.out:p.error('--out NEW_DIRECTORY is required for authorized scoring')
    out=a.out.resolve()
    if out.exists() or not out.is_relative_to(root/'analysis'):raise ValueError('Scoring output must be a new supplement analysis directory')
    seal,cfg,tasks,account,jobs=authorize(root,a.seal.resolve())
    out.mkdir(parents=True);issues=[]
    try:
        plans,caches=base.public_plans(root,tasks,account,issues,cases=CASES,expected_methods=6)
        plans=[p for p in plans if not p['view'].startswith('strict48_')]
        issues=[i for i in issues if i.get('issue')!='ARCHIVED_SELECTOR_METADATA_MISMATCH']
        audits=audit_views(plans,caches,issues)
        # Freeze all selectors and comparisons before importing/loading any private truth.
        base.write_json(out/'PUBLIC_PLANS_BEFORE_TRUTH.json',plans);base.table(out,'SELECTOR_AUDIT',audits)
        import numpy as np
        hp=root.parent/'runtime/scoring_history/src/score_results.py'
        sys.path.insert(0,str(hp.parent))
        spec=importlib.util.spec_from_file_location('observable19_original_metrics',hp);helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)
        truths={};sources=[]
        for cid in CASES:
            ts=dict(path=str(root.parent/'sealed_truth/independent24'/(cid+'.npz')),key='material_coarse',support_rule='fixed_chi_gt_1',amplitude_rule='max_truth_minus_host',host=.1,object_id=cid);sources.append(dict(case_id=cid,**ts))
            try:truths[cid]=helper.load_truth(ts)
            except Exception as exc:truths[cid]=None;issues.append(dict(case_id=cid,issue='AUTHORIZED_TRUTH_UNAVAILABLE',error=repr(exc)))
        rows,pictures=base.score_plans(plans,caches,truths,helper.material_metrics,helper.physical_counts,np,issues)
        for row in rows:
            row.update(provenance='NEW_EXPOSED24_OBSERVABLE19_FOLLOWUP',cohort_denominator=24,CPU_comparison_scope='Saved native whole-process clock; external lifetime OS CPU separate; no exact external per-candidate timestamp or stopped-rerun claim')
            if not row.get('available'):
                for metric in base.METRICS:row.setdefault(metric,None)
        costs=cost_rows(root,cfg,tasks,caches,jobs);summary=extension.summarize(rows,CASES)
        queue_costs=[];all_attempts=[]
        for q in seal['queue_receipts']:
            qr=q['snapshot'];qp=q['file']['path']
            queue_costs.append(dict(path=qp,status=qr.get('status'),queue_wall_seconds=qr.get('queue_wall_seconds'),controller_CPU_seconds=qr.get('controller_CPU_seconds'),scope='Individual queue interval; overlapping/parallel wall intervals never summed'))
            for j in qr.get('jobs',[]):
                e=j.get('external') or {};v=e.get('child_process_CPU_seconds',j.get('child_process_CPU_seconds'))
                all_attempts.append(dict(j,queue_receipt=qp,actual_OS_CPU_seconds=v,CPU_null_reason=None if base.finite(v) else 'NO_EXTERNAL_WHOLE_CHILD_CPU',scope='All queue attempts, including failed launches; repeated receipts may describe same PID, inspect before summing'))
        effects=extension.corrector_effects(rows,np,issues);diffs,members,boot=extension.pairwise(rows,CASES,METHODS,np)
        normalized=[];normalized_differences=[];normalized_boot=[]
        for view in ('common_CPU_prefix','common_physics_pool'):
            for cut in (base.PREFIXES if view=='common_CPU_prefix' else (None,)):
                group=[dict(r,view='primary') for r in rows if r['view']==view and r.get('CPU_prefix')==cut]
                d,m,b=extension.pairwise(group,CASES,METHODS,np)
                normalized.extend(dict(x,comparison_view=view,CPU_prefix=cut) for x in m)
                normalized_differences.extend(dict(x,comparison_view=view,CPU_prefix=cut) for x in d)
                normalized_boot.extend(dict(x,comparison_view=view,CPU_prefix=cut) for x in b)
        history=[dict(task_id=ident,case_id=c['task']['case_id'],solver=c['task']['method'],**h) for ident,c in caches.items() for h in c['history']]
        for name,data in [('ALL_VIEWS',rows),('SUMMARY',summary),('COST_LEDGER',costs),('QUEUE_COSTS',queue_costs),('ALL_QUEUE_ATTEMPTS',all_attempts),('PAID_HISTORY',history),('CORRECTOR_EFFECTS',effects),('PAIRED_DIFFERENCES',diffs),('RECOVERY_GAIN_LOSS',members),('PAIRED_BOOTSTRAP_DESCRIPTIVE',boot),('NORMALIZED_GAIN_LOSS',normalized),('NORMALIZED_PAIRED_DIFFERENCES',normalized_differences),('NORMALIZED_DESCRIPTIVE_BOOTSTRAP',normalized_boot),('HISTORICAL_REFERENCES',historical(root)),('TRUTH_SOURCES',sources)]:base.table(out,name,data)
        artifacts=figures(out,rows,pictures,truths,np,costs,effects);base.table(out,'ISSUES',issues)
        measured=[r['external_child_CPU_seconds'] for r in costs];total=sum(measured) if all(base.finite(v) for v in measured) else None
        receipt=dict(status='COMPLETE_WITH_RETAINED_ISSUES' if issues else 'COMPLETE',campaign=CAMPAIGN,seal=str(a.seal.resolve()),expected_tasks=144,expected_cases=24,methods=list(METHODS),issues=len(issues),truth_access='AFTER_PARENT_AUTHORIZED_SEAL',
            external_child_CPU_total=total,external_child_CPU_measured_subtotal=sum(v for v in measured if base.finite(v)),external_CPU_missing_tasks=[r['task_id'] for r in costs if not base.finite(r['external_child_CPU_seconds'])],
            CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-start,scoring_scope='Metrics/plots only; no optimizer/forward. Failed/missing objects retained in24 denominator.',
            model_scope='Exposed24 same public observations; all576 real materials; fixed19 G/S/O0 denominator bases, PGP/GP; scalar TM only',
            reference_scope='Original R14 references frozen historical; not new strict equal-cost arms',normalization='Saved completed min legal eligible originalPhi, nativeCPU prefixes and common componentwise physicalLU/RHS; full already-spent costs charged, all auxiliary/setup operations separate, no stopped rerun')
        base.write_json(out/'RECEIPT.json',receipt)
        body=['<!doctype html><meta charset="utf-8"><title>Observable19 exposed24</title><style>body{max-width:1600px;margin:auto;font-family:sans-serif}img{max-width:100%}td,th{padding:6px;border:1px solid #ddd}table{border-collapse:collapse}</style><h1>Observable19 exposed24: saved-prediction report</h1>',
            '<p>144 tasks; six arms; denominator24 per arm. Scientific interpretation reserved for parent. Historical references are separate. Recovery threshold is descriptive and never used for candidate selection.</p>',
            '<p>GT loaded only after authorized seal validation. Candidate prefixes retain full solver charges. Native process clocks differ from external OS lifetime CPU; O0 RHS are already in physical aggregate.</p>',
            '<pre>'+html.escape(json.dumps(receipt,indent=2))+'</pre>','<table><tr><th>Method</th><th>Available</th><th>Recovered/24</th><th>Mean NRMSE (available)</th><th>Mean IoU (available)</th></tr>']
        for r in summary:
            if r['view']=='primary':body.append('<tr>'+''.join('<td>'+html.escape(str(r.get(k)))+'</td>' for k in ('solver','accuracy_available_cases','recovered_posthoc_count','NRMSE_mean','IoU_mean'))+'</tr>')
        body.append('</table>');body += ['<p><a href="'+html.escape(p.name)+'">'+html.escape(p.name)+'</a></p>' for p in sorted(out.glob('*.csv'))];body += ['<img src="'+html.escape(p)+'" loading="lazy">' for p in artifacts]
        (out/'index.html').write_text('\n'.join(body),encoding='utf-8');print(json.dumps(receipt))
    except BaseException as exc:
        base.write_json(out/'FAILED_RECEIPT.json',dict(status='ANALYSIS_FAILED',error=repr(exc),CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-start,seal=str(a.seal.resolve()),no_new_physics=True));raise

if __name__=='__main__':main()
