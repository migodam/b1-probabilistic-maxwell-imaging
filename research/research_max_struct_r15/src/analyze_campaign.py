#!/usr/bin/env python3
"""R15 sealed exposed24 campaign analysis. No truth access at import/--make-seal.

First create an unauthorized seal; parent separately sets authorized=true and
 authorized_scoring=true after all264 predictions and terminal queue accounting.
Only --seal APPROVED_PATH can score. No physics, optimization or rank selection.
"""
from __future__ import annotations
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
import argparse,csv,html,importlib.util,json,sys,time
from datetime import datetime,timezone
from pathlib import Path
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1]
R14=ROOT.parent/'research_max_struct_r14'
sys.path.insert(0,str(R14/'src'))
import score_dualmetric as base
import score_extensions as extension
spec=importlib.util.spec_from_file_location('r15_inherited_view_audit',R14/'supplement_observable19/src/analyze_results.py')
view_audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(view_audit)
CAMPAIGN='twofold_exposed24_v1'
METHODS=('g19','g19_s19','twofold','twofold_hop','twofold_dr','twofold_ap','twofold_full_dr','twofold_reset','cie','dual_cie','dr')
CASES=tuple(f'r14_new_{i:04d}' for i in range(24))
SCHEMA='r15_twofold_exposed24_prediction_seal_v1'


def configuration(root,path):
    path=base.checked_path(root,path);cfg=base.read(path);tasks=cfg['tasks']
    wanted={(c,m) for c in CASES for m in METHODS}
    if cfg.get('campaign')!=CAMPAIGN or len(tasks)!=264 or len({t['task_id'] for t in tasks})!=264 or {(t['case_id'],t['method']) for t in tasks}!=wanted:raise ValueError('Exactly frozen24 x11 arms required, including failures')
    if set(cfg.get('methods',[]))!=set(METHODS) or cfg.get('exposure_status')!='EXPOSED24_FOLLOWUP':raise ValueError('Frozen methods/exposure changed')
    if cfg['protocol']['CPU_cap']!=120 or cfg['protocol']['proposal_fraction']!=.6:raise ValueError('Original120/72 protocol required')
    for t in tasks:
        if t['task_id']!=t['case_id']+'__'+t['method'] or t['campaign']!=CAMPAIGN:raise ValueError('Task identity changed')
        if base.checked_path(root,t['input_file'])!=root/'data_public/independent24'/(t['case_id']+'.npz'):raise ValueError('Unexpected public input')
    return cfg


def queues(root,paths,expected):
    paths=[base.checked_path(root,p) for p in paths] if paths else [root/'queues'/CAMPAIGN/'RECEIPT.json']
    snapshots=[];jobs={}
    for path in paths:
        if not path.is_relative_to(root/'queues'/CAMPAIGN):raise ValueError('Queue outside current campaign')
        q=base.read(path)
        if q.get('status') in (None,'RUNNING') or q.get('running') or q.get('pending_tasks',0) or q.get('pending_cases',0):raise ValueError('Terminal queue required')
        for j in q.get('jobs',[]):
            if j.get('task_id') not in expected:raise ValueError('Unknown queue task')
            jobs[j['task_id']]=dict(j,queue_receipt_path=str(path))
        snapshots.append(dict(file=base.stamp(path),snapshot=q))
    if set(jobs)!=expected:raise ValueError('Queue must account for all264 task IDs, including failed launches')
    return snapshots,jobs


def inventory(root,task):
    item=extension.task_inventory(root,task)
    folder=root/'RESULTS'/CAMPAIGN/task['task_id']
    item['candidate_pool_path']=str(folder/'CANDIDATES.npz');item['candidate_pool_file']=base.stamp(folder/'CANDIDATES.npz')
    item['all_result_files']=[base.stamp(p) for p in sorted(folder.rglob('*')) if p.is_file()] if folder.is_dir() else []
    return item


def make_seal(root,manifest,target,queue_paths,note):
    target=base.checked_path(root,target)
    if target.exists() or any(target.is_relative_to(root/n) for n in ('RESULTS','queues','configs','src')):raise ValueError('New audit/predictions seal path required')
    cfg=configuration(root,manifest);snapshots,jobs=queues(root,queue_paths,{t['task_id'] for t in cfg['tasks']})
    seal=dict(schema_version=SCHEMA,root=str(root),campaign=CAMPAIGN,authorized=False,authorized_scoring=False,
        manifest_file=base.stamp(manifest),manifest_snapshot=cfg,queue_receipts=snapshots,
        expected_tasks=264,expected_cases=24,fixed_predictions=list(base.VIEWS),tasks=[inventory(root,t) for t in cfg['tasks']],
        created_utc=datetime.now(timezone.utc).isoformat(),truth_access='NONE',source_note=note,
        authorization_note='Parent alone creates a separately approved copy with authorized=true AND authorized_scoring=true. This CLI never authorizes.',
        provenance='Existence/path/bytes and public JSON snapshots; no hash/mtime. Same-size result replacement is not detectable from size stamps alone.')
    target.parent.mkdir(parents=True,exist_ok=True);base.write_json(target,seal)
    return seal


def authorize(root,path):
    seal=base.read(path)
    if seal.get('schema_version')!=SCHEMA or seal.get('authorized') is not True or seal.get('authorized_scoring') is not True:raise PermissionError('New parent-authorized264-task R15 prediction seal required BEFORE truth access; old endpoint authorization is invalid here')
    if Path(seal['root']).resolve()!=root or seal.get('campaign')!=CAMPAIGN or seal.get('expected_tasks')!=264 or seal.get('expected_cases')!=24 or set(seal.get('fixed_predictions',[]))!=set(base.VIEWS):raise ValueError('Seal root/campaign/full matrix mismatch')
    manifest=base.validate_stamp(root,seal['manifest_file']);cfg=configuration(root,manifest)
    if cfg!=seal['manifest_snapshot']:raise ValueError('Frozen manifest changed')
    queue_paths=[]
    for q in seal['queue_receipts']:
        p=base.validate_stamp(root,q['file'])
        if base.read(p)!=q['snapshot']:raise ValueError('Terminal queue changed after seal')
        queue_paths.append(p)
    _,jobs=queues(root,queue_paths,{t['task_id'] for t in cfg['tasks']})
    originals={t['task_id']:t for t in cfg['tasks']};account={}
    for item in seal['tasks']:
        ident=item['task_id']
        if ident not in originals or ident in account:raise ValueError('Unknown/duplicate sealed task')
        task=originals[ident]
        if inventory(root,task)!=item:raise ValueError('Sealed result files/status/bytes changed: '+ident)
        if item['receipt_file']['exists'] and item['receipt_error'] is None:
            r=base.read(item['receipt_path'])
            if r.get('truth_access')!='NONE' or r.get('task',{}).get('task_id')!=ident:raise ValueError('Online truth declaration/task mismatch')
        account[ident]=item
    if set(account)!=set(originals):raise ValueError('All264 tasks including missing/failed must remain sealed')
    return seal,cfg,[dict(t,solver_label=t['method'],pair=t['method']) for t in cfg['tasks']],account,jobs


def pairwise(rows):
    lookup={(r['case_id'],r['solver']):r for r in rows if r['view']=='primary'}
    differences=[];members=[]
    for method in METHODS:
        for ref in METHODS:
            if method==ref:continue
            groups={k:[] for k in ('gain','loss','both','neither','unavailable')}
            for cid in CASES:
                a,b=lookup[cid,method],lookup[cid,ref];available=base.accuracy_available(a) and base.accuracy_available(b)
                row=dict(case_id=cid,solver=method,reference=ref,paired_available=available)
                for metric in ('NRMSE','IoU','Phi','D','R','relative_residual'):
                    av,bv=a.get(metric),b.get(metric);row['delta_'+metric]=av-bv if available and base.finite(av) and base.finite(bv) else None
                row['candidate_recovered']=base.recovered(a);row['reference_recovered']=base.recovered(b)
                differences.append(row)
                label='unavailable' if not available else 'gain' if base.recovered(a) and not base.recovered(b) else 'loss' if base.recovered(b) and not base.recovered(a) else 'both' if base.recovered(a) else 'neither'
                groups[label].append(cid)
            members.append(dict(solver=method,reference=ref,denominator=24,**{k+'_cases':v for k,v in groups.items()},scope='Posthoc descriptive original thresholds; primary paid minPhi selection; no GT selection'))
    return differences,members


def componentwise_prefix_plans(caches):
    """Saved-history LU/RHS frontiers at fixed CPU cuts; no stopped-run claim."""
    plans=[]
    for cid in CASES:
        group=[c for c in caches.values() if c['task']['case_id']==cid]
        for cutoff in base.PREFIXES:
            boundaries=[]
            for c in group:
                paid=[h for h in c['history'] if base.finite(h.get('CPU_seconds')) and h['CPU_seconds']<=cutoff]
                last=max(paid,key=lambda h:(h['CPU_seconds'],h['candidate'])) if paid else None
                boundaries.append(base.cumulative_physics(last.get('physics_cost')) if last else base.cumulative_physics(None))
            known=len(group)==11 and all(b['LU'] is not None and b['complex_RHS'] is not None for b in boundaries)
            lu=min(b['LU'] for b in boundaries) if known else None
            rhs=min(b['complex_RHS'] for b in boundaries) if known else None
            for c in group:
                eligible=[]
                for h in c['history']:
                    count=base.cumulative_physics(h.get('physics_cost'))
                    if known and base.eligible(h,cutoff) and count['LU'] is not None and count['complex_RHS'] is not None and count['LU']<=lu and count['complex_RHS']<=rhs:eligible.append(h)
                chosen=base.pick(eligible,cutoff);t=c['task']
                plans.append(dict(task_id=t['task_id'],case_id=cid,pair=t['method'],solver=t['method'],campaign=CAMPAIGN,
                    view='common_physics_prefix',selection_kind='pool',CPU_prefix=cutoff,cutoff_CPU_seconds=cutoff,
                    metadata=chosen,candidate=chosen['candidate'] if chosen else None,
                    source_file=c['item']['candidate_pool_path'],source_exists=c['item']['candidate_pool_file']['exists'],
                    common_physical_LU=lu,common_physical_complex_RHS=rhs,physics_counter_limits_available=known,
                    normalization_scope='Componentwise minima across all11 arms of last recorded paid physical LU/RHS at the same fixed native CPU cut; eligible candidates must satisfy CPU and both physical limits. Unrecorded work between audits is unknown; full actual cost retained, not a stopped rerun.'))
    return plans


def costs_and_cycles(root,tasks,caches,jobs,issues):
    costs=[];cycle_rows=[];dynamics=[]
    for task in tasks:
        ident=task['task_id'];folder=root/'RESULTS'/CAMPAIGN/ident;c=caches[ident];r=c['receipt'] or {};job=jobs.get(ident,{})
        external=job.get('external') or {};external_cpu=external.get('child_process_CPU_seconds',job.get('child_process_CPU_seconds'))
        costs.append(dict(task_id=ident,case_id=task['case_id'],solver=task['method'],run_status=c['item']['status'],queue_status=job.get('status'),returncode=job.get('returncode'),
            external_child_CPU_seconds=external_cpu,external_CPU_null_reason=None if base.finite(external_cpu) else 'MISSING_EXTERNAL_WHOLE_CHILD_CPU',
            internal_final_process_CPU_seconds=r.get('whole_process_CPU_seconds',r.get('CPU_seconds')),internal_wall_seconds=r.get('wall_seconds'),external_wall_seconds=job.get('wall_seconds'),
            physical_cost_raw=r.get('physical_cost'),physical_cost=base.cumulative_physics(r.get('physical_cost')),auxiliary_cost=r.get('auxiliary_cost',r.get('auxiliary_receipt')),fixed_projection_cost=r.get('fixed_projection_cost',r.get('projector_cost')),
            termination_reason=r.get('termination_reason',r.get('termination')),cpu=job.get('cpu'),core_type=job.get('core_type'),queue_receipt=job.get('queue_receipt_path'),
            accounting_scope='Final child CPU already includes setup/auxiliary/checkpoints. External child CPU is not added to internal CPU; physical role counts not added to aggregate RHS. Failures charged.'))
        paid={}
        for h in c['history']:
            meta=(h.get('extra') or {}).get('cycle_metadata')
            if isinstance(meta,dict) and meta.get('cycle') is not None:
                paid.setdefault(meta['cycle'],[]).append(h)
        path=folder/'CYCLES.jsonl';events=[]
        if path.is_file():
            with path.open() as f:
                for line_number,line in enumerate(f,1):
                    if not line.strip():continue
                    try:events.append(json.loads(line))
                    except Exception as exc:issues.append(dict(task_id=ident,issue='INVALID_CYCLE_LINE',line=line_number,error=repr(exc)))
        for event in events:
            audit=paid.get(event.get('cycle'),[])
            entry=dict(task_id=ident,case_id=task['case_id'],solver=task['method'],**event,paid_physical_audits=[{k:h.get(k) for k in ('candidate','CPU_seconds','Phi','D','R','eligible','phase','physics_cost')} for h in audit])
            # Keep all audits, including overbudget diagnostics; not a new selector.
            entry['paid_Phi']=audit[0].get('Phi') if audit else None
            entry['paid_Phi_null_reason']=None if audit else 'NO_MATCHING_ARCHIVED_PAID_CYCLE_AUDIT'
            cycle_rows.append(entry)
        dynamics.append(dict(task_id=ident,case_id=task['case_id'],solver=task['method'],run_status=c['item']['status'],cycle_file_exists=path.is_file(),
            archived_cycle_records=len(events) if path.is_file() else None,cycle_count_null_reason=None if path.is_file() else 'NO_ARCHIVED_CYCLE_FILE; native baselines do not share stateful cycle instrumentation',
            complete_cycles=sum(e.get('status')=='COMPLETE' for e in events) if path.is_file() else None,
            interrupted_cycles=sum(e.get('status')!='COMPLETE' for e in events) if path.is_file() else None,
            paid_cycle_audits=sum(len(v) for v in paid.values()),history_candidates=len(c['history']),last_cycle=events[-1] if events else None,
            scope='Recorded cycles only; absent gap is unknown/not instrumented, never zero; Phi is original paid audit, no causal/recovery claim'))
    return costs,cycle_rows,dynamics


def stratified_and_fixed_primary_descriptions(cfg,rows,np):
    """Frozen acquisition labels and seven fixed posthoc object-cluster summaries.

    These physical axes are coupled, not independent randomized factors. This
    function consumes scored rows only and never changes a candidate/threshold.
    """
    scope=('Already exposed same-cohort posthoc description. Manifest acquisition axes '
           'are coupled; neither marginal strata nor joint strata identify independent causal effects. '
           'No threshold search, p-values, significance or algorithm-win declaration.')
    bycase={}
    for t in cfg['tasks']:
        contract=t.get('observation_contract',{})
        aperture=contract.get('receiver_aperture')
        receiver='full' if aperture=='full_circle' else 'half' if aperture=='half_circle' else str(aperture or 'UNKNOWN')
        sources=contract.get('sources')
        noise=contract.get('relative_noise_RMS')
        mismatch=contract.get('model_mismatch_status')
        model='matched' if mismatch=='MATCHED_DISCRETE_MODEL' else 'crossgrid' if isinstance(mismatch,str) and 'CROSSGRID' in mismatch.upper() else str(mismatch or 'UNKNOWN')
        labels=dict(receiver_aperture=receiver,source_count='source'+str(sources) if sources is not None else 'UNKNOWN',
                    noise='noise'+format(float(noise),'.6g') if base.finite(noise) else 'UNKNOWN',model_match=model)
        actual=dict(receiver_aperture=aperture,sources=sources,relative_noise_RMS=noise,model_mismatch_status=mismatch)
        cid=t['case_id']
        if cid in bycase and bycase[cid]['labels']!=labels:raise ValueError('Same-case method acquisition contracts differ: '+cid)
        bycase[cid]=dict(case_id=cid,labels=labels,observation_contract_values=actual,scope=scope)
    if set(bycase)!=set(CASES):raise ValueError('Full24 stratum membership required')
    membership=[dict(case_id=c,**bycase[c]['labels'],observation_contract_values=bycase[c]['observation_contract_values'],scope=scope) for c in CASES]
    groups={}
    for c in CASES:
        labels=bycase[c]['labels']
        for axis,value in labels.items():groups.setdefault((axis,value),[]).append(c)
        joint=' | '.join(k+'='+v for k,v in labels.items())
        groups.setdefault(('joint_observation_contract',joint),[]).append(c)
    views=('primary','raw','corrected')
    lookup={(r['case_id'],r['solver'],r['view']):r for r in rows if r['view'] in views}
    summaries=[];paired=[]
    candidates=('twofold_hop','twofold_dr','twofold_reset')
    references=('g19','g19_s19','cie','dual_cie','dr')
    for (axis,value),ids in sorted(groups.items()):
        subset=[lookup[c,m,v] for c in ids for m in METHODS for v in views]
        summaries.extend(dict(r,stratum_axis=axis,stratum=value,stratum_case_ids=ids,stratum_cases=len(ids),stratification_scope=scope)
                         for r in extension.summarize(subset,ids))
        for view in views:
            for candidate in candidates:
                for reference in references:
                    buckets={k:[] for k in ('gain','loss','both','neither','unavailable')}
                    nrmse_pairs=[];nrmse_ids=[]
                    for cid in ids:
                        a,b=lookup[cid,candidate,view],lookup[cid,reference,view]
                        available=base.accuracy_available(a) and base.accuracy_available(b)
                        label='unavailable' if not available else 'gain' if base.recovered(a) and not base.recovered(b) else 'loss' if base.recovered(b) and not base.recovered(a) else 'both' if base.recovered(a) else 'neither'
                        buckets[label].append(cid)
                        if a.get('available') is True and b.get('available') is True and base.finite(a.get('NRMSE')) and base.finite(b.get('NRMSE')):
                            nrmse_pairs.append(a['NRMSE']-b['NRMSE']);nrmse_ids.append(cid)
                    paired.append(dict(stratum_axis=axis,stratum=value,stratum_case_ids=ids,view=view,solver=candidate,reference=reference,
                        denominator=len(ids),paired_recovery_available=len(ids)-len(buckets['unavailable']),
                        **{k+'_cases':v for k,v in buckets.items()},
                        candidate_recovered_count=sum(base.recovered(lookup[c,candidate,view]) for c in ids),
                        reference_recovered_count=sum(base.recovered(lookup[c,reference,view]) for c in ids),
                        paired_NRMSE_case_ids=nrmse_ids,paired_NRMSE_count=len(nrmse_pairs),
                        paired_mean_delta_NRMSE=float(np.mean(nrmse_pairs)) if nrmse_pairs else None,
                        paired_median_delta_NRMSE=float(np.median(nrmse_pairs)) if nrmse_pairs else None,
                        continuous_null_reason=None if nrmse_pairs else 'NO_AVAILABLE_PAIRED_NRMSE',scope=scope))
    fixed=(('g19_s19','g19'),('twofold','g19'),('twofold_hop','twofold'),
           ('twofold_dr','twofold_hop'),('twofold_ap','twofold_hop'),
           ('twofold_full_dr','twofold_hop'),('twofold_reset','twofold_hop'))
    seed=20261116;draws=4000;rng=np.random.default_rng(seed)
    # One material object is one paired cluster. Identical24-ID resamples are
    # shared across all seven comparisons; no treating views/methods as samples.
    indices=rng.integers(0,len(CASES),size=(draws,len(CASES)))
    bootstrap=[];pair_values=[]
    for candidate,reference in fixed:
        delta=np.full(len(CASES),np.nan);missing=[]
        for i,cid in enumerate(CASES):
            a,b=lookup[cid,candidate,'primary'],lookup[cid,reference,'primary']
            available=a.get('available') is True and b.get('available') is True and base.finite(a.get('NRMSE')) and base.finite(b.get('NRMSE'))
            if available:delta[i]=a['NRMSE']-b['NRMSE']
            else:missing.append(cid)
            pair_values.append(dict(case_id=cid,solver=candidate,reference=reference,view='primary',paired_NRMSE_available=available,
                candidate_NRMSE=a.get('NRMSE'),reference_NRMSE=b.get('NRMSE'),delta_NRMSE=float(delta[i]) if available else None,
                null_reason=None if available else 'MISSING_OR_UNAVAILABLE_PAIRED_PRIMARY_NRMSE'))
        observed=delta[np.isfinite(delta)]
        sample=delta[indices];counts=np.isfinite(sample).sum(axis=1);valid=counts>0
        if np.any(valid):
            means=np.nansum(sample[valid],axis=1)/counts[valid]
            medians=np.nanmedian(sample[valid],axis=1)
            mean_low,mean_high=map(float,np.quantile(means,[.025,.975]))
            median_low,median_high=map(float,np.quantile(medians,[.025,.975]))
        else:mean_low=mean_high=median_low=median_high=None
        bootstrap.append(dict(solver=candidate,reference=reference,view='primary',metric='delta_NRMSE_candidate_minus_reference',
            cohort_cluster_denominator=24,all_case_ids=list(CASES),paired_available=int(len(observed)),unavailable_case_ids=missing,
            mean_delta_NRMSE=float(np.mean(observed)) if len(observed) else None,
            median_delta_NRMSE=float(np.median(observed)) if len(observed) else None,
            bootstrap95_mean_low=mean_low,bootstrap95_mean_high=mean_high,
            bootstrap95_median_low=median_low,bootstrap95_median_high=median_high,
            bootstrap_seed=seed,bootstrap_draws=draws,bootstrap_cluster_draw_size=24,
            bootstrap_valid_replicates=int(valid.sum()),bootstrap_no_available_pair_replicates=int((~valid).sum()),
            bootstrap_missing_policy='Resample all24 material IDs including unavailable clusters. No imputation. Each statistic uses available paired values in its resample; empty resamples excluded and counted. Original recovery denominator remains24.',
            bootstrap_resamples='Same4000 object-index resamples shared across the seven prelisted comparisons',
            interval='Percentile2.5/97.5% descriptive intervals for paired object mean and median; already exposed same cohort, unadjusted, no population/causal or significance-win claim',
            statistic_null_reason=None if len(observed) else 'NO_AVAILABLE_PAIRED_PRIMARY_NRMSE',scope=scope))
    return [('STRATA_MEMBERSHIP',membership),('STRATIFIED_SUMMARY',summaries),
            ('STRATIFIED_PAIRED_GAIN_LOSS',paired),('FIXED_PRIMARY_PAIRED_VALUES',pair_values),
            ('FIXED_PRIMARY_CLUSTER_BOOTSTRAP',bootstrap)]


def figures(out,rows,pictures,truths,np,costs,effects,cycles):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    dest=out/'FIGURES';dest.mkdir();artifacts=[]
    primary={(r['case_id'],r['solver']):r for r in rows if r['view']=='primary'}
    def save(fig,name):
        for suffix in ('png','pdf'):
            p=dest/(name+'.'+suffix);fig.savefig(p,dpi=135,bbox_inches='tight');artifacts.append(dict(path=str(p.relative_to(out)),bytes=p.stat().st_size))
        plt.close(fig)
    # Two physical pages, every case and all eleven methods, rowwise shared scale.
    for page,start in enumerate((0,12),1):
        subset=CASES[start:start+12];fig=plt.figure(figsize=(25,27));grid=fig.add_gridspec(12,13,width_ratios=[1]*12+[.06],hspace=.15,wspace=.07)
        for i,cid in enumerate(subset):
            arrays=[truths[cid]['m'] if truths.get(cid) else None]+[pictures.get((cid,m)) for m in METHODS]
            vmax=max([.02]+[float(np.max(a)) for a in arrays if a is not None]);image=None
            for j,(label,array) in enumerate(zip(('GT',)+METHODS,arrays)):
                ax=fig.add_subplot(grid[i,j]);ax.set_xticks([]);ax.set_yticks([])
                if i==0:ax.set_title(label,rotation=35,ha='left',fontsize=8)
                if j==0:ax.set_ylabel(cid,fontsize=9)
                if array is None:ax.text(.5,.5,'MISSING',ha='center',va='center',transform=ax.transAxes,fontsize=7)
                else:image=ax.imshow(array.reshape(24,24),origin='lower',vmin=.02,vmax=vmax,cmap='viridis')
            cb=fig.add_subplot(grid[i,12])
            if image is not None:fig.colorbar(image,cax=cb)
            else:cb.axis('off')
        fig.suptitle(f'Exposed24 primary paid minPhi predictions — page {page}/2; shared scale within case',y=.995)
        save(fig,f'materials_page{page}')
    fig,axes=plt.subplots(1,2,figsize=(19,12));recovery=np.full((24,11),np.nan);delta=np.full((24,11),np.nan)
    for i,c in enumerate(CASES):
        ref=primary[c,'g19']
        for j,m in enumerate(METHODS):
            r=primary[c,m]
            if base.accuracy_available(r):recovery[i,j]=float(base.recovered(r))
            if base.accuracy_available(r) and base.accuracy_available(ref):delta[i,j]=r['NRMSE']-ref['NRMSE']
    for ax,matrix,title in [(axes[0],recovery,'Original posthoc joint recovery; blank=missing (denominator24)'),(axes[1],delta,'Paired continuous NRMSE difference versus G19 primary')]:
        bound=float(np.nanmax(abs(matrix))) if np.isfinite(matrix).any() else 1.
        im=ax.imshow(np.ma.masked_invalid(matrix),aspect='auto',cmap='Greys' if ax is axes[0] else 'coolwarm',vmin=0 if ax is axes[0] else -(bound or 1),vmax=1 if ax is axes[0] else (bound or 1))
        ax.set_xticks(range(11),METHODS,rotation=55,ha='right',fontsize=8);ax.set_yticks(range(24),CASES,fontsize=8);ax.set_title(title,fontsize=10);fig.colorbar(im,ax=ax)
    fig.tight_layout();save(fig,'paired_recovery_and_nrmse')
    fig,axes=plt.subplots(2,2,figsize=(19,13))
    for ax,metric in [(axes[0,0],'NRMSE'),(axes[0,1],'Phi')]:
        for m in METHODS:
            means=[]
            for cut in base.PREFIXES:
                rr=[r[metric] for r in rows if r['solver']==m and r['view']=='common_CPU_prefix' and r.get('CPU_prefix')==cut and base.finite(r.get(metric))]
                means.append(sum(rr)/len(rr) if rr else np.nan)
            ax.plot(base.PREFIXES,means,'.-',label=m)
        ax.set_title('Common recorded paid CPU prefix: '+metric+' (available means)');ax.set_xlabel('Native whole-process CPU seconds');ax.grid(alpha=.2)
    ax=axes[1,0]
    for j,m in enumerate(METHODS):
        vals=[r['NRMSE'] for r in rows if r['solver']==m and r['view']=='common_physics_pool' and base.finite(r.get('NRMSE'))]
        if vals:ax.bar(j,sum(vals)/len(vals));ax.text(j,sum(vals)/len(vals),f'n={len(vals)}',ha='center',va='bottom',fontsize=7)
    ax.set_title('Componentwise common final LU/RHS paid pool (case-specific limits)');ax.set_xticks(range(11),METHODS,rotation=55,ha='right',fontsize=8);ax.set_ylabel('Mean NRMSE among available')
    ax=axes[1,1]
    for j,m in enumerate(METHODS):
        vals=[c['external_child_CPU_seconds'] for c in costs if c['solver']==m and base.finite(c.get('external_child_CPU_seconds'))]
        ax.scatter([j]*len(vals),vals,s=12,alpha=.7)
    ax.set_title('All actual external child CPU; includes failed/setup/export work');ax.set_ylabel('OS child lifetime CPU seconds');ax.set_xticks(range(11),METHODS,rotation=55,ha='right',fontsize=8)
    axes[0,0].legend(fontsize=7,ncol=2);fig.tight_layout();save(fig,'cpu_and_componentwise_physics_prefixes')
    fig,axes=plt.subplots(1,2,figsize=(19,8))
    for ax,metric in zip(axes,('NRMSE','Phi')):
        for j,m in enumerate(METHODS):
            vals=[r['delta_'+metric] for r in effects if r['solver']==m and base.finite(r.get('delta_'+metric))]
            ax.scatter([j]*len(vals),vals,s=16,alpha=.7)
        ax.axhline(0,color='gray',lw=.8);ax.set_xticks(range(11),METHODS,rotation=50,ha='right',fontsize=8);ax.set_title('Corrected − raw '+metric+' (all available cases)')
    fig.tight_layout();save(fig,'corrector_help_and_harm')
    fig,axes=plt.subplots(1,3,figsize=(21,7))
    for ax,key,title in zip(axes,('full_shadow_gap','routed_shadow_gap','paid_Phi'),('Full shadow gap','Routed shadow gap','Original paid cycle physical Phi')):
        count=0;labelled=set()
        grouped={}
        for event in cycles:
            if base.finite(event.get(key)):grouped.setdefault((event['solver'],event['case_id']),[]).append(event)
        for j,m in enumerate(METHODS):
            for cid in CASES:
                rr=grouped.get((m,cid),[])
                rr.sort(key=lambda x:x.get('cycle',0))
                if rr:
                    ax.plot([x['cycle'] for x in rr],[x[key] for x in rr],alpha=.3,lw=.8,color=f'C{j%10}',label=m if m not in labelled else None)
                    labelled.add(m);count+=len(rr)
        if count:ax.set_yscale('symlog',linthresh=1e-8)
        else:ax.text(.5,.5,'No archived values (not zero)',ha='center',transform=ax.transAxes)
        ax.set_title(title);ax.set_xlabel('Recorded cycle');ax.legend(fontsize=7)
    fig.suptitle('Archived cycle diagnostics — all cases; no recovery-predictor or common-baseline instrumentation claim')
    fig.tight_layout();save(fig,'cycle_shadow_gaps_and_paid_phi')
    return artifacts


def score(root,seal_path,out):
    cpu,wall=time.process_time(),time.perf_counter()
    # This is the only path to truth loading; old endpoint authorization cannot pass.
    seal,cfg,tasks,account,jobs=authorize(root,seal_path)
    if out.exists() or out!=root/'analysis'/CAMPAIGN:raise ValueError('New canonical analysis/twofold_exposed24_v1 output required')
    out.mkdir(parents=True);issues=[]
    try:
        plans,caches=base.public_plans(root,tasks,account,issues,cases=CASES,expected_methods=11)
        plans=[p for p in plans if not p['view'].startswith('strict48_')]
        plans.extend(componentwise_prefix_plans(caches))
        issues=[i for i in issues if i.get('issue')!='ARCHIVED_SELECTOR_METADATA_MISMATCH']
        audits=view_audit.audit_views(plans,caches,issues)
        # Freeze all public selectors before even importing the original truth helper.
        base.write_json(out/'PUBLIC_PLANS_BEFORE_TRUTH.json',plans);base.table(out,'SELECTOR_AUDIT',audits)
        import numpy as np
        helper_path=R14/'runtime/scoring_history/src/score_results.py'
        sys.path.insert(0,str(helper_path.parent))
        spec=importlib.util.spec_from_file_location('r15_original_material_metrics',helper_path);helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)
        truths={};sources=[]
        for cid in CASES:
            source=dict(path=str(R14/'sealed_truth/independent24'/(cid+'.npz')),key='material_coarse',support_rule='fixed_chi_gt_1',host=.1,object_id=cid)
            sources.append(dict(case_id=cid,**source))
            try:truths[cid]=helper.load_truth(source)
            except Exception as exc:truths[cid]=None;issues.append(dict(case_id=cid,issue='AUTHORIZED_TRUTH_UNAVAILABLE',error=repr(exc)))
        rows,pictures=base.score_plans(plans,caches,truths,helper.material_metrics,helper.physical_counts,np,issues)
        for row in rows:
            row.update(provenance='R15_EXPOSED24_FROZEN264',cohort_denominator=24,selection_scope='Original paid legal/eligible minPhi only; GT never selects candidates')
            for metric in base.METRICS:row.setdefault(metric,None)
        costs,cycles,dynamics=costs_and_cycles(root,tasks,caches,jobs,issues)
        effects=extension.corrector_effects(rows,np,issues)
        for r in effects:
            a=next(x for x in rows if x['case_id']==r['case_id'] and x['solver']==r['solver'] and x['view']=='raw')
            b=next(x for x in rows if x['case_id']==r['case_id'] and x['solver']==r['solver'] and x['view']=='corrected')
            r.update(raw_recovered=base.recovered(a),corrected_recovered=base.recovered(b),recovery_gained=base.accuracy_available(a) and base.accuracy_available(b) and not base.recovered(a) and base.recovered(b),recovery_lost=base.accuracy_available(a) and base.accuracy_available(b) and base.recovered(a) and not base.recovered(b))
        differences,members=pairwise(rows);normalized_diffs=[];normalized_members=[]
        for view in ('common_CPU_prefix','common_physics_prefix','common_physics_pool'):
            for cut in (base.PREFIXES if view!='common_physics_pool' else (None,)):
                group=[dict(r,view='primary') for r in rows if r['view']==view and r.get('CPU_prefix')==cut]
                d,m=pairwise(group);normalized_diffs.extend(dict(r,comparison_view=view,CPU_prefix=cut) for r in d);normalized_members.extend(dict(r,comparison_view=view,CPU_prefix=cut) for r in m)
        primary=[r for r in rows if r['view']=='primary']
        matrix=[]
        for cid in CASES:
            record=dict(case_id=cid)
            for m in METHODS:
                r=next(x for x in primary if x['case_id']==cid and x['solver']==m)
                for k in ('run_status','available','accuracy_available','NRMSE','IoU','Phi','recovered_posthoc','secondary_fallback'):record[m+'__'+k]=r.get(k)
            matrix.append(record)
        queue_costs=[];attempts=[]
        for q in seal['queue_receipts']:
            snapshot=q['snapshot'];queue_costs.append(dict(path=q['file']['path'],status=snapshot.get('status'),queue_wall_seconds=snapshot.get('queue_wall_seconds'),controller_CPU_seconds=snapshot.get('controller_CPU_seconds'),scope='Controller CPU disjoint from children; wall intervals never added to child CPU'))
            attempts.extend(dict(j,queue_receipt=q['file']['path']) for j in snapshot.get('jobs',[]))
        datasets=[('ALL_VIEWS',rows),('SUMMARY',extension.summarize(rows,CASES)),('PRIMARY_MATRIX24',matrix),('PRIMARY_TASKS264',primary),('PAIRED_DIFFERENCES',differences),('RECOVERY_GAIN_LOSS',members),('NORMALIZED_PAIRED_DIFFERENCES',normalized_diffs),('NORMALIZED_GAIN_LOSS',normalized_members),('CORRECTOR_EFFECTS',effects),('COST_LEDGER',costs),('QUEUE_COSTS',queue_costs),('ALL_QUEUE_ATTEMPTS',attempts),('CYCLE_DIAGNOSTICS',cycles),('DYNAMICS_SUMMARY',dynamics),('TRUTH_SOURCES',sources)]
        datasets.extend(stratified_and_fixed_primary_descriptions(cfg,rows,np))
        for name,view in [('PAID_CPU_PREFIXES','paid_prefix'),('COMMON_CPU_PREFIXES','common_CPU_prefix'),('COMMON_LU_RHS_PREFIXES','common_physics_prefix'),('COMMON_FINAL_LU_RHS_POOL','common_physics_pool')]:datasets.append((name,[r for r in rows if r['view']==view]))
        for name,data in datasets:base.table(out,name,data)
        artifacts=figures(out,rows,pictures,truths,np,costs,effects,cycles);base.table(out,'ISSUES',issues)
        external=[r['external_child_CPU_seconds'] for r in costs]
        receipt=dict(status='COMPLETE_WITH_RETAINED_ISSUES' if issues else 'COMPLETE',campaign=CAMPAIGN,seal=str(seal_path),expected_cases=24,expected_tasks=264,methods=list(METHODS),
            primary_rows=len(primary),full24_matrix_rows=len(matrix),issues=len(issues),truth_access='ONLY_AFTER_NEW_PARENT_AUTHORIZED_PREDICTION_SEAL',
            material_metrics_source=str(helper_path),success_rule='Inherited posthoc NRMSE<=.15 and fixed chi>1 IoU>=.70; not candidate selection',
            selection_rule='Inherited completed paid legal/eligible original Phi; secondary minR in originalCN99, primary fallback explicitly retained',
            normalization='Saved native CPU prefixes, shared CPU cutoff, componentwise minima of final physical LU/RHS across all11 same-case arms. Full actual costs remain charged; no reruns or free auxiliary operations.',
            external_final_task_attempt_CPU_total=sum(external) if all(base.finite(x) for x in external) else None,external_final_task_attempt_CPU_measured_subtotal=sum(x for x in external if base.finite(x)),external_CPU_missing_tasks=[r['task_id'] for r in costs if not base.finite(r['external_child_CPU_seconds'])],
            stratification='Frozen manifest observation_contract labels: receiver full/half, sources6/12, noise.01/.02, matched/crossgrid, plus joint labels. Coupled physical axes; exposed-cohort posthoc descriptions, not independent causal effects.',
            fixed_primary_bootstrap='Seven prelisted candidate-reference pairs; material-object clusters24, seed20261116,4000 shared draws; paired mean/median NRMSE delta percentile95 descriptive intervals, no p-values/significance-win claims.',
            all_queue_attempts_cost_scope='All launch/process attempts retained in ALL_QUEUE_ATTEMPTS; repeated queue snapshots may overlap. Final-task totals do not assert an all-attempt total when multiple receipts supplied.',
            figures=artifacts,science_interpretation='Parent only; exposed cohort, no new blind confirmation or population inference',new_physics=0,new_optimizations=0,new_SVD=0,
            CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall)
        base.write_json(out/'RECEIPT.json',receipt)
        links=''.join('<p><a href="'+html.escape(p.name)+'">'+html.escape(p.name)+'</a></p>' for p in sorted(out.glob('*.csv')))
        imgs=''.join('<figure><img loading="lazy" style="max-width:100%" src="'+html.escape(a['path'])+'"><figcaption>'+html.escape(a['path'])+'</figcaption></figure>' for a in artifacts if a['path'].endswith('.png'))
        (out/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>R15 exposed24 sealed campaign evidence</title><main style="max-width:1800px;margin:auto;font-family:sans-serif"><h1>R15 exposed24: 24 ×11 fixed methods</h1><p>All methods and failures retained. Original descriptive thresholds do not select candidates. Scientific conclusions reserved for parent.</p><pre>'+html.escape(json.dumps(receipt,indent=2))+'</pre>'+links+imgs+'</main>',encoding='utf8')
        print(json.dumps(receipt))
    except BaseException as exc:
        base.write_json(out/'FAILED_RECEIPT.json',dict(status='FAILED',error=repr(exc),seal=str(seal_path),CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,new_physics=0));raise


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=ROOT);p.add_argument('--manifest',type=Path);p.add_argument('--queue',type=Path,action='append',default=[])
    group=p.add_mutually_exclusive_group(required=True);group.add_argument('--make-seal',type=Path);group.add_argument('--seal',type=Path)
    p.add_argument('--out',type=Path);p.add_argument('--source-note',default='')
    a=p.parse_args();root=a.root.resolve()
    if a.make_seal:
        cpu,wall=time.process_time(),time.perf_counter();target=a.make_seal.resolve();manifest=(a.manifest or root/'configs'/f'{CAMPAIGN}.json').resolve()
        make_seal(root,manifest,target,a.queue,a.source_note)
        receipt=dict(status='UNAUTHORIZED_SEAL_CREATED',seal=str(target),authorized=False,authorized_scoring=False,truth_access='NONE',CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall)
        base.write_json(target.with_suffix('.receipt.json'),receipt);print(json.dumps(receipt))
    else:score(root,a.seal.resolve(),(a.out or root/'analysis'/CAMPAIGN).resolve())


if __name__=='__main__':main()
