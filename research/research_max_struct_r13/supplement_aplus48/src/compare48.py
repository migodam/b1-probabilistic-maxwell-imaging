"""Authorized FOLLOWUP_EXISTING48 comparison, no forwards or method changes.

Import NumPy/plotting and open scored metrics/materials/truth only AFTER the
parent comparison authorization and all three original/extension seals pass.
"""
import os,time
START_CPU=time.process_time();START_WALL=time.perf_counter()
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
import argparse,csv,itertools,json,math,platform,sys
from pathlib import Path
from datetime import datetime,timezone
METHODS=('cie','dual_cie','g19','dr','a_plus')
VIEWS=('initial','raw','corrected','primary','secondary')
PREFIXES=(15,30,60,90,120)
METRICS=('RMSE','NRMSE','boxrange_NRMSE','relativeL2','IoU','mean_amplitude_bias','mean_amplitude_bias_normalized','D','R','Phi','residual_over_delta')

def read(path):return json.loads(path.read_text())
IO_ISSUES=[]
def optional(path,default):
    if not path.is_file():return default
    try:
        value=read(path)
        if not isinstance(value,type(default)):raise ValueError('Unexpected JSON structure')
        return value
    except (OSError,ValueError,TypeError) as ex:
        IO_ISSUES.append(dict(path=str(path),error=repr(ex),retained_as_missing=True));return default

def table(out,name,rows):
    (out/(name+'.json')).write_text(json.dumps(rows,indent=2,allow_nan=False)+'\n')
    fields=list(dict.fromkeys(k for row in rows for k in row))
    with (out/(name+'.csv')).open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=fields,lineterminator='\n');writer.writeheader()
        for row in rows:writer.writerow({k:json.dumps(v,separators=(',',':')) if isinstance(v,(dict,list)) else v for k,v in row.items()})

def authorize(root,sup,path):
    auth=read(path)
    if auth.get('authorized_comparison_scoring') is not True:raise PermissionError('Explicit parent authorized_comparison_scoring=true required')
    if auth.get('scope')!='FOLLOWUP_EXISTING48':raise PermissionError('Parent must label FOLLOWUP_EXISTING48, not fresh blind')
    if auth.get('root') and Path(auth['root']).resolve()!=root:raise ValueError('Authorization root mismatch')
    if auth.get('supplement_root') and Path(auth['supplement_root']).resolve()!=sup:raise ValueError('Authorization supplement root mismatch')
    sys.path.insert(0,str(root/'src'));from seal_results import validate_seal
    defaults={'confirmation_core48':root/'audit/confirmation_core48_seal.json','confirmation_mechanism16':root/'audit/confirmation_mechanism16_seal.json','aplus48_extension32':sup/'audit/aplus48_extension32_seal.json'};seals={};locations={}
    for campaign,default in defaults.items():
        value=auth.get('seals',{}).get(campaign);location=Path(value) if value else default
        if not location.is_absolute():location=root/location
        s=read(location);validate_seal(s,require_authorization=True)
        expected=sup if campaign=='aplus48_extension32' else root
        if Path(s['root']).resolve()!=expected or s['campaign']!=campaign:raise ValueError('Seal campaign/root mismatch: '+campaign)
        seals[campaign]=s;locations[campaign]=str(location.resolve())
    # Exact original public cohort and32-complement; source-selection rules frozen.
    core=seals['confirmation_core48']['manifest_snapshot'];old=seals['confirmation_mechanism16']['manifest_snapshot'];ext=seals['aplus48_extension32']['manifest_snapshot']
    cases=[f'r13_fresh_{i:04d}' for i in range(48)];oldtasks=[t for t in old['tasks'] if t['method']=='a_plus'];newtasks=ext['tasks']
    if sorted(t['case_id'] for t in oldtasks)!=cases[::3] or sorted(t['case_id'] for t in newtasks)!=sorted(set(cases)-set(cases[::3])):raise ValueError('Old16/new32 cohort differs from exact frozen complement')
    if any(t['method']!='a_plus' or t.get('algorithm','a_plus')!='a_plus' or t['task_id']!=t['case_id']+'__a_plus' for t in newtasks+oldtasks):raise ValueError('Unexpected A+ arm/task ID')
    if old['protocol']!=ext['protocol'] or old['protocol']!=core['protocol']:raise ValueError('Numerical/selection/budget protocol changed')
    bycase={t['case_id']:t for t in core['tasks'] if t['method']=='g19'}
    for t in oldtasks+newtasks:
        if t['input_file'].replace('\\','/')!=bycase[t['case_id']]['input_file'].replace('\\','/'):raise ValueError('Original observation reference changed')
    return auth,seals,locations,cases


def summarize(rows,np):
    result=[]
    for split in ['all48','reused16','new32']:
        for method,view in itertools.product(METHODS,sorted({r['view'] for r in rows})):
            rr=[r for r in rows if r['method']==method and r['view']==view and (split=='all48' or r['aplus_source_split']==split)]
            if not rr:continue
            item=dict(split=split,method=method,view=view,expected_cases=len(rr),available_cases=sum(bool(r.get('available')) for r in rr),failed_runs=sum(r.get('run_status')=='FAILED' for r in rr),ball_known_cases=sum(isinstance(r.get('noise_feasible'),bool) for r in rr),ball_feasible_cases=sum(r.get('noise_feasible') is True for r in rr))
            for metric in METRICS:
                vals=[r[metric] for r in rr if isinstance(r.get(metric),(int,float)) and math.isfinite(r[metric])];item[metric+'_n']=len(vals);item[metric+'_mean']=float(np.mean(vals)) if vals else None;item[metric+'_median']=float(np.median(vals)) if vals else None
            result.append(item)
    return result


def subgroup_summaries(rows,np):
    result=[]
    axes=('source_count','receiver_aperture','relative_noise_RMS','forward_side','model_mismatch_status','acquisition_noise_condition')
    for axis in axes:
        for value in sorted({str(r.get(axis)) for r in rows}):
            subset=[r for r in rows if str(r.get(axis))==value]
            for r in summarize(subset,np):result.append(dict(r,public_axis=axis,public_value=value,qualification='Descriptive conditional strata; acquisition/noise/source-count/grid axes can be confounded, not causal effects or new independent cohorts'))
    return result


def paired(rows,np,views=VIEWS):
    changes=[];stats=[];lookup={(r['case_id'],r['method'],r['view']):r for r in rows};cases=sorted({r['case_id'] for r in rows});rng=np.random.default_rng(131313)
    for split in ['all48','reused16','new32']:
        cc=[cid for cid in cases if split=='all48' or (cid in cases[::3])==(split=='reused16')]
        for view,reference in itertools.product(views,METHODS[:-1]):
            diffs={k:[] for k in ['NRMSE','RMSE','IoU']}
            for cid in cc:
                a=lookup.get((cid,'a_plus',view),{});b=lookup.get((cid,reference,view),{});row=dict(split=split,case_id=cid,method='a_plus',reference=reference,view=view,paired_available=bool(a.get('available') and b.get('available')))
                for metric in diffs:
                    x,y=a.get(metric),b.get(metric);delta=x-y if row['paired_available'] and x is not None and y is not None else None;row['a_plus_'+metric]=x;row['reference_'+metric]=y;row['delta_'+metric]=delta
                    if delta is not None:diffs[metric].append(delta)
                changes.append(row)
            item=dict(split=split,method='a_plus',reference=reference,view=view,expected_pairs=len(cc),paired_available=sum(x['paired_available'] for x in changes[-len(cc):]),bootstrap_seed=131313,bootstrap_replicates=4000,unit='material object; descriptive unadjusted; no success gate')
            for metric,values in diffs.items():
                d=np.asarray(values,float);boots=np.mean(d[rng.integers(0,len(d),size=(4000,len(d)))],axis=1) if len(d) else[];item[metric+'_n']=len(d);item[metric+'_mean_delta']=float(np.mean(d)) if len(d) else None;item[metric+'_median_delta']=float(np.median(d)) if len(d) else None;item[metric+'_bootstrap95_low']=float(np.quantile(boots,.025)) if len(d) else None;item[metric+'_bootstrap95_high']=float(np.quantile(boots,.975)) if len(d) else None
            stats.append(item)
    return changes,stats


def figures(out,cases,rows,prefixrows,costrows,truths,pictures,np):
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    dest=out/'FIGURES';dest.mkdir();artifacts=[];labels=['GT','CIE','Dual-CIE','G19','DR','A+'];lookup={(r['case_id'],r['method'],r['view']):r for r in rows};plt.rcParams.update({'font.size':10,'pdf.fonttype':42})
    def save(fig,name):
        for extension in ['png','pdf']:
            f=dest/(name+'.'+extension);fig.savefig(f,dpi=170,bbox_inches='tight');artifacts.append(str(f.relative_to(out)))
        plt.close(fig)
    for diagnostic in [False,True]:
        for page in range(12):
            subset=cases[page*4:(page+1)*4];fig,axes=plt.subplots(4,6,figsize=(17,12),layout='constrained')
            for i,cid in enumerate(subset):
                gt=truths.get(cid);truth=gt['m'] if gt else None;vmax=float(np.max(truth)) if diagnostic and truth is not None else 20.;vmax=max(vmax,1e-12)
                arrays=[truth]+[pictures.get((cid,m)) for m in METHODS]
                for j,array in enumerate(arrays):
                    ax=axes[i,j];ax.set_xticks([]);ax.set_yticks([])
                    if i==0:ax.set_title(labels[j])
                    if j==0:ax.set_ylabel(cid+'\n'+('row maxχ='+f'{vmax:.3f}' if diagnostic else 'physical χ[0,20]'))
                    if array is None:ax.text(.5,.5,'UNAVAILABLE',ha='center',va='center');continue
                    im=ax.imshow(np.asarray(array).reshape(24,24),origin='lower',cmap='viridis',norm=Normalize(0,vmax,clip=True),interpolation='nearest')
                    if j:
                        r=lookup.get((cid,METHODS[j-1],'primary'),{});metric=r.get('NRMSE');ax.set_xlabel(('NRMSE='+f'{metric:.4f}' if metric is not None else 'No available metric')+('\nabove-row-max='+str(np.sum(array>vmax)) if diagnostic else''),fontsize=9)
                if diagnostic and truth is not None:fig.colorbar(im,ax=axes[i,:].tolist(),fraction=.014,pad=.01,label='Physical χ, shared row scale')
            if not diagnostic:fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0,20),cmap='viridis'),ax=axes.ravel().tolist(),fraction=.015,pad=.01,label='Physical χ [0,20]')
            fig.suptitle('FOLLOWUP_EXISTING48 — all48 cases, primary saved estimates\n'+('Diagnostic GT-derived common row scale; all method saturation counts shown' if diagnostic else 'Absolute common physical scale [0,20]; unavailable panels retained'),fontsize=13);save(fig,('diagnostic_row_scale' if diagnostic else 'absolute_scale')+f'_page{page+1:02d}')
    fig,axes=plt.subplots(1,2,figsize=(13,5),layout='constrained')
    for method in METHODS:
        rr=[r for r in rows if r['method']==method and r['view']=='primary'];vals=np.sort([r['NRMSE'] for r in rr if r.get('NRMSE') is not None]);axes[0].step(vals,np.arange(1,len(vals)+1)/48,where='post',label=f'{method} ({len(vals)}/48)');points=[r for r in rr if r.get('NRMSE') is not None and r.get('external_child_CPU_seconds') is not None];axes[1].scatter([r['external_child_CPU_seconds'] for r in points],[r['NRMSE'] for r in points],label=method,s=18,alpha=.7)
    axes[0].set(xlabel='Primary amplitude NRMSE, full range',ylabel='Fraction of all48');axes[0].legend();axes[1].set(xlabel='Actual whole-child i9 CPU (s)',ylabel='Primary amplitude NRMSE');axes[1].legend();fig.suptitle('Existing48 follow-up: actual paid runs, reused16 A+ +new32 A+');save(fig,'primary_ecdf_actual_cpu_scatter')
    fig,axes=plt.subplots(1,2,figsize=(13,5),layout='constrained')
    for method in METHODS:
        rr=[r for r in costrows if r['method']==method];vals=[r['external_child_CPU_seconds'] for r in rr if r.get('external_child_CPU_seconds') is not None];axes[0].bar(method,float(np.mean(vals)) if vals else np.nan)
        yy=[]
        for prefix in PREFIXES:
            v=[r['NRMSE'] for r in prefixrows if r['method']==method and r['CPU_prefix']==prefix and r.get('NRMSE') is not None];yy.append(float(np.mean(v)) if v else np.nan)
        axes[1].plot(PREFIXES,yy,marker='o',label=method)
    axes[0].set(ylabel='Mean actual whole-child i9 CPU(s)',xlabel='48 tasks/method; missing costs not imputed');axes[1].set(xlabel='Original process CPU prefix (s)',ylabel='Mean amplitude NRMSE of public min-paid-Phi incumbent');axes[1].legend();fig.suptitle('Costs distinct from paid-prefix selection; full denominator tables retained');save(fig,'actual_cpu_and_paid_prefix_curves')
    return artifacts


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--supplement-root',type=Path);p.add_argument('--authorization',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--no-figures',action='store_true');a=p.parse_args();root=a.root.resolve();sup=(a.supplement_root or root/'supplement_aplus48').resolve();auth,seals,seal_locations,cases=authorize(root,sup,a.authorization)
    # All guards above are stdlib-only and precede score/material/truth reads.
    import numpy as np
    from score_results import material_metrics,load_truth,truth_spec,physical_counts
    out=a.out.resolve()
    if out==root or out==sup or any(out.is_relative_to(x) for x in [root/'RESULTS',root/'sealed_truth',root/'queues',sup/'RESULTS',sup/'queues']):raise ValueError('New analysis output directory required')
    out.mkdir(parents=True,exist_ok=False);runindex=[];fairrows=[];rows=[];prefixrows=[];costrows=[];historyrows=[];failures=[];pictures={};truths={};sources=[str(a.authorization.resolve())]+list(seal_locations.values());truthsources=[]
    corecfg=seals['confirmation_core48']['manifest_snapshot'];oldcfg=seals['confirmation_mechanism16']['manifest_snapshot'];extcfg=seals['aplus48_extension32']['manifest_snapshot'];oldtasks={t['case_id']:t for t in oldcfg['tasks'] if t['method']=='a_plus'};newtasks={t['case_id']:t for t in extcfg['tasks']};coretasks={(t['case_id'],t['method']):t for t in corecfg['tasks']};caps=float(corecfg['protocol']['CPU_cap']);lower,upper=map(float,corecfg['protocol']['material_box'])
    queues={};queueindices={};queuelocations={}
    for campaign,base in [('confirmation_core48',root),('confirmation_mechanism16',root),('aplus48_extension32',sup)]:
        f=base/'queues'/campaign/'RECEIPT.json';q=optional(f,{});queues[campaign]={j['task_id']:j for j in q.get('jobs',[])};queueindices[campaign]={j['task_id']:i for i,j in enumerate(q.get('jobs',[]))};queuelocations[campaign]=str(f.resolve());sources.append(str(f))
    # Reuse baseline metric outputs only if their scoring receipts identify the
    # same original seal and unchanged original candidate selection is checked.
    scored=root/'analysis/confirmation_core48';sr=read(scored/'SCORING_RECEIPT.json')
    if sr.get('campaign')!='confirmation_core48' or Path(sr['authorized_seal']).resolve()!=Path(seal_locations['confirmation_core48']):raise ValueError('Baseline score source uses another prediction seal')
    rulespath=sup/'configs/comparison_rules.json';rules=read(rulespath)
    if rules.get('schema')!='aplus48_public_comparison_rules_v1' or rules.get('methods')!=list(METHODS) or rules.get('case_ids')!=cases:raise ValueError('Registered public comparison rule cohort/method mismatch')
    sources.append(str(rulespath));rule_stat=rulespath.stat();rule_snapshot=dict(path=str(rulespath),bytes=rule_stat.st_size,mtime_ns=rule_stat.st_mtime_ns)
    publicpath=root/'data_public/fresh48/MANIFEST.json';sources.append(str(publicpath));public_cases={r['case_id']:r for r in read(publicpath)['case_files']}
    # First metadata pass: ALL FIVE actual job costs determine each limit;
    # no metric/truth file has been opened for the new comparison at this point.
    publiclimits={};budgetlimits=[];preselected_budgets={};public_budget_selections=[];public_histories={}
    def nonnegative(value):return isinstance(value,(int,float)) and not isinstance(value,bool) and math.isfinite(value) and value>=0
    def actual_lu(cost):
        cost=cost or{}
        return cost.get('physical_LU') if 'physical_LU' in cost else cost.get('LU')
    for cid in cases:
        cc=[]
        for method in METHODS:
            reused=cid in oldtasks;campaign='confirmation_core48' if method!='a_plus' else'confirmation_mechanism16' if reused else'aplus48_extension32';base=sup if campaign=='aplus48_extension32' else root;task=coretasks[(cid,method)] if method!='a_plus' else(oldtasks[cid] if reused else newtasks[cid]);folder=base/'RESULTS'/campaign/task['task_id'];receipt=optional(folder/'RECEIPT.json',{});job=queues[campaign].get(task['task_id'],{});external=(job.get('external') or{}).get('child_process_CPU_seconds',job.get('child_process_CPU_seconds'));lu=actual_lu(receipt.get('physical_cost',{}));public_histories[(cid,method)]=optional(folder/'HISTORY.json',[]);cc.append(dict(method=method,source_campaign=campaign,task_id=task['task_id'],external_whole_child_CPU_seconds=external,whole_job_physical_LU=lu))
        cpuvals=[r['external_whole_child_CPU_seconds'] for r in cc];luvals=[r['whole_job_physical_LU'] for r in cc];cpu_limit=min([120.]+cpuvals) if all(nonnegative(x) for x in cpuvals) else None;lu_limit=min(luvals) if all(nonnegative(x) for x in luvals) else None;publiclimits[cid]=dict(common_actual_CPU=cpu_limit,common_physical_LU=lu_limit)
        for item in cc:
            method=item['method'];history=public_histories[(cid,method)];valid=[h for h in history if h.get('legal') and h.get('eligible') and h['CPU_seconds']<=caps]
            pools={'common_actual_CPU':[h for h in valid if cpu_limit is not None and h['CPU_seconds']<=cpu_limit],'common_physical_LU':[h for h in valid if lu_limit is not None and nonnegative(actual_lu((h.get('physics_cost') or{}))) and actual_lu((h.get('physics_cost') or{}))<=lu_limit]}
            choices={policy:min(pool,key=lambda h:h['Phi']) if pool else None for policy,pool in pools.items()};preselected_budgets[(cid,method)]={policy:h['candidate'] if h else None for policy,h in choices.items()}
            for policy,h in choices.items():public_budget_selections.append(dict(case_id=cid,method=method,policy=policy,limit=publiclimits[cid][policy],candidate=h['candidate'] if h else None,Phi=h['Phi'] if h else None,CPU_seconds=h['CPU_seconds'] if h else None,physical_LU=actual_lu((h.get('physics_cost') or{})) if h else None,legal=h.get('legal') if h else None,eligible=h.get('eligible') if h else None,GT_fields_read_for_selection=False))
        for r in cc:budgetlimits.append(dict(case_id=cid,aplus_source_split='reused16' if cid in oldtasks else'new32',**r,common_actual_CPU=cpu_limit,common_physical_LU=lu_limit,CPU_limit_available=cpu_limit is not None,physical_LU_limit_available=lu_limit is not None,interpretation='Retrospective paid prefix estimates, not actual stopped wall; LU is physics-only diagnostic and does not equalize projector/FFT/setup'))
    (out/'PUBLIC_COMMON_BUDGET_SELECTIONS_BEFORE_METRICS.json').write_text(json.dumps(dict(scope='FOLLOWUP_EXISTING48',comparison_rules_snapshot=rule_snapshot,GT_fields_read=False,selections=public_budget_selections),indent=2,allow_nan=False)+'\n')
    retained_fields=set(METRICS)|{'available','availability_reason','material_min','material_max','material_legal','IoU_rule','IoU_threshold','amplitude','candidate','selection_candidate','eligible','legal','CPU_seconds','phase','noise_feasible','noise_certificate_eligible','noise_qualification','crossgrid','fallback_to_primary'}
    oldrows=read(scored/'SCORES.json');oldprefix=read(scored/'PAID_PREFIX_SCORES.json');sources +=[str(scored/n) for n in ['SCORING_RECEIPT.json','SCORES.json','PAID_PREFIX_SCORES.json']];baseline={(r['case_id'],r['method'],r['view']):r for r in oldrows};baselineprefix={(r['case_id'],r['method'],float(r['CPU_prefix'])):r for r in oldprefix}
    for cid in cases:
        spec=truth_spec(root,coretasks[(cid,'g19')])
        try:truths[cid]=load_truth(spec)
        except Exception as ex:truths[cid]=None;failures.append(dict(case_id=cid,method='GT',issue='TRUTH_UNAVAILABLE',error=repr(ex)))
        truthsour=dict(case_id=cid,path=spec.get('path'),scope='authorized offline coarse truth, never selection');truthsources.append(truthsour)
        for method in METHODS:
            reused=cid in oldtasks;campaign='confirmation_core48' if method!='a_plus' else'confirmation_mechanism16' if reused else'aplus48_extension32';base=sup if campaign=='aplus48_extension32' else root;task=coretasks[(cid,method)] if method!='a_plus' else(oldtasks[cid] if reused else newtasks[cid]);folder=base/'RESULTS'/campaign/task['task_id'];rp=folder/'RECEIPT.json';sources.extend(str(folder/n) for n in ['RECEIPT.json','HISTORY.json','PREFIXES.json','CANDIDATES.npz','initial.json','corrected.json','AUXILIARY_EXECUTION.json']);receipt=optional(rp,{});history=optional(folder/'HISTORY.json',[]);job=queues[campaign].get(task['task_id'],{});ext=(job.get('external') or{}).get('child_process_CPU_seconds',job.get('child_process_CPU_seconds'));common=dict(case_id=cid,object_id=cid,method=method,source_task_id=task['task_id'],parent_campaign=campaign,source_result_path=str(folder),aplus_source_split='reused16' if reused else'new32',scope='FOLLOWUP_EXISTING48',host='i9_remote',run_status=receipt.get('status','MISSING_RECEIPT'),partial_failed_run=receipt.get('status')=='FAILED',external_child_CPU_seconds=ext)
            pc=public_cases[cid];common.update(source_count=pc['sources'],receiver_count=pc['receivers'],receiver_aperture=pc['receiver_aperture'],relative_noise_RMS=pc['relative_noise_RMS'],forward_side=pc['forward_side'],model_mismatch_status=pc['model_mismatch_status'],acquisition_noise_condition=f"{pc['sources']}sources/{pc['receiver_aperture']}/{pc['relative_noise_RMS']}")
            if receipt and (receipt.get('task')!=task or receipt.get('protocol')!=corecfg['protocol']):raise ValueError('Receipt task/protocol mismatch: '+str(rp))
            if not receipt or receipt.get('status')=='FAILED':failures.append(dict(**common,issue='FAILED_OR_MISSING_RUN',error=receipt.get('error')))
            counts=physical_counts(receipt.get('physical_cost') or{});sourcekeys=counts.pop('physical_cost_source_keys');counts['physical_LU_source_key']=sourcekeys.get('LU');counts['physical_RHS_source_keys']='|'.join(sourcekeys.get('complex_RHS',[]));aux=receipt.get('auxiliary_receipt') or optional(folder/'AUXILIARY_EXECUTION.json',{});ops=aux.get('operation_counts',{});cost=dict(common,**counts,internal_CPU_seconds_reference=receipt.get('whole_process_CPU_seconds'),auxiliary_FFT_reference=aux.get('FFT'),Green_FFT_2D_transforms=ops.get('Green_FFT_2D_transforms'),Green_IFFT_2D_transforms=ops.get('Green_IFFT_2D_transforms'),coefficient_FFT_2D_transforms=ops.get('coefficient_FFT_2D_transforms'),coefficient_IFFT_2D_transforms=ops.get('coefficient_IFFT_2D_transforms'),cost_role='new_extension_job' if campaign=='aplus48_extension32' else'original_paid_run_reference_only');costrows.append(cost)
            eligible=[h for h in history if h.get('eligible') and h.get('legal') and float(h['CPU_seconds'])<=caps];byindex={int(h['candidate']):h for h in eligible}
            def idx(name):return optional(folder/(name+'.json'),{}).get('candidate')
            indices={'initial':idx('initial'),'raw':receipt.get('raw_candidate'),'corrected':idx('corrected'),'primary':receipt.get('primary_candidate'),'secondary':receipt.get('secondary_candidate')};secondary_fallback=indices['secondary'] not in byindex
            if secondary_fallback:indices['secondary']=indices['primary']
            chosenprefix={point:min([h for h in eligible if h['CPU_seconds']<=point],key=lambda h:h['Phi'])['candidate'] if any(h['CPU_seconds']<=point for h in eligible) else None for point in PREFIXES}
            limits=publiclimits[cid];cpu_limit=limits['common_actual_CPU'];lu_limit=limits['common_physical_LU']
            budgetchoices=preselected_budgets[(cid,method)]
            # These are exact frozen view IDs and the same public paid-Phi
            # prefix IDs used below. No physics/GT-based re-selection here.
            allhistory={int(h['candidate']):h for h in history}
            def index_metadata(index):
                h=allhistory.get(index,{})
                return dict(candidate=index,CPU_seconds=h.get('CPU_seconds'),phase=h.get('phase'),Phi=h.get('Phi'),D=h.get('D'),R=h.get('R'),residual_over_delta=h.get('residual_over_delta'),noise_feasible=h.get('noise_feasible'),legal=h.get('legal'),eligible=h.get('eligible'),physical_LU=actual_lu((h.get('physics_cost') or{})),physical_LU_source_key='physical_LU' if 'physical_LU' in (h.get('physics_cost') or{}) else 'LU' if 'LU' in (h.get('physics_cost') or{}) else None)
            runindex.append(dict(case_id=cid,method=method,source_campaign=campaign,source_task_id=task['task_id'],source_root=str(base),run_folder=str(folder.resolve()),original_public_input_path=str((root/coretasks[(cid,'g19')]['input_file']).resolve()),execution_public_input_path=str((base/task['input_file'].replace('\\','/')).resolve()),source_history_path=str((folder/'HISTORY.json').resolve()),source_candidates_path=str((folder/'CANDIDATES.npz').resolve()),source_receipt_path=str(rp.resolve()),source_prefixes_path=str((folder/'PREFIXES.json').resolve()),aplus_source_split=common['aplus_source_split'],run_status=common['run_status'],view_candidate_ids=dict(indices),view_metadata={name:index_metadata(index) for name,index in indices.items()},original_secondary_candidate_id=receipt.get('secondary_candidate'),secondary_fallback_to_primary=secondary_fallback,prefix_candidate_ids={str(point):index for point,index in chosenprefix.items()},prefix_metadata=[dict(CPU_prefix=point,**index_metadata(index)) for point,index in chosenprefix.items()],external_cost_source=dict(queue_receipt_path=queuelocations[campaign],job_json_pointer='jobs['+str(queueindices[campaign].get(task['task_id'],'MISSING'))+']',job_index=queueindices[campaign].get(task['task_id']),external=dict(job.get('external') or{}),external_child_CPU_seconds=ext,child_wall_seconds_reference=job.get('wall_seconds'),returncode=job.get('returncode'),child_status=job.get('child_status')),selection_policy='Frozen view IDs; prefixes min original Phi over old legal/eligible audits with CPU<=prefix and totalcap. No GT or new selection criterion.',new_physics=0))
            for policy,index in budgetchoices.items():runindex[-1][policy]=dict(limit=limits[policy],**index_metadata(index),interpretation='Retrospective common actual whole-child CPU' if policy=='common_actual_CPU' else 'Physics-only LU diagnostic; projector/FFT/setup are not free',selection_rule='min original paid Phi at declared public same-object limit; no GT')
            needed=set(indices.values())|set(chosenprefix.values())|set(budgetchoices.values());needed.discard(None);materials={}
            # Baseline arrays only for actual primary galleries; no pool metrics.
            if method!='a_plus':needed={index for index in [indices['primary']]+list(budgetchoices.values()) if index is not None}
            if needed:
                try:
                    with np.load(folder/'CANDIDATES.npz',allow_pickle=False) as z:
                        saved=z['m'];materials={i:saved[int(i)].copy().ravel() for i in needed if 0<=int(i)<len(saved)}
                except Exception as ex:failures.append(dict(**common,issue='MATERIAL_FILE_UNAVAILABLE',error=repr(ex)))
            pictures[(cid,method)]=materials.get(indices['primary'])
            def scored_row(index,view):
                h=byindex.get(index);material=materials.get(index);metric=material_metrics(material,truths[cid],lower,upper) if h is not None else material_metrics(None,truths[cid],lower,upper);record=dict(common,**metric,view=view,selection_candidate=index,eligible=h is not None,legal=h.get('legal') if h else None,CPU_seconds=h.get('CPU_seconds') if h else None,phase=h.get('phase') if h else None,noise_certificate_eligible=coretasks[(cid,'g19')].get('observation_contract',{}).get('model_mismatch_status')=='MATCHED_DISCRETE_MODEL',noise_qualification=receipt.get('noise_qualification'),crossgrid=coretasks[(cid,'g19')].get('observation_contract',{}).get('forward_side')!=24)
                for key in ['D','R','Phi','residual_over_delta','noise_feasible']:record[key]=h.get(key) if h else None
                record['physical_LU_at_candidate']=actual_lu((h.get('physics_cost') or{})) if h else None
                record['physical_complex_RHS_at_candidate']=physical_counts(h.get('physics_cost') or{})['physical_complex_RHS'] if h else None
                return record
            for view in VIEWS:
                index=indices[view]
                if method=='a_plus':r=scored_row(index,view)
                else:
                    saved=baseline.get((cid,method,view));r=dict(common,**{k:v for k,v in (saved or{}).items() if k not in common and k in retained_fields},view=view)
                    if not saved:r.update(available=False,availability_reason='BASELINE_SCORE_ROW_MISSING')
                    elif saved.get('selection_candidate')!=index:raise ValueError('Baseline frozen view index mismatch')
                if view=='secondary':r['fallback_to_primary']=secondary_fallback
                rows.append(r)
            stored=optional(folder/'PREFIXES.json',[]);storedmap={float(x['CPU_cap']):x.get('candidate') for x in stored}
            for point,index in chosenprefix.items():
                if method=='a_plus':r=scored_row(index,'paid_prefix')
                else:
                    saved=baselineprefix.get((cid,method,float(point)));r=dict(common,**{k:v for k,v in (saved or{}).items() if k not in common and k in retained_fields},view='paid_prefix')
                    if not saved:r.update(available=False,availability_reason='BASELINE_PREFIX_ROW_MISSING')
                    elif saved.get('selection_candidate')!=index:raise ValueError('Baseline min-paid-Phi prefix index mismatch')
                r.update(CPU_prefix=point,stored_prefix_matches_recomputed=storedmap.get(float(point))==index,prefix_rule='min original Phi over legal old-eligible audits with CPU<=prefix<=120; no GT');prefixrows.append(r)
            for policy,index in budgetchoices.items():
                r=scored_row(index,policy);r.update(public_budget_limit=limits[policy],public_budget_limit_available=limits[policy] is not None,public_budget_unit='external whole-child CPU seconds' if policy=='common_actual_CPU' else'cumulative physical LU only',qualification='Retrospective already-paid public prefix; not a measured stopped-deployment wall; projector/FFT/setup still charged separately');fairrows.append(r)
            if method=='a_plus':
                incumbent=None
                for h in history:
                    paid=h.get('eligible') and h.get('legal') and h['CPU_seconds']<=caps
                    if paid and (incumbent is None or h['Phi']<incumbent['Phi']):incumbent=h
                    selected={k:h.get(k) for k in ['candidate','CPU_seconds','wall_seconds','phase','Phi','D','R','relative_residual','residual_over_delta','noise_feasible','eligible','legal','deadline','lower_active','upper_active']};selected.update(common,incumbent_candidate=incumbent['candidate'] if incumbent else None,incumbent_Phi=incumbent['Phi'] if incumbent else None,acceptance_status='unknown optimizer acceptance: paid physical audit or corrector line-search trial; never infer accepted iterate from HISTORY',no_unsaved_states_interpolated=True)
                    for key,value in (h.get('auxiliary') or{}).items():
                        if value is None or isinstance(value,(str,int,float,bool)):selected['audit_aux_'+key]=value
                    historyrows.append(selected)
    (out/'RUN_INDEX.json').write_text(json.dumps(dict(schema='Aplus48_comparison_run_index_v1',scope='FOLLOWUP_EXISTING48',root=str(root),supplement_root=str(sup),authorization=str(a.authorization.resolve()),seals=seal_locations,expected_runs=240,actual_runs=len(runindex),methods=list(METHODS),prefixes=list(PREFIXES),path_policy='All file/root fields are absolute local source references; original_public_input_path always points to unchanged R13 original observation, execution_public_input_path records the actual original/extension run root. Absolute references need explicit restore remapping.',selection_policy='Index of the exact same frozen views and recomputed original min-paid-Phi prefixes emitted by compare48; no additional choice or physics.',GT_fields_in_index=False,comparison_rules_snapshot=rule_snapshot,common_public_budget_rules='common_actual_CPU=min(120,all5externalwholejobCPU);common_physical_LU=min(all5wholejobLU). Missing any value makes limit null; original minpaidPhi selection unchanged.',runs=runindex),indent=2,allow_nan=False)+'\n')
    changes,stats=paired(rows,np);summary=summarize(rows,np);fairchanges,fairstats=paired(fairrows,np,views=('common_actual_CPU','common_physical_LU'));fairsummary=summarize(fairrows,np)
    prefixsummary=[]
    for point in PREFIXES:
        for record in summarize([dict(r,view='prefix_'+str(point)) for r in prefixrows if r['CPU_prefix']==point],np):prefixsummary.append(dict(record,CPU_prefix=point))
    for name,data in [('ALL48_SCORES',rows),('PAID_PREFIX_SCORES',prefixrows),('METHOD_VIEW_SUMMARY',summary),('PREFIX_SUMMARY',prefixsummary),('PAIRED_OBJECT_CHANGES',changes),('PAIRED_BOOTSTRAP',stats),('COST_LEDGER',costrows),('APLUS_PAID_AUDIT_HISTORY',historyrows),('FAILURES',failures),('IO_ISSUES',IO_ISSUES),('TRUTH_SOURCES',truthsources)]:table(out,name,data)
    for name,data in [('PUBLIC_COMMON_BUDGET_LIMITS',budgetlimits),('COMMON_PUBLIC_BUDGET_SCORES',fairrows),('COMMON_PUBLIC_BUDGET_SUMMARY',fairsummary),('COMMON_PUBLIC_BUDGET_PAIRED_CHANGES',fairchanges),('COMMON_PUBLIC_BUDGET_BOOTSTRAP',fairstats)]:table(out,name,data)
    table(out,'PUBLIC_SUBGROUP_VIEW_SUMMARY',subgroup_summaries(rows+fairrows,np))
    subgroup_prefix=[]
    for point in PREFIXES:
        for r in subgroup_summaries([dict(x,view='prefix_'+str(point)) for x in prefixrows if x['CPU_prefix']==point],np):subgroup_prefix.append(dict(r,CPU_prefix=point))
    table(out,'PUBLIC_SUBGROUP_PREFIX_SUMMARY',subgroup_prefix)
    table(out,'PUBLIC_COHORT_FIELDS',[dict(case_id=cid,aplus_source_split='reused16' if cid in oldtasks else'new32',**{k:public_cases[cid].get(k) for k in ['sources','receivers','receiver_aperture','relative_noise_RMS','forward_side','model_mismatch_status']}) for cid in cases])
    (out/'REGISTERED_COMPARISON_RULES.json').write_text(json.dumps(rules,indent=2)+'\n')
    artifacts=figures(out,cases,rows,prefixrows,costrows,truths,pictures,np) if not a.no_figures else[]
    def total(rr,key):return sum(r[key] for r in rr if r.get(key) is not None)
    table(out,'COST_TOTALS',[dict(method=method,aplus_source_split=split,tasks=len(rr),external_CPU_known_sum=total(rr,'external_child_CPU_seconds'),external_CPU_unknown=sum(r.get('external_child_CPU_seconds') is None for r in rr),physical_LU_known_sum=total(rr,'physical_LU'),physical_complex_RHS_known_sum=total(rr,'physical_complex_RHS'),physics_unknown_tasks=sum(r.get('physical_LU') is None or r.get('physical_complex_RHS') is None for r in rr),cost_policy='Original jobs are references; only extension32 are new solver consumption. No cumulative HISTORY/auxiliary deltas/internal CPU added to whole-job totals') for method in METHODS for split in ['reused16','new32'] for rr in [[r for r in costrows if r['method']==method and r['aplus_source_split']==split]]])
    receipt=dict(status='FOLLOWUP_EXISTING48_AUTHORIZED_COMPARISON_COMPLETE',scope='FOLLOWUP_EXISTING48_NOT_FRESH_BLIND',host='Apple_local',platform=platform.platform(),CPU_seconds=time.process_time()-START_CPU,wall_seconds=time.perf_counter()-START_WALL,authorization=str(a.authorization.resolve()),seals=seal_locations,expected_material_objects=48,expected_methods=5,expected_view_rows=1200,actual_view_rows=len(rows),expected_prefix_rows=1200,actual_prefix_rows=len(prefixrows),run_index_rows=len(runindex),common_public_budget_rows=len(fairrows),comparison_rules_snapshot=rule_snapshot,new_Aplus32=32,reused_Aplus16=16,new_Maxwell_LU=0,new_RHS=0,new_optimization=0,baseline_metrics_source='Unchanged authorized core48 score tables; candidate indices checked against frozen receipts/public paid history',Aplus_metrics_source='Existing material_metrics/load_truth/truth_spec; only requested view/prefix candidates scored',artifacts=artifacts,input_source_files=sorted(set(sources)),failure_records=len(failures),online_selections_changed=False,timer_scope='Whole comparison process including authorization, imports, metric attachment, CSV/figure writes; final receipt/exit excluded. Original online CPU and new extension CPU remain separate')
    (out/'COMPARISON_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({k:receipt[k] for k in ['status','actual_view_rows','actual_prefix_rows','CPU_seconds','new_Maxwell_LU']}))
if __name__=='__main__':main()
