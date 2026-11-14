"""Authorized existing-object A+ version/NN factorial; selected states only, zero forwards.

Do not execute before the parent authorizes all four prediction seals and sets
``authorized_stack_scoring: true`` in the independent comparison authorization.
"""
import os,time
START_CPU=time.process_time();START_WALL=time.perf_counter()
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
import argparse,csv,itertools,json,math,platform,sys
from pathlib import Path
from datetime import datetime,timezone
from compare48 import optional,read,table,IO_ISSUES,METRICS
SOLVERS=('cie','dual_cie','g19','dr','a_plus_mm_r13','a_plus_arc_r7')
VIEWS=('initial','raw','corrected','primary','secondary')
PREFIXES=(15,30,60,90,120)
SCOPE='FIXED_EXISTING_VARIANTS_AND_FROZEN_NN_FOLLOWUP'
OPTIONS=dict(aplus_mode='mobius',aplus_arc_candidates=True,aplus_inherited_gtol=True,aplus_sweeps=2,aplus_step_tol=1e-8)

def finite(x):return isinstance(x,(int,float)) and not isinstance(x,bool) and math.isfinite(x) and x>=0

def authorize(root,sup,stack,path):
    """No NumPy, scored metrics, candidate arrays or truth reads in this function."""
    auth=read(path)
    if auth.get('authorized_stack_scoring') is not True:raise PermissionError('Parent authorized_stack_scoring=true required')
    if auth.get('scope')!=SCOPE:raise PermissionError('Existing-object version/NN follow-up scope required; not fresh blind')
    for key,expected in [('root',root),('supplement_root',sup),('stack_root',stack)]:
        if auth.get(key) and Path(auth[key]).resolve()!=expected:raise ValueError('Authorization '+key+' mismatch')
    sys.path.insert(0,str(root/'src'));from seal_results import validate_seal
    defaults={'confirmation_core48':root/'audit/confirmation_core48_seal.json','confirmation_mechanism16':root/'audit/confirmation_mechanism16_seal.json','aplus48_extension32':sup/'audit/aplus48_extension32_seal.json','aplus_stack112':sup/'audit/aplus_stack112_seal.json'}
    bases={'confirmation_core48':root,'confirmation_mechanism16':root,'aplus48_extension32':sup,'aplus_stack112':stack};seals={};locations={}
    for campaign,default in defaults.items():
        p=Path(auth.get('seals',{}).get(campaign,default));p=p if p.is_absolute() else root/p
        seal=read(p);validate_seal(seal,require_authorization=True)
        if seal['campaign']!=campaign or Path(seal['root']).resolve()!=bases[campaign]:raise ValueError('Seal root/campaign mismatch '+campaign)
        seals[campaign]=seal;locations[campaign]=str(p.resolve())
    cases=[f'r13_fresh_{i:04d}' for i in range(48)];nncases=cases[::3]
    configs={k:v['manifest_snapshot'] for k,v in seals.items()};protocol=configs['confirmation_core48']['protocol']
    if any(c['protocol']!=protocol for c in configs.values()):raise ValueError('Original numerical/budget/selector protocol changed')
    cfg=read(sup/'configs/aplus_stack112.json')
    if cfg!=configs['aplus_stack112']:raise ValueError('Stack seal differs from registered config')
    if cfg.get('scope')!=SCOPE or cfg.get('frozen_readout_options')!=OPTIONS:raise ValueError('Exact R7 final readout flags/options required')
    r7=read(root.parent/'research_max_struct_r7/configs/aplus_final.json')
    if r7.get('aplus_options')!=OPTIONS:raise ValueError('Registered readout options differ from R7 final configuration')
    expected={(cid,'a_plus_arc_r7') for cid in cases}|{(cid,m) for cid in nncases for m in ('nn_cie','nn_dual_cie','nn_a_plus','nn_a_plus_arc_r7')}
    ts=cfg['tasks']
    if len(ts)!=112 or {(t['case_id'],t['method']) for t in ts}!=expected:raise ValueError('Stack112 factorial membership changed')
    for t in ts:
        alg=t['method'].removeprefix('nn_');alg='a_plus' if alg=='a_plus_arc_r7' else alg
        if t.get('algorithm')!=alg or t['task_id']!=t['case_id']+'__'+t['method']:raise ValueError('Stack112 algorithm/task mapping changed')
    old=configs['confirmation_mechanism16']['tasks'];ext=configs['aplus48_extension32']['tasks']
    if sorted(t['case_id'] for t in old if t['method']=='a_plus')!=nncases or sorted(t['case_id'] for t in ext)!=sorted(set(cases)-set(nncases)):raise ValueError('Original MM16+32 membership changed')
    if any(t['method']!='a_plus' for t in ext):raise ValueError('Unexpected extension arm')
    corebycase={t['case_id']:t for t in configs['confirmation_core48']['tasks'] if t['method']=='g19'}
    for t in ts:
        expected_directory='fresh16_nn' if t['method'].startswith('nn_') else'fresh48'
        normalized=t['input_file'].replace('\\','/')
        if not normalized.endswith('/data_public/'+expected_directory+'/'+t['case_id']+'.npz'):raise ValueError('Original observation/NN input reference changed')
        if t.get('observation_contract') and t['observation_contract']!=corebycase[t['case_id']].get('observation_contract'):raise ValueError('Public observation contract changed')
    return auth,configs,locations,bases,cases,nncases


def summary(rows,np):
    output=[]
    groups={}
    for r in rows:
        # Full48 uniform and common16 factorial are deliberately separate.
        cohorts=['uniform48'] if r['initialization']=='uniform' and r.get('budget_group')!='factorial_twelve' else[]
        if r['case_id'] in r['_nncases']:cohorts.append('factorial16')
        for cohort in cohorts:
            key=(cohort,r['solver'],r['initialization'],r['view'],r.get('CPU_prefix'),r.get('budget_group'))
            groups.setdefault(key,[]).append(r)
    for (cohort,solver,init,view,prefix,budget),rr in groups.items():
        z=dict(cohort=cohort,solver=solver,initialization=init,view=view,CPU_prefix=prefix,budget_group=budget,expected_cases=len(rr),available_cases=sum(bool(r.get('available')) for r in rr),failed_runs=sum(r['run_status']=='FAILED' for r in rr),ball_known=sum(isinstance(r.get('noise_feasible'),bool) for r in rr),ball_inside=sum(r.get('noise_feasible') is True for r in rr),exposure='EXISTING_EXPOSED_OBJECTS_VERSION_AND_NN_CONTROLS_NOT_NEW_BLIND')
        for metric in METRICS:
            v=[r[metric] for r in rr if isinstance(r.get(metric),(int,float)) and math.isfinite(r[metric])]
            z.update({metric+'_n':len(v),metric+'_mean':float(np.mean(v)) if v else None,metric+'_median':float(np.median(v)) if v else None})
        output.append(z)
    return output


def paired_deltas(rows,np,kind,reference_solver="a_plus_mm_r13"):
    lookup={(r['case_id'],r['solver'],r['initialization'],r['view'],r.get('CPU_prefix'),r.get('budget_group')):r for r in rows};changes=[]
    for key,a in lookup.items():
        cid,solver,init,view,prefix,budget=key
        if kind=='NN_minus_uniform':
            if init!='NN':continue
            b=lookup.get((cid,solver,'uniform',view,prefix,budget));label=kind
        elif kind=='corrected_minus_raw':
            if view!='corrected':continue
            b=lookup.get((cid,solver,init,'raw',prefix,budget));label=kind
        else:
            if solver!='a_plus_arc_r7':continue
            b=lookup.get((cid,reference_solver,init,view,prefix,budget));label='R7_arc_guard_minus_'+reference_solver
        z=dict(case_id=cid,solver=solver,initialization=init,view=view,CPU_prefix=prefix,budget_group=budget,contrast=label,paired_available=bool(a.get('available') and b and b.get('available')),reference_solver=b.get('solver') if b else None,reference_initialization=b.get('initialization') if b else None)
        for metric in METRICS+('CPU_seconds','external_child_CPU_seconds','physical_LU','physical_complex_RHS','physical_LU_at_candidate','physical_complex_RHS_at_candidate'):
            x,y=a.get(metric),(b or{}).get(metric);z['candidate_'+metric]=x;z['reference_'+metric]=y;z['delta_'+metric]=x-y if z['paired_available'] and x is not None and y is not None else None
        changes.append(z)
    groups={}
    for r in changes:
        key=tuple(r.get(k) for k in ('contrast','solver','initialization','view','CPU_prefix','budget_group'));groups.setdefault(key,[]).append(r)
    stats=[];rng=np.random.default_rng(131313)
    for key,rr in groups.items():
        z=dict(zip(('contrast','solver','initialization','view','CPU_prefix','budget_group'),key));z.update(expected_pairs=len(rr),paired_available=sum(r['paired_available'] for r in rr),bootstrap_seed=131313,bootstrap_replicates=4000,unit='one material object; descriptive unadjusted bootstrap; missing pairs not imputed; full expected denominator retained')
        for metric in ('NRMSE','RMSE','IoU','Phi','D','R'):
            v=np.asarray([r['delta_'+metric] for r in rr if r.get('delta_'+metric) is not None],float);bs=np.mean(v[rng.integers(0,len(v),(4000,len(v)))],axis=1) if len(v) else[]
            z.update({metric+'_n':len(v),metric+'_mean_delta':float(v.mean()) if len(v) else None,metric+'_median_delta':float(np.median(v)) if len(v) else None,metric+'_bootstrap95_low':float(np.quantile(bs,.025)) if len(v) else None,metric+'_bootstrap95_high':float(np.quantile(bs,.975)) if len(v) else None})
        stats.append(z)
    return changes,stats


def plots(out,cases,nncases,rows,prefixes,truths,pictures,np,families):
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    dest=out/'FIGURES';dest.mkdir();artifacts=[];lookup={(r['case_id'],r['solver'],r['initialization'],r['view']):r for r in rows}
    def save(fig,name):
        for ext in ('png','pdf'):
            p=dest/(name+'.'+ext);fig.savefig(p,dpi=160,bbox_inches='tight');artifacts.append(str(p.relative_to(out)))
        plt.close(fig)
    ordered=sorted(cases,key=lambda c:(families.get(c,''),c))
    for start in range(0,48,6):
        chunk=ordered[start:start+6];fig,axes=plt.subplots(len(chunk),3,figsize=(10,2.9*len(chunk)),squeeze=False,layout='constrained')
        for i,cid in enumerate(chunk):
            arrays=[(truths.get(cid) or{}).get('m'),pictures.get((cid,'a_plus_mm_r13','uniform')),pictures.get((cid,'a_plus_arc_r7','uniform'))]
            for j,array in enumerate(arrays):
                ax=axes[i,j];ax.set_xticks([]);ax.set_yticks([])
                if i==0:ax.set_title(('GT','R13 MM / no guard','R7 final readout flags')[j])
                if j==0:ax.set_ylabel(cid+'\n'+families.get(cid,'family unavailable'))
                if array is None:ax.text(.5,.5,'UNAVAILABLE',ha='center',va='center');continue
                ax.imshow(np.asarray(array).reshape(24,24),origin='lower',cmap='viridis',norm=Normalize(0,20),interpolation='nearest')
                if j:
                    r=lookup.get((cid,('a_plus_mm_r13','a_plus_arc_r7')[j-1],'uniform','primary'),{});v=r.get('NRMSE');ax.set_xlabel(f'NRMSE {v:.4f}' if v is not None else 'Metric unavailable')
        fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0,20),cmap='viridis'),ax=axes.ravel().tolist(),fraction=.02,pad=.01,label='Physical chi [0,20]');fig.suptitle('Existing48 version controls; full IDs and all failures retained');save(fig,f'arc_vs_MM_all48_page{start//6+1:02d}')
    fig,axes=plt.subplots(2,3,figsize=(14,8),layout='constrained')
    for ax,solver in zip(axes.ravel(),SOLVERS):
        for init in ('uniform','NN'):
            yy=[]
            for point in PREFIXES:
                v=[r['NRMSE'] for r in prefixes if r['solver']==solver and r['initialization']==init and r['case_id'] in nncases and r['CPU_prefix']==point and r.get('NRMSE') is not None];yy.append(float(np.mean(v)) if v else np.nan)
            ax.plot(PREFIXES,yy,marker='o',label=init)
        ax.set(title=solver,xlabel='Paid process CPU prefix (s)',ylabel='Mean amplitude NRMSE');ax.legend()
    fig.suptitle('Same predefined16: frozen NN versus uniform; no truth selection');save(fig,'NN_paired16_paid_CPU_curves')
    fig,axes=plt.subplots(2,3,figsize=(14,8),layout='constrained')
    for ax,solver in zip(axes.ravel(),SOLVERS):
        xx=[];yy=[]
        for cid in nncases:
            a=lookup.get((cid,solver,'uniform','primary'),{}).get('NRMSE');b=lookup.get((cid,solver,'NN','primary'),{}).get('NRMSE')
            if a is not None and b is not None:xx.append(a);yy.append(b)
        ax.scatter(xx,yy);lim=max(xx+yy+[.01]);ax.plot([0,lim],[0,lim],color='gray',ls=':');ax.set(title=f'{solver} ({len(xx)}/16)',xlabel='uniform primary NRMSE',ylabel='NN primary NRMSE')
    fig.suptitle('Within-object paired outcomes; training prior additional');save(fig,'NN_paired16_primary_scatter')
    for start in range(0,len(nncases),2):
        chunk=nncases[start:start+2];fig,axes=plt.subplots(2*len(chunk),7,figsize=(20,5.4*len(chunk)),squeeze=False,layout='constrained')
        for ci,cid in enumerate(chunk):
            for ii,init in enumerate(('uniform','NN')):
                i=2*ci+ii;first=(truths.get(cid) or{}).get('m') if init=='uniform' else pictures.get((cid,'initial','NN'))
                arrays=[first]+[pictures.get((cid,solver,init)) for solver in SOLVERS]
                for j,array in enumerate(arrays):
                    ax=axes[i,j];ax.set_xticks([]);ax.set_yticks([])
                    if j==0:ax.set_ylabel(cid+'\n'+init);ax.set_title('GT' if init=='uniform' else 'Frozen NN initial material')
                    elif i==0:ax.set_title(SOLVERS[j-1])
                    if array is None:ax.text(.5,.5,'UNAVAILABLE',ha='center',va='center');continue
                    ax.imshow(np.asarray(array).reshape(24,24),origin='lower',cmap='viridis',norm=Normalize(0,20),interpolation='nearest')
                    if j:
                        rr=lookup.get((cid,SOLVERS[j-1],init,'primary'),{});v=rr.get('NRMSE');ax.set_xlabel(f'NRMSE {v:.4f}' if v is not None else 'Metric unavailable')
        fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0,20),cmap='viridis'),ax=axes.ravel().tolist(),fraction=.012,pad=.005,label='Physical chi [0,20]')
        fig.suptitle('All predefined16 NN interfaces; same object/data, paid min-Phi primary; failures retained',fontsize=12)
        save(fig,f'NN_interfaces_all16_page{start//2+1:02d}')
    return artifacts


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--supplement-root',type=Path);p.add_argument('--stack-root',type=Path);p.add_argument('--authorization',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--baseline-analysis',type=Path);p.add_argument('--family-labels',type=Path);p.add_argument('--include-direct-reference',action='store_true');p.add_argument('--no-figures',action='store_true');a=p.parse_args()
    root=a.root.resolve();sup=(a.supplement_root or root/'supplement_aplus48').resolve();stack=(a.stack_root or sup/'stack_run').resolve();auth,configs,sealloc,bases,cases,nncases=authorize(root,sup,stack,a.authorization)
    # All four seals and exact factorial/options validated before scientific imports.
    import numpy as np
    from score_results import material_metrics,load_truth,truth_spec,physical_counts
    out=a.out.resolve()
    if not out.is_relative_to(sup/'analysis') or any(out.is_relative_to(x) for x in (root/'RESULTS',sup/'RESULTS',stack/'RESULTS')):raise ValueError('Use a new directory below supplement_root/analysis')
    out.mkdir(parents=True,exist_ok=False);cfg=configs['confirmation_core48'];cap=float(cfg['protocol']['CPU_cap']);lower,upper=map(float,cfg['protocol']['material_box']);failures=[];runs=[];publicselections=[];sources=[str(a.authorization.resolve())]+list(sealloc.values());cache={};truths={};pictures={}
    public={r['case_id']:r for r in read(root/'data_public/fresh48/MANIFEST.json')['case_files']}
    tasks={campaign:{(t['case_id'],t['method']):t for t in config['tasks']} for campaign,config in configs.items()};queues={};jobs={}
    for campaign,base in bases.items():
        qp=base/'queues'/campaign/'RECEIPT.json';queues[campaign]=optional(qp,{});jobs[campaign]={j['task_id']:j for j in queues[campaign].get('jobs',[])};sources.append(str(qp))
    def register(cid,solver,init):
        if init=='NN':
            method='nn_'+('a_plus' if solver=='a_plus_mm_r13' else solver);campaign='confirmation_mechanism16' if solver in ('g19','dr','direct') else'aplus_stack112'
        elif solver=='a_plus_arc_r7':method=solver;campaign='aplus_stack112'
        elif solver=='a_plus_mm_r13':method='a_plus';campaign='confirmation_mechanism16' if cid in nncases else'aplus48_extension32'
        else:method=solver;campaign='confirmation_core48'
        task=tasks[campaign].get((cid,method))
        if task is None:raise ValueError('Missing registered task '+cid+'/'+method)
        base=bases[campaign];folder=base/'RESULTS'/campaign/task['task_id'];rp=folder/'RECEIPT.json';receipt=optional(rp,{});history=optional(folder/'HISTORY.json',[]);job=jobs[campaign].get(task['task_id'],{});ext=(job.get('external') or{}).get('child_process_CPU_seconds',job.get('child_process_CPU_seconds'))
        if receipt and (receipt.get('task')!=task or receipt.get('protocol')!=cfg['protocol']):raise ValueError('Receipt task/protocol mismatch '+str(rp))
        eligible=[h for h in history if h.get('eligible') and h.get('legal') and finite(h.get('CPU_seconds')) and h['CPU_seconds']<=cap and finite(h.get('Phi'))];byindex={int(h['candidate']):h for h in eligible}
        indices={'initial':optional(folder/'initial.json',{}).get('candidate'),'raw':receipt.get('raw_candidate'),'corrected':optional(folder/'corrected.json',{}).get('candidate'),'primary':receipt.get('primary_candidate'),'secondary':receipt.get('secondary_candidate')};fallback=indices['secondary'] not in byindex
        if fallback:indices['secondary']=indices['primary']
        def select(limit):
            pool=[h for h in eligible if limit is not None and h['CPU_seconds']<=limit];return min(pool,key=lambda h:(h['Phi'],int(h['candidate'])))['candidate'] if pool else None
        selections=[dict(view=v,candidate=i,CPU_prefix=None,budget_group=None,public_budget_limit=None) for v,i in indices.items()]+[dict(view='paid_prefix',candidate=select(point),CPU_prefix=point,budget_group=None,public_budget_limit=point) for point in PREFIXES]
        pc=public[cid];record=dict(case_id=cid,solver=solver,initialization=init,source_method=method,source_campaign=campaign,source_task_id=task['task_id'],source_folder=str(folder.resolve()),source_receipt=str(rp.resolve()),queue_receipt=str((base/'queues'/campaign/'RECEIPT.json').resolve()),run_status=receipt.get('status','MISSING_RECEIPT'),external_child_CPU_seconds=ext if finite(ext) else None,secondary_fallback_to_primary=fallback,version_label='R7 final readout options copied exactly; R13 protocol/current wrapper' if solver=='a_plus_arc_r7' else'R13 frozen quintic MM, no inherited gtol guard' if solver=='a_plus_mm_r13' else'original frozen R13',exposure='EXISTING_EXPOSED_OBJECTS_NOT_FRESH_BLIND',source_count=pc['sources'],receiver_count=pc['receivers'],receiver_aperture=pc['receiver_aperture'],relative_noise_RMS=pc['relative_noise_RMS'],forward_side=pc['forward_side'],model_mismatch_status=pc['model_mismatch_status'],_nncases=nncases)
        aux=receipt.get('auxiliary_receipt') or optional(folder/'AUXILIARY_EXECUTION.json',{});ops=aux.get('operation_counts') or{}
        record.update(auxiliary_sourcewise_cg=aux.get('sourcewise_cg'),initialize_physical_current=aux.get('initialize_physical_current'),aplus_arc_stationary_roots=ops.get('aplus_arc_stationary_roots'),aplus_conditional_polynomial_roots=ops.get('aplus_conditional_polynomial_roots'),aplus_inherited_gtol_checks=ops.get('aplus_inherited_gtol_checks'),auxiliary_FFT_reference=aux.get('FFT'),auxiliary_CPU_seconds_reference=aux.get('cpu_seconds'),auxiliary_status=aux.get('status'))
        counts=physical_counts(receipt.get('physical_cost') or{});record.update({k:v for k,v in counts.items() if k!='physical_cost_source_keys'});record['physical_cost_source_keys']=counts.get('physical_cost_source_keys');record['internal_CPU_seconds_reference']=receipt.get('whole_process_CPU_seconds');record['cost_role']='new_stack112_run' if campaign=='aplus_stack112' else'reused_original_run_reference_only'
        if not receipt or receipt.get('status')=='FAILED':failures.append(dict(case_id=cid,solver=solver,initialization=init,issue='FAILED_OR_MISSING_RUN',error=receipt.get('error'),partial_eligible_candidates_retained=True))
        run=dict(record=record,folder=folder,task=task,receipt=receipt,history=history,byindex=byindex,selections=selections,select=select);runs.append(run);sources.append(str(rp))
    for cid in cases:
        for solver in SOLVERS:register(cid,solver,'uniform')
    for cid in nncases:
        for solver in SOLVERS:register(cid,solver,'NN')
        if a.include_direct_reference:register(cid,'direct','NN')
    # Per-object actual CPU limits use ALL SIX, never overwrite prior-five table.
    limits=[]
    for cid in cases:
        groups=[('uniform_six', [r for r in runs if r['record']['case_id']==cid and r['record']['initialization']=='uniform'])]
        if cid in nncases:
            groups += [('NN_six',[r for r in runs if r['record']['case_id']==cid and r['record']['initialization']=='NN' and r['record']['solver'] in SOLVERS]),('factorial_twelve',[r for r in runs if r['record']['case_id']==cid and r['record']['solver'] in SOLVERS])]
        for label,rr in groups:
            values=[r['record']['external_child_CPU_seconds'] for r in rr];limit=min([120.]+values) if all(finite(x) for x in values) else None
            limits.append(dict(case_id=cid,budget_group=label,method_runs=len(rr),common_actual_CPU=limit,unknown_cost_runs=sum(x is None for x in values),rule='min(120,actual whole-child CPU of every compared run); missing any cost=>null'))
            for run in rr:run['selections'].append(dict(view='common_actual_CPU_six' if label!='factorial_twelve' else'common_actual_CPU_twelve',candidate=run['select'](limit),CPU_prefix=None,budget_group=label,public_budget_limit=limit))
    for run in runs:
        r=run['record'];allhistory={int(h['candidate']):h for h in run['history']}
        for selection in run['selections']:
            h=allhistory.get(selection['candidate'],{});publicselections.append({k:v for k,v in r.items() if not k.startswith('_')}|selection|{k:h.get(k) for k in ('Phi','D','R','CPU_seconds','eligible','legal','noise_feasible')}|dict(GT_fields_read=False))
    # Persist ALL selected view/prefix/common-budget IDs before cache metrics or GT.
    table(out,'PUBLIC_SELECTED_INDICES_BEFORE_METRICS',publicselections);table(out,'PUBLIC_SIX_AND_TWELVE_CPU_LIMITS',limits)
    (out/'PUBLIC_SELECTION_RECEIPT.json').write_text(json.dumps(dict(scope=SCOPE,GT_fields_read=False,selected_rows=len(publicselections),seals=sealloc,selection='Frozen IDs and original min-paid-Phi; all indices serialized before metric attachment'),indent=2)+'\n')
    (out/'RUN_INDEX.json').write_text(json.dumps(dict(scope=SCOPE,GT_fields_read=False,seals=sealloc,runs=[{k:v for k,v in r['record'].items() if not k.startswith('_')}|dict(selected_indices=r['selections']) for r in runs],selection_index='PUBLIC_SELECTED_INDICES_BEFORE_METRICS.json'),indent=2,allow_nan=False)+'\n')
    # Exact-index cache reuse. Never reuse a metric with another frozen selection.
    cachedir=(a.baseline_analysis or sup/'analysis/comparison48').resolve();cachepaths=[cachedir/'ALL48_SCORES.json',cachedir/'PAID_PREFIX_SCORES.json',root/'analysis/confirmation_mechanism16/SCORES.json',root/'analysis/confirmation_mechanism16/PAID_PREFIX_SCORES.json']
    for cp in cachepaths:
        for r in optional(cp,[]):
            index=r.get('selection_candidate',r.get('candidate'));method=r.get('method');cache[(r.get('case_id'),method,index)]=r
        if cp.exists():sources.append(str(cp))
    families={}
    if a.family_labels:
        if not auth.get('authorized_family_labels') or Path(auth['authorized_family_labels']).resolve()!=a.family_labels.resolve():raise PermissionError('Parent must explicitly name authorized_family_labels in comparison authorization')
        families=read(a.family_labels);sources.append(str(a.family_labels.resolve()))
        if not isinstance(families,dict) or not set(families).issubset(cases):raise ValueError('Family labels must map known case IDs to labels')
    rows=[];prefixrows=[];fairrows=[];cachehits=0
    for run in runs:
        record=run['record'];cid=record['case_id'];folder=run['folder'];byindex=run['byindex']
        if cid not in truths:
            try:truths[cid]=load_truth(truth_spec(root,tasks['confirmation_core48'][(cid,'g19')]))
            except Exception as ex:truths[cid]=None;failures.append(dict(case_id=cid,issue='TRUTH_UNAVAILABLE',error=repr(ex)))
        requested={s['candidate'] for s in run['selections'] if s['candidate'] in byindex};materials={}
        # Only indexed materials are copied; no candidate-pool GT metric search.
        if requested:
            try:
                with np.load(folder/'CANDIDATES.npz',allow_pickle=False) as z:
                    saved=z['m'];materials={i:saved[int(i)].copy().ravel() for i in requested if 0<=int(i)<len(saved)}
            except Exception as ex:failures.append(dict(case_id=cid,solver=record['solver'],initialization=record['initialization'],issue='MATERIAL_UNAVAILABLE',error=repr(ex)))
        for sel in run['selections']:
            index=sel['candidate'];h=byindex.get(index);old=cache.get((cid,record['source_method'],index));material=materials.get(index)
            if old is not None and old.get('available') and h is not None and material is not None:
                # Material-only metric fields are reused. Physics metadata remains public history.
                keys=set(METRICS)-{'D','R','Phi','residual_over_delta'};keys|={'available','availability_reason','material_min','material_max','material_legal','IoU_rule','IoU_threshold','amplitude'};metric={k:v for k,v in old.items() if k in keys};cachehits+=1;metric_source='exact original candidate-index cache'
            else:metric=material_metrics(material if h else None,truths[cid],lower,upper);metric_source='existing material_metrics on requested indexed state'
            z=dict(record,**metric,**sel,selection_candidate=index,metric_source=metric_source,CPU_seconds=h.get('CPU_seconds') if h else None,phase=h.get('phase') if h else None,eligible=h is not None,legal=h.get('legal') if h else None,family=families.get(cid),noise_certificate_eligible=record['model_mismatch_status']=='MATCHED_DISCRETE_MODEL',noise_qualification=run['receipt'].get('noise_qualification'))
            for key in ('D','R','Phi','residual_over_delta','noise_feasible'):z[key]=h.get(key) if h else None
            hc=physical_counts(h.get('physics_cost') or{}) if h else{};z['physical_LU_at_candidate']=hc.get('physical_LU');z['physical_complex_RHS_at_candidate']=hc.get('physical_complex_RHS')
            if sel['view']=='paid_prefix':prefixrows.append(z)
            elif sel['view'].startswith('common_actual_CPU'):fairrows.append(z)
            else:rows.append(z)
            if sel['view']=='primary':pictures[(cid,record['solver'],record['initialization'])]=material
            if sel['view']=='initial' and record['initialization']=='NN':pictures[(cid,'initial','NN')]=material
    def clean(rr):return [{k:v for k,v in r.items() if not k.startswith('_')} for r in rr]
    for name,rr in [('FACTORIAL_SCORES',rows),('PAID_PREFIX_SCORES',prefixrows),('COMMON_ACTUAL_CPU_SIX_AND_TWELVE_SCORES',fairrows),('METHOD_VIEW_SUMMARY',summary(rows,np)),('PREFIX_SUMMARY',summary(prefixrows,np)),('COMMON_ACTUAL_CPU_SIX_AND_TWELVE_SUMMARY',summary(fairrows,np)),('COST_LEDGER',[r['record'] for r in runs]),('FAILURES',failures),('IO_ISSUES',IO_ISSUES)]:table(out,name,clean(rr))
    for kind,reference in [('NN_minus_uniform',None),('corrected_minus_raw',None)]+[('arc_minus_'+ref,ref) for ref in SOLVERS if ref!='a_plus_arc_r7']:
        delta,stats=paired_deltas(rows+prefixrows+fairrows,np,kind,reference_solver=reference);table(out,kind.upper()+'_OBJECT_CHANGES',delta);table(out,kind.upper()+'_BOOTSTRAP',stats)
    subgroups=[]
    for axis in ('relative_noise_RMS','source_count','receiver_aperture','forward_side','model_mismatch_status','family'):
        for value in sorted({str(r.get(axis)) for r in rows}):
            for z in summary([r for r in rows if str(r.get(axis))==value],np):subgroups.append(dict(z,public_axis=axis,public_value=value,qualification='Descriptive bundled conditions; not isolated causal factor; family null if unavailable'))
    table(out,'SUBGROUP_SUMMARY',subgroups)
    artifacts=plots(out,cases,nncases,rows,prefixrows,truths,pictures,np,families) if not a.no_figures else[]
    totals=[]
    for role in ('new_stack112_run','reused_original_run_reference_only'):
        rr=[r['record'] for r in runs if r['record']['cost_role']==role];z=dict(cost_role=role,tasks=len(rr),policy='Whole-job costs counted once per role; do not add nested internal/history/auxiliary timers')
        for key in ('external_child_CPU_seconds','physical_LU','physical_complex_RHS'):
            values=[r[key] for r in rr if r.get(key) is not None];z[key+'_known_sum']=sum(values) if values else None;z[key+'_unknown_tasks']=len(rr)-len(values)
        totals.append(z)
    table(out,'COST_TOTALS',totals)
    receipt=dict(status='AUTHORIZED_EXISTING_STACK_COMPARISON_COMPLETE',scope=SCOPE,exposure='EXISTING_EXPOSED_OBJECTS_NOT_FRESH_BLIND',authorization=str(a.authorization.resolve()),seals=sealloc,root=str(root),supplement_root=str(sup),stack_root=str(stack),objects48=48,factorial_objects16=16,uniform_solver_runs=288,NN_solver_runs=96,direct_reference_runs=16 if a.include_direct_reference else 0,new_stack_runs=112,reused_original_runs=len(runs)-112,view_rows=len(rows),prefix_rows=len(prefixrows),common_budget_rows=len(fairrows),metric_cache_hits=cachehits,failures=len(failures),new_Maxwell_forward=0,new_Maxwell_LU=0,new_RHS=0,new_optimization=0,NN_training_performed=False,NN_inference_performed=False,historical_R8_physics_training_and_generation_Mac_CPU_seconds=402.558889,historical_cost_source='research/research_max_struct_r8/NEURAL_WARM_STARTS.md:55; additional Born/domain costs separate; not repeated consumption',NN_inference='Already paid/reused preprocessing; its receipt remains separate and must be charged in independent deployment',family_labels=str(a.family_labels.resolve()) if a.family_labels else None,online_selections_changed=False,prior_five_method_comparison_replaced=False,all_six_common_actual_CPU='min(all six whole-child actual CPU,120) per object; new table; all12 matched factorial also separate',readout_options=OPTIONS,version_claim='Exact R7 final READOUT options within R13 protocol, not full historical R7 experiment',artifacts=artifacts,input_sources=sorted(set(sources)),platform=platform.platform(),CPU_seconds=time.process_time()-START_CPU,wall_seconds=time.perf_counter()-START_WALL,timer_scope='Whole analysis including guards/imports/score attachment/plots/CSV; final receipt write/exit excluded; missing costs null')
    (out/'STACK_COMPARISON_RECEIPT.json').write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n');print(json.dumps({k:receipt[k] for k in ('status','view_rows','prefix_rows','common_budget_rows','failures','new_Maxwell_LU')}))

if __name__=='__main__':main()
