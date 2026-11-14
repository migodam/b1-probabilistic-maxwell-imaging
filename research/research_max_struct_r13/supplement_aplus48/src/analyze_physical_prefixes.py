"""Post-authorized original-Phi diagnostics for the fixed Aplus48 RUN_INDEX.

OFFLINE_DIAGNOSTIC: same finite discretization, no inversion, candidate change,
GT array read, Hessian/SVD sweep, or native-step reconstruction. Numerical
imports and ALL public/material arrays follow the explicit comparison guard.
"""
from __future__ import annotations
import os,time
START_CPU=time.process_time();START_WALL=time.perf_counter()
for _key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_key]='1'
import argparse,csv,importlib.util,json,sys
from pathlib import Path
from collections import defaultdict
sys.dont_write_bytecode=True
METHODS=('cie','dual_cie','g19','dr','a_plus')
VIEWS=('raw','corrected','primary')
PREFIXES=(15,30,60,90,120)


def read(path):return json.loads(Path(path).read_text())


def path(root,value):
    item=Path(str(value).replace('\\','/'))
    return item.resolve() if item.is_absolute() else(root/item).resolve()


def dump(target,value):
    def encode(item):
        if isinstance(item,Path):return str(item)
        if hasattr(item,'tolist'):return item.tolist()
        if hasattr(item,'item'):return item.item()
        raise TypeError(type(item).__name__)
    target.write_text(json.dumps(value,indent=2,default=encode,allow_nan=False)+'\n')


def table(out,name,rows):
    dump(out/(name+'.json'),rows)
    fields=list(dict.fromkeys(key for row in rows for key in row))
    with (out/(name+'.csv')).open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields,lineterminator='\n');writer.writeheader()
        for row in rows:writer.writerow({k:json.dumps(v) if isinstance(v,(dict,list)) else v for k,v in row.items()})


def authorize(root,sup,authorization,index_path):
    # compare48.authorize and seal_results are stdlib-only. Never call main,
    # scoring helpers, stationarity.fixed_requests, or a truth loader here.
    source=sup/'src/compare48.py'
    spec=importlib.util.spec_from_file_location('fixed_comparison_authorization',source)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    auth,seals,locations,cases=module.authorize(root,sup,authorization)
    index=read(index_path)
    if index.get('schema')!='Aplus48_comparison_run_index_v1' or index.get('scope')!='FOLLOWUP_EXISTING48':
        raise ValueError('Exact coordinated RUN_INDEX schema/scope required')
    if path(root,index['root'])!=root or path(root,index['supplement_root'])!=sup:
        raise ValueError('RUN_INDEX root/supplement mismatch')
    if path(root,index['authorization'])!=authorization.resolve():raise ValueError('RUN_INDEX authorization differs')
    if set(index.get('seals',{}))!=set(locations):raise ValueError('RUN_INDEX seal set differs')
    for name,location in locations.items():
        if path(root,index['seals'][name])!=Path(location).resolve():raise ValueError('RUN_INDEX seal location differs')
    runs=index.get('runs',[])
    if index.get('expected_runs')!=240 or index.get('actual_runs')!=240 or len(runs)!=240:
        raise ValueError('Retain exactly48cases x5methods, including missing evidence rows')
    if {(r['case_id'],r['method']) for r in runs}!={(c,m) for c in cases for m in METHODS}:
        raise ValueError('RUN_INDEX case/method set differs from frozen scope')
    if len({(r['case_id'],r['method']) for r in runs})!=len(runs):raise ValueError('Duplicate indexed run')
    taskmaps={name:{t['task_id']:t for t in seal['manifest_snapshot']['tasks']} for name,seal in seals.items()}
    original={t['case_id']:t for t in seals['confirmation_core48']['manifest_snapshot']['tasks'] if t['method']=='g19'}
    for run in runs:
        cid,method=run['case_id'],run['method'];campaign=run['source_campaign']
        expected_campaign=('confirmation_core48' if method!='a_plus' else
            'confirmation_mechanism16' if cid in cases[::3] else'aplus48_extension32')
        if campaign!=expected_campaign:raise ValueError('Indexed campaign does not match registered source split')
        task=taskmaps[campaign].get(run['source_task_id'])
        if task is None or task['case_id']!=cid or task['method']!=method:raise ValueError('Indexed task differs from frozen source manifest')
        base=sup if campaign=='aplus48_extension32' else root
        folder=(base/'RESULTS'/campaign/task['task_id']).resolve()
        checks={'source_root':base,'run_folder':folder,'source_history_path':folder/'HISTORY.json',
            'source_candidates_path':folder/'CANDIDATES.npz','source_receipt_path':folder/'RECEIPT.json',
            'source_prefixes_path':folder/'PREFIXES.json',
            'original_public_input_path':path(root,original[cid]['input_file']),
            'execution_public_input_path':path(base,task['input_file'])}
        for key,expected in checks.items():
            if not Path(run[key]).is_absolute() or path(root,run[key])!=Path(expected).resolve():
                raise ValueError('Indexed absolute source path differs: '+key)
        if set(str(k) for k in run['prefix_candidate_ids'])!={str(t) for t in PREFIXES}:
            raise ValueError('Exactly five fixed paid CPU prefixes required')
    return auth,seals,locations,cases,runs,index


def requests_for(run):
    result=[]
    for view in VIEWS:
        result.append(dict(view=view,CPU_prefix=None,candidate_id=run['view_candidate_ids'].get(view),
            saved_metadata=run.get('view_metadata',{}).get(view,{})))
    meta={float(v['CPU_prefix']):v for v in run.get('prefix_metadata',[])}
    for point in PREFIXES:
        result.append(dict(view='paid_prefix',CPU_prefix=point,
            candidate_id=run['prefix_candidate_ids'].get(str(point)),saved_metadata=meta.get(float(point),{})))
    return result


def load_public(np,source):
    with np.load(source,allow_pickle=False) as z:
        if 'components'in z and int(z['components'])!=1:raise ValueError('Existing48 fixed scope is scalar TM only')
        public={k:z[k].copy() for k in ('G','S','Y','grid_shape')}
        public['E0']=z['E0' if 'E0'in z else'E'].copy()
        for k in ('lower','upper'):
            if k in z:public[k]=float(z[k])
    if tuple(map(int,public['grid_shape']))!=(24,24):raise ValueError('Original576-material grid required')
    return public


def same_operator(np,a,b):
    return all(a[k].shape==b[k].shape and np.array_equal(a[k],b[k]) for k in ('G','E0','S','Y','grid_shape'))


def evaluate_material(np,model,m,Y,lower,upper):
    """Adapted from R13 analyze_stationarity.evaluate_material; no fixed400 call."""
    value,g=model.data_value_gradient(m,Y) # beta=None is the ORIGINAL Phi.
    state=model.state(m)                 # matching retained LU/state, no new solve.
    prior,gr,_=model.reg.evaluate(m);gd=g-gr
    residual=state['Y']-Y;normY=float(np.linalg.norm(Y));rr=float(np.linalg.norm(residual))
    data=.5*rr*rr/(normY*normY);mapping=m-np.clip(m-g,lower,upper)
    nd,nr=float(np.linalg.norm(gd)),float(np.linalg.norm(gr))
    metrics=dict(Phi=float(value),D=data,R=float(prior),Phi_minus_sum_D_R=float(value-data-prior),
        original_data_relative_residual=rr/normY,original_gradient_L2=float(np.linalg.norm(g)),
        original_gradient_inf=float(np.linalg.norm(g,np.inf)),data_gradient_L2=nd,prior_gradient_L2=nr,
        data_gradient_inf=float(np.linalg.norm(gd,np.inf)),prior_gradient_inf=float(np.linalg.norm(gr,np.inf)),
        data_prior_gradient_cosine=float(np.dot(gd,gr)/(nd*nr)) if nd>0 and nr>0 else None,
        data_prior_cosine_undefined_reason=None if nd>0 and nr>0 else'zero data/prior gradient',
        unit_step_box_gradient_mapping_inf=float(np.linalg.norm(mapping,np.inf)),
        unit_step_box_gradient_mapping_L2=float(np.linalg.norm(mapping)),unit_step=1.,
        lower_active_exact=int(np.count_nonzero(m==lower)),upper_active_exact=int(np.count_nonzero(m==upper)),
        below_lower=int(np.count_nonzero(m<lower)),above_upper=int(np.count_nonzero(m>upper)),
        box_violation_max=float(max(np.max(lower-m,initial=0.),np.max(m-upper,initial=0.))),
        material_min=float(m.min()),material_max=float(m.max()),
        physical_state_relative_residual=float(state['state_relative_residual']),
        gradient_scale='original unscaled physical Phi; optimizer1e6 scaling not used')
    if not all(np.isfinite(v).all() for v in (g,gd,gr)) or not np.isfinite(metrics['Phi']):
        raise FloatingPointError('Nonfinite physical objective/gradient')
    return metrics,dict(gPhi=np.array(g,copy=True),gD=np.array(gd,copy=True),gR=np.array(gr,copy=True))


def transition(np,a,b,entries,kind):
    row=dict(case_id=a['case_id'],method=a['method'],transition_kind=kind,
        from_view=a['view'],to_view=b['view'],from_CPU_prefix=a['CPU_prefix'],to_CPU_prefix=b['CPU_prefix'],
        from_candidate_id=a['candidate_id'],to_candidate_id=b['candidate_id'],
        from_evaluation_id=a.get('unique_evaluation_id'),to_evaluation_id=b.get('unique_evaluation_id'),
        interpretation='paid-incumbent estimate change, NOT a native material/optimizer step',status='UNAVAILABLE')
    ea,eb=entries.get(a.get('unique_evaluation_id')),entries.get(b.get('unique_evaluation_id'))
    if not ea or not eb or ea['status']!='EVALUATED' or eb['status']!='EVALUATED':return row
    if ea['operator_group_id']!=eb['operator_group_id']:raise ValueError('Transition crosses different physical operators/priors')
    d=eb['m']-ea['m'];g=ea['vectors']['gPhi'];dn,gn=float(np.linalg.norm(d)),float(np.linalg.norm(g))
    first=float(np.dot(g,d));actual=eb['metrics']['Phi']-ea['metrics']['Phi']
    row.update(status='EVALUATED',material_change_L2=dn,material_change_inf=float(np.linalg.norm(d,np.inf)),
        same_saved_material=bool(np.array_equal(ea['m'],eb['m'])),
        step_cosine_with_minus_gPhi=float(np.dot(d,-g)/(dn*gn)) if dn>0 and gn>0 else None,
        step_cosine_undefined_reason=None if dn>0 and gn>0 else'zero incumbent change or zero physical gradient',
        actual_delta_Phi=actual,first_order_prediction_gPhi_dot_delta_m=first,
        first_order_predicted_end_Phi=ea['metrics']['Phi']+first,
        actual_minus_first_order=actual-first,
        first_order_data_part=float(np.dot(ea['vectors']['gD'],d)),
        first_order_prior_part=float(np.dot(ea['vectors']['gR'],d)),
        actual_delta_D=eb['metrics']['D']-ea['metrics']['D'],actual_delta_R=eb['metrics']['R']-ea['metrics']['R'],
        comparison_NRMSE_delta=(b.get('comparison_NRMSE')-a.get('comparison_NRMSE'))
            if b.get('comparison_NRMSE') is not None and a.get('comparison_NRMSE') is not None else None)
    return row


def summaries(np,rows,changes):
    keys=('comparison_NRMSE','comparison_RMSE','Phi','D','R','original_gradient_L2',
        'unit_step_box_gradient_mapping_inf','data_gradient_L2','prior_gradient_L2',
        'data_prior_gradient_cosine','lower_active_exact','upper_active_exact')
    def summarize(rr,extra):
        good=[r for r in rr if r['status']=='EVALUATED'];item=dict(extra,expected_objects=48,
            retained_objects=len(rr),evaluated_objects=len(good),interpretation='descriptive object summaries, not a success threshold')
        for key in keys:
            values=[r[key]for r in good if r.get(key)is not None]
            item[key+'_n']=len(values);item[key+'_mean']=float(np.mean(values))if values else None
            item[key+'_median']=float(np.median(values))if values else None
        return item
    prefix=[summarize([r for r in rows if r['method']==m and r['CPU_prefix']==t],dict(method=m,CPU_prefix=t))
        for m in METHODS for t in PREFIXES]
    views=[summarize([r for r in rows if r['method']==m and r['view']==v],dict(method=m,view=v))
        for m in METHODS for v in VIEWS]
    steps=[]
    for m in METHODS:
        for kind,left,right in [('adjacent_paid_prefix',a,b)for a,b in zip(PREFIXES[:-1],PREFIXES[1:])]+[('raw_to_corrected_selected_estimate',None,None)]:
            rr=[r for r in changes if r['method']==m and r['transition_kind']==kind and r['from_CPU_prefix']==left and r['to_CPU_prefix']==right]
            good=[r for r in rr if r['status']=='EVALUATED'];item=dict(method=m,transition_kind=kind,from_CPU_prefix=left,to_CPU_prefix=right,
                expected_objects=48,retained_objects=len(rr),evaluated_objects=len(good),zero_material_changes=sum(r['same_saved_material']for r in good),
                interpretation='paid-incumbent changes, NOT native material steps')
            for key in('actual_delta_Phi','first_order_prediction_gPhi_dot_delta_m','actual_minus_first_order',
                    'material_change_L2','step_cosine_with_minus_gPhi','comparison_NRMSE_delta'):
                values=[r[key]for r in good if r.get(key)is not None];item[key+'_n']=len(values)
                item[key+'_mean']=float(np.mean(values))if values else None;item[key+'_median']=float(np.median(values))if values else None
            steps.append(item)
    return prefix,views,steps


def plot_compact(np,out,rows,transitions):
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    dest=out/'FIGURES';dest.mkdir();artifacts=[]
    def save(fig,name):
        for suffix in ('png','pdf'):
            target=dest/(name+'.'+suffix);fig.savefig(target,dpi=160,bbox_inches='tight');artifacts.append(str(target.relative_to(out)))
        plt.close(fig)
    fig,axes=plt.subplots(1,4,figsize=(18,4.4),layout='constrained')
    for method in METHODS:
        for ax,key,label in zip(axes,('comparison_NRMSE','Phi','original_gradient_L2','unit_step_box_gradient_mapping_inf'),
                ('Existing comparison NRMSE','Original physical Phi','Original gradient L2 norm','Unit-step box mapping infinity norm')):
            values=[]
            for point in PREFIXES:
                xx=[r[key] for r in rows if r['method']==method and r['CPU_prefix']==point and r.get('status')=='EVALUATED' and r.get(key) is not None]
                values.append(float(np.mean(xx)) if xx else np.nan)
            ax.plot(PREFIXES,values,marker='o',label=method)
            ax.set(xlabel='Paid incumbent available by original process CPU (s)',ylabel=label)
    axes[0].legend(fontsize=8);fig.suptitle('OFFLINE_DIAGNOSTIC: existing48, frozen paid-min-Phi prefixes; descriptive means, missing counts in table')
    save(fig,'paid_prefix_error_Phi_mapping')
    lookup={(r['case_id'],r['method'],r['view']):r for r in rows if r['view']in VIEWS}
    fig,axes=plt.subplots(1,4,figsize=(18,4.4),layout='constrained')
    for method in METHODS:
        for ax,key,label in zip(axes,('comparison_NRMSE','Phi','original_gradient_L2','unit_step_box_gradient_mapping_inf'),
                ('Existing comparison NRMSE','Physical Phi','Original gradient L2 norm','Unit-step mapping infinity norm')):
            pairs=[(lookup[c,method,'raw'],lookup[c,method,'corrected']) for c in sorted({r['case_id']for r in rows})
                if (c,method,'raw')in lookup and(c,method,'corrected')in lookup]
            xy=[(a.get(key),b.get(key))for a,b in pairs if a.get(key)is not None and b.get(key)is not None and a['status']=='EVALUATED' and b['status']=='EVALUATED']
            ax.scatter([v[0]for v in xy],[v[1]for v in xy],s=16,alpha=.65,label=f'{method} ({len(xy)}/48)')
            ax.set(xlabel='Raw '+label,ylabel='Corrected '+label)
    for ax in axes:
        limits=[*ax.get_xlim(),*ax.get_ylim()];lo,hi=min(limits),max(limits);ax.plot([lo,hi],[lo,hi],color='.4',lw=.7,ls='--')
    axes[0].legend(fontsize=8);fig.suptitle('Frozen raw/corrected estimates: no new correction, no stationarity/success gate')
    save(fig,'raw_corrected_error_Phi_mapping')
    fig,axes=plt.subplots(1,2,figsize=(11,4.3),layout='constrained')
    for method in METHODS:
        rr=[r for r in transitions if r['method']==method and r['transition_kind']=='adjacent_paid_prefix' and r['status']=='EVALUATED']
        cos=[r['step_cosine_with_minus_gPhi']for r in rr if r['step_cosine_with_minus_gPhi']is not None]
        axes[0].hist(cos,bins=np.linspace(-1,1,21),histtype='step',label=f'{method} ({len(cos)})')
        axes[1].scatter([r['first_order_prediction_gPhi_dot_delta_m']for r in rr],
            [r['actual_delta_Phi']for r in rr],s=12,alpha=.5,label=method)
    axes[0].set(xlabel='Cosine: incumbent change versus minus original gradient',ylabel='Adjacent-prefix changes');axes[0].legend(fontsize=8)
    axes[1].set(xlabel='gPhi at old incumbent dot material change',ylabel='Actual selected-estimate deltaPhi')
    fig.suptitle('Paid-incumbent changes are NOT native material steps; zero changes have undefined cosine')
    save(fig,'paid_incumbent_direction_and_first_order')
    return artifacts


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True);p.add_argument('--supplement-root',type=Path)
    p.add_argument('--authorization',type=Path,required=True);p.add_argument('--run-index',type=Path,required=True)
    p.add_argument('--legacy-repo',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--no-figures',action='store_true');a=p.parse_args()
    root=a.root.resolve();sup=(a.supplement_root or root/'supplement_aplus48').resolve()
    auth,seals,locations,cases,runs,index=authorize(root,sup,a.authorization.resolve(),a.run_index.resolve())
    stat=a.run_index.stat();index_snapshot=dict(path=str(a.run_index.resolve()),bytes=stat.st_size,mtime_ns=stat.st_mtime_ns)
    metric_dir=a.run_index.resolve().parent;annotations={};annotation_sources=[]
    comparison_receipt=read(metric_dir/'COMPARISON_RECEIPT.json')
    if comparison_receipt.get('status')!='FOLLOWUP_EXISTING48_AUTHORIZED_COMPARISON_COMPLETE' or comparison_receipt.get('scope')!='FOLLOWUP_EXISTING48_NOT_FRESH_BLIND' or path(root,comparison_receipt['authorization'])!=a.authorization.resolve():
        raise ValueError('Existing comparison metrics must have the same completed authorized source')
    for campaign,location in locations.items():
        if path(root,comparison_receipt['seals'][campaign])!=Path(location).resolve():raise ValueError('Comparison metric seal differs')
    for name in ('ALL48_SCORES.json','PAID_PREFIX_SCORES.json'):
        source=metric_dir/name
        if not source.is_file():continue
        annotation_sources.append(str(source))
        for item in read(source):
            key=(item['case_id'],item['method'],item['view'],item.get('CPU_prefix'))
            annotations[key]=item
    out=a.out.resolve()
    banned=[root/'RESULTS',root/'queues',root/'sealed_truth',root/'src',root/'configs',
        sup/'RESULTS',sup/'queues',sup/'sealed_truth',sup/'src',sup/'configs',root/'data_public',sup/'data_public']
    if out in(root,sup)or any(out.is_relative_to(folder)for folder in banned):raise ValueError('New diagnostic output outside frozen inputs/source required')
    for seal in seals.values():
        if any(out.is_relative_to(Path(d['path']).resolve())for d in seal['directories']):raise ValueError('Output overlaps a sealed prediction directory')
    out.mkdir(parents=True,exist_ok=False);(out/'GRADIENTS').mkdir()
    # All authorization/index source validation above precedes numerical arrays.
    import numpy as np
    legacy=a.legacy_repo.resolve();sys.path.insert(0,str(legacy))
    import research.research_fbg_r3_experiments.src.physics as physics
    if not Path(physics.__file__).resolve().is_relative_to(legacy):raise RuntimeError('Physical module is outside requested qualified legacy repo')
    rows=[];unique=[];setups=[];changes=[];case_runs=defaultdict(list);uid=0
    for run in runs:case_runs[run['case_id']].append(run)
    for cid in cases:
        groups=[];public_cache={};entries={};case_refs=[];gradient_arrays={}
        for run in case_runs[cid]:
            campaign=run['source_campaign'];protocol=seals[campaign]['manifest_snapshot']['protocol']
            prior=protocol.get('prior')
            if not isinstance(prior,dict)or not all(k in prior for k in ('smooth','tv','epsilon')):
                raise ValueError('Exact shared frozen scalar prior required')
            prior_key=tuple(float(prior[k])for k in ('smooth','tv','epsilon'))
            folder=Path(run['run_folder']);requests=requests_for(run);history=None;pool=None;run_error=None
            try:
                history=read(Path(run['source_history_path']));byid={int(h['candidate']):h for h in history}
                if len(byid)!=len(history):raise ValueError('Duplicate candidate IDs in immutable HISTORY')
                receipt=read(Path(run['source_receipt_path']))
                if receipt.get('protocol')!=protocol or receipt.get('truth_access')!='NONE':raise ValueError('Frozen receipt protocol/truth declaration changed')
                original_task=next(t for t in seals[campaign]['manifest_snapshot']['tasks']if t['task_id']==run['source_task_id'])
                if receipt.get('task')!=original_task:raise ValueError('Frozen receipt task differs from indexed source manifest')
                cap=float(protocol['CPU_cap']);eligible=[h for h in history if h.get('eligible')and h.get('legal')and h['CPU_seconds']<=cap]
                expected_views={'raw':receipt.get('raw_candidate'),'primary':receipt.get('primary_candidate'),
                    'corrected':read(folder/'corrected.json').get('candidate')if(folder/'corrected.json').is_file()else None}
                if any(run['view_candidate_ids'].get(k)!=v for k,v in expected_views.items()):raise ValueError('Indexed saved view IDs differ from immutable receipt')
                for point in PREFIXES:
                    available=[h for h in eligible if h['CPU_seconds']<=point]
                    expected=min(available,key=lambda h:h['Phi'])['candidate']if available else None
                    if run['prefix_candidate_ids'].get(str(point))!=expected:raise ValueError('Indexed paid-min-Phi prefix differs from immutable public rule')
                with np.load(run['source_candidates_path'],allow_pickle=False)as z:
                    saved=z['m']
                    if len(saved)!=len(history):raise ValueError('CANDIDATES/HISTORY row count differs')
                    needed={int(r['candidate_id'])for r in requests if r['candidate_id']is not None}
                    pool={i:saved[i].copy()for i in needed if 0<=i<len(saved)}
                originals=Path(run['original_public_input_path']);execution=Path(run['execution_public_input_path'])
                for source in (originals,execution):
                    if source not in public_cache:public_cache[source]=load_public(np,source)
                public=public_cache[originals]
                if not same_operator(np,public,public_cache[execution]):raise ValueError('Execution input differs from original exact observation/operators')
                box=protocol['material_box'];lower=float(public.get('lower',box[0]));upper=float(public.get('upper',box[1]))
                if (lower,upper)!=(float(public_cache[execution].get('lower',box[0])),float(public_cache[execution].get('upper',box[1]))):raise ValueError('Execution box differs')
                group=next((g for g in groups if g['prior_key']==prior_key and g['box']==(lower,upper)and same_operator(np,g['public'],public)),None)
                if group is None:
                    c,w=time.process_time(),time.perf_counter();model=physics.MaxwellModel(public)
                    model.reg=physics.Regularizer(model.shape,smooth=prior_key[0],tv=prior_key[1],epsilon=prior_key[2])
                    group=dict(id=cid+'__operator'+str(len(groups)+1),public=public,prior_key=prior_key,box=(lower,upper),model=model,materials=[])
                    groups.append(group);setups.append(dict(case_id=cid,operator_group_id=group['id'],prior=dict(zip(('smooth','tv','epsilon'),prior_key)),box=[lower,upper],
                        CPU_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-w,setup_counts=dict(model.cost)))
            except Exception as exc:run_error=repr(exc)
            for request in requests:
                row=dict(case_id=cid,method=run['method'],source_campaign=campaign,source_task_id=run['source_task_id'],
                    aplus_source_split=run.get('aplus_source_split'),source_run_folder=run['run_folder'],source_candidates_path=run['source_candidates_path'],
                    source_public_input_path=run['original_public_input_path'],view=request['view'],CPU_prefix=request['CPU_prefix'],
                    candidate_id=request['candidate_id'],saved_metadata=request['saved_metadata'],status='UNAVAILABLE',error=run_error,
                    unique_evaluation_id=None,truth_array_access='NONE',scope='OFFLINE_DIAGNOSTIC_SAME_DISCRETIZATION')
                rows.append(row);case_refs.append(row)
                if run_error or request['candidate_id']is None:continue
                ident=int(request['candidate_id']);m=pool.get(ident);h=byid.get(ident)
                if m is None or h is None:row['error']='INDEXED_CANDIDATE_UNAVAILABLE';continue
                if np.iscomplexobj(m)or m.ndim!=1 or not np.isfinite(m).all():row['error']='INVALID_FROZEN_REAL_MATERIAL';continue
                row.update(saved_candidate_CPU_seconds=h.get('CPU_seconds'),saved_phase=h.get('phase'),saved_Phi=h.get('Phi'),saved_D=h.get('D'),saved_R=h.get('R'),
                    original_eligible=h.get('eligible'),original_legal=h.get('legal'),box=list(group['box']))
                entry=next((item for item in group['materials']if np.array_equal(item['m'],m)),None)
                reference=dict(case_id=cid,method=run['method'],view=row['view'],CPU_prefix=row['CPU_prefix'],candidate_id=ident)
                if entry is None:
                    uid+=1;eid=cid+'__eval%05d'%uid;model=group['model'];before=dict(model.cost);c,w=time.process_time(),time.perf_counter()
                    entry=dict(unique_evaluation_id=eid,operator_group_id=group['id'],case_id=cid,m=m.copy(),status='FAILED',metrics=None,vectors=None,references=[reference],error=None,
                        gradient_archive='GRADIENTS/'+cid+'.npz',gradient_key_prefix=eid)
                    try:
                        if m.shape!=(model.n,):raise ValueError('Frozen material shape differs from exact physical model')
                        entry['metrics'],entry['vectors']=evaluate_material(np,model,m,public['Y'],*group['box']);entry['status']='EVALUATED'
                    except Exception as exc:entry['error']=repr(exc)
                    finally:
                        counts={k:model.cost.get(k,0)-before.get(k,0)for k in set(model.cost)|set(before)}
                        entry['cost']=dict(physical_model_delta=counts,LU=counts.get('LU',0),forward_RHS=counts.get('source_RHS',0),
                            adjoint_RHS=counts.get('adjoint_RHS',0),actual_complex_RHS=counts.get('actual_complex_RHS',0),
                            CPU_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-w,
                            additional_prior_gradient_evaluations=1 if entry['status']=='EVALUATED'else None)
                    if entry['vectors'] is None:entry['gradient_archive']=None;entry['gradient_key_prefix']=None
                    group['materials'].append(entry);entries[eid]=entry;unique.append({k:v for k,v in entry.items()if k not in('m','vectors')})
                    if entry['vectors']:
                        for name,array in entry['vectors'].items():gradient_arrays[eid+'__'+name]=array
                else:entry['references'].append(reference)
                row.update(status=entry['status'],error=entry['error'],unique_evaluation_id=entry['unique_evaluation_id'],cost_reference=entry['unique_evaluation_id'])
                if entry['metrics']:
                    row.update(entry['metrics'])
                    row['RUN_INDEX_saved_Phi_minus_history']=(request['saved_metadata'].get('Phi')-h['Phi'])if request['saved_metadata'].get('Phi')is not None and h.get('Phi')is not None else None
                    for name in('Phi','D','R'):row[name+'_minus_saved']=row[name]-h[name]if h.get(name)is not None else None
                annkey=(cid,run['method'],row['view'],row['CPU_prefix']);ann=annotations.get(annkey)
                if ann and ann.get('selection_candidate',ann.get('candidate'))==ident:
                    for name in('RMSE','NRMSE','IoU','relativeL2'):row['comparison_'+name]=ann.get(name)
                    for name in('crossgrid','noise_certificate_eligible','noise_qualification'):row['comparison_'+name]=ann.get(name)
                    row['comparison_metric_source']='existing authorized comparison; no GT arrays read here'
                elif ann:row['comparison_annotation_issue']='CANDIDATE_ID_MISMATCH; annotation not attached'
            # Material pool is limited to this indexed run; selected unique
            # vectors survive only for this case's transitions and archive.
            pool=None
        lookup={(r['method'],r['view'],r['CPU_prefix']):r for r in case_refs}
        for method in METHODS:
            for left,right in zip(PREFIXES[:-1],PREFIXES[1:]):
                changes.append(transition(np,lookup[method,'paid_prefix',left],lookup[method,'paid_prefix',right],entries,'adjacent_paid_prefix'))
            changes.append(transition(np,lookup[method,'raw',None],lookup[method,'corrected',None],entries,'raw_to_corrected_selected_estimate'))
        if gradient_arrays:np.savez_compressed(out/'GRADIENTS'/(cid+'.npz'),**gradient_arrays)
    # References point to actual unique work, never charge duplicates as solves.
    byid={u['unique_evaluation_id']:u for u in unique}
    for row in rows:
        item=byid.get(row['unique_evaluation_id'])
        if item:
            n=len(item['references']);row.update(shared_reference_count=n,
                allocated_offline_CPU_seconds=item['cost']['CPU_seconds']/n,
                cost_policy='one actual unique-material evaluation; equal CPU shares are allocations, not extra solves')
    stat2=a.run_index.stat()
    if(stat2.st_size,stat2.st_mtime_ns)!=(index_snapshot['bytes'],index_snapshot['mtime_ns']):raise RuntimeError('RUN_INDEX changed during diagnostic')
    sys.path.insert(0,str(root/'src'));from seal_results import validate_seal
    for seal in seals.values():validate_seal(seal,require_authorization=True)
    prefix_summary,view_summary,transition_summary=summaries(np,rows,changes)
    table(out,'PHYSICAL_REFERENCES',rows);table(out,'UNIQUE_PHYSICAL_EVALUATIONS',unique)
    table(out,'PREFIX_SUMMARY',prefix_summary);table(out,'VIEW_SUMMARY',view_summary);table(out,'TRANSITION_SUMMARY',transition_summary)
    table(out,'INCUMBENT_TRANSITIONS',changes);table(out,'MODEL_SETUP_COSTS',setups)
    plotting_cpu,plotting_wall=time.process_time(),time.perf_counter();artifacts=plot_compact(np,out,rows,changes)if not a.no_figures else[]
    plot_cost=dict(CPU_seconds=time.process_time()-plotting_cpu,wall_seconds=time.perf_counter()-plotting_wall)
    for seal in seals.values():validate_seal(seal,require_authorization=True)
    receipt=dict(status='OFFLINE_DIAGNOSTIC_COMPLETE',scope='FOLLOWUP_EXISTING48_NOT_FRESH_BLIND',
        physics_scope='exact original same finite discretization/prior, NOT independent validation',
        parent_authorization=str(a.authorization.resolve()),validated_seals=locations,run_index=index_snapshot,
        expected_runs=240,expected_references=1920,actual_references=len(rows),expected_transitions=1200,actual_transitions=len(changes),
        evaluated_references=sum(r['status']=='EVALUATED'for r in rows),unique_material_evaluations=len(unique),
        failed_unique_evaluations=sum(u['status']!='EVALUATED'for u in unique),
        actual_LU=sum(u['cost']['LU']for u in unique),forward_RHS=sum(u['cost']['forward_RHS']for u in unique),
        adjoint_RHS=sum(u['cost']['adjoint_RHS']for u in unique),actual_complex_RHS=sum(u['cost']['actual_complex_RHS']for u in unique),
        unique_evaluation_CPU_seconds=sum(u['cost']['CPU_seconds']for u in unique),
        unique_evaluation_wall_seconds=sum(u['cost']['wall_seconds']for u in unique),
        CPU_seconds=time.process_time()-START_CPU,wall_seconds=time.perf_counter()-START_WALL,
        plotting_cost=plot_cost,annotation_sources=annotation_sources,GT_array_reads=0,
        new_inversions=0,new_candidates=0,new_Hessian_or_SVD=0,online_selection_changed=False,
        gradients='original unscaled physical Phi/data/prior; unit-step box mapping; full vectors retained in per-case GRADIENTS',
        transition_interpretation='successive paid-min-Phi incumbents, NOT native material steps or a reconstructed optimizer trajectory',
        selection='exact indexed saved views/paid prefixes; no new selection by Phi gradient, truth error, or stationarity',
        cost_scope='offline setup/unique gradients/serialization/plots charged separately from old online runs; timing components overlap script total',
        timer_scope='entry/guards/imports/reads/evaluations/gradient archives/tables/plots included; final receipt write and interpreter exit excluded',
        extra_RUN_INDEX_fair_budget_IDs='ignored: common_actual_CPU/common_physical_LU do not extend the fixed1920 physical-reference scope',
        artifacts=artifacts,no_success_gate=True,scientific_conclusion_owner='parent')
    dump(out/'RECEIPT.json',receipt)
    print(json.dumps({k:receipt[k]for k in('status','actual_references','unique_material_evaluations','actual_LU','actual_complex_RHS','GT_array_reads')}))


if __name__=='__main__':main()
