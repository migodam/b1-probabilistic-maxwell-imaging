"""Post-seal offline metrics/figures. Never opens truth before parent authorization.

Usage: score_results.py --root R13 --seal reviewed_seal.json --out NEW_OUTPUT
The parent alone authorizes the seal. No method, NN, optimizer or forward solve
is called by this module. Missing/failed expected tasks remain table rows.
"""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS'):os.environ[k]='1'
import argparse,csv,json,math,itertools
from pathlib import Path
from datetime import datetime,timezone
try:
    from .seal_results import read,validate_seal
except ImportError:
    from seal_results import read,validate_seal


def json_default(v):
    import numpy as np
    if isinstance(v,np.generic):return v.item()
    if isinstance(v,np.ndarray):return v.tolist()
    if isinstance(v,Path):return str(v)
    raise TypeError(type(v).__name__)


def table(out,name,rows):
    (out/(name+'.json')).write_text(json.dumps(rows,indent=2,default=json_default,allow_nan=False)+'\n')
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with (out/(name+'.csv')).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader()
        for r in rows:w.writerow({k:json.dumps(v,default=json_default) if isinstance(v,(dict,list)) else v for k,v in r.items()})


def otsu(m):
    import numpy as np
    if np.ptp(m)<=1e-12:return np.zeros(m.size,bool),float(m[0])
    hist,edges=np.histogram(m,bins=256);centers=(edges[:-1]+edges[1:])/2;weight=np.cumsum(hist);total=np.cumsum(hist*centers);den=weight*(weight[-1]-weight)
    score=np.full(256,-np.inf);good=den>0;score[good]=(total[-1]*weight[good]-total[good]*weight[-1])**2/den[good];threshold=float(centers[np.argmax(score)])
    return m>threshold,threshold


def truth_spec(root,task):
    """Resolve locations only; no private content is read here."""
    cid=task['case_id'];spec=dict(task.get('truth_spec',{}))
    if spec:
        if spec.get('path'):spec['path']=str((root/spec['path']).resolve())
        return spec
    if cid.startswith('r13_fresh_'):
        return dict(path=str(root/'sealed_truth/fresh48'/f'{cid}.npz'),key='material_coarse',support_rule='fixed_chi_gt_1',amplitude_rule='max_truth_minus_host',host=.1,object_id=cid)
    if '__' in cid and cid.split('__')[0] in ['A','B','C']:
        cohort,ident=cid.split('__',1)
        return dict(path=str(root.parents[1]/'research/research_max_struct_r12/supplement_paired/offline_truth'/cohort/f'{ident}.npz'),key='m',support_key='support',amplitude_key='amplitude',support_rule='prediction_Otsu256_vs_saved_support' if cohort=='C' else'fixed_chi_gt_1',object_id=cid)
    vector=root/'delegated/physics/new3d/sealed_truth'/f'{cid}.npz'
    if vector.is_file():return dict(path=str(vector),key='chi',support_rule='fixed_chi_gt_1',amplitude_rule='max_truth_minus_host',host=.1,object_id=cid.split('_crossgrid')[0])
    return dict(path=None,material_truth_available=False,support_rule='NO_PIXEL_GT',object_id=cid)


def load_truth(spec):
    """Called only after validated authorized seal, never at import."""
    import numpy as np
    if not spec.get('path'):return None
    with np.load(spec['path'],allow_pickle=False) as z:
        truth=z[spec.get('key','m')].copy().ravel()
        support=z[spec['support_key']].copy().ravel().astype(bool) if spec.get('support_key') else truth>float(spec.get('support_threshold',1.))
        amp=float(z[spec['amplitude_key']]) if spec.get('amplitude_key') else float(spec.get('amplitude',truth.max()-spec.get('host',.1)))
    return dict(m=truth,support=support,amplitude=amp,rule=spec['support_rule'],path=spec['path'])


def material_metrics(m,truth,lower,upper):
    import numpy as np
    if m is None:return dict(available=False,availability_reason='NO_PAID_CANDIDATE')
    m=np.asarray(m).ravel()
    if np.iscomplexobj(m) or not np.isfinite(m).all():return dict(available=False,availability_reason='INVALID_MATERIAL')
    result=dict(available=True,availability_reason='AVAILABLE',material_min=float(m.min()),material_max=float(m.max()),material_legal=bool(m.min()>=lower-1e-6 and m.max()<=upper+1e-6))
    if truth is None:return dict(**result,RMSE=None,NRMSE=None,boxrange_NRMSE=None,relativeL2=None,IoU=None,mean_amplitude_bias=None,mean_amplitude_bias_normalized=None,IoU_rule='NO_PIXEL_GT')
    t=truth['m'];amp=truth['amplitude']
    if t.shape!=m.shape:return dict(available=False,availability_reason='MATERIAL_TRUTH_SHAPE_MISMATCH')
    e=m-t;rmse=float(np.linalg.norm(e)/np.sqrt(t.size));seg,threshold=otsu(m) if truth['rule'].startswith('prediction_Otsu') else(m>1.,1.);sup=truth['support'];union=np.count_nonzero(seg|sup)
    bias=float(np.mean(e[sup])) if np.any(sup) else float(np.mean(e))
    return dict(**result,RMSE=rmse,NRMSE=rmse/amp if amp>0 else None,boxrange_NRMSE=rmse/(upper-lower) if upper>lower else None,relativeL2=float(np.linalg.norm(e)/np.linalg.norm(t)) if np.linalg.norm(t)>0 else None,IoU=float(np.count_nonzero(seg&sup)/union) if union else 1.,IoU_rule=truth['rule'],IoU_threshold=threshold,amplitude=amp,mean_amplitude_bias=bias,mean_amplitude_bias_normalized=bias/amp if amp>0 else None)


def physical_counts(cost):
    """Canonical disjoint counters; aggregate RHS takes precedence over roles.

    Scalar FrozenModel actual_complex_RHS already includes source/adjoint/tangent
    and receiver-transpose solves. VectorPublicModel instead increments disjoint
    source_RHS and adjoint_RHS (and explicitly tracked derivative roles, if any).
    No setup Gram, propagation columns, history totals or nested deltas are added.
    """
    lu_key='physical_LU' if 'physical_LU' in cost else 'LU' if 'LU' in cost else None
    if 'actual_complex_RHS' in cost:
        rhs_keys=['actual_complex_RHS'];rhs=cost['actual_complex_RHS'];scope='aggregate actual_complex_RHS; role counters not added again'
    else:
        rhs_keys=[k for k in ['source_RHS','adjoint_RHS','tangent_RHS','receiver_transpose_RHS'] if k in cost]
        rhs=sum(cost[k] for k in rhs_keys) if rhs_keys else None
        scope='sum of explicitly disjoint source/adjoint/tangent/receiver-transpose solve RHS; operator columns and Gram RHS excluded'
    return dict(physical_LU=cost[lu_key] if lu_key else None,physical_complex_RHS=rhs,
                physical_source_RHS=cost.get('source_RHS'),physical_adjoint_RHS=cost.get('adjoint_RHS'),
                physical_cost_source_keys={'LU':lu_key,'complex_RHS':rhs_keys},physical_RHS_scope=scope)


def paired_stats(rows,reference):
    import numpy as np
    from scipy.stats import wilcoxon
    primary=[r for r in rows if r['view']=='primary'];methods=sorted({r['method'] for r in primary});lookup={(r['case_id'],r['method']):r for r in primary};cases=sorted({r['case_id'] for r in primary});out=[]
    for method in methods:
        if method==reference:continue
        groups={}
        for cid in cases:
            a=lookup.get((cid,method),{});b=lookup.get((cid,reference),{})
            if a.get('available') and b.get('available') and a.get('NRMSE') is not None and b.get('NRMSE') is not None:groups.setdefault(a['object_id'],[]).append(a['NRMSE']-b['NRMSE'])
        d=np.array([np.mean(v) for v in groups.values()]);rng=np.random.default_rng(131313)
        boots=np.mean(d[rng.integers(0,len(d),size=(4000,len(d)))],axis=1) if len(d) else None
        try:p=float(wilcoxon(d,alternative='two-sided').pvalue) if len(d) and np.any(d!=0) else 1. if len(d) else None
        except ValueError:p=None
        out.append(dict(method=method,reference=reference,expected_case_pairs=len(cases),available_case_pairs=sum(len(x) for x in groups.values()),material_objects=len(d),mean_object_paired_NRMSE_difference=float(np.mean(d)) if len(d) else None,bootstrap95_low=float(np.quantile(boots,.025)) if len(d) else None,bootstrap95_high=float(np.quantile(boots,.975)) if len(d) else None,wilcoxon_descriptive_uncorrected_p=p,bootstrap_seed=131313,bootstrap_replicates=4000,unit='material object; average repeated conditions within object',interpretation='descriptive unadjusted interval/p-value; no new success gate'))
    return out


def figures(out,rows,candidates,pictures,truths,costs,methods,structure_gallery=False):
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    figdir=out/'FIGURES';figdir.mkdir()
    def save(fig,name):
        fig.savefig(figdir/(name+'.png'),dpi=160,bbox_inches='tight');fig.savefig(figdir/(name+'.pdf'),bbox_inches='tight');plt.close(fig)
    lookup={(r['case_id'],r['method'],r['view']):r for r in rows};cases=sorted(truths)
    def show(ax,array,title,lo=0,hi=20):
        ax.set_title(title,fontsize=8);ax.set_xticks([]);ax.set_yticks([])
        if array is None:ax.text(.5,.5,'UNAVAILABLE',ha='center');return
        n=array.size;side=int(round(n**.5))
        if side*side==n:image=array.reshape(side,side)
        else:
            side=int(round(n**(1/3)))
            if side**3!=n:ax.text(.5,.5,'Unplottable material shape',ha='center');return
            vol=array.reshape(side,side,side);mid=side//2;image=np.concatenate([vol[mid],vol[:,mid,:],vol[:,:,mid]],axis=1);title+='\n3D central XY/XZ/YZ slices';ax.set_title(title,fontsize=8)
        ax.imshow(image,origin='lower',cmap='viridis',vmin=lo,vmax=hi,interpolation='nearest')
    for mg in range(0,len(methods),5):
        subset=methods[mg:mg+5]
        for page,start in enumerate(range(0,len(cases),4),1):
            batch=cases[start:start+4];fig,axes=plt.subplots(len(batch),len(subset)+2,figsize=(2.5*(len(subset)+2),2.9*len(batch)),squeeze=False)
            for i,cid in enumerate(batch):
                truth=truths[cid];show(axes[i,0],truth['m'] if truth else None,cid+'\nOffline truth')
                initial=next((pictures[(cid,m,'initial')] for m in methods if (cid,m,'initial') in pictures),None)
                show(axes[i,1],initial,cid+'\nCommon initial material')
                for j,method in enumerate(subset,2):
                    r=lookup.get((cid,method,'primary'),{});metric=r.get('NRMSE');suffix=f'NRMSE={metric:.4f}' if metric is not None else'No material GT' if r.get('available') else'UNAVAILABLE'
                    show(axes[i,j],pictures.get((cid,method,'primary')),cid+'\n'+method+'\n'+suffix)
            fig.suptitle('Full case IDs; shared contrast color scale [0,20]. Primary paid selection.',fontsize=11);fig.tight_layout();save(fig,f'materials_methods{mg//5+1}_page{page:02d}')
    if structure_gallery:
        eligible_cases=[cid for cid in cases if truths[cid] is not None]
        for mg in range(0,len(methods),5):
            subset=methods[mg:mg+5]
            for page,start in enumerate(range(0,len(eligible_cases),4),1):
                batch=eligible_cases[start:start+4];fig,axes=plt.subplots(len(batch),len(subset)+2,figsize=(2.6*(len(subset)+2),3.1*len(batch)),squeeze=False)
                for i,cid in enumerate(batch):
                    gt=truths[cid]['m'];vmax=float(np.max(gt));vmax=vmax if vmax>0 else 1.
                    initial=next((pictures[(cid,m,'initial')] for m in methods if (cid,m,'initial') in pictures),None)
                    arrays=[('Offline truth',gt),('Common initial',initial)]+[(method,pictures.get((cid,method,'primary'))) for method in subset]
                    for j,(label,array) in enumerate(arrays):
                        clipped=int(np.count_nonzero(array>vmax)) if array is not None else 0
                        title=(cid+'\n' if j==0 else '')+label+f'\nrow vmax={vmax:.4g}\nsaturated cells={clipped}'
                        show(axes[i,j],array,title,lo=0,hi=vmax)
                fig.suptitle('DIAGNOSTIC VIEW ONLY: GT-derived scale shared by every method in each row. Absolute [0,20] gallery retained.',fontsize=10)
                fig.subplots_adjust(left=.025,right=.935,top=.90,bottom=.025,hspace=.62,wspace=.14)
                for i in range(len(batch)):
                    pos=axes[i,-1].get_position();bar=fig.add_axes([.95,pos.y0,.008,pos.height])
                    fig.colorbar(axes[i,0].images[0],cax=bar,label='physical contrast chi; row-common')
                save(fig,f'diagnostic_structure_methods{mg//5+1}_page{page:02d}')
    # Per-case paid Phi incumbent trajectories; CPU and physical LU are separate axes.
    grouped={};costlookup={(c['case_id'],c['method']):c for c in costs}
    for r in candidates:
        if r.get('eligible') and r.get('legal') and r.get('NRMSE') is not None:grouped.setdefault((r['case_id'],r['method']),[]).append(r)
    for page,start in enumerate(range(0,len(cases),4),1):
        batch=cases[start:start+4];fig,axes=plt.subplots(len(batch),2,figsize=(12,3.2*len(batch)),squeeze=False)
        for i,cid in enumerate(batch):
            for method in methods:
                seq=sorted(grouped.get((cid,method),[]),key=lambda r:r['CPU_seconds']);best=None;xs=[];ls=[];ys=[]
                for r in seq:
                    if best is None or r['Phi']<best['Phi']:best=r
                    xs.append(r['CPU_seconds']);ls.append(r.get('physical_LU_at_candidate'));ys.append(best['NRMSE'])
                if seq:
                    axes[i,0].plot(xs,ys,label=method,lw=1)
                    external=costlookup.get((cid,method),{}).get('external_child_CPU_seconds')
                    if external is not None:axes[i,0].scatter([external],[ys[-1]],s=16,marker='x')
                    valid=[j for j,v in enumerate(ls) if v is not None];axes[i,1].plot([ls[j] for j in valid],[ys[j] for j in valid],label=method,lw=1)
            for j,xlabel in enumerate(['Actual elapsed process CPU at paid audit (s)','Cumulative physical LU at paid audit (projector work separate)']):axes[i,j].set_title(cid);axes[i,j].set_xlabel(xlabel);axes[i,j].set_ylabel('NRMSE of min-Phi paid incumbent');axes[i,j].grid(alpha=.2)
        handles,labels=axes[0,0].get_legend_handles_labels()
        if handles:axes[0,0].legend(fontsize=7,ncol=3)
        else:
            for ax in axes.ravel():ax.text(.5,.5,'Material NRMSE unavailable: no pixel GT / paid candidate',ha='center',transform=ax.transAxes,fontsize=8)
        fig.tight_layout();save(fig,f'paid_cost_curves_page{page:02d}')
    for method in methods:
        selected=[r for r in rows if r['method']==method and r['view'] in ['raw','corrected'] and r.get('available')];pairs={}
        for r in selected:pairs.setdefault(r['case_id'],{})[r['view']]=r
        if not pairs:continue
        fig,axes=plt.subplots(1,3,figsize=(14,4))
        for cid,pair in pairs.items():
            if set(pair)!=set(['raw','corrected']):continue
            a,b=pair['raw'],pair['corrected']
            for ax,key in zip(axes,['D','R','Phi']):
                if a.get(key) is not None and b.get(key) is not None:ax.scatter(a[key],b[key],s=20);ax.annotate(cid,(a[key],b[key]),fontsize=5)
            if a.get('NRMSE') is not None and b.get('NRMSE') is not None:
                pass
        for ax,key in zip(axes,['D','R','Phi']):
            lim=ax.get_xlim();ylim=ax.get_ylim();lo=min(lim[0],ylim[0]);hi=max(lim[1],ylim[1]);ax.plot([lo,hi],[lo,hi],':',color='gray');ax.set_xlabel('raw '+key);ax.set_ylabel('corrected '+key)
        fig.suptitle(method+' corrector decomposition; lower objective need not improve material');fig.tight_layout();save(fig,'corrector_components_'+method)
        fig,ax=plt.subplots(figsize=(7,5))
        for cid,pair in pairs.items():
            if set(pair)!=set(['raw','corrected']):continue
            aa,bb=pair['raw'],pair['corrected']
            if any(r.get('NRMSE') is None or r.get('residual_over_delta') is None for r in [aa,bb]):continue
            x,y=aa['residual_over_delta'],aa['NRMSE'];xx,yy=bb['residual_over_delta'],bb['NRMSE']
            ax.scatter([x],[y],marker='o',s=18);ax.scatter([xx],[yy],marker='x',s=22);ax.annotate('',xy=(xx,yy),xytext=(x,y),arrowprops=dict(arrowstyle='->',lw=.8));ax.annotate(cid,(xx,yy),fontsize=6)
        ax.axvline(1.,color='gray',ls=':');ax.set_xlabel('Physical residual / declared delta');ax.set_ylabel('Amplitude NRMSE');ax.set_title(method+' raw to corrected: arrows retain both candidates');ax.grid(alpha=.2);fig.tight_layout();save(fig,'corrector_material_arrows_'+method)
    material=[r for r in rows if r['view']=='primary' and r.get('NRMSE') is not None]
    if material:
        fig,ax=plt.subplots(figsize=(max(8,len(methods)*.8),5));lookup0={(r['case_id'],r['method']):r for r in material}
        for cid in cases:
            vals=[lookup0.get((cid,m),{}).get('NRMSE',np.nan) for m in methods];ax.plot(range(len(methods)),vals,alpha=.3,lw=.7)
        ax.set_xticks(range(len(methods)),methods,rotation=45,ha='right');ax.set_ylabel('Primary amplitude NRMSE');ax.set_title('Paired material cases; missing runs remain absent points');ax.grid(axis='y',alpha=.2);fig.tight_layout();save(fig,'paired_primary_NRMSE')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--seal',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--reference',default='dr');p.add_argument('--no-figures',action='store_true');p.add_argument('--structure-gallery',action='store_true',help='Additional offline diagnostic gallery with one shared GT-derived physical scale per row; never used for selection');a=p.parse_args()
    # This guard precedes NumPy loading, every result material and ALL truth access.
    root=a.root.resolve();seal=read(a.seal);validate_seal(seal,require_authorization=True)
    if Path(seal['root']).resolve()!=root:raise RuntimeError('Seal root differs from requested study root')
    out=a.out.resolve()
    if any(out.is_relative_to(Path(d['path'])) for d in seal['directories']):
        raise ValueError('Scoring output must be outside sealed result/queue directories to preserve the prediction snapshot.')
    out.mkdir(parents=True,exist_ok=False)
    import numpy as np
    cfg=seal['manifest_snapshot'];campaign=cfg['campaign'];tasks=cfg['tasks'];methods=list(dict.fromkeys(t['method'] for t in tasks));now=datetime.now(timezone.utc).isoformat()
    queue_path=root/'queues'/campaign/'RECEIPT.json';queue=read(queue_path) if queue_path.is_file() else{};job_list=queue.get('jobs',[]);jobs={r['task_id']:r for r in job_list}
    attempt_rows=[]
    for key,values in queue.items():
        if not isinstance(values,list) or not(key=='jobs' or 'attempt' in key.lower() or 'interrupt' in key.lower()):continue
        for index,value in enumerate(values):
            if not isinstance(value,dict):value={'raw_value':value}
            external=value.get('external',{}) or{};cpu=external.get('child_process_CPU_seconds',value.get('child_process_CPU_seconds'))
            attempt_rows.append(dict(source_field=key,source_index=index,task_id=value.get('task_id'),external_child_CPU_seconds=cpu,record=value,count_scope='final authoritative jobs list' if key=='jobs' else 'separate history; possible overlap with final jobs, not added without explicit attempt identity'))
    rows=[];candidate_rows=[];prefixrows=[];costs=[];failures=[];truths={};pictures={};truthrecords=[];arrows=[];eventrows=[]
    for task in tasks:
        cid,method,tid=task['case_id'],task['method'],task['task_id'];folder=root/'RESULTS'/campaign/tid;spec=truth_spec(root,task)
        if cid not in truths:
            try:
                truth=load_truth(spec);truth_error=None
            except Exception as ex:
                truth=None;truth_error=repr(ex);failures.append(dict(campaign=campaign,case_id=cid,task_id=tid,method=method,status='TRUTH_LOAD_FAILED',issues=[truth_error]))
            truths[cid]=truth;truthrecords.append(dict(case_id=cid,object_id=spec.get('object_id',cid),path=spec.get('path'),support_rule=spec.get('support_rule'),truth_accessed_at_utc=now,qualification='No material GT' if truth is None else'Sealed truth; offline scoring only',truth_load_error=truth_error))
        truth=truths[cid];hist=[];mats=None;receipt={};issue=[]
        public_kind='';components=1;known_noise=True;crossgrid=False;public_noise_role=None
        try:
            with np.load(root/task['input_file'],allow_pickle=False) as public_info:
                public_kind=str(public_info['kind']) if 'kind' in public_info else ''
                components=int(public_info['components']) if 'components' in public_info else 1
                known_noise=bool(public_info['noise_sigma_known']) if 'noise_sigma_known' in public_info else True
                public_noise_role=str(public_info['sigma_role']) if 'sigma_role' in public_info else None
                crossgrid=int(public_info['forward_side'])!=int(public_info['inverse_side']) if 'forward_side' in public_info and 'inverse_side' in public_info else False
        except Exception as ex:issue.append('PUBLIC_QUALIFICATION_UNAVAILABLE:'+repr(ex))
        inferred_kind=public_kind or ('vector3d' if components==3 else 'scalar_tm_2d')
        kind=str(task.get('physics_kind',inferred_kind));cohort=str(task.get('cohort','fresh48' if cid.startswith('r13_fresh_') else 'development' if cid.split('__')[0] in ['A','B','C'] else kind))
        diagnostic=not known_noise or any(x in (kind+' '+cohort+' '+cid+' '+str(public_noise_role)).lower() for x in ['fresnel','austria','classic','diagnostic'])
        domain=kind+'/'+cohort
        common=dict(campaign=campaign,case_id=cid,object_id=spec.get('object_id',cid),task_id=tid,method=method,physics_kind=kind,cohort=cohort,domain=domain,noise_certificate_eligible=not diagnostic,data_model_noise_certificate=not diagnostic and not crossgrid,crossgrid=crossgrid,noise_qualification=public_noise_role)
        try:
            receipt=read(folder/'RECEIPT.json');hist=read(folder/'HISTORY.json')
            with np.load(folder/'CANDIDATES.npz',allow_pickle=False) as z:mats=z['m'].copy()
            if len(mats)!=len(hist):raise ValueError('History/candidate length mismatch')
            if receipt.get('truth_access')!='NONE':raise ValueError('Online truth access declaration inconsistent')
        except Exception as ex:issue.append(repr(ex));hist=[];mats=None
        lower=float(cfg['protocol'].get('material_box',[.02,20])[0]);upper=float(cfg['protocol'].get('material_box',[.02,20])[1])
        try:
            with np.load(root/task['input_file'],allow_pickle=False) as public:
                lower=float(public['lower']) if 'lower' in public else lower
                upper=float(public['upper']) if 'upper' in public else upper
        except Exception as ex:issue.append('PUBLIC_INPUT_UNAVAILABLE:'+repr(ex))
        cap=float(task.get('CPU_cap',cfg['protocol']['CPU_cap']))
        status=receipt.get('status','MISSING');scored={};valid={}
        common['noise_qualification']=receipt.get('noise_qualification') or public_noise_role or ('Diagnostic discrepancy; not a CN noise certificate' if diagnostic else 'iid circular complex RMS; CN99 noise-only ball; cross-grid mismatch excluded' if crossgrid else 'iid circular complex RMS; CN99 raw noise ball')
        if any(x in common['noise_qualification'].lower() for x in ['diagnostic','not known noise','not_noise']):common['noise_certificate_eligible']=False;common['data_model_noise_certificate']=False
        for h,m in zip(hist,mats if mats is not None else[]):
            index=int(h['candidate']);metric=material_metrics(m,truth,lower,upper);cost=h.get('physics_cost',{});canonical=physical_counts(cost)
            row=dict(common,**metric,candidate=index,run_status=status,partial_failed_run=status=='FAILED',eligible=bool(h.get('eligible',False)),legal=bool(h.get('legal',metric.get('material_legal',False))),CPU_seconds=h['CPU_seconds'],phase=h['phase'],Phi=h['Phi'],D=h['D'],R=h['R'],residual_over_delta=h.get('residual_over_delta'),noise_feasible=h.get('noise_feasible'),discrepancy_feasible=h.get('noise_feasible'),physical_LU_at_candidate=canonical['physical_LU'],physical_complex_RHS_at_candidate=canonical['physical_complex_RHS'],physical_cost_source_keys=canonical['physical_cost_source_keys'],metadata=h)
            candidate_rows.append(row);scored[index]=row
            if row['eligible'] and row['legal'] and metric['available'] and h['CPU_seconds']<=cap:valid[index]=row
        def index_from_json(name):
            path=folder/(name+'.json')
            return read(path).get('candidate') if path.is_file() else None
        indices={'initial':index_from_json('initial'),'raw':receipt.get('raw_candidate'),'corrected':index_from_json('corrected'),'primary':receipt.get('primary_candidate'),'secondary_only':receipt.get('secondary_candidate')}
        indices['secondary']=indices['secondary_only'] if indices['secondary_only'] in valid else indices['primary']
        raw=valid.get(indices['raw']);indices['feasible_raw_retain']=indices['raw'] if raw is not None and raw.get('noise_feasible') else indices['primary']
        for view,index in indices.items():
            selected=valid.get(index);fallback=view=='secondary' and indices['secondary_only'] not in valid or view=='feasible_raw_retain' and index!=indices['raw']
            if selected:
                r=dict(selected,view=view,selection_candidate=index,fallback_to_primary=bool(fallback),requested_view=view,selection_rule='minR among paid noise-feasible, fallback primary' if view=='secondary' else'raw retained if noise-feasible, otherwise primary' if view=='feasible_raw_retain' else'runner frozen '+view)
                pictures[(cid,method,view)]=mats[index].copy()
            else:r=dict(common,view=view,available=False,availability_reason='NO_LEGAL_ELIGIBLE_PAID_'+view.upper(),run_status=status,partial_failed_run=status=='FAILED',selection_candidate=index,fallback_to_primary=bool(fallback),NRMSE=None)
            rows.append(r)
        stored=read(folder/'PREFIXES.json') if(folder/'PREFIXES.json').is_file() else[];stored_map={float(r['CPU_cap']):r.get('candidate') for r in stored}
        for point in cfg['protocol'].get('prefixes',[]):
            pool=[r for r in valid.values() if r['CPU_seconds']<=point];selected=min(pool,key=lambda r:r['Phi']) if pool else None;index=selected['candidate'] if selected else None
            match=stored_map.get(float(point))==index
            if not match:issue.append(f'PREFIX_MISMATCH_AT_{point}')
            prefixrows.append(dict(selected or common,view='paid_prefix',CPU_prefix=point,available=selected is not None,selection_candidate=index,stored_prefix_matches_recomputed=match,eligibility='legal, runner eligible, audit process CPU <= prefix and total cap',NRMSE=selected.get('NRMSE') if selected else None))
        raw=valid.get(indices['raw']);corrected=valid.get(indices['corrected'])
        if raw and corrected:
            ar=dict(common,raw_candidate=raw['candidate'],corrected_candidate=corrected['candidate'])
            for k in ['D','R','Phi','NRMSE','RMSE','IoU','residual_over_delta']:
                ar['raw_'+k]=raw.get(k);ar['corrected_'+k]=corrected.get(k);ar['delta_'+k]=corrected[k]-raw[k] if raw.get(k) is not None and corrected.get(k) is not None else None
            ar['raw_noise_feasible']=raw.get('noise_feasible');ar['corrected_noise_feasible']=corrected.get('noise_feasible');arrows.append(ar)
        if(folder/'EVENTS.json').is_file():
            for e in read(folder/'EVENTS.json'):eventrows.append(dict(common,event=e))
        job=jobs.get(tid,{});external=job.get('external',{});cpu=external.get('child_process_CPU_seconds',job.get('child_process_CPU_seconds'));phys=receipt.get('physical_cost',{});proj=receipt.get('projector_cost') or{}
        costs.append(dict(common,status=status,returncode=job.get('returncode'),external_child_CPU_seconds=cpu,external_child_CPU_scope='external queue child-process total' if cpu is not None else'UNKNOWN_NOT_IMPUTED',internal_process_CPU_seconds=receipt.get('whole_process_CPU_seconds'),queue_task_wall_seconds=job.get('wall_seconds'),internal_wall_seconds=receipt.get('wall_seconds'),**physical_counts(phys),physical_cost=phys,projector_cost=proj,auxiliary_receipt_not_added_again=receipt.get('auxiliary_receipt'),physical_count_scope='physical_cost already includes model calls; operation categories/LU/auxiliary deltas not added together',endpoint_RSS_bytes=receipt.get('endpoint_RSS_bytes'),peak_memory_bytes=external.get('peak_RSS_bytes'),peak_memory_scope=external.get('memory_scope','NOT_MEASURED; end-of-process resident snapshot is not peak memory')))
        if issue or status=='FAILED' or not valid:failures.append(dict(common,status=status,issues=issue,runner_error=receipt.get('error'),valid_paid_candidates=len(valid),policy='Retained in full expected-task denominator; valid failed-run candidates separately labeled partial.'))
    # Full rows are retained; availability denominators accompany every mean.
    summary=[];domains=list(dict.fromkeys(r['domain'] for r in rows))
    for domain,method,view in itertools.product(domains,methods,list(dict.fromkeys(r['view'] for r in rows))):
        selected=[r for r in rows if r['domain']==domain and r['method']==method and r['view']==view]
        if not selected:
            continue
        available=[r for r in selected if r.get('available')]
        record=dict(domain=domain,physics_kind=selected[0]['physics_kind'],cohort=selected[0]['cohort'],method=method,view=view,expected_tasks=len(selected),available_tasks=len(available),failed_tasks=sum(r.get('run_status')=='FAILED' for r in selected),fallbacks_to_primary=sum(r.get('fallback_to_primary',False) for r in selected))
        for key in ['RMSE','NRMSE','boxrange_NRMSE','relativeL2','IoU','mean_amplitude_bias','D','R','Phi','residual_over_delta']:
            vals=[r[key] for r in available if r.get(key) is not None];record[key+'_mean_available']=float(np.mean(vals)) if vals else None;record[key+'_available_count']=len(vals)
        summary.append(record)
    totals={'expected_tasks':len(tasks),'external_queue_reported_child_CPU_seconds':queue.get('finished_child_CPU_seconds',queue.get('completed_child_CPU_seconds')),'external_final_jobs_CPU_known_sum':sum(r['external_child_CPU_seconds'] for r in attempt_rows if r['source_field']=='jobs' and r['external_child_CPU_seconds'] is not None),'queue_history_attempt_rows':sum(r['source_field']!='jobs' for r in attempt_rows),'queue_history_unknown_CPU_rows':sum(r['source_field']!='jobs' and r['external_child_CPU_seconds'] is None for r in attempt_rows),'queue_history_cost_policy':'All explicit partial/failed/interrupted attempt records retained separately. History may overlap final jobs; never sum overlapping snapshots or invent missing attempt CPU.','external_child_CPU_seconds_known_sum':sum(c['external_child_CPU_seconds'] for c in costs if c['external_child_CPU_seconds'] is not None),'external_child_CPU_unknown_tasks':sum(c['external_child_CPU_seconds'] is None for c in costs),'controller_CPU_seconds_separate':queue.get('controller_CPU_seconds'),'queue_wall_seconds':queue.get('wall_seconds'),'physical_LU_known_sum':sum(c['physical_LU'] for c in costs if c['physical_LU'] is not None),'physical_complex_RHS_known_sum':sum(c['physical_complex_RHS'] for c in costs if c['physical_complex_RHS'] is not None),'physical_count_unknown_tasks':sum(c['physical_LU'] is None for c in costs),'projector_count_totals':{k:sum(c['projector_cost'].get(k,0) for c in costs) for k in sorted({k for c in costs for k in c['projector_cost']})},'domain_physical_totals':{domain:{'physical_LU_known_sum':sum(c['physical_LU'] for c in costs if c['domain']==domain and c['physical_LU'] is not None),'physical_complex_RHS_known_sum':sum(c['physical_complex_RHS'] for c in costs if c['domain']==domain and c['physical_complex_RHS'] is not None),'tasks':sum(c['domain']==domain for c in costs)} for domain in domains},'no_nested_double_count':True,'notes':'Do not add internal CPU to external CPU; do not add cumulative history/auxiliary physical deltas to final physical_cost; projector RHS are separate from physical solves; queue elapsed is not sum of child walls.'}
    for name,data in [('SCORES',rows),('ALL_PAID_AND_OVERRUN_CANDIDATES',candidate_rows),('PAID_PREFIX_SCORES',prefixrows),('METHOD_VIEW_SUMMARY',summary),('PAIRED_STATS',[dict(stat,domain=domain) for domain in domains for stat in paired_stats([r for r in rows if r['domain']==domain],a.reference)]),('CORRECTOR_DECOMPOSITION',arrows),('EVENTS',eventrows),('COST_LEDGER',costs),('QUEUE_ATTEMPT_HISTORY',attempt_rows),('FAILURES',failures),('TRUTH_SOURCES',truthrecords)]:table(out,name,data)
    (out/'COST_TOTALS.json').write_text(json.dumps(totals,indent=2,default=json_default,allow_nan=False)+'\n')
    if not a.no_figures:
        import re
        for domain in domains:
            subset=[r for r in rows if r['domain']==domain];ids={r['case_id'] for r in subset};dm=list(dict.fromkeys(r['method'] for r in subset));dout=out/('DOMAIN_'+re.sub('[^a-zA-Z0-9_-]+','_',domain));dout.mkdir()
            figures(dout,subset,[r for r in candidate_rows if r['domain']==domain],{k:v for k,v in pictures.items() if k[0] in ids},{k:v for k,v in truths.items() if k in ids},[c for c in costs if c['domain']==domain],dm,structure_gallery=a.structure_gallery)
    report=['# R13 post-seal offline report','',f'Campaign: `{campaign}`. Expected tasks: {len(tasks)}; task issues/failures: {len(failures)}.',f'Parent-authorized seal: `{a.seal.resolve()}`. Truth opened only after authorization and byte/mtime validation.', '', '## Estimators and denominators','', '`primary` is frozen paid minimum-Phi selection. `secondary` uses a paid noise-feasible minimum-R candidate and falls back to primary; `secondary_only` exposes actual availability. `feasible_raw_retain` preserves an already noise-feasible raw candidate and otherwise uses primary. Raw/corrected stages are retained independently. Missing tasks are never removed from expected denominators; valid paid evidence from failed runs is labeled partial.', '', '## Material metric definitions','', 'RMSE uses all material cells; amplitude NRMSE divides by saved inherited amplitude for development cases, or max(cell-average truth)-host for new synthetic data. Boxrange NRMSE divides RMSE by the public allowed material box width (upper-lower); relative L2 is reported separately. Development C uses prediction Otsu256 versus saved support; other synthetic rows use fixed chi>1 support. Mean amplitude bias is mean signed material error on truth support, alongside its amplitude-normalized value. Continuous-field support IoU is descriptive. No new joint success threshold is introduced. Fresnel without trustworthy pixel GT has null material accuracy metrics.', '', '## Paid time and costs','', 'Prefix candidates must be legal, eligible, already audited by the requested cumulative process CPU and within total budget. Cost curves use paid physical min-Phi incumbents, never GT-selected roots. Whole external child CPU, internal audit CPU, physical LU/RHS and projector work remain separate. Missing cost counters remain null. Explicit queue partial/interrupted/failed attempt history is retained separately in QUEUE_ATTEMPT_HISTORY; unknown costs are not imputed and possible overlapping history is not added again. Complete cost records are in COST_LEDGER/COST_TOTALS. Scalar aggregate RHS and vector disjoint source/adjoint RHS use documented physical_cost_source_keys; nested deltas are never added again.', '', '## Physics domains and statistical scope','', 'METHOD_VIEW_SUMMARY, PAIRED_STATS and figure directories are separated by task physics_kind/cohort domain. No single pooled material accuracy mean combines classic, measured Fresnel and synthetic 3D/2D experiments. Noise/discrepancy eligibility retains the runner qualification; Fresnel/classic diagnostics are not CN noise certificates.', '', 'Paired effects resample material objects, averaging repeated conditions within an object. Fixed bootstrap seed131313,4000replicates; confidence intervals and Wilcoxon p-values are descriptive and unadjusted. Available paired counts are displayed against full case denominators. Development objects are exposed; cross-grid mismatch and prior selection can affect material error independently of data-ball feasibility.', '', '## Artifacts and review','', 'JSON/CSV retain all estimator, candidate, prefix, cost, failure and source rows. Paginated figures include full case IDs and shared color scale, paid CPU/LU curves, corrector objective decomposition and paired material errors. Three-dimensional pictures are central orthogonal slices, not complete volume validation. Optional --structure-gallery adds an explicitly diagnostic row-shared GT-derived physical scale, a colorbar and counts of saturated prediction cells. It never changes selection or the absolute gallery. Fake verification fixtures under delegated/evidence/generation/scorer_fixture are excluded from these campaign output roots and must never enter publication. Parent visual/scientific review and publication remain outstanding.', '', '## Failures','']
    report.extend([f"- `{r['task_id']}`: {r['status']}; {', '.join(r['issues']) or 'runner failure / no valid candidate'}" for r in failures] or ['No missing or failed task detected by this scorer.'])
    (out/'REPORT.md').write_text('\n'.join(report)+'\n');(out/'SCORING_RECEIPT.json').write_text(json.dumps({'campaign':campaign,'authorized_seal':str(a.seal.resolve()),'truth_accessed_at_utc':now,'expected_tasks':len(tasks),'score_rows':len(rows),'candidate_rows':len(candidate_rows),'new_physical_forward':0,'new_adjoint':0,'new_optimization':0,'parent_scientific_review':'PENDING'},indent=2)+'\n')
    print(json.dumps({'campaign':campaign,'expected_tasks':len(tasks),'score_rows':len(rows),'failures':len(failures),'output':str(out)}))
if __name__=='__main__':main()
