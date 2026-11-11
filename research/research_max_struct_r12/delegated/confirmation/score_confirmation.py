"""Root-authorized postseal scoring/figures; zero forward/optimizer calls.

Primary and secondary retain every object; unavailable is never imputed.
Bootstrap samples whole material IDs, 10000 draws seeded20261112.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key]='1'
from pathlib import Path
import argparse,json,csv,itertools,time,shutil,warnings
import numpy as np
from scipy.stats import chi2
ROOT=Path(__file__).resolve().parent
REPO=ROOT.parents[3]
PUBLIC=ROOT/'selected_public_outputs'
METHODS=('cie_single','cie_dual','global_dual','a_plus','fullG19','dr','multistart','nn_direct')
LABELS=dict(zip(METHODS,('CIE','Dual','Global','A+','G19','DR','Multi','NN')))
CAPS=(30,60,120,240)
SEED=20261112
REPLICATES=10000

def read(path):return json.loads(Path(path).read_text())
def dump(path,obj):
    def convert(x):
        if isinstance(x,np.ndarray):return x.tolist()
        if isinstance(x,np.generic):return x.item()
        if isinstance(x,Path):return str(x)
        raise TypeError(type(x).__name__)
    Path(path).write_text(json.dumps(obj,indent=2,allow_nan=False,default=convert)+'\n')
def csvwrite(path,rows):
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with Path(path).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader()
        for r in rows:w.writerow({k:json.dumps(v) if isinstance(v,(list,dict)) else v for k,v in r.items()})
def table(out,name,rows):dump(out/(name+'.json'),rows);csvwrite(out/(name+'.csv'),rows)
def relpath(path):return Path(str(path).replace('\\','/'))
def loadm(path,key='m'):
    with np.load(path,allow_pickle=False) as z:return z[key].copy()
def normalized_stage(row):
    # Windows compact paths must be split before pathlib on POSIX. The old
    # companion's stage field accidentally contains the complete Windows path.
    return str(row['compact_primary']).replace('\\','/').split('/')[0]

def material_metrics(m,t,meta):
    valid=bool(m.shape==(576,) and not np.iscomplexobj(m) and np.isfinite(m).all())
    if not valid:return dict(available=False,availability_reason='INVALID_SEALED_MATERIAL')
    d=m-t;rmse=float(np.linalg.norm(d)/24);amplitude=float(t.max()-.1)
    tp=t>1.;mp=m>1.;union=np.count_nonzero(tp|mp)
    iou=float(np.count_nonzero(tp&mp)/union) if union else 1.
    legal=bool(m.min()>=.02-1e-12 and m.max()<=20.+1e-12)
    ans=dict(available=True,availability_reason='AVAILABLE',RMSE=rmse,NRMSE=rmse/amplitude,
        relativeL2=float(np.linalg.norm(d)/np.linalg.norm(t)),IoU=iou,amplitude=amplitude,
        legal_inherited=legal,legal_exact=bool(m.min()>=.02 and m.max()<=20.),
        material_min=float(m.min()),material_max=float(m.max()),
        descriptive_success=bool(legal and rmse/amplitude<=.15 and iou>=.7))
    masks={'background':np.isclose(t,.1,rtol=0,atol=1e-12)}
    # Continuous fields are not relabeled as two-phase material. Their peak
    # pixel alone is not a high-phase region.
    for name,value in [('low',2.5),('high',float(meta['high_chi']))]:
        masks[name]=np.isclose(t,value,rtol=0,atol=1e-12) if meta['family']!='continuous_random' else np.zeros(576,bool)
    for name,mask in masks.items():
        ans[name+'_true_cells']=int(mask.sum())
        ans[name+'_RMSE']=float(np.sqrt(np.mean(d[mask]**2))) if mask.any() else None
        ans[name+'_MAE']=float(np.mean(abs(d[mask]))) if mask.any() else None
        ans[name+'_bias']=float(np.mean(d[mask])) if mask.any() else None
    return ans

def bootstats(values,indices):
    v=np.asarray([np.nan if x is None else float(x) for x in values])
    finite=np.isfinite(v);n=int(finite.sum())
    if not n:return dict(n=0,mean=None,median=None,mean_low=None,mean_high=None,median_low=None,median_high=None,valid_bootstraps=0)
    x=v[indices];den=np.isfinite(x).sum(axis=1);ok=den>0
    means=np.nansum(x[ok],axis=1)/den[ok]
    medians=np.nanmedian(x[ok],axis=1)
    return dict(n=n,mean=float(np.mean(v[finite])),median=float(np.median(v[finite])),
        mean_low=float(np.quantile(means,.025)),mean_high=float(np.quantile(means,.975)),
        median_low=float(np.quantile(medians,.025)),median_high=float(np.quantile(medians,.975)),valid_bootstraps=int(ok.sum()))

def classify(row):
    if not row.get('available'):return 'UNAVAILABLE'
    if row.get('noise_feasible') is None:return 'NO_SAVED_PHYSICAL_METRIC'
    if not row['noise_feasible']:return 'OUTSIDE_CN99'
    return 'INSIDE_CN99_DESCRIPTIVE_SUCCESS' if row['descriptive_success'] else 'INSIDE_CN99_MATERIAL_ERROR'

def score(out):
    out.mkdir(parents=True,exist_ok=False)
    audit=out/'frozen_inputs';audit.mkdir()
    for name in ('ROOT_CORE_FREEZE.json','SECONDARY_SELECTION_FREEZE.json'):
        shutil.copy2(ROOT/name,audit/name)
    shutil.copy2(Path(__file__),audit/'score_confirmation_source.py')
    import sys
    sys.path.insert(0,str(REPO))
    from research.research_z_alg_r4_diversity.src.physics import Regularizer
    reg=Regularizer((24,24))  # pure original R; no model/forward instantiation.
    source=read(PUBLIC/'COMPACT_PUBLIC_MANIFEST.json')
    truth={};metadata={};public={};cases={}
    for stage in ('dev6','confirm24'):
        seal=read(ROOT/('SEALED_'+stage.upper()+'_PUBLIC_COHORT.json'))
        if seal['status']!='COMPLETE_PUBLIC_COHORT_SEALED_BEFORE_GT':raise ValueError('Cohort seal absent')
        metas=read(ROOT/'sealed_truth'/stage/'generator_metadata.json')['records']
        cases[stage]=sorted(r['case_id'] for r in metas)
        metadata.update({r['case_id']:r for r in metas})
        for case in cases[stage]:
            truth[case]=loadm(ROOT/'sealed_truth'/stage/(case+'.npz'))
            with np.load(ROOT/'data_public'/stage/(case+'.npz'),allow_pickle=False) as z:
                y=z['Y'];sigma=float(z['sigma'])
                public[case]=dict(Ynorm=float(np.linalg.norm(y)),sigma=sigma,
                    delta=float(sigma*np.sqrt(chi2.ppf(.99,2*y.size)/2)))
    endpoints=[];prefixes=[];matrices={};index={}
    for src in source['rows']:
        stage=normalized_stage(src);case=src['case_id'];method=src['method'];meta=metadata[case]
        directory=PUBLIC/relpath(src['compact_primary']).parent
        receipt=read(directory/'RECEIPT.json');counts=receipt.get('physical_model_cost',receipt.get('physics_cost',{}))
        base=dict(stage=stage,case_id=case,method=method,family=meta['family'],regime=meta['regime'],
            high_chi=meta['high_chi'],external_CPU_seconds=src['completed_external_child_CPU_seconds'],
            physical_576_LU=counts.get('LU',0),complex_RHS=counts.get('actual_complex_RHS',0),
            execution_status=receipt['status'],auxiliary_iterations=src.get('auxiliary_iterations_from_paid_history'),
            source_compact_directory=str(directory),truth_file=str(ROOT/'sealed_truth'/stage/(case+'.npz')))
        primary_final_meta=read(directory/'R11_CAP_PREDICTION.json') if (directory/'R11_CAP_PREDICTION.json').exists() else src.get('primary_paid_audit',{})
        for selector in ('primary','secondary'):
            row=dict(base,selector=selector,view='endpoint',CPU_cap=240,available=False,
                availability_reason=('NOT_REGISTERED_FOR_DEV6' if stage=='dev6' else 'NO_NOISE_FEASIBLE_PAID_CANDIDATE'),descriptive_success=None)
            path=src['compact_primary'] if selector=='primary' else src.get('compact_secondary')
            if path:
                m=loadm(PUBLIC/relpath(path));row.update(material_metrics(m,truth[case],meta))
                R=float(reg.evaluate(m)[0]);row['R']=R
                if selector=='primary':
                    row['Phi']=src['Phi'];row['selected_CPU_seconds']=primary_final_meta.get('paid_whole_child_CPU_seconds',primary_final_meta.get('CPU_seconds'))
                    rel=primary_final_meta.get('physical_relative_residual',primary_final_meta.get('relative_data_residual'))
                    if rel is None:rel=float(np.sqrt(max(0.,2*(row['Phi']-R))))
                    row['relative_data_residual']=float(rel)
                    row['physical_metric_source']='Saved paid physical residual; otherwise derived sealedPhi minus sameR (noFW)'
                else:
                    sel=read(directory/'NOISE_FEASIBLE_MIN_PRIOR.json')['selection']
                    row.update(Phi=sel['Phi'],selected_CPU_seconds=sel['CPU_seconds'],relative_data_residual=sel['raw_residual_norm']/public[case]['Ynorm'],
                        physical_metric_source='Saved secondary paid physical residual',saved_secondary_R=sel['R'])
                row.update(delta_relative=public[case]['delta']/public[case]['Ynorm'],
                    noise_feasible=bool(row['relative_data_residual']<=public[case]['delta']/public[case]['Ynorm']))
                row['physical_material_class']=classify(row)
                matrices[(stage,selector,case,method,'endpoint')]=m
            else:row['physical_material_class']='UNAVAILABLE'
            endpoints.append(row);index[(stage,selector,case,method,'endpoint')]=row
        for selector in ('primary','secondary'):
            if selector=='primary' and (directory/'R11_CPU_CHECKPOINTS.json').exists():
                specs={int(r['CPU_cap_seconds']):r for r in read(directory/'R11_CPU_CHECKPOINTS.json')}
                mats=loadm(directory/'R11_CHECKPOINT_MATERIALS.npz');kind='baseline'
            elif selector=='primary':
                specs={int(r['cap']):r for r in receipt.get('prefixes',[])};kind='search'
                mats=None
            else:
                specs={int(r['cap']):r for r in read(directory/'NOISE_FEASIBLE_MIN_PRIOR_TIMEPREFIX.json')} if (directory/'NOISE_FEASIBLE_MIN_PRIOR_TIMEPREFIX.json').exists() else {};kind='secondary';mats=None
            archive=directory/('paid_prefix.npz' if kind=='search' else 'NOISE_FEASIBLE_MIN_PRIOR_TIMEPREFIX.npz')
            for cap in CAPS:
                spec=specs.get(cap,{});m=None
                if kind=='baseline' and 'material_index' in spec:m=mats[spec['material_index']]
                elif kind!='baseline' and archive.exists():
                    with np.load(archive,allow_pickle=False) as z:
                        if 'm_'+str(cap) in z.files:m=z['m_'+str(cap)].copy()
                row=dict(base,selector=selector,view='prefix',CPU_cap=cap,available=False,
                    availability_reason='NO_RETAINED_CAP_ELIGIBLE_STATE' if stage=='confirm24' or selector=='primary' else 'NOT_REGISTERED_FOR_DEV6',descriptive_success=None)
                if m is not None:
                    row.update(material_metrics(m,truth[case],meta));R=float(reg.evaluate(m)[0]);phi=spec.get('original_Phi',spec.get('Phi'))
                    row.update(R=R,Phi=phi,selected_CPU_seconds=spec.get('state_paid_CPU_seconds',spec.get('CPU_seconds',spec.get('CPU'))))
                    rel=spec.get('raw_residual_norm')
                    row['relative_data_residual']=float(rel/public[case]['Ynorm']) if rel is not None else float(np.sqrt(max(0.,2*(phi-R))))
                    row['physical_metric_source']='Saved residual or sealedPhi minus sameR; noFW'
                    row.update(delta_relative=public[case]['delta']/public[case]['Ynorm'],noise_feasible=bool(row['relative_data_residual']<=public[case]['delta']/public[case]['Ynorm']))
                    matrices[(stage,selector,case,method,cap)]=m
                row['physical_material_class']=classify(row);prefixes.append(row)
                index[(stage,selector,case,method,cap)]=row
    table(out,'per_object',endpoints);table(out,'prefix_per_object',prefixes)
    aggregates=[];pairs=[];objectpairs=[]
    rng_indices={stage:np.random.default_rng(SEED).integers(0,len(ids),(REPLICATES,len(ids))) for stage,ids in cases.items()}
    for stage in ('dev6','confirm24'):
        ids=cases[stage];draws=rng_indices[stage]
        for selector in ('primary','secondary'):
            for view in ('endpoint',)+CAPS:
                for method in METHODS:
                    rows=[index[(stage,selector,c,method,view)] for c in ids]
                    a=dict(stage=stage,selector=selector,view=view,method=method,objects=len(ids),
                        available=sum(r['available'] for r in rows),unavailable=sum(not r['available'] for r in rows),
                        descriptive_success=sum(r['descriptive_success'] is True for r in rows),success_denominator=len(ids),
                        inside_CN99_material_error=sum(r['physical_material_class']=='INSIDE_CN99_MATERIAL_ERROR' for r in rows),
                        outside_CN99=sum(r['physical_material_class']=='OUTSIDE_CN99' for r in rows),
                        mean_external_whole_run_CPU_seconds=float(np.mean([r['external_CPU_seconds'] for r in rows])),
                        sum_physical_576_LU=sum(r['physical_576_LU'] for r in rows),sum_complex_RHS=sum(r['complex_RHS'] for r in rows))
                    for metric in ('RMSE','NRMSE','relativeL2','IoU','relative_data_residual','R'):
                        values=[r.get(metric) for r in rows]
                        stats=bootstats(values,draws)
                        for k,v in stats.items():a[metric+'_'+k]=v
                    aggregates.append(a)
                for left,right in itertools.combinations(METHODS,2):
                    ls=[index[(stage,selector,c,left,view)] for c in ids];rs=[index[(stage,selector,c,right,view)] for c in ids]
                    p=dict(stage=stage,selector=selector,view=view,left=left,right=right,sign='left_minus_right',objects=len(ids),
                        paired_available=sum(l['available'] and r['available'] for l,r in zip(ls,rs)),
                        left_available=sum(r['available'] for r in ls),right_available=sum(r['available'] for r in rs),
                        bootstrap_replicates=REPLICATES,bootstrap_seed=SEED,
                        resampling='Whole original material roster; missing values retained asNA. Continuous statistics conditional on bothavailable; no imputation.',
                        mean_external_CPU_difference=float(np.mean([l['external_CPU_seconds']-r['external_CPU_seconds'] for l,r in zip(ls,rs)])))
                    for metric in ('NRMSE','relativeL2','RMSE','IoU'):
                        diffs=[l[metric]-r[metric] if l['available'] and r['available'] else None for l,r in zip(ls,rs)]
                        for k,v in bootstats(diffs,draws).items():p[metric+'_difference_'+k]=v
                    coverage=[float(l['available'])-float(r['available']) for l,r in zip(ls,rs)]
                    success=[float(l['descriptive_success'] is True)-float(r['descriptive_success'] is True) for l,r in zip(ls,rs)]
                    p['availability_difference']=bootstats(coverage,draws);p['coverage_aware_descriptive_success_difference']=bootstats(success,draws)
                    pairs.append(p)
                    if view=='endpoint':
                        for c,l,r in zip(ids,ls,rs):
                            o=dict(stage=stage,selector=selector,case_id=c,left=left,right=right,left_available=l['available'],right_available=r['available'])
                            for metric in ('NRMSE','relativeL2','RMSE','IoU','Phi','R','relative_data_residual'):
                                o[metric+'_left_minus_right']=l.get(metric)-r.get(metric) if l['available'] and r['available'] else None
                            o.update(left_class=l['physical_material_class'],right_class=r['physical_material_class']);objectpairs.append(o)
    table(out,'aggregate',aggregates);table(out,'pairs',pairs);table(out,'per_object_pairs',objectpairs)
    # Exact object-level DR comparisons: all objects and signed differences,
    # no GT-selected image shortlist or new binary material-success threshold.
    drdiff=[]
    for stage in cases:
        for selector in ('primary','secondary'):
            for c in cases[stage]:
                d=index[(stage,selector,c,'dr','endpoint')]
                for control in ('cie_dual','fullG19'):
                    r=index[(stage,selector,c,control,'endpoint')]
                    item=dict(stage=stage,selector=selector,case_id=c,control=control,DR_available=d['available'],control_available=r['available'],
                        DR_class=d['physical_material_class'],control_class=r['physical_material_class'],
                        DR_NRMSE=d.get('NRMSE'),control_NRMSE=r.get('NRMSE'),DR_relativeL2=d.get('relativeL2'),control_relativeL2=r.get('relativeL2'))
                    if d['available'] and r['available']:
                        delta=d['NRMSE']-r['NRMSE'];item.update(DR_minus_control_NRMSE=delta,direction='IMPROVEMENT' if delta<0 else 'DEGRADATION' if delta>0 else 'TIE',
                            relative_NRMSE_change=delta/r['NRMSE'] if r['NRMSE'] else None,DR_minus_control_IoU=d['IoU']-r['IoU'])
                    drdiff.append(item)
    table(out,'DR_vs_Dual_G19_per_object',drdiff)
    # Compare saved raw DR proposal shadows without new projection/FW and
    # without choosing a shadow by truth. All audit errors are exported.
    raw=[];rawsummary=[];rawroot=ROOT/'raw_dr_public_audits'
    if (rawroot/'RAW_DR_PUBLIC_INDEX.json').exists():
        for entry in read(rawroot/'RAW_DR_PUBLIC_INDEX.json')['rows']:
            stage=entry['stage'];case=entry['case_id'];meta=read(rawroot/relpath(entry['metadata']))
            with np.load(rawroot/relpath(entry['npz']),allow_pickle=False) as z:ms=z['candidate_materials']
            candidate_rows=[]
            for i,(m,physical) in enumerate(zip(ms,meta['records'])):
                r=dict(stage=stage,case_id=case,method='dr',view='raw_audited_shadow',audit_index=i,**physical)
                r.update(material_metrics(m,truth[case],metadata[case]));r['physical_material_class']=classify(r)
                raw.append(r);candidate_rows.append(r)
            eligible=[r for r in candidate_rows if r['eligible_proposal144']]
            if eligible:
                choices={'raw_proposal_minPhi':min(eligible,key=lambda r:r['Phi']),
                    'raw_last_cap_eligible_audited_shadow':eligible[-1]}
                for view,r in choices.items():
                    summary=dict(r,view=view,selection='Saved proposal physicalminPhi or lastcapeligibleaudit; neverGTselected')
                    rawsummary.append(summary);matrices[(stage,view,case,'dr','endpoint')]=ms[r['audit_index']]
    table(out,'raw_DR_all_audits',raw);table(out,'raw_DR_fixed_selections',rawsummary)
    # Endpoint within-object primary→secondary change never omits unavailable.
    changes=[]
    for stage in cases:
        for case in cases[stage]:
            for method in METHODS:
                p=index[(stage,'primary',case,method,'endpoint')];s=index[(stage,'secondary',case,method,'endpoint')]
                r=dict(stage=stage,case_id=case,method=method,secondary_available=s['available'],primary_class=p['physical_material_class'],secondary_class=s['physical_material_class'])
                for metric in ('NRMSE','relativeL2','IoU','Phi','R','relative_data_residual'):
                    r[metric+'_secondary_minus_primary']=s.get(metric)-p.get(metric) if s['available'] else None
                changes.append(r)
    table(out,'primary_secondary_per_object',changes)
    status=dict(status='OFFLINE_SCORING_COMPLETE',objects={k:len(v) for k,v in cases.items()},endpoints=len(endpoints),prefix_rows=len(prefixes),
        bootstrap_seed=SEED,bootstrap_replicates=REPLICATES,source_stage_field='Ignored malformedcompanion; derivefirstnormalizedcompactpathcomponent',
        legal_metric='InheritedR11 descriptive legality uses1e-12 tolerance; exactlegality alsoexported; secondaryonlinelegalunchanged',
        raw_DR='Savedalreadyphysicalauditedshadowmaterials; un-auditedfinalnextshadow not reconstructed' if raw else 'UNAVAILABLE_IN_COMPACT_NO_RECONSTRUCTION',
        GT_access='ROOT_EXPLICIT_AUTHORIZED_POSTSEAL_OFFLINE_ONLY',new_FW=0,new_algorithms=0,scoring_CPU_seconds_near_write=time.process_time())
    dump(out/'SCORING_RECEIPT.json',status)
    print(json.dumps(dict(status=status['status'],endpoints=len(endpoints),prefixes=len(prefixes),raw_audits=len(raw))),flush=True)
    return endpoints,prefixes,aggregates,pairs,drdiff,rawsummary,truth,metadata,matrices

def figures(out,fig,data):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    fig.mkdir(parents=True,exist_ok=False)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'svg.fonttype':'none'})
    rows,prefix,aggregate,pairs,drdiff,raw,truth,meta,materials=data
    index={(r['stage'],r['selector'],r['case_id'],r['method']):r for r in rows}
    colors={m:plt.get_cmap('tab10')(i) for i,m in enumerate(METHODS)}
    def save(f,name):
        f.savefig(fig/(name+'.png'),dpi=170,bbox_inches='tight',facecolor='white')
        f.savefig(fig/(name+'.pdf'),bbox_inches='tight',facecolor='white');plt.close(f)
    def grid(stage,selector,cases,name):
        f,axes=plt.subplots(len(cases),9,figsize=(14.5,1.52*len(cases)+.7),squeeze=False)
        f.subplots_adjust(left=.06,right=.965,bottom=.025,top=.95,wspace=.055,hspace=.35)
        for i,case in enumerate(cases):
            for j,method in enumerate(('GT',)+METHODS):
                ax=axes[i,j];ax.set_xticks([]);ax.set_yticks([])
                if method=='GT':arr=truth[case];caption=case.rsplit('_',1)[-1]+' '+meta[case]['family']+'\n'+meta[case]['regime']+' '+str(int(meta[case]['high_chi']))
                else:
                    r=index[(stage,selector,case,method)];arr=materials.get((stage,selector,case,method,'endpoint'))
                    caption=f"N={r['NRMSE']:.3f}  I={r['IoU']:.2f}" if r['available'] else 'unavailable'
                if arr is None:
                    ax.set_facecolor('#e7e9ec');ax.text(.5,.5,'unavailable',ha='center',va='center',transform=ax.transAxes,fontsize=7)
                else:ax.imshow(arr.reshape(24,24),vmin=0,vmax=20,cmap='viridis',origin='lower',interpolation='nearest')
                ax.set_xlabel(caption,fontsize=6.3 if method=='GT' else 7,color=('firebrick' if method!='GT' and arr is not None and not r['descriptive_success'] else 'black'))
                if i==0:ax.set_title('GT' if method=='GT' else LABELS[method],fontsize=10)
        f.suptitle(stage+' / '+selector+' / all listed objects; N=amplitude NRMSE, I=chi>1 IoU',fontsize=11,y=.993)
        bar=f.add_axes([.976,.2,.008,.58]);f.colorbar(plt.cm.ScalarMappable(norm=Normalize(0,20),cmap='viridis'),cax=bar,label='contrast chi')
        save(f,name)
    for stage in ('dev6','confirm24'):
        cases=sorted(c for c in truth if c.startswith('r12_'+stage+'_'))
        for selector in ('primary','secondary'):
            grid(stage,selector,cases,stage+'_'+selector+'_ALL_MATERIALS')
            if stage=='confirm24':
                for family in ('levelsets','branches','mixed_islands','continuous_random'):
                    subset=[c for c in cases if meta[c]['family']==family]
                    grid(stage,selector,subset,stage+'_'+selector+'_'+family)
            f,axes=plt.subplots(1,2,figsize=(11,4.1))
            for j,m in enumerate(METHODS):
                a=next(r for r in aggregate if r['stage']==stage and r['selector']==selector and r['view']=='endpoint' and r['method']==m)
                values=[r['NRMSE'] for r in rows if r['stage']==stage and r['selector']==selector and r['method']==m and r['available']]
                if values:
                    jitter=np.linspace(-.17,.17,len(values));axes[0].scatter(j+jitter,values,s=15,alpha=.7,color=colors[m])
                    axes[0].errorbar(j,a['NRMSE_mean'],yerr=[[a['NRMSE_mean']-a['NRMSE_mean_low']],[a['NRMSE_mean_high']-a['NRMSE_mean']]],fmt='D',color='black',ms=4,capsize=3)
                    axes[1].scatter(a['mean_external_whole_run_CPU_seconds'],a['NRMSE_mean'],color=colors[m],s=45)
                    axes[1].annotate(LABELS[m]+f" ({a['available']}/{len(cases)})",(a['mean_external_whole_run_CPU_seconds'],a['NRMSE_mean']),xytext=(4,4),textcoords='offset points',fontsize=8)
                axes[0].text(j,-.12,f"{a['available']}/{len(cases)}",ha='center',transform=axes[0].get_xaxis_transform(),fontsize=8)
            axes[0].set_xticks(range(8),[LABELS[m] for m in METHODS],rotation=30);axes[0].set_ylabel('Amplitude NRMSE');axes[0].set_title('All objects shown; black=mean with95% object bootstrap CI')
            axes[1].set_xlabel('Mean whole child CPU seconds (all runs)');axes[1].set_ylabel('Mean NRMSE, available only');axes[1].set_title('Observed cost-quality; no equal-spend claim')
            for ax in axes:ax.grid(alpha=.2)
            f.suptitle(stage+' / '+selector+'; availability denominators preserved');f.tight_layout();save(f,stage+'_'+selector+'_endpoint_cost_quality')
            f,axes=plt.subplots(2,1,figsize=(8,6),sharex=True)
            for m in METHODS:
                a=[next(r for r in aggregate if r['stage']==stage and r['selector']==selector and r['view']==cap and r['method']==m) for cap in CAPS]
                axes[0].plot(CAPS,[np.nan if r['NRMSE_mean'] is None else r['NRMSE_mean'] for r in a],'-o',ms=3,color=colors[m],label=LABELS[m])
                axes[1].plot(CAPS,[r['available']/len(cases) for r in a],'-o',ms=3,color=colors[m])
            axes[0].set_ylabel('Mean NRMSE, available only');axes[0].legend(ncol=4,fontsize=8)
            axes[1].set_ylabel('Availability / complete roster');axes[1].set_ylim(-.03,1.03);axes[1].set_xlabel('Allowed paid CPU prefix seconds (not new runs / not additive)');axes[1].set_xticks(CAPS)
            for ax in axes:ax.grid(alpha=.25)
            f.suptitle(stage+' / '+selector+' frozen CPU prefixes');f.tight_layout();save(f,stage+'_'+selector+'_CPU_PREFIXES')
    f,axes=plt.subplots(1,2,figsize=(11,4.5))
    for ax,selector in zip(axes,('primary','secondary')):
        for y,control in enumerate(m for m in METHODS if m!='dr'):
            p=next(r for r in pairs if r['stage']=='confirm24' and r['selector']==selector and r['view']=='endpoint' and {r['left'],r['right']}=={'dr',control})
            sign=1 if p['left']=='dr' else -1;mid=p['NRMSE_difference_mean'];lo=p['NRMSE_difference_mean_low'];hi=p['NRMSE_difference_mean_high']
            if mid is not None:
                mid,lo,hi=(mid,lo,hi) if sign==1 else (-mid,-hi,-lo)
                ax.errorbar(mid,y,xerr=[[mid-lo],[hi-mid]],fmt='o',color=colors[control],capsize=3)
            ax.text(.98,y,f"n={p['paired_available']}/24",transform=ax.get_yaxis_transform(),ha='right',va='center',fontsize=8)
        ax.set_yticks(range(7),[LABELS[m] for m in METHODS if m!='dr']);ax.axvline(0,color='gray',ls='--');ax.set_title(selector);ax.set_xlabel('Mean NRMSE difference: DR minus control (negative favors DR)');ax.grid(alpha=.2)
    f.suptitle('confirm24 paired material bootstrap10000 / seed20261112 / 95% percentile CI\nSecondary continuous comparison conditional on both available; full24 roster retained');f.tight_layout();save(f,'confirm24_DR_paired_bootstrap_forest')
    # Raw proposal vs corrector/secondary, using fixed physical selectors only.
    if raw:
        cases=sorted(c for c in truth if c.startswith('r12_confirm24_'))
        views=('raw_proposal_minPhi','raw_last_cap_eligible_audited_shadow','primary','secondary')
        f,axes=plt.subplots(24,5,figsize=(9,33));f.subplots_adjust(hspace=.45,wspace=.04,top=.97,bottom=.02)
        for i,c in enumerate(cases):
            for j,v in enumerate(('GT',)+views):
                ax=axes[i,j];ax.set_xticks([]);ax.set_yticks([])
                arr=truth[c] if v=='GT' else materials.get(('confirm24',v,c,'dr','endpoint'))
                if arr is not None:ax.imshow(arr.reshape(24,24),vmin=0,vmax=20,cmap='viridis',origin='lower',interpolation='nearest')
                else:ax.text(.5,.5,'unavailable',ha='center',va='center',transform=ax.transAxes,fontsize=6)
                if j==0:ax.set_ylabel(c.rsplit('_',1)[-1],rotation=0,labelpad=15,fontsize=7)
                if i==0:ax.set_title(['GT','raw minPhi','raw last audited','primary','secondary'][j],fontsize=9)
        f.suptitle('All24 DR objects: saved raw physical-audited shadows and frozen final selections',fontsize=11);save(f,'confirm24_DR_RAW_PRIMARY_SECONDARY_ALL24')

def report(out,fig,data):
    rows,prefix,aggregates,pairs,drdiff,raw,truth,meta,materials=data
    lines=['# R12 独立材料：封存后离线评分','',
        '本报告只给已冻结实验的数值、图像与范围；科学结论由主线程审查。确认24与开发6分开，主规则与新增次规则分开。所有材料均保留；unavailable不删对象、不补主规则结果。0次新正演、0次新优化/反演。','',
        '确认主规则在24生成前冻结；次规则在24数据生成后、首个确认逆解/GT评分前登记。次规则是新的候选选择估计器，不是纯DR机制增益。开发6未实施次规则，全部明确标NA。','',
        '指标：RMSE=||m−GT||₂/24；amplitude NRMSE=RMSE/(maxGT−.1)；relativeL2=||m−GT||₂/||GT||₂；IoU使用固定χ>1超水平集，连续材料中它不代表真实物理support。旧描述成功为继承合法box(1e−12数值容差)、NRMSE≤.15且IoU≥.7，不能替代连续指标。高/低相误差只在结构化二相生成区计算，连续随机族不强行分相。','',
        '配对bootstrap：10000次、seed20261112，以完整独立材料ID作cluster、无分层重抽；所有方法使用同一重抽索引。次规则NA留在全对象重抽册中；连续差异仅在该次重抽的双侧可用对象上计算，必须与paired_n/availability一起阅读；没有把缺失误差填0。95%为percentile区间，未作多重比较校正。','',
        '路径修正：compact manifest中stage字段含完整Windows相对路径；评分忽略该字段，以compact_primary先将反斜杠换成斜杠再取第一段。原manifest及冻结源未改。','']
    for stage in ('confirm24','dev6'):
        lines+=['## '+stage,'','|规则|方法|可用/全体|旧描述成功/全体|均值NRMSE(可用)|中位NRMSE|均值relativeL2|球内但材料差|球外|均值整child CPU秒|','|---|---|---:|---:|---:|---:|---:|---:|---:|---:|']
        for selector in ('primary','secondary'):
            for m in METHODS:
                a=next(r for r in aggregates if r['stage']==stage and r['selector']==selector and r['view']=='endpoint' and r['method']==m)
                def v(x):return 'NA' if x is None else f'{x:.6f}'
                lines.append(f"|{selector}|{LABELS[m]}|{a['available']}/{a['objects']}|{a['descriptive_success']}/{a['objects']}|{v(a['NRMSE_mean'])}|{v(a['NRMSE_median'])}|{v(a['relativeL2_mean'])}|{a['inside_CN99_material_error']}|{a['outside_CN99']}|{a['mean_external_whole_run_CPU_seconds']:.3f}|")
        lines+=['','次规则均值条件于可用子集，不能拿它与主规则全体均值直接声称改进；不同方法的可用集合也可能不同。完整缺失册、逐对象同物理数据判别和共同可用配对差异均在RESULTS中。','']
    lines+=['## DR逐对象物理诊断','',
        '物理分类区分：CN99球内且旧描述材料标准未通过、球外、球内且旧标准通过、unavailable。球内“材料差”是描述标签，连续NRMSE/relativeL2/IoU原值都保留；球外也不自动等于未成像。','']
    for selector in ('primary','secondary'):
        dr=[r for r in rows if r['stage']=='confirm24' and r['selector']==selector and r['method']=='dr']
        failures=[r for r in dr if r['available'] and not r['descriptive_success']]
        inside=[r for r in failures if r['noise_feasible']];outside=[r for r in failures if not r['noise_feasible']]
        lines.append(f"- {selector} DR：可用{sum(r['available'] for r in dr)}/24；旧描述未通过{len(failures)}个中，球内{len(inside)}、球外{len(outside)}；另有NA{sum(not r['available'] for r in dr)}个。球内未过对象：{', '.join(r['case_id'] for r in inside) or '无'}。球外对象：{', '.join(r['case_id'] for r in outside) or '无'}。")
    lines+=['','DR相对Dual/G19的全部对象差异在DR_vs_Dual_G19_per_object；没有另设事后“实质改善”门槛。下面按绝对NRMSE差排序列出双方前三个最大改善/退化，完整24册仍全部保留。','']
    for selector in ('primary','secondary'):
        for control in ('cie_dual','fullG19'):
            d=[r for r in drdiff if r['stage']=='confirm24' and r['selector']==selector and r['control']==control and r.get('DR_minus_control_NRMSE') is not None]
            improve=sorted([r for r in d if r['DR_minus_control_NRMSE']<0],key=lambda r:r['DR_minus_control_NRMSE'])[:3]
            worsen=sorted([r for r in d if r['DR_minus_control_NRMSE']>0],key=lambda r:-r['DR_minus_control_NRMSE'])[:3]
            fmt=lambda xs:'；'.join(f"{r['case_id']} ΔN={r['DR_minus_control_NRMSE']:+.6f}" for r in xs) or '无'
            lines.append(f"- {selector} DR−{LABELS[control]}：最大改善 {fmt(improve)}；最大退化 {fmt(worsen)}。")
    lines+=['','原始DR报告采用两种非GT选择：proposal阶段cap内已物理审计材料的最小原Phi，以及最后一个cap内已审计C-shadow；不是再跑一次DR，也未重建未审计的next-step shadow。源projection audit材料已从远端原记录抽取，全部audit评分在raw_DR_all_audits；与primary/secondary分开的固定raw选择在raw_DR_fixed_selections。','']
    for view in ('raw_proposal_minPhi','raw_last_cap_eligible_audited_shadow'):
        rs=[r for r in raw if r['stage']=='confirm24' and r['view']==view]
        if rs:
            lines.append(f"- {view}：{len(rs)}/24可用；旧描述成功{sum(r['descriptive_success'] for r in rs)}/24；均值NRMSE={np.mean([r['NRMSE'] for r in rs]):.6f}；球内材料差{sum(r['physical_material_class']=='INSIDE_CN99_MATERIAL_ERROR' for r in rs)}、球外{sum(r['physical_material_class']=='OUTSIDE_CN99' for r in rs)}。")
    lines+=['','## 成本、环境和文件','',
        '确认192个有效child外部GetProcessTimes CPU为28147.796875秒；开发48为7189.6875秒。13次零物理启动失败另计.90625秒，失败源与receipt保留。主/次规则共享同一次已付费运行，不应将二者成本相加。576LU与rank19smallLU、FFT等不同维度操作不能合并；外部CPU已经包含这些工作。','',
        '所有实验保留240CPU总上限、proposal144、corrector96、单BLAS线程和8个固定P核；不同方法/自然停止的实际消耗不同。30/60/120/240图为已付费前缀，不是四次新实验，也不是可相加的耗时。原minPhi主规则没有被本次评分改变。','',
        'Windows五个baseline wrapper源码显式设自身BelowNormal；search继承controller。EXECUTION_PRIORITY_ENVIRONMENT.json只是一张中途actual nice快照（当时baseline BelowNormal、所见controller/search Normal），不能证明整队列priority历史。开始/近结束/空闲结束频率快照同样只是成本环境采样，不证明某个硬件限频原因或严格同吞吐。没有改电源、频率、并发或线程。','',
        'NN使用冻结R8 physics权重，含既有监督训练信息；训练成本不属于本次inference CPU。multistart源码先用uniform.1到proposal deadline(或自然停止)，之后Gaussian扰动RMS1/4/10循环、每次最多24CPU/1000iter；最后再同corrector。NN-direct proposal只有feature/inference，随后最多96CPU corrector，未用闲余proposal补新优化。实际算法定义、耗时与费用须共同阅读。','',
        f'详细表：{out.relative_to(ROOT)}；科学图：{fig.relative_to(ROOT)}。每个JSON对应CSV；per_object含每个失败与NA，prefix_per_object含四个固定前缀，aggregate/pairs包含全体分母及条件bootstrap说明，per_object_pairs与DR专表给全部逐对象差异。','',
        '全24网格：confirm24_primary_ALL_MATERIALS与confirm24_secondary_ALL_MATERIALS；另有四家族逐页大图。灰格为NA、红caption为旧描述标准未过；所有对象都画出。DR_RAW_PRIMARY_SECONDARY_ALL24并列GT、两个raw定义及两个最终选择。','',
        '此处不作新颖性、普适优越性或连续Maxwell精度判断；主线程据公式、完整数据、可用性与科学边界审查。']
    (ROOT/'CONFIRMATION_RESULTS_ZH.md').write_text('\n'.join(lines)+'\n')

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,default=ROOT/'RESULTS/offline_scoring_v1')
    ap.add_argument('--figures',type=Path,default=ROOT/'FIGURES/offline_scoring_v1');a=ap.parse_args()
    start=time.perf_counter();data=score(a.output);figures(a.output,a.figures,data);report(a.output,a.figures,data)
    dump(a.output/'SCORING_AND_FIGURES_COST.json',dict(status='COMPLETE',own_process_CPU_seconds_near_write=time.process_time(),wall_seconds=time.perf_counter()-start,
        new_FW=0,new_optimizer_calls=0,GT='EXPLICIT_ROOT_POSTSEAL_OFFLINE_AUTHORIZATION',monetary_cost=None))
    print(json.dumps(dict(status='COMPLETE_SCORE_FIGURES_REPORT',results=str(a.output),figures=str(a.figures),report=str(ROOT/'CONFIRMATION_RESULTS_ZH.md'))),flush=True)
if __name__=='__main__':main()
