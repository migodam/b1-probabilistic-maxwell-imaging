"""Offline fixed physical-complement anchor analysis; never runs an inverse solve.

pilot: existing scored validation pilot and baseline only, no labels read.
development: score 24 extended tasks only after all are terminal/public-sealed.
seal-confirmation: public-only 126-task seal, independent of original750 scoring.
confirmation: requires both seals/scoring and explicit private scoring authorization.
"""
from pathlib import Path
from collections import defaultdict
import argparse,json,time,sys
PROCESS_BEGIN_CPU=time.process_time();PROCESS_BEGIN_WALL=time.perf_counter()
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
sys.path.insert(0,str(REPO));sys.path.insert(0,str(Path(__file__).resolve().parent))
import analyse_warm_results as shared
from research.research_max_struct_r7.experiments.score_primary import material_metrics,noise_metrics
OUT=ROOT/'delegated/warm_anchor_analysis';FIG=ROOT/'FIGURES'
ARMS=('cie_dual','a_plus','a_plus_sourcewise')
TERMINAL=('SEALED_PUBLIC_PREDICTION','FAILED')
load=shared.load;save=shared.save;finite=shared.finite


def resolve(p):
    p=Path(p);return p if p.is_absolute() else REPO/p


def table(name,rows):
    old=shared.OUT
    try:shared.OUT=OUT;shared.table(name,rows)
    finally:shared.OUT=old


def stamp(p):
    p=resolve(p);s=p.stat()
    return dict(path=str(p.relative_to(REPO)),bytes=s.st_size,mtime_ns=s.st_mtime_ns)


def verify_stamps(seal):
    for st in seal.get('evidence',[])+[x for e in seal['entries'] for x in e['files']]:
        if stamp(st['path'])!=st:raise RuntimeError('Changed sealed public input: '+st['path'])


def tasks(name):return load(ROOT/f'configs/{name}.json')['tasks']


def public_snapshot(tt):
    entries=[]
    for t in tt:
        d=resolve(t['output_dir']);p=d/'RECEIPT.json'
        if not p.exists():raise RuntimeError('PENDING: missing receipt '+t['task_id'])
        rec=load(p)
        if rec['status'] not in TERMINAL:raise RuntimeError('PENDING: nonterminal '+t['task_id'])
        if (rec['task_id'],rec['case_id'],rec['arm'])!=(t['task_id'],t['case_id'],t['arm']):raise ValueError('Identity mismatch')
        names=['RECEIPT.json','PUBLIC_TASK.json']
        if rec['status']=='SEALED_PUBLIC_PREDICTION':names+=['prediction.npz','prediction.json','BUDGET_STATES.npz']
        files=[stamp(d/n) for n in names]
        files += [stamp(p) for p in t['public_files']]
        if t.get('initialization_file'):files.append(stamp(t['initialization_file']))
        entries.append(dict(task=t,status=rec['status'],files=files))
    return dict(entries=entries,tasks=len(tt),objects=len({t['case_id'] for t in tt}),truth_access='NONE',hash_checks='NONE')


def normalized_baselines():
    rr=load(ROOT/'RESULTS/warm_development_analysis/SCORES.json')['rows'];tt=tasks('warm_development_manifest')
    features,infer=development_prefixes()
    by={(t['case_id'],t['arm'],t.get('initialization','uniform')):t for t in tt}
    out=[]
    for r in rr:
        t=by[(r['case_id'],r['arm'],r['initialization'])];rec=load(resolve(t['output_dir'])/'RECEIPT.json');c=rec.get('physical_model_cost') or {}
        out.append(dict(r,task_id=t['task_id'],method=r['arm']+' / '+r['initialization'],cohort='VALIDATION_DEVELOPMENT',
                        hardware='Mac',status=rec['status'],original_Phi=r['Phi'],solver_cpu_seconds=r['cpu_seconds'],
                        LU=c.get('LU'),complex_RHS=c.get('actual_complex_RHS')))
        o=out[-1];nn=r['initialization'].startswith('nn');prefix=features[r['case_id']]+infer[r['case_id']] if nn else 0.
        o.update(deployment_preprocessing_hardware='Mac',deployment_prefix_cpu_seconds=prefix,warm_deployment_cpu_seconds=r['cpu_seconds']+prefix,
                 feature_cpu_seconds=features[r['case_id']] if nn else 0.,inference_cpu_seconds=infer[r['case_id']] if nn else 0.)
    return out


def development_prefixes():
    gen=load(ROOT/'audit/generation_validation.json')
    features={r['case_id']:f['cpu_seconds'] for r,f in zip(gen['per_object'],gen['feature_cost'])}
    ir=load(ROOT/'data_public/initializers/development_physics/RECEIPT.json')
    infer={r['case_id']:r['inference_cpu_seconds'] for r in ir['records']}
    return features,infer


def score_entries(seal,truths,caseinfo,baseline_costs=None):
    rows=[];budgets=[]
    features,infer=development_prefixes() if baseline_costs is None else ({},{})
    for e in seal['entries']:
        t=e['task'];ident=t['case_id'];d=resolve(t['output_dir']);rec=load(d/'RECEIPT.json');cost=rec.get('physical_model_cost') or {}
        init=t.get('initialization','uniform');isnn=init.startswith('nn');host='XINAN' if (rec.get('execution_platform') or {}).get('system')=='Windows' else 'Mac'
        r=dict(task_id=t['task_id'],case_id=ident,arm=t['arm'],initialization=init,method=t['arm']+' / '+init,
            cohort=t['cohort'],hardware=host,status=rec['status'],material_nrmse=None,otsu_iou=None,
            material_image_success=False,original_Phi=rec.get('physical_objective'),relative_residual=rec.get('physical_relative_residual'),
            solver_cpu_seconds=rec.get('cpu_seconds'),solver_wall_seconds=rec.get('wall_seconds'),
            LU=cost.get('LU'),complex_RHS=cost.get('actual_complex_RHS'),corrector_stop=(rec.get('corrector') or {}).get('stop'),
            auxiliary_iterations=rec.get('auxiliary_iterations'),selected_auxiliary_stage=rec.get('auxiliary_selected_stage'),
            current_anchor='physical_complement',new_scoring_physical_calls=0)
        if baseline_costs is None:
            prefix=features[ident]+infer[ident] if isnn else 0.
            r.update(feature_cpu_seconds=features[ident] if isnn else 0.,inference_cpu_seconds=infer[ident] if isnn else 0.,
                     deployment_prefix_cpu_seconds=prefix,deployment_preprocessing_hardware='Mac',actual_preprocessing_hardware='Mac',
                     warm_deployment_cpu_seconds=rec.get('cpu_seconds',0)+prefix,actual_workflow_cpu_seconds=rec.get('cpu_seconds',0)+prefix,
                     prefix_scope='Development original actual Mac generation/inference timings; no new timing measurement')
        else:
            b=baseline_costs[(ident,t['arm'])]
            if b['hardware']!=host:raise ValueError('Frozen same-hardware anchor/baseline mismatch')
            prefix=b['deployment_prefix_cpu_seconds']
            for k in ('feature_cpu_seconds','inference_cpu_seconds','deployment_prefix_cpu_seconds','deployment_preprocessing_hardware','actual_preprocessing_hardware'):
                r[k]=b[k]
            r.update(warm_deployment_cpu_seconds=rec.get('cpu_seconds',0)+prefix,
                     actual_workflow_cpu_seconds=rec.get('cpu_seconds',0)+b['feature_cpu_seconds']+b['inference_cpu_seconds'])
            r['cold_deployment_cpu_seconds']=rec.get('cpu_seconds',0)+b['cold_deployment_cpu_seconds']-b['solver_cpu_seconds']
        if rec['status']=='SEALED_PUBLIC_PREDICTION':
            truth,support,amp=truths[ident]
            with np.load(d/'prediction.npz',allow_pickle=False) as z:
                m=z['m'].ravel();yp=z['Y_prediction'].copy() if 'Y_prediction' in z.files else None
            r.update(material_metrics(m,truth,support,amp))
            if baseline_costs is not None and yp is not None:
                public={}
                for file in t['public_files']:
                    with np.load(resolve(file),allow_pickle=False) as z:
                        for k in ('Y','sigma','noise_radius'):
                            if k in z.files:public[k]=z[k].copy()
                r.update(noise_metrics(yp,public,caseinfo[ident]['condition']['generation_grid']==24))
            with np.load(d/'BUDGET_STATES.npz',allow_pickle=False) as z:cs=z['costs'].copy();ms=z['m'].copy()
            for index,(c,m) in enumerate(zip(cs,ms)):
                met=material_metrics(m,truth,support,amp)
                budgets.append(dict(task_id=t['task_id'],case_id=ident,cohort=t['cohort'],hardware=host,method=r['method'],index=index,
                    CPU_seconds=float(c[0])+prefix,LU=float(c[1]),complex_RHS=float(c[2]),original_Phi=float(c[3]),
                    material_nrmse=met['material_nrmse'],otsu_iou=met['otsu_iou']))
        rows.append(r)
    verify_stamps(seal)
    return rows,budgets


def comparison_pairs(rows):
    pairs=[]
    for arm in ARMS:
        pairs += [(arm+' / nn_affine_anchor',arm+' / nn'),(arm+' / nn_affine_anchor','physical_lbfgs / nn')]
    pairs += [('a_plus / uniform_affine_anchor','a_plus / uniform'),
              ('a_plus / nn_affine_anchor','a_plus / nn_joint')]
    return pairs


def summaries(rows,label):
    groups=defaultdict(list)
    for r in rows:groups[(r['cohort'],r['method'])].append(r)
    table(label+'_cohort_material.csv',[dict(cohort=c,method=m,**shared.summary(rr,False)) for (c,m),rr in sorted(groups.items())])
    groups=defaultdict(list)
    for r in rows:groups[(r['cohort'],r['hardware'],r['method'])].append(r)
    table(label+'_hardware_costs.csv',[dict(cohort=c,hardware=h,method=m,**shared.summary(rr)) for (c,h,m),rr in sorted(groups.items())])
    effects=[]
    for a,b in comparison_pairs(rows):
        for host in sorted({r['hardware'] for r in rows}):
            selected=[r for r in rows if r['method'] in (a,b) and r['hardware']==host]
            effects.extend(dict(x,hardware=host) for x in shared.bootstrap_pairs(selected,b))
    table(label+'_paired_effects.csv',effects);table(label+'_all_objects.csv',rows)


def plots(rows,label):
    cases=sorted({r['case_id'] for r in rows});methods=sorted({r['method'] for r in rows})
    with PdfPages(FIG/f'warm_anchor_{label}_objects.pdf') as pdf:
        for cohort in sorted({r['cohort'] for r in rows}):
            cc=[c for c in cases if any(r['case_id']==c and r['cohort']==cohort for r in rows)]
            fig,axes=plt.subplots(2,1,figsize=(15,9))
            for method in methods:
                by={r['case_id']:r for r in rows if r['method']==method and r['cohort']==cohort}
                vals=[by.get(c,{}).get('material_nrmse') for c in cc]
                axes[0].plot(range(len(cc)),[v if finite(v) else np.nan for v in vals],'o-',label=method,markersize=3)
            axes[0].set(yscale='log',ylabel='Amplitude normalized NRMSE',title=cohort+'; all objects including failures/missing metrics')
            axes[0].legend(fontsize=7,ncol=3)
            # Absolute solver CPU: separate line per hardware; never pool machines.
            for method in methods:
                for host in sorted({r['hardware'] for r in rows}):
                    by={r['case_id']:r for r in rows if r['method']==method and r['hardware']==host}
                    vals=[by.get(c,{}).get('solver_cpu_seconds') for c in cc]
                    if any(finite(v) for v in vals):axes[1].plot(range(len(cc)),[v if finite(v) else np.nan for v in vals],'o-',label=method+' / '+host,markersize=3)
            axes[1].set(ylabel='Actual solver CPU seconds (hardware-labelled)',xlabel='Fixed object identity')
            axes[1].legend(fontsize=6,ncol=3)
            for ax in axes:ax.set_xticks(range(len(cc)),[c[-4:] for c in cc],rotation=45)
            fig.tight_layout();pdf.savefig(fig);fig.savefig(FIG/f'warm_anchor_{label}_{cohort}.png',dpi=150);plt.close(fig)


def paired_budgets(rows,budgets,label):
    by={(r['case_id'],r['method']):r for r in rows};meta={r['task_id']:r for r in rows};tracks=defaultdict(list)
    for b in budgets:tracks[b['task_id']].append(b)
    out=[]
    def coord(b,axis):
        if axis not in ('CPU_seconds','solver_CPU_seconds'):return b.get(axis)
        r=meta[b['task_id']]
        if r.get('deployment_preprocessing_hardware')!=r['hardware']:return None
        return b['CPU_seconds']-(r.get('deployment_prefix_cpu_seconds',0) if axis=='solver_CPU_seconds' else 0)
    for a,b in comparison_pairs(rows):
        for case in sorted({r['case_id'] for r in rows}):
            ra,rb=by.get((case,a)),by.get((case,b))
            if not ra or not rb:continue
            for axis,final in [('CPU_seconds','warm_deployment_cpu_seconds'),('solver_CPU_seconds','solver_cpu_seconds'),('LU','LU'),('complex_RHS','complex_RHS')]:
                row=dict(case_id=case,cohort=ra['cohort'],hardware=ra['hardware'],method=a,reference=b,axis=axis)
                if ra['hardware']!=rb['hardware'] or not finite(ra.get(final)) or not finite(rb.get(final)):
                    out.append(dict(row,status='UNAVAILABLE_SAME_HARDWARE_FINAL_COST'));continue
                B=min(ra[final],rb[final]);row['allowance']=B;ok=True
                for prefix,r in [('method',ra),('reference',rb)]:
                    bs=tracks[r['task_id']]
                    paid=[x for x in bs if finite(coord(x,axis)) and coord(x,axis)<=B]
                    if not paid:ok=False;break
                    last=max(paid,key=lambda x:x['index']);inc=min(paid,key=lambda x:(x['original_Phi'],-x['index']))
                    above=[coord(x,axis) for x in bs if finite(coord(x,axis)) and coord(x,axis)>B]
                    row[prefix+'_paid_cost']=coord(last,axis);row[prefix+'_unused_allowance']=B-coord(last,axis)
                    row[prefix+'_checkpoint_gap']=min(above)-coord(last,axis) if above else None
                    row[prefix+'_incumbent_index']=inc['index']
                    for metric in ('original_Phi','material_nrmse','otsu_iou'):row[prefix+'_'+metric]=inc.get(metric)
                row['status']='PAIRED_PAID_INCUMBENTS' if ok else 'MISSING_PAID_CHECKPOINT'
                out.append(row)
    table(label+'_budget_matched_endpoints.csv',out)
    groups=defaultdict(list)
    for r in out:groups[(r['cohort'],r['hardware'],r['method'],r['reference'],r['axis'])].append(r)
    effects=[];rng=np.random.default_rng(2026110611)
    for (cohort,hardware,method,reference,axis),rr in sorted(groups.items()):
        for metric in ('material_nrmse','otsu_iou','original_Phi'):
            valid=[r for r in rr if finite(r.get('method_'+metric)) and finite(r.get('reference_'+metric))]
            delta=np.array([r['method_'+metric]-r['reference_'+metric] for r in valid]);ci=[None,None]
            if len(delta):ci=np.quantile(np.mean(delta[rng.integers(0,len(delta),(10000,len(delta)))],axis=1),[.025,.975]).tolist()
            effects.append(dict(cohort=cohort,hardware=hardware,method=method,reference=reference,axis=axis,metric=metric,
                objects=len(rr),paired_objects=len(delta),missing_pairs=len(rr)-len(delta),mean_paired_difference=float(np.mean(delta)) if len(delta) else None,
                percentile95_low=ci[0],percentile95_high=ci[1],sampling='10000 independent-object paired descriptive bootstrap; no significance claim'))
    table(label+'_budget_paired_effects.csv',effects)
    # Absolute at-most curves: stopped runs retain the last paid incumbent, no new compute.
    curves=[]
    for cohort,host in sorted({(r['cohort'],r['hardware']) for r in rows}):
        rr=[r for r in rows if (r['cohort'],r['hardware'])==(cohort,host)]
        for axis in ('CPU_seconds','solver_CPU_seconds','LU','complex_RHS'):
            valid=[coord(b,axis) for r in rr for b in tracks[r['task_id']] if finite(coord(b,axis))]
            if not valid:continue
            for B in np.linspace(min(valid),max(valid),60):
                for method in sorted({r['method'] for r in rr}):
                    selected=[];unused=[]
                    for r in rr:
                        if r['method']!=method:continue
                        paid=[b for b in tracks[r['task_id']] if finite(coord(b,axis)) and coord(b,axis)<=B]
                        if not paid:continue
                        selected.append(min(paid,key=lambda b:(b['original_Phi'],-b['index'])))
                        unused.append(B-coord(max(paid,key=lambda b:b['index']),axis))
                    if selected:
                        curves.append(dict(cohort=cohort,hardware=host,axis=axis,allowance=float(B),method=method,
                            objects=len(selected),expected_objects=sum(r['method']==method for r in rr),
                            mean_material_nrmse=float(np.mean([b['material_nrmse'] for b in selected])),
                            mean_original_Phi=float(np.mean([b['original_Phi'] for b in selected])),mean_unused_allowance=float(np.mean(unused)),
                            semantics='At most allowance; stopped method holds paid terminal incumbent; not equal spent cost; no interpolation/GT selection'))
    table(label+'_at_most_budget_curves.csv',curves)
    with PdfPages(FIG/f'warm_anchor_{label}_budget.pdf') as pdf:
        for cohort,host,axis in sorted({(r['cohort'],r['hardware'],r['axis']) for r in curves}):
            fig,ax=plt.subplots(figsize=(12,6))
            for method in sorted({r['method'] for r in curves}):
                rr=[r for r in curves if (r['cohort'],r['hardware'],r['axis'],r['method'])==(cohort,host,axis,method)]
                if rr:ax.plot([r['allowance'] for r in rr],[r['mean_material_nrmse'] for r in rr],label=method)
            ax.set(xlabel='At most '+axis,ylabel='Observed-object mean NRMSE',title=cohort+' / '+host+'; stopped runs retain terminal state; availability in CSV')
            ax.legend(fontsize=7,ncol=2);fig.tight_layout();pdf.savefig(fig);plt.close(fig)


def saved_baseline_budgets(truths):
    features,infer=development_prefixes();out=[]
    for t in tasks('warm_development_manifest'):
        d=resolve(t['output_dir']);ident=t['case_id'];truth,support,amp=truths[ident]
        prefix=features[ident]+infer[ident] if t.get('initialization','uniform').startswith('nn') else 0.
        with np.load(d/'BUDGET_STATES.npz',allow_pickle=False) as z:cs=z['costs'].copy();ms=z['m'].copy()
        for i,(c,m) in enumerate(zip(cs,ms)):
            met=material_metrics(m,truth,support,amp)
            out.append(dict(task_id=t['task_id'],case_id=ident,index=i,CPU_seconds=float(c[0])+prefix,LU=float(c[1]),complex_RHS=float(c[2]),
                original_Phi=float(c[3]),material_nrmse=met['material_nrmse'],otsu_iou=met['otsu_iou']))
    return out


def all_materials(rows,truths,taskmap,label):
    methods=['physical_lbfgs / nn','a_plus / nn_joint']+[arm+' / nn' for arm in ARMS]+[arm+' / nn_affine_anchor' for arm in ARMS]
    selected={(r['case_id'],r['method']):r for r in rows};fixed={f'r8_test_{i:04d}' for i in (0,3,24,27,36,39)}
    with PdfPages(FIG/f'warm_anchor_{label}_all_materials.pdf') as pdf:
        for case in sorted(truths):
            truth=truths[case][0].reshape(24,24);images={}
            for method in methods:
                r=selected.get((case,method))
                if r and r['status']=='SEALED_PUBLIC_PREDICTION':
                    with np.load(resolve(taskmap[r['task_id']]['output_dir'])/'prediction.npz',allow_pickle=False) as z:images[method]=z['m'].reshape(24,24).copy()
            mats=[truth]+list(images.values());vmin=min(.02,*(float(m.min()) for m in mats));vmax=max(.02,*(float(m.max()) for m in mats))
            fig,axes=plt.subplots(3,3,figsize=(13,13));axes=axes.ravel()
            for ax in axes:ax.axis('off')
            im=axes[0].imshow(truth,cmap='viridis',vmin=vmin,vmax=vmax);axes[0].set_title(f'GT (offline)\nmin={truth.min():.3g}; max={truth.max():.3g}')
            for ax,method in zip(axes[1:],methods):
                r=selected.get((case,method));m=images.get(method)
                if m is None:
                    reason='NOT_AVAILABLE: control not run' if r is None else ('FAILED' if r['status']=='FAILED' else 'NO_SAVED_PREDICTION')
                    ax.set_title(method+'\n'+reason,fontsize=8);continue
                ax.imshow(m,cmap='viridis',vmin=vmin,vmax=vmax)
                ax.set_title(method+f"\nNRMSE={r['material_nrmse']:.3g}; min={m.min():.3g}; max={m.max():.3g}",fontsize=8)
            fig.suptitle(f'{case}; full shared contrast range [{vmin:g}, {vmax:g}]; all predefined objects retained')
            fig.subplots_adjust(top=.88,left=.02,right=.92,bottom=.03,hspace=.25,wspace=.08)
            cax=fig.add_axes([.94,.12,.012,.68]);fig.colorbar(im,cax=cax,label='Contrast, full common range')
            pdf.savefig(fig)
            if case in fixed:fig.savefig(FIG/f'warm_anchor_confirmation_fixed_{case}.png',dpi=160)
            plt.close(fig)


def report(label,pending=False,rows=None):
    (OUT/(label+'_REPORT.md')).write_text(
        '# Fixed physical-complement anchor: '+label+'\n\n'
        'The fixed anchor retains JD+(I-PD)jphys(m0) and optimizes the existing Fourier coefficients. '
        'The experiment has no extra optimized degree of freedom; computational projection/storage are charged by solver receipts. '
        'Interpretation, scientific claims and novelty judgment remain with the root.\n\n'
        +('Only existing scored pilot8 and baseline48 rows are plotted; extended24 scoring is pending until all24 are public terminal. No labels were read for this pilot plot.\n\n' if pending else '')+
        'All predefined development/confirmation objects are retained, including failures. Original NN Dual/global A+ have matched unanchored development baselines; '
        'sourcewise A+ has NO identical unanchored development baseline. Its comparison to NN direct is descriptive, not a matched anchor ablation. '
        'The separate42-object confirmation has all three same-arm anchor/material-NN pairs and the critical global A+ '
        'anchor vs nn_joint comparison. nn_joint projects the physical initial current into the old truncated affine family; '
        'the anchor retains the full physical APCS offset. Dual/sourcewise lack an identical truncated physical-current '
        'projection control, so their anchor comparisons combine current initialization and affine-space changes and '
        'cannot isolate a pure high-mode effect. Development observations cannot establish confirmation performance.\n\n'
        'Absolute solver/warm/cold CPU are stratified by hardware. feature_cpu_seconds and inference_cpu_seconds are original Mac timings; '
        'deployment_prefix_cpu_seconds is the same-solver-hardware independent measurement for confirmation. '
        'Validation prefixes reuse actual Mac preprocessing/inference records and are labelled separately. Training remains a separate cost.\n\n'
        'Paid original-Phi incumbents determine budget states; material scores never select an iterate. Matched endpoints use per-object '
        'min(finalA,finalB), report unused allowance/checkpoint gaps, and are same allowance rather than equal actual spend. '
        'At-most-budget curves retain a stopped method’s paid terminal state without continuing computation. Missing/failed checkpoint availability is explicit. '
        'Bootstrap material effects resample independent objects; no significance claim. No new inverse solve is performed by this analysis.\n')


    if rows is not None:
        selected={(r['case_id'],r['method']):r for r in rows}
        columns=[('NNdirect','physical_lbfgs / nn'),('old_global_NN','a_plus / nn'),
                 ('old_global_NN_joint','a_plus / nn_joint'),('anchor_Dual','cie_dual / nn_affine_anchor'),
                 ('anchor_global_A+','a_plus / nn_affine_anchor'),('anchor_source_A+','a_plus_sourcewise / nn_affine_anchor'),
                 ('uniform_anchor_global','a_plus / uniform_affine_anchor')]
        text='\n## All-object final amplitude NRMSE\n\nThese are saved final paid-Phi selections, not equal-CPU claims; IoU, failures, actual costs and matched-budget results are in the accompanying tables.\n\n'
        text+='|Object|Cohort|'+'|'.join(k for k,_ in columns)+'|\n'
        text+='|---|---|'+'|'.join('---:' for _ in columns)+'|\n'
        for case in sorted({r['case_id'] for r in rows}):
            cohort=next(r['cohort'] for r in rows if r['case_id']==case)
            values=[]
            for _,method in columns:
                r=selected.get((case,method));v=r.get('material_nrmse') if r else None
                values.append(f'{v:.6f}' if finite(v) else 'FAILED / NOT_AVAILABLE')
            text+='|'+case+'|'+cohort+'|'+'|'.join(values)+'|\n'
        with (OUT/(label+'_REPORT.md')).open('a') as f:f.write(text)


def pilot():
    rows=normalized_baselines()
    for r in load(ROOT/'RESULTS/warm_anchor_development_analysis/SCORES.json')['rows']:
        c=r.get('physical_model_cost') or {}
        rows.append(dict(r,cohort='VALIDATION_DEVELOPMENT',hardware='Mac',status='SEALED_PUBLIC_PREDICTION',
             method=r['arm']+' / '+r['initialization'],original_Phi=r['Phi'],solver_cpu_seconds=r['cpu_seconds'],
             LU=c.get('LU'),complex_RHS=c.get('actual_complex_RHS')))
    summaries(rows,'PILOT');plots(rows,'pilot');report('PILOT',True)
    save(OUT/'PILOT_RECEIPT.json',dict(status='EXISTING_PILOT_METRICS_PLOTTED',baseline_rows=48,pilot_rows=8,
         truth_access='NONE: existing validation score rows only',new_physical_calls=0,extended24_scoring='PENDING'))


def development():
    tt=tasks('warm_anchor_development_manifest')+tasks('warm_anchor_development_extended_manifest')
    if len(tt)!=32:raise ValueError('Expected pilot8 plus extended24')
    seal=public_snapshot(tt);seal.update(status='ALL_ANCHOR_DEVELOPMENT_PUBLIC_TERMINAL')
    save(OUT/'DEVELOPMENT_PUBLIC_SNAPSHOT.json',seal)
    # Validation labels are accessed only after all32 public predictions are terminal.
    cases={r['case_id']:r for r in load(ROOT/'data_public/warm_development/cases.json')['records']}
    with np.load(ROOT/'data_training/validation.npz',allow_pickle=False) as z:materials=z['m'][:8].reshape(8,576).copy()
    truths={c:(materials[int(c[-4:])],materials[int(c[-4:])]>.1+.5*r['condition']['amplitude'],r['condition']['amplitude']) for c,r in cases.items()}
    rows,budgets=score_entries(seal,truths,cases);baseline=normalized_baselines()
    save(OUT/'DEVELOPMENT_SCORES.json',dict(rows=rows,truth_access='EXISTING_VALIDATION_LABELS_ONLY; post public snapshot',new_physical_calls=0))
    save(OUT/'DEVELOPMENT_BUDGET_SCORES.json',dict(rows=budgets,selection='Paid Phi incumbent only',new_physical_calls=0))
    combined=rows+baseline;summaries(combined,'DEVELOPMENT');plots(combined,'development');report('DEVELOPMENT',rows=combined)
    paired_budgets(combined,budgets+saved_baseline_budgets(truths),'DEVELOPMENT')
    taskmap={t['task_id']:t for t in tt+tasks('warm_development_manifest')}
    all_materials(combined,truths,taskmap,'development')
    save(OUT/'DEVELOPMENT_RECEIPT.json',dict(status='ALL32_VALIDATION_ANCHOR_SCORED',objects=8,anchor_rows=32,baseline_rows=48,
        sourcewise_identical_unanchored_development_baseline='NOT_AVAILABLE',truth_access='VALIDATION_ONLY',new_physical_calls=0))


def original_gate():
    seal=load(ROOT/'audit/WARM_ALL_PUBLIC_SEAL.json');rec=load(ROOT/'RESULTS/warm_scoring/RECEIPT.json')
    if seal.get('status')!='ALL_FROZEN_WARM_TASKS_SEALED' or seal.get('tasks')!=750 or rec.get('status')!='SCORED_ALL_FROZEN_PREDICTIONS':
        raise RuntimeError('Original750 complete public seal and scoring required')
    verify_stamps(seal)
    rr=load(ROOT/'RESULTS/warm_scoring/SCORES.json')['rows']
    if len(rr)!=750 or len({r['task_id'] for r in rr})!=750 or {r['task_id'] for r in rr}!={t['task_id'] for t in tasks('warm_all_manifest')}:raise ValueError('Original750 scoring frozen identities')
    return rr


def seal_confirmation():
    tt=tasks('warm_anchor_confirmation_manifest')
    if (len(tt)!=126 or len({t['task_id'] for t in tt})!=126 or len({t['case_id'] for t in tt})!=42
        or any({t['arm'] for t in tt if t['case_id']==c}!=set(ARMS) for c in {t['case_id'] for t in tt})):
        raise ValueError('Expected separate126=42x3 exact anchor arms')
    seal=public_snapshot(tt)
    evidence=[ROOT/'configs/warm_anchor_confirmation_manifest.json',ROOT/'configs/warm_anchor_extension.json',
              ROOT/'RESULTS/training/physics/RECEIPT.json',ROOT/'RESULTS/training/physics/weights.pt',ROOT/'audit/MAC_FEATURE_INFERENCE.json',ROOT/'audit/XINAN_FEATURE_INFERENCE.json',
              ROOT/'src/cie_fft.py',ROOT/'experiments/run_case.py',Path(__file__)]
    seal.update(status='ALL_ANCHOR_TASKS_SEALED',evidence=[stamp(p) for p in evidence])
    path=ROOT/'audit/WARM_ANCHOR_PUBLIC_SEAL.json'
    if path.exists():raise FileExistsError(path)
    save(path,seal)


def confirmation(authorized):
    if not authorized:raise PermissionError('Explicit root authorization required: --authorized-confirmation-scoring')
    original=original_gate();seal=load(ROOT/'audit/WARM_ANCHOR_PUBLIC_SEAL.json')
    if seal.get('status')!='ALL_ANCHOR_TASKS_SEALED' or seal.get('tasks')!=126:raise RuntimeError('Separate126 anchor public seal required')
    verify_stamps(seal)
    frozen=tasks('warm_anchor_confirmation_manifest')
    if {e['task']['task_id'] for e in seal['entries']}!={t['task_id'] for t in frozen}:raise ValueError('Anchor frozen identity mismatch')
    baselines=[r for r in original if (r['initialization']=='nn' and r['arm'] in ARMS+('physical_lbfgs',))
               or (r['arm']=='a_plus' and r['initialization']=='nn_joint')]
    costmap={(r['case_id'],r['arm']):r for r in baselines if r['initialization']=='nn'}
    cases={r['case_id']:r for r in load(ROOT/'data_public/warm_test/test_manifest.json')['records']}
    # First test truth access: both complete seals, original scoring, and authorization verified above.
    truths={}
    for ident in sorted({t['case_id'] for t in frozen}):
        amp=costmap[(ident,'a_plus')]['amplitude']
        with np.load(ROOT/f'sealed_truth/warm_test/{ident}.npz',allow_pickle=False) as z:
            truths[ident]=(z['m'].ravel().copy(),z['mask'].ravel().astype(bool),amp)
    rows,budgets=score_entries(seal,truths,cases,costmap)
    save(OUT/'CONFIRMATION_SCORES.json',dict(rows=rows,objects=42,tasks=126,new_physical_calls=0,truth_access='AUTHORIZED_POST_BOTH_PUBLIC_SEALS_AND_ORIGINAL_SCORING'))
    save(OUT/'CONFIRMATION_BUDGET_SCORES.json',dict(rows=budgets,selection='Paid Phi incumbent only',new_physical_calls=0))
    combined=baselines+rows;oldbudgets=load(ROOT/'RESULTS/warm_scoring/BUDGET_SCORES.json')['rows']
    selected={r['task_id'] for r in baselines};allbudgets=budgets+[b for b in oldbudgets if b['task_id'] in selected]
    summaries(combined,'CONFIRMATION');plots(combined,'confirmation');paired_budgets(combined,allbudgets,'CONFIRMATION');report('CONFIRMATION',rows=combined)
    taskmap={t['task_id']:t for t in frozen+tasks('warm_primary_manifest')}
    all_materials(combined,truths,taskmap,'confirmation')
    save(OUT/'CONFIRMATION_RECEIPT.json',dict(status='SEPARATE126_ANCHOR_OFFLINE_SCORED',tasks=126,objects=42,
        original750_unchanged=True,cohorts='ID24 / OOD12 / cross-grid6 reported separately',new_physical_calls=0,
        truth_access='AUTHORIZED_POST_BOTH_PUBLIC_SEALS_AND_ORIGINAL_SCORING'))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=('pilot','development','seal-confirmation','confirmation'),default='pilot',nargs='?')
    p.add_argument('--authorized-confirmation-scoring',action='store_true');a=p.parse_args()
    OUT.mkdir(parents=True,exist_ok=True);FIG.mkdir(exist_ok=True)
    receipt_name={'pilot':'PILOT_RECEIPT.json','development':'DEVELOPMENT_RECEIPT.json',
                  'seal-confirmation':'PUBLIC_SEAL_RECEIPT.json','confirmation':'CONFIRMATION_RECEIPT.json'}[a.mode]
    target=OUT/receipt_name
    if target.exists():
        previous=load(target)
        if 'cpu_seconds' not in previous:
            previous.update(cpu_seconds=None,wall_seconds=None,cost_measurement='PAST_EXECUTION_NOT_ISOLATED; no inferred cost')
            archive=OUT/(target.stem+'_prior_cost_unknown.json')
            if not archive.exists():save(archive,previous)
    error=None
    try:
        if a.mode=='pilot':pilot()
        elif a.mode=='development':development()
        elif a.mode=='seal-confirmation':seal_confirmation()
        else:confirmation(a.authorized_confirmation_scoring)
    except Exception as exc:
        error=type(exc).__name__+': '+str(exc)
        raise
    finally:
        elapsed=dict(mode=a.mode,cpu_seconds=time.process_time()-PROCESS_BEGIN_CPU,
                     wall_seconds=time.perf_counter()-PROCESS_BEGIN_WALL,
                     cost_measurement='Actual outer process elapsed including library imports after standard-library startup',
                     status='FAILED' if error else 'COMPLETE',error=error,new_inverse_solves=0)
        with (OUT/'EXECUTION_COSTS.jsonl').open('a') as f:f.write(json.dumps(elapsed)+'\n')
        if error is None:
            receipt=load(target) if target.exists() else dict(status='PUBLIC_ANCHOR_SEAL_CREATED')
            receipt.update({k:v for k,v in elapsed.items() if k not in ('status','error')});save(target,receipt)
