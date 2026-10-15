"""Independent Research G scorer and gallery; no Maxwell/search imports.

Metrics are frozen in configs/OFFLINE_SCORING_FROZEN.json. The input adapter uses
the real online runner's explicitly supplied frozen manifest contract.
The online program must not import this truth-bearing scorer.
"""
from __future__ import annotations
import argparse
import csv
import json
import os
from pathlib import Path
import sys
import time
import numpy as np
from scipy.ndimage import label

ROOT=Path(__file__).resolve().parents[3]
STAGE=ROOT/'research/research_g_20261015'


def dump(path, value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def support_mask(field, background):
    return (abs(np.asarray(field)-background)/1.2>.25).reshape(24,24)


def hole_masks(field, background):
    foreground=support_mask(field,background)
    components,n=label(~foreground,structure=np.ones((3,3),int))
    boundary=set(np.r_[components[0],components[-1],components[:,0],components[:,-1]].tolist())
    return [components==i for i in range(1,n+1) if i not in boundary and np.count_nonzero(components==i)>=6]


def field_metrics(field, truth, background):
    field,truth,background=map(np.asarray,[field,truth,background])
    if any(x.size!=576 for x in [field,truth,background]):raise ValueError('A30_GRID_SIZE_MISMATCH')
    if not all(np.isfinite(x).all() for x in [field,truth,background]):raise ValueError('NONFINITE_MATERIAL_FIELD')
    field,truth,background=[x.reshape(576) for x in [field,truth,background]]
    error=field-truth
    historical_distance=float(np.sqrt(np.mean(abs(error)**2))/1.2)
    primary_distance=float(np.sqrt(np.mean((error.real/4.5)**2+(error.imag/.5)**2)))
    support,truth_support=support_mask(field,background),support_mask(truth,background)
    union=np.count_nonzero(support|truth_support)
    support_iou=float(np.count_nonzero(support&truth_support)/union) if union else 1.
    predicted_holes,true_holes=hole_masks(field,background),hole_masks(truth,background)
    covered=[any(np.count_nonzero(t&p)/np.count_nonzero(t)>=.5 for p in predicted_holes) for t in true_holes]
    return dict(real_MAE=float(np.mean(abs(error.real))),imag_MAE=float(np.mean(abs(error.imag))),complex_MAE=float(np.mean(abs(error))),
        historical_material_distance=historical_distance,G_material_distance=primary_distance,
        historical_material_coverage=bool(historical_distance<=.1),G_material_coverage=bool(primary_distance<=.1),
        support_IoU=support_iou,support_coverage=bool(support_iou>=.8),
        predicted_digital_holes=len(predicted_holes),truth_digital_holes=len(true_holes),
        truth_hole_recall=float(np.mean(covered)) if covered else None,
        all_truth_holes_covered=bool(all(covered)) if covered else None,
        hole_recall_status='DEFINED' if covered else 'NOT_APPLICABLE_NO_TRUTH_HOLES',
        hole_count_exact_match=bool(len(predicted_holes)==len(true_holes)))


def prefix_metrics(online_selected_ids, metrics_by_id, k):
    ids=list(online_selected_ids[:k])
    if len(set(ids))!=len(ids):raise ValueError('DUPLICATE_ONLINE_TOPK_ID')
    missing=[i for i in ids if i not in metrics_by_id]
    if missing:raise ValueError('TOPK_REFERS_TO_UNVALIDATED_OR_MISSING_CANDIDATE: '+str(missing))
    values=[metrics_by_id[i] for i in ids]
    return dict(requested_K=k,available_K=len(ids),archive_shortfall=bool(len(ids)<k),candidate_ids=ids,
        Best_real_MAE=min([x['real_MAE'] for x in values],default=None),
        Best_imag_MAE=min([x['imag_MAE'] for x in values],default=None),
        Best_complex_MAE=min([x['complex_MAE'] for x in values],default=None),
        historical_material_coverage=any(x['historical_material_coverage'] for x in values),
        G_material_coverage=any(x['G_material_coverage'] for x in values),
        support_coverage=any(x['support_coverage'] for x in values),
        truth_hole_coverage=any(x['all_truth_holes_covered'] is True for x in values) if any(x['all_truth_holes_covered'] is not None for x in values) else None,
        no_valid_candidate=not bool(values))


def load_exposed_truth(a30_root, case_ids):
    """Only call in this independent offline process after online output closes."""
    source=Path(a30_root).resolve()/'inputs/a29_release/results/a29/OFFLINE_TARGETS.npz'
    if any('locked' in p.lower() or 'test'==p.lower() for p in source.parts):raise ValueError('FORBIDDEN_TRUTH_SOURCE')
    with np.load(source,allow_pickle=False) as z:
        ids=z['case_id'].astype(str).tolist()
        missing=[c for c in case_ids if c not in ids]
        if missing:raise ValueError('MISSING_EXPOSED_CASE_TRUTH: '+str(missing))
        truth_array=z['truth']
        truths={c:truth_array[ids.index(c)].copy() for c in case_ids}
    return truths,str(source)


def _path(value):
    value=Path(value)
    if value.is_absolute():
        return value.resolve()
    # Online locations are phase-relative; permit repository-prefixed paths
    # only when their prefix is explicitly present, without probing fallbacks.
    base = ROOT if value.parts and value.parts[0] == 'research' else STAGE
    return (base/value).resolve()


def _flatten_csv(path,rows):
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with Path(path).open('w',newline='') as f:
        writer=csv.DictWriter(f,keys);writer.writeheader()
        for row in rows:writer.writerow({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(dict,list)) else v for k,v in row.items()})


def G_material_distance(a,b):
    d=np.asarray(a)-b
    return float(np.sqrt(np.mean((d.real/4.5)**2+(d.imag/.5)**2)))


def G_structure_distance(a,b,background):
    from scipy.ndimage import gaussian_filter
    aa,bb=[(abs(np.asarray(x)-background)/1.2).reshape(24,24) for x in [a,b]]
    masks_xor=int(np.count_nonzero((aa>.25)!=(bb>.25)))
    views=[]
    for amplitude in [aa,bb]:
        views.append(np.array([np.clip(((gaussian_filter(amplitude,s,mode='constant',cval=0.) if s else amplitude)-.20)/.10,0,1) for s in [0.,1.,2.]]))
    d=float(np.sqrt(np.mean((views[0]-views[1])**2)))
    return dict(dS=d,mask_XOR_cells=masks_xor,structure_distinct=bool(d>=.08 and masks_xor>=6))


def _finite_json(value):
    if isinstance(value,float) and not np.isfinite(value):raise ValueError('NONFINITE_OUTPUT_METRIC')
    if isinstance(value,dict):
        for v in value.values():_finite_json(v)
    if isinstance(value,list):
        for v in value:_finite_json(v)


def score(manifest,a30_root,output):
    manifest,output=Path(manifest).resolve(),Path(output).resolve()
    original_text=manifest.read_text()
    online=json.loads(original_text)
    if 'COMPLETE' not in str(online.get('status','')) and 'FROZEN' not in str(online.get('status','')):
        raise ValueError('ONLINE_MANIFEST_NOT_COMPLETE_OR_FROZEN')
    cases=online['cases']
    if not all(isinstance(c,str) for c in cases):raise ValueError('CASE_IDS_MUST_BE_FROZEN_STRINGS')
    allowed={'p250069_weak_a1_snr20','p250073_weak_a1_snr20','p250077_weak_a1_snr20'}
    if not set(cases)<=allowed or len(set(cases))!=len(cases):raise ValueError('UNREGISTERED_OR_DUPLICATE_CASE')
    if set(online['rounds'])!={'first','quality'}:raise ValueError('UNREGISTERED_ROUNDS')
    for method in online['methods']:
        if method['case_id'] not in cases:raise ValueError('METHOD_CASE_NOT_REGISTERED')
    output.mkdir(parents=True,exist_ok=True)
    config=json.loads((STAGE/'configs/OFFLINE_SCORING_FROZEN.json').read_text())
    dump(output/'ONLINE_INPUT_SNAPSHOT.json',online)
    truths,truth_source=load_exposed_truth(a30_root,cases)
    cache={};array_access=[]
    def archive(path):
        path=_path(path)
        if path not in cache:
            cache[path]=np.load(path,allow_pickle=False)
            array_access.append(str(path))
        return cache[path]
    def field(record,prefix=None):
        z=archive(record['array_path']);key=(prefix or record['array_prefix'])+'__chi'
        x=z[key]
        if x.shape!=(576,) or not np.iscomplexobj(x) or not np.isfinite(x).all():raise ValueError('ILLEGAL_SAVED_FIELD_FORMAT: '+key)
        return x
    base_by_case={b['case_id']:b for b in online['bases']}
    backgrounds={}
    for case in cases:
        z=archive(base_by_case[case]['array_path'])
        if '__background' in z:backgrounds[case]=z['__background']
        else:
            ep=next((e for m in online['methods'] if m['case_id']==case for e in m['episodes']),None)
            if ep is None:raise ValueError('MISSING_KNOWN_BACKGROUND: '+case)
            backgrounds[case]=archive(ep['array_path'])['__background']
    methods_output=[];candidate_rows=[];metric_rows=[];episode_rows=[]
    for method in online['methods']:
        case,round_id,method_id=method['case_id'],method['round_id'],method['method']
        if round_id not in online['rounds']:raise ValueError('UNREGISTERED_ROUND')
        bg,truth=backgrounds[case],truths[case]
        success=method['success_ledger'];ids=[c['candidate_id'] for c in success]
        if len(set(ids))!=len(ids):raise ValueError('DUPLICATE_SUCCESS_LEDGER_ID')
        by_id={c['candidate_id']:c for c in success};metrics_by_id={}
        for candidate in success:
            receipt = candidate.get('validation_receipt')
            if not ((isinstance(receipt,str) and receipt.strip()) or isinstance(receipt,dict)):
                raise ValueError('MISSING_VALIDATION_RECEIPT')
            cm=field_metrics(field(candidate),truth,bg);metrics_by_id[candidate['candidate_id']]=cm
            candidate_rows.append(dict(case_id=case,parent_id=int(case.split('_')[0][1:]),round_id=round_id,method=method_id,
                candidate_id=candidate['candidate_id'],is_base=candidate['is_base'],full_residual=candidate['full_residual'],
                active_archive_member=candidate['candidate_id'] in method['active_archive_ids'],
                material_representative=candidate['candidate_id'] in method['material_representative_ids'],
                structural_representative=candidate['candidate_id'] in method['structural_representative_ids'],**cm))
        for name in ['active_archive_ids','material_representative_ids','structural_representative_ids']:
            if not set(method[name])<=set(ids):raise ValueError('ARCHIVE_ID_OUTSIDE_VALID_SUCCESS_LEDGER: '+name)
        top1=method['top1_id']
        if top1 is not None and top1 not in metrics_by_id:raise ValueError('TOP1_NOT_FULL_VALIDATED')
        km={}
        for k in [2,4,8]:
            selected=method['topk_ids'][str(k)]
            # The frozen coverage selector preserves a legal base from the
            # immutable success ledger even after active representative repair.
            allowed_topk=set(method['active_archive_ids']) | {c['candidate_id'] for c in success if c['is_base']}
            if not set(selected)<=allowed_topk:raise ValueError('TOPK_OUTSIDE_ACTIVE_ARCHIVE_AND_LEGAL_BASE')
            value=prefix_metrics(selected,metrics_by_id,k);km[str(k)]=value
            metric_rows.append(dict(case_id=case,parent_id=int(case.split('_')[0][1:]),round_id=round_id,method=method_id,K=k,
                active_archive_count=len(method['active_archive_ids']),material_representative_count=len(method['material_representative_ids']),structural_representative_count=len(method['structural_representative_ids']),
                selected_top1_id=top1,selected_top1_complex_MAE=metrics_by_id[top1]['complex_MAE'] if top1 is not None else None,**value))
        for episode in method['episodes']:
            z=archive(episode['array_path']);raw=field(episode,episode['raw_array_prefix'])
            validation_key=episode['array_prefix']+'__chi'
            end=z[validation_key] if validation_key in z else None
            base=field(base_by_case[case]);erow=dict(case_id=case,parent_id=int(case.split('_')[0][1:]),round_id=round_id,method=method_id,
                trajectory=episode['trajectory'],endpoint_id=episode['endpoint_id'],full_valid=episode['valid'],offer_status=episode['offer_status'],
                online_new_structural_world_flag=episode['new_structural_world'],online_retention=episode['retention'],
                online_bounds_fraction=episode['bounds_fraction'],online_correction_displacement=episode['correction_displacement'],
                cost_prefix=episode['cost_prefix'],cost_delta=episode['cost_delta'],
                raw_G_distance_to_base=G_material_distance(raw,base),endpoint_G_distance_to_base=G_material_distance(end,base) if end is not None else None,
                raw_to_endpoint_G_distance=G_material_distance(raw,end) if end is not None else None,
                collapse_to_base_neighborhood=bool(G_material_distance(raw,base)>=.03 and G_material_distance(end,base)<.03) if end is not None else None,
                raw_endpoint_support_retention_IoU=field_metrics(raw,end,bg)['support_IoU'] if end is not None else None,
                raw_to_endpoint_structure=G_structure_distance(raw,end,bg) if end is not None else None,
                raw_array_key=episode['raw_array_prefix']+'__chi',endpoint_array_key=validation_key,
                source_array=str(_path(episode['array_path'])))
            episode_rows.append(erow)
        methods_output.append(dict(case_id=case,round_id=round_id,method=method_id,archive_summary=method['archive_summary'],
            active_archive_ids=method['active_archive_ids'],success_ledger_count=len(success),Top1_id=top1,Top1_metrics=metrics_by_id[top1] if top1 is not None else None,
            TopK=km,candidate_metrics=metrics_by_id))
    summaries=[]
    for round_id in online['rounds']:
        for method_id in dict.fromkeys(m['method'] for m in online['methods']):
            for k in [2,4,8]:
                rows=[r for r in metric_rows if r['round_id']==round_id and r['method']==method_id and r['K']==k]
                valid=[r for r in rows if r['Best_complex_MAE'] is not None]
                summaries.append(dict(round_id=round_id,method=method_id,K=k,expected_cases=len(cases),available_cases=len(rows),
                    missing_cases=sorted(set(cases)-{r['case_id'] for r in rows}),no_valid_cases=sum(r['no_valid_candidate'] for r in rows),
                    mean_Best_complex_MAE_conditional_on_valid=float(np.mean([r['Best_complex_MAE'] for r in valid])) if valid else None,
                    mean_Best_real_MAE_conditional_on_valid=float(np.mean([r['Best_real_MAE'] for r in valid])) if valid else None,
                    mean_Best_imag_MAE_conditional_on_valid=float(np.mean([r['Best_imag_MAE'] for r in valid])) if valid else None,
                    all_available_case_historical_coverage=float(np.mean([r['historical_material_coverage'] for r in rows])) if rows else None,
                    all_available_case_G_material_coverage=float(np.mean([r['G_material_coverage'] for r in rows])) if rows else None,
                    all_available_case_support_coverage=float(np.mean([r['support_coverage'] for r in rows])) if rows else None,
                    expected_case_historical_coverage_lower_bound=sum(r['historical_material_coverage'] for r in rows)/len(cases),
                    archive_shortfall_cases=sum(r['archive_shortfall'] for r in rows)))
    paired=[]
    for row in metric_rows:
        if row['round_id']!='first':continue
        match=next((r for r in metric_rows if r['round_id']=='quality' and r['case_id']==row['case_id'] and r['method']==row['method'] and r['K']==row['K']),None)
        if match is None:continue
        paired.append(dict(case_id=row['case_id'],method=row['method'],K=row['K'],
            delta_quality_minus_first_Best_complex_MAE=match['Best_complex_MAE']-row['Best_complex_MAE'] if match['Best_complex_MAE'] is not None and row['Best_complex_MAE'] is not None else None,
            delta_historical_coverage=int(match['historical_material_coverage'])-int(row['historical_material_coverage']),
            delta_material_representatives=match['material_representative_count']-row['material_representative_count'],
            delta_structural_representatives=match['structural_representative_count']-row['structural_representative_count'],
            interpretation='Descriptive paired exposed parents; no independent confirmation or superiority inference'))
    results=dict(status='COMPLETE_INDEPENDENT_OFFLINE_SCORE',scientific_scope=config['scope'],probability='candidate_only',posterior_weights=None,
        online_manifest=str(manifest),online_status=online['status'],target=online['target'],truth_source=truth_source,config=config,
        methods=methods_output,summary=summaries,paired_differences=paired,resources_copied_not_recharged=online['resources'],
        no_new_physics=True,no_new_discovery_claim=True)
    _finite_json(results)
    dump(output/'METRICS.json',results)
    _flatten_csv(output/'PARENT_METHOD_METRICS.csv',metric_rows)
    _flatten_csv(output/'CANDIDATE_METRICS.csv',candidate_rows)
    _flatten_csv(output/'RAW_ENDPOINT_ATTRIBUTION.csv',episode_rows)
    _flatten_csv(output/'SUMMARY_METRICS.csv',summaries)
    _flatten_csv(output/'PAIRED_FIRST_QUALITY.csv',paired)
    gallery=write_gallery(online,backgrounds,truths,field,archive,STAGE/'MULTISOLUTION_GALLERY')
    for z in cache.values():z.close()
    if manifest.read_text()!=original_text:raise RuntimeError('ONLINE_MANIFEST_CHANGED_DURING_SCORING')
    receipt=dict(online_input_unchanged=True,online_array_files_read_only=True,source_truth_modified=False,truth_source=truth_source,
        online_manifest=str(manifest),cases=cases,metric_rows=len(metric_rows),candidate_metric_rows=len(candidate_rows),episodes=len(episode_rows),
        array_files_accessed=array_access,gallery=gallery,new_full_state_attempts=0,new_complex_full_solver_RHS=0)
    dump(output/'OFFLINE_ACCESS_RECEIPT.json',receipt)
    return receipt


def write_gallery(online,backgrounds,truths,field,archive,destination):
    import html
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    destination=Path(destination);destination.mkdir(parents=True,exist_ok=True)
    images=[];pages=[]
    labels={'raw':'Raw proposal','correction1':'Correction 1','correction2':'Correction 2','validation':'Independent full validation'}
    def plot_fields(columns,filename,heading):
        fig,axes=plt.subplots(2,len(columns),figsize=(3.15*len(columns),6),squeeze=False,layout='constrained')
        for j,(title,value) in enumerate(columns):
            for row,(component,vmin,vmax,cmap) in enumerate([('real',-.5,4.,'viridis'),('imag',0.,.5,'magma')]):
                shown=np.asarray(value).reshape(24,24).real if component=='real' else np.asarray(value).reshape(24,24).imag
                im=axes[row,j].imshow(shown,origin='lower',vmin=vmin,vmax=vmax,cmap=cmap,interpolation='nearest')
                axes[row,j].set_title(title+(' | Re chi' if row==0 else ' | Im chi'),fontsize=9)
                axes[row,j].set_xlabel('grid x');axes[row,j].set_ylabel('grid y')
            if j==len(columns)-1:
                for row in [0,1]:fig.colorbar(axes[row,j].images[0],ax=axes[row,:].tolist(),shrink=.8,label='chi, dimensionless')
        fig.suptitle(heading,fontsize=10)
        fig.savefig(destination/filename,dpi=110);plt.close(fig);images.append(filename)
    for case in online['cases']:
        base=next(b for b in online['bases'] if b['case_id']==case)
        plot_fields([('Known background',backgrounds[case]),('Full-validated base',field(base)),('OFFLINE truth; scoring only',truths[case])],case+'__reference.png',case+' | fixed physical color bounds')
    for method in online['methods']:
        identity=method['round_id']+'__'+method['case_id']+'__'+method['method']
        filename=identity+'.html';pages.append(filename)
        snippets=[]
        for episode in method['episodes']:
            z=archive(episode['array_path']);columns=[];missing=[]
            for prefix in ['raw','correction1','correction2','validation']:
                if prefix+'__chi' in z:columns.append((labels[prefix],z[prefix+'__chi']))
                else:missing.append(prefix)
            image_name=identity+'__t'+str(episode['trajectory'])+'.png'
            if columns:
                plot_fields(columns,image_name,episode['endpoint_id']+' | valid='+str(episode['valid'])+' | '+episode['offer_status'])
                snippets.append('<h2>'+html.escape(episode['endpoint_id'])+'</h2><p>Full-valid: '+str(episode['valid'])+'; archive status: '+html.escape(episode['offer_status'])+'; missing stages: '+html.escape(', '.join(missing) or 'none')+'</p><img src="'+image_name+'" alt="Actual raw, corrected and independently validated material fields">')
            else:snippets.append('<p>'+html.escape(episode['endpoint_id'])+': no saved material stage; MISSING.</p>')
        body=''.join(snippets)
        (destination/filename).write_text('<!doctype html><meta charset="utf-8"><title>'+html.escape(identity)+'</title><style>body{font:16px system-ui;max-width:1350px;margin:2rem auto;padding:0 1rem}img{width:100%}code{overflow-wrap:anywhere}</style><p><a href="index.html">Gallery index</a></p><h1>'+html.escape(identity)+'</h1><p>Candidate-only. Every image uses identical physical Re/Im chi color limits. Failed raw proposals remain visible. Truth never guides online selection.</p><p>Frozen active archive IDs: <code>'+html.escape(', '.join(method['active_archive_ids']))+'</code></p>'+body)
    links=''.join('<li><a href="'+f+'">'+html.escape(f[:-5])+'</a></li>' for f in pages)
    refs=''.join('<h2>'+html.escape(c)+'</h2><img src="'+c+'__reference.png" alt="Known background, base and offline truth">' for c in online['cases'])
    (destination/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>Research G multisolution gallery</title><style>body{font:16px system-ui;max-width:1200px;margin:2rem auto;padding:0 1rem}img{width:100%}</style><h1>Research G actual material gallery</h1><p>Three exposed 24x24 scalar-TM parents. Candidate-only; no posterior weights. Actual raw/corrected/full-validation arrays, including failures. Shared Re chi [-.5,4] and Im chi [0,.5] scales. Offline truth is evaluation only; these are not general 3D images.</p>'+refs+'<h2>All methods and rounds</h2><ul>'+links+'</ul>')
    dump(destination/'GALLERY_MANIFEST.json',dict(status='ACTUAL_SAVED_FIELDS',case_count=len(online['cases']),method_round_pages=len(pages),material_figure_count=len(images),pages=pages,figures=images,
        physical_color_bounds={'real':[-.5,4.],'imag':[0.,.5]},offline_truth_labeled=True,failed_proposals_included=True,source_arrays_modified=False))
    return dict(path=str(destination),pages=len(pages),figures=len(images))

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--manifest',type=Path,required=True)
    parser.add_argument('--a30-root',type=Path,required=True)
    parser.add_argument('--output',type=Path,default=STAGE/'results/E1_OFFLINE')
    args=parser.parse_args()
    started,cpu=time.perf_counter(),time.process_time();status='FAILED'
    try:
        receipt=score(args.manifest,args.a30_root,args.output)
        print(json.dumps(receipt))
        status='COMPLETE'
    finally:
        dump(args.output/'OFFLINE_EXECUTION_RECEIPT.json',dict(status=status,controlled_execution_wall_seconds=time.perf_counter()-started,process_cpu_seconds=time.process_time()-cpu,
            pid=os.getpid(),parent_pid=os.getppid(),independent_offline_process=True,new_full_state_attempts=0,new_complex_full_solver_RHS=0,GPU_seconds=0,
            physics_adapter_modules_imported=[n for n in sys.modules if n.startswith(('a30.','research_a.','b1.physics','g_physics'))],
            original_truth_source_modified=False,online_outputs_modified=False))

if __name__=='__main__':main()
