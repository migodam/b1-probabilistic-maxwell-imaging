"""Offline paired R12 scoring. No forward model, solver, or training imports.

Truth is opened only after --seal validates a root-authored complete 44-task seal.
Selection always uses the runner's public paid Phi/R metadata, never truth.
"""
import os
for _key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[_key] = '1'
import argparse
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
STUDY = ROOT.parent
WORKTREES = STUDY.parents[2]
METHODS = ('g19', 'dr', 'ap', 'direct')
VIEWS = ('initial', 'raw', 'corrected', 'primary', 'secondary')
REPORT_VIEWS = VIEWS + ('secondary_deploy',)
COLORS = dict(g19='#0072B2', dr='#D55E00', ap='#009E73', direct='#CC79A7')


def read(path):
    return json.loads(Path(path).read_text())


def dump(path, obj):
    def convert(x):
        if isinstance(x, np.generic): return x.item()
        if isinstance(x, np.ndarray): return x.tolist()
        if isinstance(x, Path): return str(x)
        raise TypeError(type(x).__name__)
    Path(path).write_text(json.dumps(obj, indent=2, default=convert, allow_nan=False)+'\n')


def table(out, name, rows):
    dump(out/(name+'.json'), rows)
    keys = list(dict.fromkeys(k for row in rows for k in row))
    with (out/(name+'.csv')).open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: json.dumps(v) if isinstance(v, (dict, list)) else v for k,v in row.items()})


def loadm(path, key='m'):
    with np.load(path, allow_pickle=False) as z:
        return z[key].copy().ravel()


def validate_seal(path, tasks):
    seal = read(path)
    if seal.get('status') != 'ALL44_PREDICTIONS_SEALED':
        raise RuntimeError('Require root seal status ALL44_PREDICTIONS_SEALED')
    expected = {t['task_id'] for t in tasks}
    if len(tasks) != 44 or set(seal.get('task_ids', [])) != expected:
        raise RuntimeError('Root seal must list exactly the frozen 44 task IDs')
    if not seal.get('root_authorized_offline_scoring', False):
        raise RuntimeError('Root must explicitly authorize offline scoring in seal')
    receipts = {}
    for task in tasks:
        directory = ROOT/'RESULTS'/task['task_id']
        receipt = read(directory/'RECEIPT.json')
        if receipt.get('status') not in ('COMPLETE', 'CPU_BUDGET_REACHED', 'FAILED'):
            raise RuntimeError('Nonterminal receipt '+task['task_id'])
        if receipt.get('truth_access') != 'NONE':
            raise RuntimeError('Unexpected online truth declaration')
        if receipt.get('task', {}).get('task_id') != task['task_id']:
            raise RuntimeError('Receipt task identity mismatch')
        # Failed tasks remain in all denominators; existing partial evidence is scored.
        for name in ('HISTORY.json', 'PREFIXES.json', 'CANDIDATES.npz'):
            if not (directory/name).is_file():
                raise RuntimeError('Incomplete terminal evidence '+str(directory/name))
        history = read(directory/'HISTORY.json')
        with np.load(directory/'CANDIDATES.npz', allow_pickle=False) as z:
            if len(z['m']) != len(history): raise RuntimeError('Candidate/history length mismatch')
        for view in VIEWS:
            if (directory/(view+'.npz')).exists() != (directory/(view+'.json')).exists():
                raise RuntimeError('Incomplete selected snapshot '+task['task_id']+'/'+view)
        receipts[task['task_id']] = receipt
    # Optional root stamps enforce immutable outputs without cryptographic hashes.
    for stamp in seal.get('files', []):
        p = Path(stamp['path'])
        if not p.is_absolute(): p = ROOT/p
        stat = p.stat()
        if stat.st_size != stamp['size'] or stat.st_mtime_ns != stamp['mtime_ns']:
            raise RuntimeError('Sealed file changed '+str(p))
    return seal, receipts


def load_truth(case):
    """Only called by main after complete seal validation; no private reads at import."""
    ident, cohort = case['case_id'], case['cohort']
    portable = ROOT/'offline_truth'/cohort/(ident+'.npz')
    if portable.is_file():
        with np.load(portable, allow_pickle=False) as z:
            return z['m'].copy().ravel(), z['support'].copy().ravel().astype(bool), float(z['amplitude']), str(portable)
    if cohort == 'A':
        candidates = [STUDY/'sealed_truth/basin'/(ident+'.npz'),
                      WORKTREES/'r10_worktree/research/research_max_struct_r10/sealed_truth/twod'/(ident+'.npz')]
        path = next((p for p in candidates if p.is_file()), None)
        if path is None: raise FileNotFoundError('No actual basin GT for '+ident)
        truth = loadm(path)
        return truth, truth>1, float(truth.max()-.1), str(path)
    if cohort == 'B':
        path = STUDY/'delegated/confirmation/sealed_truth/confirm24'/(ident+'.npz')
        truth = loadm(path)
        return truth, truth>1, float(truth.max()-.1), str(path)
    r8 = WORKTREES/'r8_worktree/research/research_max_struct_r8'
    path = r8/'sealed_truth/warm_test'/(ident+'.npz')
    metadata = {r['case_id']:r for r in read(r8/'sealed_truth/warm_test/score_metadata.json')['records']}
    with np.load(path, allow_pickle=False) as z:
        return z['m'].copy().ravel(), z['mask'].copy().ravel().astype(bool), float(metadata[ident]['amplitude']), str(path)


def otsu(m):
    if np.ptp(m) <= 1e-12: return np.zeros(m.size, bool), float(m[0])
    h,e = np.histogram(m, bins=256); c = (e[:-1]+e[1:])/2
    w = np.cumsum(h); s = np.cumsum(h*c); den = w*(w[-1]-w)
    score = np.full(256, -np.inf); good = den>0
    score[good] = (s[-1]*w[good]-s[good]*w[-1])**2/den[good]
    threshold = float(c[np.argmax(score)])
    return m>threshold, threshold


def metrics(m, truth, support, amplitude, cohort):
    if m is None: return dict(available=False, availability_reason='NO_SEALED_CANDIDATE')
    if m.shape != truth.shape or np.iscomplexobj(m) or not np.isfinite(m).all():
        return dict(available=False, availability_reason='INVALID_SEALED_MATERIAL')
    rmse = float(np.linalg.norm(m-truth)/np.sqrt(truth.size))
    seg, threshold = otsu(m) if cohort=='C' else (m>1., 1.)
    union = np.count_nonzero(seg|support)
    iou = float(np.count_nonzero(seg&support)/union) if union else 1.
    tol_low = float(np.spacing(np.float32(.02))) if m.dtype==np.float32 and cohort=='C' else 1e-12
    tol_high = float(np.spacing(np.float32(20.))) if m.dtype==np.float32 and cohort=='C' else 1e-12
    legal = bool(m.min()>=.02-tol_low and m.max()<=20.+tol_high)
    return dict(available=True, availability_reason='AVAILABLE', RMSE=rmse, NRMSE=rmse/amplitude,
                relativeL2=float(np.linalg.norm(m-truth)/np.linalg.norm(truth)), IoU=iou,
                amplitude=amplitude, IoU_threshold=threshold,
                IoU_rule='prediction_Otsu_vs_saved_mask' if cohort=='C' else 'fixed_chi_gt_1',
                legal_inherited=legal, material_min=float(m.min()), material_max=float(m.max()),
                descriptive_success=bool(legal and rmse/amplitude<=.15 and iou>=.7))


def parent_cpu(case):
    cost = case.get('parent_original_actual_cost', {})
    if case['cohort']=='C':
        # Shared import/model/setup is itemized once in ledger, never four times.
        p = cost.get('preprocessing_record', {})
        return sum(p.get(k, {}).get('cpu_seconds', 0.) for k in ('features', 'inference'))
    external=cost.get('completed_external_child_CPU_seconds')
    return external if external is not None else cost.get('whole_child_CPU_seconds', cost.get('cpu_seconds'))


def parent_cpu_scope(case):
    cost=case.get('parent_original_actual_cost', {})
    if case['cohort']=='C':return 'Warm feature/inference nested deployment view; shared total charged once separately'
    if cost.get('completed_external_child_CPU_seconds') is not None:return 'Completed external child-process CPU, including exit'
    return 'Original receipt CPU scope only; completed external child CPU unavailable'


def export_offline_truth(case, truth, support, amplitude, source):
    # Postseal reproduction export; never used or imported by the online runner.
    path=ROOT/'offline_truth'/case['cohort']/(case['case_id']+'.npz')
    path.parent.mkdir(parents=True,exist_ok=True)
    if not path.is_file():
        np.savez_compressed(path,m=truth,support=support,amplitude=np.asarray(amplitude))
    meta=path.with_suffix('.json')
    if not meta.exists():
        dump(meta,dict(case_id=case['case_id'],cohort=case['cohort'],source=source,
                       export_scope='OFFLINE_ONLY_AFTER_COMPLETE_ROOT_SEAL',online_truth_access='NONE'))
    return path


def controller_jobs(path):
    if path is None: return {}, None
    obj = read(path); rows = obj if isinstance(obj, list) else obj.get('jobs', obj.get('records', obj.get('tasks', [])))
    if isinstance(rows, dict): rows = [dict(task_id=k, **v) for k,v in rows.items()]
    result = {}
    for row in rows:
        ident = row.get('task_id', row.get('task', {}).get('task_id') if isinstance(row.get('task'), dict) else row.get('task'))
        if ident:
            if ident in result: raise RuntimeError('Duplicate controller cost '+ident)
            result[ident] = row
    return result, obj


def savefig(fig, folder, name):
    fig.tight_layout(rect=(0,0,.94,.97) if getattr(fig,'_material_bar',False) else (0,0,1,.97))
    fig.savefig(folder/(name+'.png'), dpi=160, bbox_inches='tight')
    fig.savefig(folder/(name+'.pdf'), bbox_inches='tight')
    import matplotlib.pyplot as plt
    plt.close(fig)


def material_colorbar(fig, error=False):
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize
    fig._material_bar=True
    axis=fig.add_axes([.955,.18,.012,.64])
    bar=fig.colorbar(ScalarMappable(norm=Normalize(-20 if error else .02,20),cmap='coolwarm' if error else 'viridis'),cax=axis)
    bar.set_label('signed χ error (m−GT), fixed ±20' if error else 'χ, fixed range [0.02, 20]')


def figures(cases, rows, prefixes, pictures, arrows, out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    # Full identifiers are explicit, including the 20 untriggered B cases.
    fig, ax = plt.subplots(figsize=(14, 10)); ax.axis('off')
    text = ['Cohort / complete case ID / label / rescue trigger']
    for c in cases:
        trigger = any(r['cohort']==c['cohort'] and r['case_id']==c['case_id'] and r.get('rescue_triggered') for r in rows)
        status=('TRIGGERED' if trigger else 'RETAIN ORIGINAL') if c['cohort']=='B' else 'FIXED PAIRED'
        text.append(f"{c['cohort']}  {c['case_id']}  {c.get('label', '')}  {status}")
    ax.text(0, 1, '\n'.join(text), va='top', family='monospace', fontsize=10)
    savefig(fig, out, 'cohort_mapping_full_IDs')
    index = {(r['cohort'],r['case_id'],r['method'],r['view']):r for r in rows}
    for cohort in ('A','C'):
        cohort_cases = [c for c in cases if c['cohort']==cohort]
        for view in ('raw','corrected'):
            fig, axes = plt.subplots(len(cohort_cases), 6, figsize=(18,3*len(cohort_cases)), squeeze=False)
            for i,c in enumerate(cohort_cases):
                ident = c['case_id']; keys = [('GT','truth'),('same initial','initial')]+[(m,m+'_'+view) for m in METHODS]
                for ax,(label,key) in zip(axes[i],keys):
                    arr = pictures.get((cohort,ident,key))
                    if arr is not None: ax.imshow(arr.reshape(24,24),vmin=.02,vmax=20,cmap='viridis')
                    title = ident+'\n'+label if key in ('truth','initial') else label
                    row = index.get((cohort,ident,label,view))
                    if label=='direct' and view=='raw':title+=' (initial; no auxiliary)' 
                    if row and row.get('available'): title += f"\nN={row['NRMSE']:.3f}, I={row['IoU']:.2f}"
                    if arr is None: title += '\nNA'
                    ax.set_title(title,fontsize=9); ax.set_xticks([]); ax.set_yticks([])
            fig.suptitle(f'{cohort}: all objects / {view}; inherited metrics; no GT selection')
            material_colorbar(fig)
            savefig(fig,out,f'{cohort}_{view}_all_objects')
    # B all24: GT, original, all four rescue strategies; shared material/error scales.
    for page in range(3):
        subset = [c for c in cases if c['cohort']=='B'][8*page:8*(page+1)]
        for error in (False,True):
            fig, axes = plt.subplots(len(subset),6,figsize=(16,2.3*len(subset)),squeeze=False)
            for i,c in enumerate(subset):
                ident=c['case_id']; truth=pictures[('B',ident,'truth')]
                for ax,key in zip(axes[i],('truth','initial')+tuple(m+'_primary' for m in METHODS)):
                    arr=pictures.get(('B',ident,key)); row=index.get(('B',ident,key.split('_')[0],'primary'))
                    if arr is not None:
                        ax.imshow((arr-truth if error else arr).reshape(24,24),cmap='coolwarm' if error else 'viridis',vmin=-20 if error else .02,vmax=20)
                    title=ident+'\n'+key if key=='truth' else ('original G19 primary' if key=='initial' else key)
                    if row and row.get('available'): title+=f" N={row['NRMSE']:.3f} {'TRIGGERED' if row['rescue_triggered'] else 'RETAIN'}"
                    if arr is None:title+=' NA'
                    ax.set_title(title,fontsize=7);ax.set_xticks([]);ax.set_yticks([])
            fig.suptitle('B all24 / original G19 + four frozen strategies / '+('signed error' if error else 'material'))
            material_colorbar(fig,error)
            savefig(fig,out,f'B_all24_{"errors" if error else "materials"}_{page+1}')
    for cohort in ('A','B','C'):
        cs=[c for c in cases if c['cohort']==cohort]
        pages=[cs[i:i+8] for i in range(0,len(cs),8)]
        for page,subset in enumerate(pages,1):
            fig,axes=plt.subplots(len(subset),2,figsize=(12,max(3,2*len(subset))),squeeze=False)
            for i,c in enumerate(subset):
                for method in METHODS:
                    rr=[r for r in prefixes if r['cohort']==cohort and r['case_id']==c['case_id'] and r['method']==method]
                    for ax,key in zip(axes[i],('NRMSE','Phi')):
                        values=[r.get(key,np.nan) if r['available'] else np.nan for r in rr]
                        ax.plot([r['CPU_cap'] for r in rr],values,'-o',label=method,color=COLORS[method],ms=3)
                        if key=='Phi' and any(np.isfinite(v) and v>0 for v in values):ax.set_yscale('log')
                        ax.set_title(c['case_id']+' / '+key,fontsize=8);ax.set_xlabel('Additional whole-child CPU cap (s)');ax.grid(alpha=.2)
                axes[i,0].legend(fontsize=7,ncol=4)
            fig.suptitle(f'{cohort}: every object, public minPhi paid-prefix selections; NA retained / page {page}')
            savefig(fig,out,f'{cohort}_paid_CPU_prefixes_{page}')
    for cohort in ('A','B','C'):
        fig,axes=plt.subplots(1,5,figsize=(20,5))
        for item in [r for r in arrows if r['cohort']==cohort]:
            for ax,key in zip(axes,('D','R','Phi','NRMSE','residual_over_delta')):
                a,b=item.get('raw_'+key),item.get('corrected_'+key)
                if a is not None and b is not None:
                    ax.annotate('',xy=(1,b),xytext=(0,a),arrowprops=dict(arrowstyle='->',alpha=.6,color=COLORS[item['method']]))
                    ax.scatter([0,1],[a,b],s=12,color=COLORS[item['method']])
                ax.set_xticks([0,1],['raw / direct initial','corrected']);ax.set_ylabel(key)
                ax.grid(alpha=.2)
                if key=='residual_over_delta':ax.axhline(1,color='gray',ls=':')
        from matplotlib.lines import Line2D
        fig.legend(handles=[Line2D([0],[0],color=COLORS[m],label=m) for m in METHODS],loc='lower center',ncol=4)
        fig.suptitle(cohort+' corrector: each frozen method/object; D/R/Phi/material/noise compared separately')
        savefig(fig,out,cohort+'_corrector_arrows')



def historical_figures(context, folder, cases):
    """Reuse saved decompositions only. Keep historical cohorts out of paired estimates."""
    import matplotlib.pyplot as plt
    for cohort, obj in context.items():
        if 'saved_arrays' not in obj: continue
        saved = obj['saved_arrays']
        if cohort == 'original_DR24':
            source_rows = [r for r in saved if r.get('stage')=='confirm24']
            data = []
            for r in source_rows:
                item = dict(case_id=r['case_id'], method='DR')
                for side,prefix in [('raw','raw_proposal_minPhi'),('corrected','primary')]:
                    rel=r.get(prefix+'_relative_data_residual')
                    for metric,key in [('Phi','Phi'),('R','R'),('material','relativeL2')]:
                        item[side+'_'+metric]=r.get(prefix+'_'+key)
                    item[side+'_D']=.5*rel**2 if rel is not None else None
                    case=next(c for c in cases if c['cohort']=='B' and c['case_id']==r['case_id'])
                    with np.load(ROOT/case['input_file'],allow_pickle=False) as z:ynorm=float(np.linalg.norm(z['Y']))
                    item[side+'_residual_over_delta']=rel*ynorm/case['raw_CN99_radius'] if rel is not None else None
                data.append(item)
            labels=('D','R','Phi','material','residual_over_delta')
        else:
            data=[]
            for r in saved['rows']:
                item=dict(case_id=r['case_id']+'/'+r['method'], method=r['method'])
                for side in ('raw','corrected'):
                    x=r[side]
                    for metric,key in [('D','data_term'),('R','prior'),('Phi','Phi'),('material','relative_material_l2')]:item[side+'_'+metric]=x.get(key)
                    noise=x.get('noise',{})
                    ratio=noise.get('residual_over_CN99_delta') if side=='raw' else x.get('residual_over_CN99_delta')
                    item[side+'_residual_over_delta']=ratio
                data.append(item)
            labels=('D','R','Phi','material','residual_over_delta')
        fig,axes=plt.subplots(1,5,figsize=(20,5))
        method_names=sorted({r['method'] for r in data})
        historical_colors={m:plt.get_cmap('tab10')(i%10) for i,m in enumerate(method_names)}
        historical_colors['DR']=COLORS['dr']
        for ax,key in zip(axes,labels):
            for r in data:
                a,b=r.get('raw_'+key),r.get('corrected_'+key)
                if a is not None and b is not None:
                    color=historical_colors[r['method']];is_dr=r['method']=='DR'
                    ax.annotate('',xy=(1,b),xytext=(0,a),arrowprops=dict(arrowstyle='->',color=color,alpha=.85 if is_dr else .35,lw=1.8 if is_dr else .8),zorder=4 if is_dr else 2)
                    ax.scatter([0,1],[a,b],s=18 if is_dr else 10,color=color,alpha=.85 if is_dr else .4,zorder=4 if is_dr else 2)
            ax.set_xticks([0,1],['raw','corrected']);ax.set_ylabel('relative material L2' if key=='material' else key);ax.grid(alpha=.2)
            if key=='residual_over_delta':ax.axhline(1,color='gray',ls=':')
        from matplotlib.lines import Line2D
        fig.legend(handles=[Line2D([0],[0],color=historical_colors[m],label=m,lw=2 if m=='DR' else 1) for m in method_names],loc='lower center',ncol=min(4,len(method_names)),fontsize=8)
        fig.suptitle(cohort+' / historical saved decompositions only / distinct objectives and budgets')
        savefig(fig,folder,cohort+'_context_decomposition')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seal',type=Path,required=True)
    parser.add_argument('--controller',type=Path,required=True,help='External full-exit controller job costs; never add nested CPU')
    parser.add_argument('--output-dir',type=Path,default=ROOT/'RESULTS/offline_v1',help='New offline scoring output directory; never overwrite an earlier score ledger')
    args=parser.parse_args();started=time.process_time()
    tasks=read(ROOT/'configs/TASKS.json')['tasks'];cases=read(ROOT/'inputs/manifest.json')['cases']
    seal,receipts=validate_seal(args.seal,tasks)
    jobs,controller=controller_jobs(args.controller)
    if args.controller and set(t['task_id'] for t in tasks)-set(jobs):
        raise RuntimeError('Controller missing task IDs; complete costs required')
    out=args.output_dir.resolve()
    if out.exists(): raise RuntimeError(str(out)+' already exists; preserve audit history')
    out.mkdir(parents=True);figdir=ROOT/'FIGURES';figdir.mkdir(exist_ok=True)
    dump(out/'TRUTH_ACCESS_BEGIN.json',dict(utc=datetime.now(timezone.utc).isoformat(),seal=str(args.seal),new_forward_calls=0))
    rows=[];prefixes=[];pictures={};costs=[];truth_sources=[]
    for case in cases:
        cohort,ident=case['cohort'],case['case_id']
        truth,support,amplitude,source=load_truth(case)
        portable=export_offline_truth(case,truth,support,amplitude,source)
        truth_sources.append(dict(cohort=cohort,case_id=ident,path=source,portable_offline_export=str(portable)))
        with np.load(ROOT/case['input_file'],allow_pickle=False) as z:
            m0=z['m_init'].copy().ravel(); ynorm=float(np.linalg.norm(z['Y']))
        parent_relative=case.get('parent_paid_state', {}).get('physical_relative_residual')
        parent_meta=dict(Phi=case.get('parent_primary_Phi'))
        if parent_relative is not None:
            parent_meta.update(relative_residual=parent_relative, residual_over_delta=parent_relative*ynorm/case['raw_CN99_radius'],
                               D=.5*parent_relative**2, noise_feasible=parent_relative*ynorm<=case['raw_CN99_radius'])
            if parent_meta['Phi'] is not None:parent_meta['R']=parent_meta['Phi']-parent_meta['D']
        pictures[(cohort,ident,'truth')]=truth;pictures[(cohort,ident,'initial')]=m0
        triggered=any(t['cohort']==cohort and t['case_id']==ident for t in tasks)
        init_metrics=metrics(m0,truth,support,amplitude,cohort)
        for method in METHODS:
            tid=f'{cohort}__{ident}__{method}';directory=ROOT/'RESULTS'/tid
            common=dict(cohort=cohort,case_id=ident,method=method,rescue_triggered=triggered,
                        parent_CPU_seconds=parent_cpu(case),parent_CPU_scope=parent_cpu_scope(case),scope='EXPOSED_POSTHOC_PAIRED')
            if tid not in receipts:
                if cohort!='B':raise RuntimeError('Unexpected unrun case '+tid)
                for view in VIEWS:
                    rows.append(dict(**common,view=view,**init_metrics,retained_original=True,
                                     **parent_meta))
                    pictures[(cohort,ident,method+'_'+view)]=m0
                for cap in (30,60,120):
                    prefixes.append(dict(**common,CPU_cap=cap,**init_metrics,**parent_meta,retained_original=True))
                continue
            receipt=receipts[tid]
            for view in VIEWS:
                file=directory/(view+'.npz');m=loadm(file) if file.exists() else None
                meta=read(directory/(view+'.json')) if file.exists() else {}
                row=dict(**common,view=view,**metrics(m,truth,support,amplitude,cohort),**meta,
                         task_status=receipt['status'],retained_original=False)
                if receipt['status']=='FAILED':
                    row.update(available=False,availability_reason='FAILED_TASK_PARTIAL_EVIDENCE_ONLY')
                    for metric in ('NRMSE','RMSE','IoU','relativeL2','descriptive_success'):row.pop(metric,None)
                rows.append(row)
                if m is not None:pictures[(cohort,ident,method+'_'+view)]=m
            history=read(directory/'HISTORY.json')
            with np.load(directory/'CANDIDATES.npz',allow_pickle=False) as z:materials=z['m'].copy()
            for prefix in read(directory/'PREFIXES.json'):
                cap=prefix['CPU_cap'];eligible=[r for r in history if r['eligible'] and r['CPU_seconds']<=cap]
                chosen=min(eligible,key=lambda r:r['Phi']) if eligible else None
                index=chosen['candidate'] if chosen else None
                if prefix['candidate']!=index:raise RuntimeError('Paid prefix selection mismatch '+tid)
                meta=chosen or {};m=materials[index].ravel() if index is not None else None
                score=metrics(m,truth,support,amplitude,cohort) if receipt['status']!='FAILED' else dict(available=False,availability_reason='FAILED_TASK_PARTIAL_EVIDENCE_ONLY')
                prefixes.append(dict(**common,CPU_cap=cap,**score,**meta))
            job=jobs.get(tid,{})
            external=job.get('external', {}).get('child_process_CPU_seconds')
            if external is None:
                external=next((job[k] for k in ('CPU_seconds','cpu_seconds','whole_child_CPU_seconds','process_cpu_seconds') if isinstance(job.get(k),(int,float))),None)
            costs.append(dict(task_id=tid,cohort=cohort,case_id=ident,method=method,
                              external_whole_child_CPU_seconds=external,internal_pre_exit_CPU_seconds=receipt.get('whole_child_CPU_seconds'),
                              external_controller_record=job or None,physics_cost=receipt.get('physics_cost'),
                              projection_cost=receipt.get('projection_cost'),G19_execution=receipt.get('G19_execution'),
                              accounting='Use external whole-child CPU only; internal/physics/projection are nested, never additive'))
    index={(r['cohort'],r['case_id'],r['method'],r['view']):r for r in rows}
    deployments=[]
    for case in cases:
        for method in METHODS:
            key=(case['cohort'],case['case_id'],method)
            secondary=index[key+('secondary',)];primary=index[key+('primary',)]
            # Availability is public sealed-candidate availability, never GT quality.
            use_secondary=secondary['available']
            chosen=secondary if use_secondary else primary
            deployed=dict(chosen,view='secondary_deploy',secondary_available=use_secondary,
                          fallback_to_primary=not use_secondary,
                          deployment_source_view='secondary' if use_secondary else 'primary')
            deployments.append(deployed)
            source_picture=pictures.get((key[0],key[1],method+'_'+deployed['deployment_source_view']))
            if source_picture is not None:pictures[(key[0],key[1],method+'_secondary_deploy')]=source_picture
    rows.extend(deployments)
    index={(r['cohort'],r['case_id'],r['method'],r['view']):r for r in rows}
    arrows=[];pairs=[];g19_pairs=[];summary=[]
    for case in cases:
        for method in METHODS:
            key=(case['cohort'],case['case_id'],method)
            raw=index[key+('raw',)];corrected=index[key+('corrected',)]
            item=dict(cohort=key[0],case_id=key[1],method=method)
            for metric in ('D','R','Phi','NRMSE','IoU','relativeL2','residual_over_delta'):
                item['raw_'+metric]=raw.get(metric);item['corrected_'+metric]=corrected.get(metric)
                a,b=raw.get(metric),corrected.get(metric);item['delta_'+metric]=b-a if a is not None and b is not None else None
            arrows.append(item)
        for control,comparisons,destination in [('direct',('g19','dr','ap'),pairs),('g19',('dr','ap','direct'),g19_pairs)]:
            for method in comparisons:
                for view in ('raw','corrected','primary','secondary','secondary_deploy'):
                    a=index[(case['cohort'],case['case_id'],method,view)]
                    b=index[(case['cohort'],case['case_id'],control,view)]
                    paired=a['available'] and b['available']
                    item=dict(cohort=case['cohort'],case_id=case['case_id'],method=method,control=control,
                              view=view,paired_available=paired,method_fallback_to_primary=a.get('fallback_to_primary'),
                              control_fallback_to_primary=b.get('fallback_to_primary'))
                    for metric in ('NRMSE','IoU','relativeL2','Phi','D','R','residual_over_delta'):
                        x,y=a.get(metric),b.get(metric)
                        item[metric+'_minus_'+control]=x-y if paired and x is not None and y is not None else None
                    destination.append(item)
    for cohort in ('A','B','C'):
        n=sum(c['cohort']==cohort for c in cases)
        for method in METHODS:
            for view in REPORT_VIEWS:
                rs=[r for r in rows if r['cohort']==cohort and r['method']==method and r['view']==view]
                item=dict(cohort=cohort,method=method,view=view,objects=n,available=sum(r['available'] for r in rs),
                          unavailable=sum(not r['available'] for r in rs),fallback_to_primary=sum(r.get('fallback_to_primary',False) for r in rs),
                          inherited_descriptive_success=sum(r.get('descriptive_success',False) for r in rs))
                for metric in ('NRMSE','IoU','relativeL2','Phi','D','R','residual_over_delta'):
                    values=[r[metric] for r in rs if r['available'] and r.get(metric) is not None]
                    item[metric+'_mean_available']=float(np.mean(values)) if values else None
                    item[metric+'_n']=len(values)
                summary.append(item)
    for name,values in [('SCORES',rows),('PAID_PREFIX_SCORES',prefixes),('CORRECTOR_DECOMPOSITION',arrows),('PAIRED_VS_DIRECT',pairs),('PAIRED_VS_G19',g19_pairs),('SECONDARY_DEPLOYMENT',deployments),('SUMMARY',summary),('CHILD_COST_LEDGER',costs),('TRUTH_SOURCES',truth_sources)]:table(out,name,values)
    projections=[]
    for tid,receipt in receipts.items():
        if receipt['method']=='g19':
            projections.append(dict(task_id=tid,cohort=receipt['task']['cohort'],case_id=receipt['case_id'],
                                    task_status=receipt['status'],
                                    projection_available=receipt.get('initial_current_projection') is not None,
                                    **(receipt.get('initial_current_projection') or {}),
                                    interpretation='Inherited LS projection diagnostics only; no new currents or forward calls'))
    table(out,'G19_INITIAL_PROJECTION_DEFECT',projections)
    parent_entries=[dict(cohort=c['cohort'],case_id=c['case_id'],parent_task_id=c['parent_task_id'],CPU_seconds=parent_cpu(c),CPU_scope=parent_cpu_scope(c),source_cost=c.get('parent_original_actual_cost')) for c in cases]
    table(out,'PARENT_COST_LEDGER',parent_entries)
    shared=next((c['parent_original_actual_cost'].get('shared_preprocess_summary') for c in cases if c['cohort']=='C'),None)
    dump(out/'COST_ACCOUNTING.json',dict(parent_acquisition_charged_once_per_case=parent_entries,shared_C_preprocess_charged_once=shared,
         warning='C shared preprocessing total already includes feature/inference/model load; report as one total, not sum with warm per-case costs. Parent costs are separate from additional rescue costs. Strategy totals add unique A/B parents once; no addition of nested internal/physics counters.',external_controller=controller))
    strategy_costs=[]
    for cohort in ('A','B','C'):
        parent_rows=[r for r in parent_entries if r['cohort']==cohort]
        parent_values=[r['CPU_seconds'] for r in parent_rows]
        # C shared preprocessing already includes all features/inference and imports.
        inherited=(shared.get('cpu_seconds') if shared else None) if cohort=='C' else (sum(parent_values) if all(x is not None for x in parent_values) else None)
        for method in METHODS:
            child_rows=[r for r in costs if r['cohort']==cohort and r['method']==method]
            values=[r['external_whole_child_CPU_seconds'] for r in child_rows]
            additional=sum(values) if all(x is not None for x in values) else None
            strategy_costs.append(dict(cohort=cohort,method=method,objects=len(parent_rows),rescue_children=len(child_rows),
                unique_parent_CPU_seconds=inherited,additional_external_child_CPU_seconds=additional,
                total_deployment_CPU_seconds=inherited+additional if inherited is not None and additional is not None else None,
                child_costs_available=sum(x is not None for x in values),controller_overhead='Reported separately in source controller receipt; not added once per strategy'))
    table(out,'STRATEGY_COST_TOTALS',strategy_costs)
    cost_quality=[]
    for cost in strategy_costs:
        cohort,method,n=cost['cohort'],cost['method'],cost['objects']
        for view in ('primary','secondary_deploy'):
            quality=next(r for r in summary if r['cohort']==cohort and r['method']==method and r['view']==view)
            item=dict(cohort=cohort,method=method,view=view,objects=n,available=quality['available'],
                      NRMSE_mean_available=quality['NRMSE_mean_available'],relativeL2_mean_available=quality['relativeL2_mean_available'],
                      selection='Frozen public minPhi or feasible minR with primary fallback; no GT-selected points')
            for source,target in [('unique_parent_CPU_seconds','mean_parent_CPU_seconds'),('additional_external_child_CPU_seconds','mean_added_rescue_CPU_seconds'),('total_deployment_CPU_seconds','mean_total_deployment_CPU_seconds')]:
                item[target]=cost[source]/n if cost[source] is not None else None
            cost_quality.append(item)
    table(out,'PERFORMANCE_VS_ACTUAL_MEAN_CPU',cost_quality)
    context={}
    for name,path in [('original_DR24',STUDY/'delegated/confirmation/RESULTS/offline_scoring_v1/DR_raw_primary_secondary_comparison.json'),('original_vector3D',STUDY/'delegated/vector_high/RAW_CORRECTOR_DECOMPOSITION.json')]:
        if path.is_file():context[name]=dict(source=str(path),saved_arrays=read(path),scope='Historical context only; no pooled cohort or equal-cost inference')
        else:context[name]=dict(source=str(path),available=False)
    dump(out/'HISTORICAL_CONTEXT_SEPARATE.json',context)
    figures(cases,rows,prefixes,pictures,arrows,figdir)
    historical_figures(context,figdir,cases)
    dump(out/'SCORING_RECEIPT.json',dict(status='COMPLETE_OFFLINE_ONLY',utc=datetime.now(timezone.utc).isoformat(),
         seal=str(args.seal),tasks=44,cases=31,rows=len(rows),new_forward=0,new_adjoint=0,new_training=0,
         CPU_seconds=time.process_time()-started,controller_costs_available=bool(jobs),
         scope='Exposed posthoc paired; inherited success labels descriptive; no new pass threshold or scientific verdict'))
    print(json.dumps(dict(status='COMPLETE_OFFLINE_ONLY',out=str(out),rows=len(rows))))


if __name__=='__main__':main()
