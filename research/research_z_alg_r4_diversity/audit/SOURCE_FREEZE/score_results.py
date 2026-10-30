"""Independent offline scoring. No solver import or new physical solve."""
import argparse, itertools, json, pathlib, time
import numpy as np
ROOT=pathlib.Path(__file__).resolve().parent
METHODS=('native_trf','reciprocal_trf','data_lbfgs','globaldual')

def read(path):return json.loads(pathlib.Path(path).read_text(encoding='utf-8-sig'))
def dump(path,value):
    path=pathlib.Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')

def seal(manifest_path,results,output):
    manifest=read(manifest_path);tasks=manifest['tasks'];receipts=[]
    # Parent calls only after closing the experiment. Missing planned tasks stay NOT_RUN.
    for task in tasks:
        path=pathlib.Path(results)/(task['task_id']+'.json')
        receipt=read(path) if path.exists() else None
        disposition='NOT_RUN' if receipt is None else 'complete' if receipt.get('status')=='OK' else 'failed'
        receipts.append({'task_id':task['task_id'],'task':task,'path':str(path.resolve()) if receipt else None,
                         'disposition':disposition})
    dump(output,{'status':'PLANNED_RESULTS_CLOSED','sealed_unix':time.time(),
                 'planned':len(tasks),'receipts':receipts,'manifest':str(pathlib.Path(manifest_path).resolve()),
                 'no_truth_read':True,'no_new_physics':True})

def otsu_threshold(m):
    """256 bins on reconstruction's own min/max, never inspect truth."""
    low,high=float(m.min()),float(m.max())
    if high<=low:return high
    hist,edges=np.histogram(m,bins=256,range=(low,high));centers=(edges[:-1]+edges[1:])/2
    weight=np.cumsum(hist);moment=np.cumsum(hist*centers);total=moment[-1]
    denom=weight*(len(m)-weight)
    score=np.divide((total*weight-len(m)*moment)**2,denom,out=np.zeros(256),where=denom>0)
    return float(centers[int(np.argmax(score))])

def iou(pred,mask):
    union=np.count_nonzero(pred|mask)
    return float(np.count_nonzero(pred&mask)/union) if union else 1.

def metrics(m,truth,mask,amplitude,background):
    m=np.asarray(m,float).ravel();truth=truth.ravel();mask=mask.ravel().astype(bool)
    rmse=float(np.sqrt(np.mean((m-truth)**2)));threshold=otsu_threshold(m)
    return {'RMSE':rmse,'NRMSE':rmse/amplitude,'MAE':float(np.mean(np.abs(m-truth))),
            'IoU_Otsu':iou(m>threshold,mask),'otsu_threshold_reconstruction_only':threshold,
            'IoU_known_midpoint_OFFLINE':iou(m>background+.5*amplitude,mask)}

def summarize(rows,group_fields):
    groups={}
    for row in rows:groups.setdefault(tuple(row[k] for k in group_fields),[]).append(row)
    out=[]
    for key,group in groups.items():
        valid=[r for r in group if r['finite_endpoint']]
        not_run=sum(r['transport_status']=='NOT_RUN' for r in group)
        item={**dict(zip(group_fields,key)),'n_planned':len(group),'n_finite':len(valid),
              'n_failed':len(group)-len(valid)-not_run,'n_NOT_RUN':not_run}
        for gate in ('data_gate','material_gate','joint_gate'):
            item[gate+'_counts']={state:sum(r[gate]==state for r in group) for state in ('PASS','FAIL','NOT_RUN')}
        for name in ('RMSE','NRMSE','MAE','IoU_Otsu','IoU_known_midpoint_OFFLINE'):
            item[name+'_mean']=float(np.mean([r[name] for r in valid])) if valid else None
            item[name+'_median']=float(np.median([r[name] for r in valid])) if valid else None
        for name in ('LU','actual_complex_RHS','wall_seconds','cpu_seconds'):
            values=[r[name] for r in group if r.get(name) is not None]
            item[name+'_sum']=float(np.sum(values)) if values else None
            item[name+'_mean']=float(np.mean(values)) if values else None
            item[name+'_known_count']=len(values)
        out.append(item)
    return out

def pairwise(rows):
    out=[]
    for backend in ('cpu','cuda'):
        available={ (r['observation_id'],r['method']):r for r in rows if r['backend']==backend}
        ids=sorted({k[0] for k in available})
        for left,right in itertools.combinations(METHODS,2):
            wins=ties=losses=missing=0;deltas=[];pairs=[]
            for ident in ids:
                a,b=available.get((ident,left)),available.get((ident,right))
                if not a or not b or not a['finite_endpoint'] or not b['finite_endpoint']:missing+=1;continue
                delta=a['NRMSE']-b['NRMSE'];tol=max(.02,.10*min(a['NRMSE'],b['NRMSE']))
                deltas.append(delta)
                pairs.append({'observation_id':ident,'delta_NRMSE_left_minus_right':delta,'practical_tie_tolerance':tol})
                if abs(delta)<=tol:ties+=1
                elif delta<0:wins+=1
                else:losses+=1
            out.append({'backend':backend,'left':left,'right':right,'left_wins':wins,'ties':ties,'right_wins':losses,
                'missing_pairs':missing,'paired_NRMSE_delta_mean':float(np.mean(deltas)) if deltas else None,
                'pairs':pairs,'tie_rule':'abs(delta)<=max(.02,.10*min(NRMSE_left,NRMSE_right))'})
    return out

def score(seal_path,truth_dir,output):
    frozen=read(seal_path)
    if frozen['status'] not in ('PLANNED_RESULTS_CLOSED','ALL_PLANNED_PREDICTIONS_SEALED'):
        raise RuntimeError('Parent must first close planned result dispositions')
    # Only now disclose local truth to the independent scoring process.
    truth_dir=pathlib.Path(truth_dir);mapping=read(truth_dir/'mapping_PRIVATE.json')['records']
    lookup={r['observation_id']:r for r in mapping};rows=[]
    gates=read(ROOT/'configs/CONTRACT.json')['metrics']['joint_success']
    for file in frozen['receipts']:
        receipt=read(file['path']) if file.get('path') else {'task':file['task'],'status':'NOT_RUN'}
        task=receipt['task'];result=receipt.get('result') or {}
        ident=pathlib.PurePosixPath(task['public_path'].replace('\\','/')).stem;case=lookup[ident]
        row={**case,'backend':task['backend'],'method':task['method'],'task_id':task['task_id'],
             'transport_status':receipt['status'],'solver_status':result.get('status'),
             'stop_reason':result.get('stop_reason'),'finite_endpoint':False,
             'wall_seconds':result.get('wall_seconds',receipt.get('wall_observed_s')),
             'cpu_seconds':result.get('cpu_seconds'),'LU':result.get('costcounts',{}).get('LU'),
             'actual_complex_RHS':result.get('costcounts',{}).get('actual_complex_RHS'),
             'costs_status':receipt.get('costs_status','COMPLETE' if 'costcounts' in result else 'UNKNOWN')}
        row['relative_data_residual']=(result.get('top1') or {}).get('relative_residual')
        row['oracle_bestK_RMSE_OFFLINE']=None;row['oracle_bestK_NRMSE_OFFLINE']=None
        row['oracle_candidate_semantics']='finite correction-phase optimizer candidates, not independent branches'
        row['oracle_finite_candidate_count']=0
        material=result.get('m')
        if material is not None:
            material=np.asarray(material,float)
            if material.size==576 and np.all(np.isfinite(material)):
                with np.load(truth_dir/(ident+'.npz'),allow_pickle=False) as t:
                    row.update(metrics(material,t['m'],t['mask'],case['amplitude'],case['background']))
                    candidates=[np.asarray(c['m'],float).ravel() for c in result.get('candidates',[]) if c.get('m') is not None]
                    candidates=[c for c in candidates if c.size==576 and np.all(np.isfinite(c))]
                    row['oracle_finite_candidate_count']=len(candidates)
                    if candidates:
                        best=min(float(np.sqrt(np.mean((c-t['m'].ravel())**2))) for c in candidates)
                        row['oracle_bestK_RMSE_OFFLINE']=best;row['oracle_bestK_NRMSE_OFFLINE']=best/case['amplitude']
                row['finite_endpoint']=True
        not_run=row['transport_status']=='NOT_RUN'
        residual_limit=gates['relative_data_residual_clean_max' if case['noise']=='clean' else 'relative_data_residual_noisy_max']
        data_ok=row['relative_data_residual'] is not None and row['relative_data_residual']<=residual_limit
        material_ok=row['finite_endpoint'] and row['NRMSE']<=gates['NRMSE_max'] and row['IoU_Otsu']>=gates['Otsu_IoU_min']
        for name,ok in (('data_gate',data_ok),('material_gate',material_ok),('joint_gate',data_ok and material_ok)):
            row[name]='NOT_RUN' if not_run else 'PASS' if ok else 'FAIL'
        rows.append(row)
    paired=[]
    lookup_rows={(r['backend'],r['shape'],r['amplitude'],r['noise'],r['method']):r for r in rows}
    for row in rows:
        if row['noise']!='clean':continue
        noisy=lookup_rows.get((row['backend'],row['shape'],row['amplitude'],'noisy',row['method']))
        if noisy and row['finite_endpoint'] and noisy['finite_endpoint']:
            paired.append({k:row[k] for k in ('backend','shape','amplitude','method')}|{
                'clean_NRMSE':row['NRMSE'],'noisy_NRMSE':noisy['NRMSE'],'noise_NRMSE_delta':noisy['NRMSE']-row['NRMSE']})
    result={'scope':'Finite discrete 24x24 scalar-TM full-VIE comparison; CPU/CUDA summarized separately',
        'seal':str(seal_path),'truth_disclosed_only_after_seal':True,'new_physics_LU':0,'rows':rows,
        'summary':summarize(rows,('backend','method')),
        'by_amplitude_noise':summarize(rows,('backend','amplitude','noise','method')),
        'pairwise':pairwise(rows),'noise_pairs':paired}
    dump(output,result);print(json.dumps({'rows':len(rows),'finite':sum(r['finite_endpoint'] for r in rows),'output':str(output)}))

if __name__=='__main__':
    ap=argparse.ArgumentParser();sub=ap.add_subparsers(dest='command',required=True)
    p=sub.add_parser('seal');p.add_argument('--manifest',required=True);p.add_argument('--results',required=True);p.add_argument('--output',required=True)
    p=sub.add_parser('score');p.add_argument('--seal',required=True);p.add_argument('--truth',default=ROOT/'private_truth');p.add_argument('--output',required=True)
    a=ap.parse_args()
    if a.command=='seal':seal(a.manifest,a.results,a.output)
    else:score(a.seal,a.truth,a.output)
