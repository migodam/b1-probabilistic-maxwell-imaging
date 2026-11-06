"""Offline scoring after complete prediction closure; no search/physics import."""
import argparse,json,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
HIST=ROOT.parent/'research_z_alg_r4_diversity'
def read(p):return json.loads(Path(p).read_text())
def dump(p,v):Path(p).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
def otsu(m):
    low,high=float(m.min()),float(m.max())
    if high<=low:return high
    h,e=np.histogram(m,bins=256,range=(low,high));c=(e[:-1]+e[1:])/2;w=np.cumsum(h);v=np.cumsum(h*c);d=w*(len(m)-w)
    s=np.divide((v[-1]*w-len(m)*v)**2,d,out=np.zeros(256),where=d>0)
    return float(c[np.argmax(s)])
def iou(p,q):
    u=np.count_nonzero(p|q);return float(np.count_nonzero(p&q)/u)if u else 1.
def metrics(m,truth,mask,amplitude,background):
    m=np.asarray(m).ravel();truth=truth.ravel();mask=mask.ravel().astype(bool);t=otsu(m);rmse=float(np.sqrt(np.mean((m-truth)**2)))
    return dict(RMSE=rmse,NRMSE=rmse/amplitude,MAE=float(np.mean(np.abs(m-truth))),IoU_Otsu=iou(m>t,mask),otsu_threshold_reconstruction_only=t,IoU_known_midpoint_OFFLINE=iou(m>background+.5*amplitude,mask))
def seal(include_tight=True):
    receipts=[]
    for exp in ('e2','e3'):
        directory=ROOT/'results'/exp;manifest=read(directory/'MANIFEST.json')
        planned=manifest['tasks']
        if exp=='e2'and include_tight:
            tight=directory/'STOPPING_MANIFEST.json'
            if tight.exists():planned+=read(tight)['tasks']
        for task in planned:
            p=directory/(task['task_id']+'.json');a=read(p)if p.exists()else None
            receipts.append(dict(task=task,path=str(p.relative_to(ROOT))if a else None,disposition='NOT_RUN'if a is None else a['status']))
    output=ROOT/'results/e3/ALL_SEARCH_PREDICTIONS_SEALED.json';dump(output,dict(status='ALL_PLANNED_PREDICTIONS_SEALED',sealed_unix=time.time(),receipts=receipts,no_truth_read_by_this_stage=True))
    return output

def score(seal_path):
    frozen=read(seal_path)
    if frozen['status']!='ALL_PLANNED_PREDICTIONS_SEALED':raise ValueError('Predictions must be closed before scoring')
    mapping={r['observation_id']:r for r in read(HIST/'truth_after_freeze/mapping_PRIVATE.json')['records']};rows=[]
    for item in frozen['receipts']:
        task=item['task'];case=mapping[task['observation_id']];row=dict(task,transport_status=item['disposition'],finite_endpoint=False)
        if not item['path']:rows.append(row);continue
        p=ROOT/item['path'];receipt=read(p);result=receipt.get('result',{})
        row.update(solver_status=result.get('status'),stop_reason=result.get('stop_reason'),wall_seconds=result.get('wall_seconds'),cpu_seconds=result.get('cpu_seconds'),costcounts=result.get('costcounts',{}),top1_online=result.get('top1'),reused_from=receipt.get('reused_from'),incremental_compute_seconds=receipt.get('incremental_compute_seconds',result.get('wall_seconds')))
        npz=p.parent/(task['task_id']+'.npz')
        if not npz.exists():rows.append(row);continue
        with np.load(npz,allow_pickle=False)as z,np.load(HIST/'truth_after_freeze'/(task['observation_id']+'.npz'),allow_pickle=False)as gt:
            truth,mask=gt['m'],gt['mask'];m=z['top1'];anchor=z['anchor']
            row.update(metrics(m,truth,mask,case['amplitude'],case['background']));row['anchor_metrics']=metrics(anchor,truth,mask,case['amplitude'],case['background']);row['finite_endpoint']=bool(np.all(np.isfinite(m)))
            before=z['proposal_before'].reshape(-1,576)if 'proposal_before'in z else np.empty((0,576));after=z['proposal_after'].reshape(-1,576)if 'proposal_after'in z else np.empty((0,576))
            row['candidate_before_metrics']=[metrics(x,truth,mask,case['amplitude'],case['background'])for x in before]
            row['candidate_after_metrics']=[metrics(x,truth,mask,case['amplitude'],case['background'])for x in after]
            candidates=np.vstack((anchor[None],m[None],before,after));errors=np.sqrt(np.mean((candidates-truth.reshape(1,-1))**2,axis=1))/case['amplitude']
            row['oracle_bestK_NRMSE_OFFLINE']=float(errors.min());row['oracle_candidate_count']=len(candidates);row['oracle_candidate_semantics']='anchor,online top1,finite directional seeds and per-seed correction best; not roots'
            h=z['proposal_filtered']if 'proposal_filtered'in z else np.empty((0,576));norms=np.linalg.norm(h,axis=1)
            row['attempted_seed_count']=len(h);row['nonzero_seed_count_gt_1e_12']=int(np.count_nonzero(norms>1e-12));row['seed_displacement_norms']=norms.tolist()
            row['top1_minus_anchor_norm']=float(np.linalg.norm(m-anchor));row['top1_delta_NRMSE']=row['NRMSE']-row['anchor_metrics']['NRMSE']
            if 'anchor_J_singular_values'in z:
                s=z['anchor_J_singular_values'];V=z['anchor_J_right_vectors'];eps=1e-4 if task['noise']=='clean'else 1e-3;strong=s*.03*np.sqrt(576)>eps
                e=truth.ravel()-m;ea=truth.ravel()-anchor;delta=m-anchor
                row['J_strong_dimension']=int(strong.sum());row['material_error_strong_squared']=float(np.linalg.norm(V[:,strong].T@e)**2);row['material_error_weak_squared']=float(np.linalg.norm(V[:,~strong].T@e)**2)
                row['anchor_error_strong_squared']=float(np.linalg.norm(V[:,strong].T@ea)**2);row['anchor_error_weak_squared']=float(np.linalg.norm(V[:,~strong].T@ea)**2)
                row['selected_change_strong_squared']=float(np.linalg.norm(V[:,strong].T@delta)**2);row['selected_change_weak_squared']=float(np.linalg.norm(V[:,~strong].T@delta)**2)
        residual=(result.get('top1')or{}).get('relative_residual');row['strict_data_qualification']=residual is not None and residual<=(1e-4 if task['noise']=='clean'else 1e-3)
        row['historical_r4_material_gate']=row['NRMSE']<=.15 and row['IoU_Otsu']>=.7
        rows.append(row)
    for exp in ('e2','e3'):
        subset=[r for r in rows if r['experiment']==exp];groups={}
        for r in subset:groups.setdefault((r['arm'],r['reg_scale'],r.get('stopping_profile','frozen'),r['noise']),[]).append(r)
        summary=[]
        for key,group in groups.items():
            valid=[r for r in group if r['finite_endpoint']];out=dict(zip(('arm','reg_scale','stopping_profile','noise'),key));out.update(n_planned=len(group),n_finite=len(valid))
            for metric in ('NRMSE','IoU_Otsu','oracle_bestK_NRMSE_OFFLINE','top1_delta_NRMSE','nonzero_seed_count_gt_1e_12','wall_seconds'):
                vals=[r[metric]for r in valid];out[metric+'_mean']=float(np.mean(vals))if vals else None;out[metric+'_median']=float(np.median(vals))if vals else None
            summary.append(out)
        dump(ROOT/'results'/exp/'SCORES.json',dict(seal=str(Path(seal_path).relative_to(ROOT)),historical_truth_already_disclosed=True,formal_scoring_after_prediction_seal=True,new_physics_LU=0,rows=subset,summary=summary))
    print(json.dumps(dict(rows=len(rows),finite=sum(r['finite_endpoint']for r in rows))))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('command',choices=['seal','score']);p.add_argument('--seal');a=p.parse_args()
    if a.command=='seal':print(seal())
    else:score(a.seal or ROOT/'results/e3/ALL_SEARCH_PREDICTIONS_SEALED.json')
