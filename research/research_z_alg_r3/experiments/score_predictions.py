"""Offline scorer, run only after prediction files are frozen.

The solver never imports this module. The output explicitly labels oracle
best@K separately from the online residual-selected top1.
"""
import argparse, datetime, json, time
from pathlib import Path
import numpy as np
from validate_predictions import frozen_files
ROOT=Path(__file__).resolve().parents[1]

def score(path,split,result=None,truth_root=None):
    result=json.loads(path.read_text()) if result is None else result
    case=result.get('case_id')
    if case is None:
        case=next(x.stem for x in (ROOT/'public_data'/split).glob('*.npz') if x.stem in path.name)
    truth_root=ROOT/'private_truth' if truth_root is None else Path(truth_root)
    with np.load(truth_root/split/(case+'.npz'),allow_pickle=False) as f:truth=f['m']
    candidates=result['candidates'];order=result['ranking']
    errors=[float(np.sqrt(np.mean((np.asarray(c['m'])-truth)**2))) for c in candidates]
    curve=[];best=float('inf')
    for k,i in enumerate(order,1):
        best=min(best,errors[i]);curve.append({'K':k,'best_RMS':best})
    mats=np.asarray([c['m'] for c in candidates]);pairs=[]
    for i in range(len(mats)):
        for j in range(i):pairs.append(float(np.sqrt(np.mean((mats[i]-mats[j])**2))))
    archive_time=[];top=None;best_error=float('inf')
    for i,c in sorted(enumerate(candidates),key=lambda item:item[1]['discovery_seconds']):
        if top is None or c['relative_residual'] < candidates[top]['relative_residual']:
            top=i
        best_error=min(best_error,errors[i])
        archive_time.append(dict(seconds=c['discovery_seconds'],candidate_index=i,
                                 online_archive_top1_index=top,online_archive_top1_RMS=errors[top],
                                 online_archive_top1_relative_residual=candidates[top]['relative_residual'],
                                 oracle_archive_best_RMS=best_error,
                                 full_complex_rhs=c.get('discovery_full_complex_rhs'),
                                 reduced_complex_rhs=c.get('discovery_reduced_complex_rhs')))
    anytime=[]
    for row in result.get('best_so_far',[]):
        if 'm' in row:
            anytime.append({k:v for k,v in row.items() if k!='m'} | {'top1_RMS':float(np.sqrt(np.mean((np.asarray(row['m'])-truth)**2)))})
    accepted=[i for i,c in enumerate(candidates) if c['relative_residual']<=1e-6]
    return {'file':str(path.relative_to(ROOT)),'case_id':case,'method':result['method'],'solver_revision':result.get('solver_revision'),'status':result['status'],'wall_seconds':result['wall_seconds'],'candidate_count':len(candidates),'top1_index':result['top1_index'],'top1_RMS':errors[result['top1_index']] if errors else None,'top1_relative_residual':result['top1']['relative_residual'] if errors else None,'oracle_best_RMS':min(errors) if errors else None,'oracle_best_index':int(np.argmin(errors)) if errors else None,'oracle_best_at_K':curve,'all_candidate_RMS':errors,'pair_RMS_min':min(pairs) if pairs else None,'pair_RMS_max':max(pairs) if pairs else None,'archive_available_by_time':archive_time,'online_anytime_top1':anytime,'anytime_trace_status':'FULL_STORED_RESIDUAL_IMPROVEMENTS' if anytime else 'ARCHIVE_AVAILABLE_BY_TIME_ONLY_NOT_FULL_ANYTIME','data_accepted_threshold':1e-6,'data_accepted_count':len(accepted),'data_accepted_candidate_indices':accepted,'cost':result['accounting'],'GT_values_disclosed':False,'metric_boundary':'best@K and archive-best-by-time are offline oracle coverage; top1 uses stored online all-data residual ranking; all returned roots are not necessarily data accepted; no scorer feedback to frozen predictions'}

def main():
    p=argparse.ArgumentParser();p.add_argument('--split',required=True,choices=['dev','holdout','confirm']);p.add_argument('--freeze-manifest',required=True,type=Path);p.add_argument('--output',type=Path,required=True);p.add_argument('--truth-root',type=Path,default=ROOT/'private_truth');p.add_argument('--disclosed-truth-replay',action='store_true');a=p.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    started=time.perf_counter()
    snapshots=list(frozen_files(a.freeze_manifest,a.disclosed_truth_replay));rows=[score(x,a.split,r,a.truth_root) for x,r in snapshots]
    payload={'scored_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'freeze_manifest':str(a.freeze_manifest),'prediction_files':[str(x.relative_to(ROOT)) for x,_ in snapshots],'rows':rows,'new_physics_rhs':0,'scoring_only':True,'truth_read_after_prediction_freeze':True,'truth_root':str(a.truth_root),'disclosed_truth_replay_not_new_blind':a.disclosed_truth_replay,'snapshot_identity_mode':'byte length only, allowing clone mtimes; disclosed replay cannot establish blindness' if a.disclosed_truth_replay else 'original byte length and mtime_ns','scoring_wall_seconds':time.perf_counter()-started}
    a.output.write_text(json.dumps(payload,indent=2))
    for r in rows:print(json.dumps({k:r[k] for k in ['file','top1_RMS','oracle_best_RMS','top1_relative_residual','candidate_count','wall_seconds']}))

if __name__=='__main__':main()
