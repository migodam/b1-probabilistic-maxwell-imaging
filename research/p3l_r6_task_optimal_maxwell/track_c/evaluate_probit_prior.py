"""Seal all prior-Probit64/128 actual/q Det predictions before CDEV scoring.

Additional strong-baseline DEVELOPMENT only; no new locked-D arm or statistics.
No matching-head fit, physical call, D truth, or old test truth is permitted here.
"""
import os
for _name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_name]='1'
import argparse,hashlib,io,json,sys,time
from pathlib import Path
import numpy as np
TRACK=Path(__file__).resolve().parent;OUT=TRACK.parent;sys.path.insert(0,str(OUT))
from track_c import evaluate_dev as base
from track_c.frozen_heads import load_action,predict_contexts,SEEDS
from track_c.kkt_closure.closure_action import Action as ClosureAction

CONFIGS=('prior_x0_probit64','prior_x0_probit128')
MODES=('actual','q')


def artifact(path,**fields):
    content=path.read_bytes()
    return dict(path=path.name,bytes=len(content),blake2b=hashlib.blake2b(content).hexdigest(),**fields)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=TRACK/'probit_prior_baseline/evaluation_CDEV')
    parser.add_argument('--predict-only',action='store_true')
    parser.add_argument('--score-only',action='store_true')
    args=parser.parse_args()
    if args.predict_only and args.score_only:raise ValueError('at most one phase flag')
    if not args.output.resolve().is_relative_to(TRACK.resolve()):raise ValueError('Track C outputs only')
    freeze=json.loads((TRACK/'probit_prior_baseline/DEVELOPMENT_SOURCE_FREEZE.json').read_text())
    for row in freeze['files']:
        content=Path(row['path']).read_bytes()
        if len(content)!=row['bytes'] or hashlib.blake2b(content).hexdigest()!=row['blake2b']:raise RuntimeError('development baseline source/protocol changed')
    parents=base.expected_parents();marker=json.loads((TRACK/'probit_prior_baseline/ANCHORS_COMPLETE.json').read_text())
    if marker.get('status')!='COMPLETE' or marker.get('n')!=64 or set(marker.get('configs',[]))!=set(CONFIGS):
        raise RuntimeError('all64 x64/128 physical baseline paths must complete before predictions')
    reports={config:[json.loads((TRACK/'contexts/CDEV'/p/config/'REPORT.json').read_text()) for p in parents] for config in CONFIGS}
    contexts={config:base.load_contexts(parents,config) for config in CONFIGS}
    args.output.mkdir(parents=True,exist_ok=True);seal_path=args.output/'ALL_PREDICTIONS_SEALED.json'
    if seal_path.exists():seal=json.loads(seal_path.read_text())
    else:
        if args.score_only:raise RuntimeError('existing complete seal required')
        if list(args.output.glob('*.npz')):raise RuntimeError('partial prediction archive; explicit audit before replay')
        records=[];loads=[]
        for config in CONFIGS:
            for mode in MODES:
                for seed in SEEDS:
                    tick,cpu=time.perf_counter(),time.process_time()
                    action=load_action('deterministic_matched',seed) if mode=='actual' else ClosureAction('deterministic_matched',seed,'q',1.)
                    loads.append(dict(config=config,mode=mode,seed=seed,load_wall_seconds=time.perf_counter()-tick,load_CPU_seconds=time.process_time()-cpu))
                    prediction=[];wall=[];CPU=[]
                    for parent,context in zip(parents,contexts[config]):
                        tick,cpu=time.perf_counter(),time.process_time();value=predict_contexts([context],action,[parent],1)[0]
                        wall.append(time.perf_counter()-tick);CPU.append(time.process_time()-cpu)
                        if value.shape!=(576,) or not np.isfinite(value).all():raise ValueError('finite full-grid prediction required; no silent fallback')
                        prediction.append(value)
                    path=args.output/f'{config}__{mode}_s{seed}.npz'
                    np.savez_compressed(path,parent=np.array(parents),prediction=np.array(prediction),inference_wall_seconds=np.array(wall),inference_CPU_seconds=np.array(CPU))
                    records.append(artifact(path,config=config,mode=mode,seed=seed));del action
        seal=dict(status='ALL_PRIOR64_128_ACTUAL_Q_DEVELOPMENT_PREDICTIONS_SEALED',parents=parents,configs=list(CONFIGS),modes=list(MODES),
            seeds=list(SEEDS),truth_read=False,CDEV_truth_previously_exposed=True,new_RHS=0,new_training_updates=0,
            model_load_costs=loads,artifacts=records,not_locked_D_arm=True)
        base.write(seal_path,seal)
    snapshots={}
    if seal.get('truth_read') is not False or seal['parents']!=parents:raise ValueError('matching global seal required')
    for row in seal['artifacts']:
        path=args.output/row['path'];content=path.read_bytes()
        if Path(row['path']).name!=row['path'] or len(content)!=row['bytes'] or hashlib.blake2b(content).hexdigest()!=row['blake2b']:raise ValueError('sealed content changed/escaping path')
        with np.load(io.BytesIO(content),allow_pickle=False) as z:
            if z['parent'].astype(str).tolist()!=parents or z['prediction'].shape!=(64,576):raise ValueError('prediction schema/parents')
            snapshots[(row['config'],row['mode'],row['seed'])]={key:z[key].copy() for key in ('prediction','inference_wall_seconds','inference_CPU_seconds')}
    expected={(c,m,s) for c in CONFIGS for m in MODES for s in SEEDS}
    if set(snapshots)!=expected:raise ValueError('complete12-file prediction family required before scoring')
    if args.predict_only:print(json.dumps(dict(status=seal['status'],truth_read=False)),flush=True);return
    old_dir=TRACK/'evaluation_CDEV';old_seal=json.loads((old_dir/'ALL_PREDICTIONS_SEALED.json').read_text())
    old=base.validate_seal(old_dir,old_seal,parents)
    qdir=TRACK/'kkt_closure';qseal=json.loads((qdir/'ALL_PREDICTIONS_SEALED.json').read_text());old_q={}
    for row in qseal['artifacts']:
        if row['config']!='native24_probit64' or row['mode']!='q':continue
        content=(qdir/row['path']).read_bytes()
        if len(content)!=row['bytes'] or hashlib.blake2b(content).hexdigest()!=row['blake2b']:raise ValueError('previous q24/P64 sealed bytes changed')
        with np.load(io.BytesIO(content),allow_pickle=False) as z:
            if z['parent'].astype(str).tolist()!=parents:raise ValueError('previous q parent identity')
            old_q[row['seed']]=z['prediction'].copy()
    if set(old_q)!=set(SEEDS):raise ValueError('complete previous q24/P64 reference')
    # First actual CDEV material-truth open in this evaluation. Never D or oldtest.
    truth=[]
    for parent in parents:
        with np.load(OUT/'raw/offline_truth'/f'{parent}.npz',allow_pickle=False) as z:truth.append(z['chi'].copy())
    truth=np.array(truth)
    labels=['OriginalFullDet96P64','ExistingQDet24P64']+[f'{config}__{mode}' for config in CONFIGS for mode in MODES]
    mse=np.zeros((64,len(labels),3));mae=np.zeros_like(mse);wholewall=np.zeros_like(mse);ensemble=np.zeros((64,len(labels)))
    old_reports=[json.loads((TRACK/'contexts/CDEV'/p/'native96_probit64/REPORT.json').read_text()) for p in parents]
    qreports=[json.loads((TRACK/'contexts/CDEV'/p/'native24_probit64/REPORT.json').read_text()) for p in parents]
    for j,label in enumerate(labels):
        predictions=[]
        for k,seed in enumerate(SEEDS):
            if j==0:
                row=old[('native96_probit64','deterministic_matched',seed)];pred=row['prediction'];wall=row['inference_wall_seconds'];anchor=[r['end_to_end_anchor_wall_seconds'] for r in old_reports]
            elif j==1:
                pred=old_q[seed]
                with np.load(qdir/f'native24_probit64__q_s{seed}.npz',allow_pickle=False) as z:wall=z['inference_wall_seconds'].copy()
                anchor=[r['end_to_end_anchor_wall_seconds'] for r in qreports]
            else:
                config,mode=label.split('__');row=snapshots[config,mode,seed];pred=row['prediction'];wall=row['inference_wall_seconds'];anchor=[r['end_to_end_anchor_wall_seconds'] for r in reports[config]]
            error=abs(pred-truth);mse[:,j,k]=(error**2).mean(1);mae[:,j,k]=error.mean(1);wholewall[:,j,k]=np.array(anchor)+wall;predictions.append(pred)
        ensemble[:,j]=(abs(np.mean(predictions,0)-truth)**2).mean(1)
    means=mse.mean(2);walls=wholewall.mean(2);draw=base.bootstrap(np.stack((means,walls),-1));arms=[]
    for j,label in enumerate(labels):
        ratio=draw[:,j,0]/draw[:,0,0];saving=1-draw[:,j,1]/draw[:,0,1]
        upper=float(np.quantile(ratio,1-base.ALPHA));lower=float(np.quantile(saving,base.ALPHA))
        difference=means[:,j]-means[:,1];dd=draw[:,j,0]-draw[:,1,0]
        arms.append(dict(arm=label,MSE_meanseed=float(means[:,j].mean()),MAE_meanseed=float(mae[:,j].mean()),ensemble_MSE_descriptive=float(ensemble[:,j].mean()),
            ratio_vs_OriginalFullDet=float(means[:,j].mean()/means[:,0].mean()),upper98_333_ratio_vs_OriginalFullDet=upper,
            NI1pct=upper<=1.01,mean_whole_wall_seconds=float(walls[:,j].mean()),wall_saving_lower98_333=lower,wall20pct=lower>=.20,
            development_joint_gate=bool(j>1 and upper<=1.01 and lower>=.20),
            MSE_minus_existing_q24P64=float(difference.mean()),MSE_difference_CI98_333_vs_existing_q24P64=np.quantile(dd,[base.ALPHA/2,1-base.ALPHA/2]).tolist()))
    np.savez_compressed(args.output/'PARENT_LOSSES.npz',parent=np.array(parents),arms=np.array(labels),seed=np.array(SEEDS),MSE_seedwise=mse,MAE_seedwise=mae,
        MSE_meanseed=means,MAE_meanseed=mae.mean(2),MSE_ensemble=ensemble,whole_wall_meanseed=walls)
    result=dict(status='ADDITIONAL_PRIOR_PROBIT_STRONG_BASELINE_DEVELOPMENT_COMPLETE',parents=64,arms=arms,bootstrap=20000,
        confidence=1-base.ALPHA,development_only=True,CDEV_previously_exposed=True,D_GT_read=False,oldtest_GT_read=False,
        Native_optimizer_ran=False,source_freeze=str(TRACK/'probit_prior_baseline/DEVELOPMENT_SOURCE_FREEZE.json'),
        whole_wall_scope='actual standalone prior-point+Probit path report plus measured current-load CPU1 frozen head chain; final D stats unchanged',
        actual_physics_cost_by_config={c:dict(RHS=sum(r['path_resource_cost']['RHS'] for r in rows),wall_seconds=sum(r['end_to_end_anchor_wall_seconds'] for r in rows)) for c,rows in reports.items()},
        stronger_baseline_claim_limit='If direct prior Probit competes with or dominates the proposed allocation, final TAP/novelty claims must disclose it; no unannounced locked-D arm/statistic change',
        new_prediction_RHS=0,new_head_training_updates=0,seal=str(seal_path))
    base.write(args.output/'RESULTS.json',result);print(json.dumps(result),flush=True)


if __name__=='__main__':main()
