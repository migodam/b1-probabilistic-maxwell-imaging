"""CDEV selection gates for matched512 heads versus newFull512 AND OriginalR5.

All requested frozen head predictions are globally BLAKE2b sealed before this
evaluation's truth opens. CDEV has already selected checkpoints, so these gates
are development diagnostics, never independent confirmation. Native cold phase
cost remains charged; actual standalone D timing is required for final wall claims.
"""
import os
for _name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_name]='1'
import argparse,hashlib,io,json,sys,time
from pathlib import Path
import numpy as np
TRACK=Path(__file__).resolve().parent;OUT=TRACK.parent;sys.path.insert(0,str(OUT))
from track_c.matched_heads import MatchedAction,context_batch
from track_c.frozen_heads import predict_contexts,SEEDS
from track_c import evaluate_dev as base


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--systems',nargs='+',choices=['full','native','low'],default=['full','native'])
    parser.add_argument('--heads',type=Path,default=TRACK/'matched_heads_TRAIN512')
    parser.add_argument('--output',type=Path,default=TRACK/'evaluation_matched512_full_native')
    args=parser.parse_args()
    if 'full' not in args.systems or len(set(args.systems))!=len(args.systems):raise ValueError('paired newFull512 required')
    if not args.output.resolve().is_relative_to(TRACK.resolve()):raise ValueError('Track C output only')
    systems=['full']+[arm for arm in args.systems if arm!='full']
    parents=base.expected_parents();args.output.mkdir(parents=True,exist_ok=True)
    contexts={};anchor_wall={};anchor_rhs={}
    for arm in systems:
        config={'full':'native96_probit64','native':'native96_only','low':'native48_probit32'}[arm]
        root=TRACK/'native96_only_contexts/CDEV' if arm=='native' else TRACK/'contexts/CDEV'
        contexts[arm]=[];anchor_wall[arm]=[];anchor_rhs[arm]=[]
        for parent in parents:
            folder=root/parent/config
            with np.load(folder/'CONTEXT.npz',allow_pickle=False) as z:
                contexts[arm].append({key:z[key].copy() for key in (*base.ONLINE_KEYS,'parent')})
            report=json.loads((folder/'REPORT.json').read_text())
            if arm=='native':
                anchor_wall[arm].append(report['cold_native_wall_seconds']+report['build_wall_seconds'])
                anchor_rhs[arm].append(report['cold_native_RHS'])
            else:
                anchor_wall[arm].append(report['end_to_end_anchor_wall_seconds'])
                anchor_rhs[arm].append(report['path_resource_cost']['RHS'])
    seal_path=args.output/'ALL_MATCHED_PREDICTIONS_SEALED.json'
    if seal_path.exists():seal=json.loads(seal_path.read_text())
    else:
        if list(args.output.glob('*_s*.npz')):raise RuntimeError('partial matched inference requires audit')
        records=[]
        for arm in systems:
            for seed in SEEDS:
                checkpoint=args.heads/f'{arm}_s{seed}/BEST.pt'
                if not checkpoint.exists() or not checkpoint.with_name('RESULT.json').exists():raise RuntimeError('all requested head fits complete before inference')
                action=MatchedAction(checkpoint);predictions=[];wall=[];cpu=[]
                for parent,context in zip(parents,contexts[arm]):
                    tick,cpu0=time.perf_counter(),time.process_time()
                    predictions.append(predict_contexts([context],action,[parent],1)[0])
                    wall.append(time.perf_counter()-tick);cpu.append(time.process_time()-cpu0)
                path=args.output/f'{arm}_s{seed}.npz'
                np.savez_compressed(path,parent=np.array(parents),prediction=np.array(predictions),inference_wall=np.array(wall),inference_CPU=np.array(cpu))
                content=path.read_bytes();records.append(dict(path=path.name,arm=arm,seed=seed,bytes=len(content),
                    blake2b=hashlib.blake2b(content).hexdigest(),checkpoint=str(checkpoint),
                    checkpoint_blake2b=hashlib.blake2b(checkpoint.read_bytes()).hexdigest()))
        seal=dict(status='ALL_MATCHED_CDEV_PREDICTIONS_SEALED',parents=parents,systems=systems,truth_read=False,
                  CDEV_previously_used_for_checkpoint_selection=True,new_RHS=0,artifacts=records)
        base.write(seal_path,seal)
    if seal['parents']!=parents or seal['systems']!=systems or seal['truth_read'] is not False:raise ValueError('matching sealed family')
    snapshots={}
    for row in seal['artifacts']:
        content=(args.output/row['path']).read_bytes()
        if len(content)!=row['bytes'] or hashlib.blake2b(content).hexdigest()!=row['blake2b']:raise ValueError('matched prediction content changed')
        with np.load(io.BytesIO(content),allow_pickle=False) as z:
            if z['parent'].astype(str).tolist()!=parents or z['prediction'].shape!=(64,576):raise ValueError('matched prediction schema')
            snapshots[(row['arm'],row['seed'])]={key:z[key].copy() for key in ('prediction','inference_wall','inference_CPU')}
    if len(snapshots)!=len(systems)*3:raise ValueError('complete three-seed family required')
    old_destination=TRACK/'evaluation_CDEV';old_seal=json.loads((old_destination/'ALL_PREDICTIONS_SEALED.json').read_text())
    original=base.validate_seal(old_destination,old_seal,parents)
    # First actual truth open in this frozen prediction evaluation.
    truth=[]
    for parent in parents:
        with np.load(OUT/'raw/offline_truth'/f'{parent}.npz',allow_pickle=False) as z:truth.append(z['chi'].copy())
    truth=np.array(truth)
    shape=(64,len(systems)+1,3);mse=np.zeros(shape);mae=np.zeros(shape);wall=np.zeros(shape)
    ensemble_mse=np.zeros(shape[:2]);ensemble_mae=np.zeros(shape[:2])
    for j,arm in enumerate([*systems,'original_R5_full1024']):
        predictions=[]
        for k,seed in enumerate(SEEDS):
            row=(snapshots[(arm,seed)] if arm in systems else original[('native96_probit64','distilled_w128_n32768',seed)])
            prediction=row['prediction'];error=abs(prediction-truth);predictions.append(prediction)
            mse[:,j,k]=(error**2).mean(1);mae[:,j,k]=error.mean(1)
            timing=row['inference_wall'] if arm in systems else row['inference_wall_seconds']
            wall[:,j,k]=np.array(anchor_wall['full' if arm=='original_R5_full1024' else arm])+timing
        error=abs(np.mean(predictions,0)-truth);ensemble_mse[:,j]=(error**2).mean(1);ensemble_mae[:,j]=error.mean(1)
    meanseed=mse.mean(2);wholewall=wall.mean(2);draws=base.bootstrap(np.stack((meanseed,wholewall),-1))
    results=[]
    for j,arm in enumerate(systems):
        ratio_new=draws[:,j,0]/draws[:,0,0];ratio_original=draws[:,j,0]/draws[:,-1,0]
        saving=1-draws[:,j,1]/draws[:,-1,1]
        upper_new=float(np.quantile(ratio_new,1-base.ALPHA));upper_old=float(np.quantile(ratio_original,1-base.ALPHA));lower_wall=float(np.quantile(saving,base.ALPHA))
        results.append(dict(system=arm,MSE_meanseed=float(meanseed[:,j].mean()),MAE_meanseed=float(mae[:,j].mean()),
            ensemble_MSE_descriptive=float(ensemble_mse[:,j].mean()),ensemble_MAE_descriptive=float(ensemble_mae[:,j].mean()),
            risk_ratio_vs_newFull512=float(meanseed[:,j].mean()/meanseed[:,0].mean()),risk_ratio_vs_originalFull1024=float(meanseed[:,j].mean()/meanseed[:,-1].mean()),
            upper98_333_vs_newFull512=upper_new,upper98_333_vs_originalFull1024=upper_old,
            same_N_noninferiority_gate=upper_new<=1.01,original_R5_noninferiority_gate=upper_old<=1.01,
            composed_cold_whole_wall_mean=float(wholewall[:,j].mean()),wall_saving_lower98_333_vs_originalFull1024=lower_wall,
            composed_wall_reduction20_gate=lower_wall>=.20,
            development_risk_candidate=bool(arm!='full' and upper_new<=1.01 and upper_old<=1.01),
            final_performance_claim='HOLD_UNTIL_ACTUAL_STANDALONE_NATIVE_D_TIMING_WITH_MATCHED_LOGGING'))
    np.savez_compressed(args.output/'PARENT_MATCHED_LOSSES.npz',parent=np.array(parents),systems=np.array([*systems,'original_R5_full1024']),
        seeds=np.array(SEEDS),MSE_seedwise=mse,MAE_seedwise=mae,MSE_meanseed=meanseed,MAE_meanseed=mae.mean(2),MSE_ensemble=ensemble_mse,MAE_ensemble=ensemble_mae,
        composed_whole_wall=wholewall)
    result=dict(status='MATCHED_CDEV_DEVELOPMENT_EVALUATED',systems=results,original_R5_full1024_MSE=float(meanseed[:,-1].mean()),
        matched_newFull512_MSE=float(meanseed[:,0].mean()),bootstrap=20000,one_sided_confidence=1-base.ALPHA,
        checkpoint_selection_set='CDEV64, same parents as this development evaluation; independent D still required',new_RHS=0,D_GT_read=False,
        timing_qualification='MEASURED_CURRENT_LOAD inference; Native cold phase+cache validation/preparation composition; not standalone end-to-end timing',
        original_reference_training_n=1024,new_matched_training_n=512,new_provider_updates=0,
        phase_definition='Original full96/64 anchor common masks; phase stratification retained in original CDEV artifacts',
        training_ledger=str(TRACK/'MATCHED_TRAINING_LEDGER.jsonl'),seal=str(seal_path))
    base.write(args.output/'MATCHED_CDEV_RESULTS.json',result)
    print(json.dumps(result),flush=True)


if __name__=='__main__':main()
