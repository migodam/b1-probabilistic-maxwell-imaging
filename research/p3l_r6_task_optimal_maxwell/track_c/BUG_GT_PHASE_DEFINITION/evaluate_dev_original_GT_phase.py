"""CDEV64 online prediction -> global BLAKE2b seal -> offline paired scoring.

No paid Target/Book, solver, training module, or locked-D evaluator is imported.
First truth open is structurally after all 50 original frozen prediction files
have been written and revalidated against the global seal. Existing seals permit
re-scoring identical byte snapshots, never prediction replacement after exposure.
"""
from __future__ import annotations
import os
for _key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS',
             'VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_key] = '1'
import argparse
import csv
import hashlib
import io
import json
import sys
import time
from pathlib import Path
import numpy as np

OUT = Path(__file__).resolve().parents[1]
REPO = OUT.parents[1]
sys.path.insert(0, str(OUT))
from track_c.anchors import CONFIGS
from track_c.frozen_heads import load_action, predict_contexts, ONLINE_KEYS, SEEDS

CONFIG_NAMES = tuple(CONFIGS)
METHODS = ('physics','strongCNN','distilled_w128_n32768','deterministic_matched')
PHASES = ('outside_roi','host_plateau','transition','inclusion_plateau')
BASELINE = 'native96_probit64'
ALPHA = 0.05/3


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')
    temp.replace(path)


def expected_parents():
    manifest = json.loads((OUT/'configs/parents_CDEV.json').read_text())
    parents = [str(row['parent']) for row in manifest['parents']]
    if len(parents) != 64 or len(set(parents)) != 64 or any(not p.startswith('R6_CDEV_') for p in parents):
        raise ValueError('exact new CDEV64 namespace required')
    return parents


def require_complete(parents):
    marker_path = OUT/'track_c/ANCHORS_COMPLETE_CDEV.json'
    if not marker_path.exists():
        raise RuntimeError('NOT_READY: CDEV64 all-config completion marker missing')
    marker = json.loads(marker_path.read_text())
    if marker.get('status') != 'COMPLETE' or marker.get('n') != 64 or set(marker.get('configs',[])) != set(CONFIG_NAMES):
        raise RuntimeError('NOT_READY: exact CDEV64 x five configurations required')
    reports = {}
    for config in CONFIG_NAMES:
        reports[config] = []
        for parent in parents:
            folder = OUT/'track_c/contexts/CDEV'/parent/config
            report = json.loads((folder/'REPORT.json').read_text())
            if report.get('parent') != parent or report.get('GT_read') is not False or report.get('path_complete') is not True:
                raise ValueError('complete truth-free matching anchor report required')
            if not (folder/'CONTEXT.npz').exists():
                raise RuntimeError('NOT_READY: missing context')
            reports[config].append(report)
    return reports


def load_contexts(parents, config):
    contexts = []
    for parent in parents:
        path = OUT/'track_c/contexts/CDEV'/parent/config/'CONTEXT.npz'
        with np.load(path, allow_pickle=False) as z:
            c = {key:z[key].copy() for key in (*ONLINE_KEYS,'parent')}
        if str(c['parent']) != parent or c['observation'].shape != (32,8) or c['J'].shape != (512,66):
            raise ValueError('original full-acquisition context contract')
        if any(not np.isfinite(c[key]).all() for key in ONLINE_KEYS):
            raise ValueError('nonfinite online context')
        contexts.append(c)
    return contexts


def seal_predictions(destination, parents):
    """Online phase: this function has no truth path or truth loader."""
    destination.mkdir(parents=True, exist_ok=True)
    seal_path = destination/'ALL_PREDICTIONS_SEALED.json'
    if seal_path.exists():
        return json.loads(seal_path.read_text())
    if list(destination.glob('*.npz')):
        raise RuntimeError('partial predictions: explicit audit required before resuming; no overwrite')
    contexts = {config:load_contexts(parents,config) for config in CONFIG_NAMES}
    records, model_loads = [], []
    for method in METHODS:
        for seed in (SEEDS[:1] if method == 'physics' else SEEDS):
            wall0, cpu0 = time.perf_counter(), time.process_time()
            action = load_action(method, seed)
            model_loads.append(dict(method=method,seed=seed,wall_seconds=time.perf_counter()-wall0,
                                    CPU_seconds=time.process_time()-cpu0))
            for config in CONFIG_NAMES:
                predictions, wall, cpu = [], [], []
                # True single-parent CPU1 inference. A seed-mean loss represents
                # a randomly chosen frozen seed, not a three-seed ensemble job.
                for parent, context in zip(parents, contexts[config]):
                    wall0, cpu0 = time.perf_counter(), time.process_time()
                    prediction = predict_contexts([context],action,[parent],batch_size=1)[0]
                    wall.append(time.perf_counter()-wall0);cpu.append(time.process_time()-cpu0)
                    if prediction.shape != (576,) or not np.isfinite(prediction).all():
                        raise ValueError('finite full-grid material prediction required')
                    predictions.append(prediction)
                path = destination/f'{config}__{method}_s{seed}.npz'
                np.savez_compressed(path,parent=np.array(parents),prediction=np.asarray(predictions),
                                    inference_wall_seconds=np.array(wall),inference_CPU_seconds=np.array(cpu))
                payload = path.read_bytes()
                records.append(dict(path=path.name,config=config,method=method,seed=seed,
                                    bytes=len(payload),blake2b=hashlib.blake2b(payload).hexdigest()))
                print(json.dumps(dict(event='CDEV_PREDICTION_SAVED',config=config,method=method,seed=seed)),flush=True)
            del action
    seal = dict(status='ALL_PREDICTIONS_SEALED',parents=parents,configurations=list(CONFIG_NAMES),
                methods=list(METHODS),seeds=list(SEEDS),truth_read=False,new_physics_RHS=0,
                inference_backend='CPU1_SINGLE_PARENT',model_loads=model_loads,artifacts=records,
                created_unix=time.time(),hash_algorithm='BLAKE2b-512',SHA256_used=False)
    write(seal_path,seal)
    return seal


def validate_seal(destination, seal, parents):
    if seal.get('truth_read') is not False or seal.get('parents') != parents:
        raise ValueError('truth-free matching global seal required')
    expected = {(c,m,s) for c in CONFIG_NAMES for m in METHODS
                for s in (SEEDS[:1] if m == 'physics' else SEEDS)}
    if len(seal['artifacts']) != len(expected):
        raise ValueError('exact complete prediction family required')
    snapshots = {}
    for row in seal['artifacts']:
        key = (row['config'],row['method'],row['seed'])
        if key not in expected or key in snapshots or Path(row['path']).name != row['path']:
            raise ValueError('unexpected/duplicate/escaping sealed artifact')
        content = (destination/row['path']).read_bytes()
        if len(content) != row['bytes'] or hashlib.blake2b(content).hexdigest() != row['blake2b']:
            raise ValueError('prediction bytes changed after sealing')
        with np.load(io.BytesIO(content),allow_pickle=False) as z:
            value = {name:z[name].copy() for name in ('parent','prediction','inference_wall_seconds','inference_CPU_seconds')}
        if value['parent'].astype(str).tolist() != parents or value['prediction'].shape != (64,576):
            raise ValueError('sealed prediction parent/shape mismatch')
        if any(not np.isfinite(value[name]).all() for name in ('prediction','inference_wall_seconds','inference_CPU_seconds')):
            raise ValueError('nonfinite sealed predictions/timing')
        snapshots[key] = value
    if set(snapshots) != expected:
        raise ValueError('incomplete sealed family')
    return snapshots


def bootstrap(values, replicates=20000, seed=202610270608):
    """Identical independent-parent draws across every config/method/metric."""
    if replicates != 20000:
        raise ValueError('frozen 20000 bootstrap draws required')
    values = np.asarray(values,float)
    result = np.empty((replicates,*values.shape[1:]))
    rng = np.random.default_rng(seed)
    for start in range(0,replicates,256):
        size = min(256,replicates-start)
        indices = rng.integers(0,len(values),(size,len(values)))
        result[start:start+size] = values[indices].mean(1)
    return result


def make_prior():
    r4 = REPO.parent/'b1-p3l-r4-conditional-transport'
    sys.path.insert(0,str(r4/'src'))
    from p1_material.prior import MaterialPrior
    cfg = json.loads((r4/'research/p3l_r4_conditional_transport/configs/TARGET.json').read_text())
    points = json.loads((r4/'research/p1_material_target_20261018/inputs/geometry.json').read_text())['points']
    return MaterialPrior(cfg,points)


def frozen_learning_cost():
    path = REPO.parent/'b1-p3l-r5-parallel-ab/research/p3l_r5_parallel_ab/stage_b_material_risk/COST_LEDGER.json'
    ledger = json.loads(path.read_text())
    relevant = [r for r in ledger['records'] if r['operation'].startswith(('experiments/deterministic_matched_',
                  'learned_context/distilled_w128_n32768_')) or r['operation']=='TEACHER_GENERATION_32768.json']
    return dict(source=str(path),historical_records=relevant,
                selected_head_and_provider_recorded_wall_seconds=sum(r['wall_seconds'] for r in relevant),
                R4_provider_and_strongCNN_training='NOT_RECONSTRUCTED_FROM_THIS_SUBLEDGER',
                new_learning_CPU_seconds=0,new_training_RHS=0,amortized_training_per_prediction=None,
                amortization_status='NO_DEPLOYMENT_COUNT_ASSUMED; historical MPS and prior physics costs kept separate',
                matched_lowcost_head_training='NOT_RUN')


def evaluate_paths(parents, prior, truth, phases, destination):
    """All paid native/probit endpoint material risks by exact deterministic decode.

    This is offline descriptive analysis. It cannot feed stopping or selection;
    no counterfactual continuation and no new Maxwell solve occurs.
    """
    counts = dict(native=0,probit=0,accepted_iterations=0)
    path = destination/'ALL_EVALUATED_PATH_RISKS.jsonl'
    with path.open('w') as stream:
        for n,parent in enumerate(parents):
            for config in CONFIG_NAMES:
                folder = OUT/'track_c/contexts/CDEV'/parent/config
                files = (('native',folder/f'original_optimizer/inversion/{parent}/full_TRF/trace.jsonl','EVALUATION'),
                         ('probit',folder/'original_optimizer/probit/MODE_HISTORY.jsonl','EVALUATED'))
                seen = set()
                for phase,source,event in files:
                    for line in source.read_text().splitlines():
                        row = json.loads(line)
                        if row['event'] != event or (phase,row['endpoint']) in seen:
                            continue
                        seen.add((phase,row['endpoint']))
                        chi = prior.decode(np.array(row['q']))[0]
                        error = abs(chi-truth[n])
                        out = dict(parent=parent,config=config,optimizer_phase=phase,endpoint=row['endpoint'],
                                   material_complex_MSE=float(np.mean(error**2)),material_complex_MAE=float(np.mean(error)),
                                   source_trace=str(source),GT_offline_only=True,new_RHS=0,
                                   phase_risk={name:dict(cells=int(mask.sum()),
                                     MSE=float(np.mean(error[mask]**2)) if mask.any() else None,
                                     MAE=float(np.mean(error[mask])) if mask.any() else None)
                                     for name,mask in zip(PHASES,phases[n])})
                        for key in ('fresh_forward_states','fresh_state_attempts','wall_seconds','objective','Phi','native_objective','transformed_energy'):
                            if key in row:out[key]=row[key]
                        stream.write(json.dumps(out,allow_nan=False)+'\n');counts[phase]+=1
                # Accepted-step trace already contains actual prefix costs and
                # outputs, including any optional frozen material action output.
                for line in (folder/'path/PATH.jsonl').read_text().splitlines():
                    row = json.loads(line)
                    if row.get('output_file') is None:continue
                    with np.load(folder/'path'/row['output_file'],allow_pickle=False) as z:
                        error = abs(z['current_output']-truth[n])
                    stream.write(json.dumps(dict(parent=parent,config=config,event='ACCEPTED_CURRENT_OUTPUT_RISK',
                        optimizer_phase=row['phase'],iteration=row['iteration'],
                        source=row['current_output_source'],MSE=float(np.mean(error**2)),MAE=float(np.mean(error)),
                        prefix_resource_cost=row['path_resource_cost'],GT_offline_only=True,new_RHS=0),allow_nan=False)+'\n')
                    counts['accepted_iterations']+=1
    return dict(output=str(path),counts=counts,role='DESCRIPTIVE_OFFLINE_ONLY_NO_RULE_RETUNING')


def score(destination, parents, reports, snapshots, seal):
    # THE FIRST TRUTH OPEN. Every prediction has already been revalidated into
    # immutable in-memory snapshots. No old R5B path is constructed anywhere.
    prior = make_prior()
    truth, phases = [], []
    for parent in parents:
        with np.load(OUT/'raw/offline_truth'/f'{parent}.npz',allow_pickle=False) as z:
            chi,q = z['chi'].copy(),z['q'].copy()
            if str(z['parent']) != parent:raise ValueError('offline truth parent mismatch')
        decoded,meta = prior.decode(q)
        if not np.array_equal(decoded,chi):raise ValueError('truth material/declared decoder mismatch')
        s = meta['field']/prior.temperature
        truth.append(chi)
        phases.append(np.stack((~prior.mask,prior.mask&(s<=-1),prior.mask&(s>-1)&(s<1),prior.mask&(s>=1))))
    truth,phases = np.array(truth),np.array(phases)
    shape = (64,len(CONFIG_NAMES),len(METHODS),3)
    mse,mae = np.zeros(shape),np.zeros(shape)
    ensemble_mse,ensemble_mae = np.zeros(shape[:3]),np.zeros(shape[:3])
    inference_wall,inference_cpu = np.zeros(shape),np.zeros(shape)
    phase_sse,phase_sae = np.zeros((*shape,len(PHASES))),np.zeros((*shape,len(PHASES)))
    ensemble_phase_sse,ensemble_phase_sae = np.zeros((*shape[:3],len(PHASES))),np.zeros((*shape[:3],len(PHASES)))
    for c,config in enumerate(CONFIG_NAMES):
        for m,method in enumerate(METHODS):
            predictions = []
            for k,seed in enumerate(SEEDS):
                row = snapshots[(config,method,SEEDS[0] if method=='physics' else seed)]
                prediction = row['prediction'];predictions.append(prediction)
                error = abs(prediction-truth);mse[:,c,m,k] = (error**2).mean(1);mae[:,c,m,k] = error.mean(1)
                inference_wall[:,c,m,k] = row['inference_wall_seconds']
                inference_cpu[:,c,m,k] = row['inference_CPU_seconds']
                for p in range(len(PHASES)):
                    phase_sse[:,c,m,k,p] = (error**2*phases[:,p]).sum(1)
                    phase_sae[:,c,m,k,p] = (error*phases[:,p]).sum(1)
            error = abs(np.mean(predictions,0)-truth)
            ensemble_mse[:,c,m] = (error**2).mean(1);ensemble_mae[:,c,m] = error.mean(1)
            for p in range(len(PHASES)):
                ensemble_phase_sse[:,c,m,p] = (error**2*phases[:,p]).sum(1)
                ensemble_phase_sae[:,c,m,p] = (error*phases[:,p]).sum(1)
    anchor_wall = np.array([[reports[c][n]['end_to_end_anchor_wall_seconds'] for c in CONFIG_NAMES] for n in range(64)])
    whole_wall = anchor_wall[:,:,None]+inference_wall.mean(3)
    ensemble_wall = anchor_wall[:,:,None]+inference_wall.sum(3)
    # The deterministic physics action is one action, not three repeated calls.
    ensemble_wall[:,:,0] = whole_wall[:,:,0]
    meanseed_mse,meanseed_mae = mse.mean(3),mae.mean(3)
    boot = bootstrap(np.stack((meanseed_mse,meanseed_mae,whole_wall),-1))
    phase_counts = phases.sum(2)
    summaries,passing = [],[]
    for c,config in enumerate(CONFIG_NAMES):
        for m,method in enumerate(METHODS):
            risk_denominator = boot[:,0,m,0]
            risk_ratio = boot[:,c,m,0]/risk_denominator
            saving = 1-boot[:,c,m,2]/boot[:,0,m,2]
            risk_upper = float(np.quantile(risk_ratio,1-ALPHA))
            saving_lower = float(np.quantile(saving,ALPHA))
            risk_gate = risk_upper <= 1.01
            cost_gate = saving_lower >= .20
            row = dict(config=config,method=method,full_complex_MSE=float(meanseed_mse[:,c,m].mean()),
                full_complex_MAE=float(meanseed_mae[:,c,m].mean()),
                ensemble_MSE_descriptive=float(ensemble_mse[:,c,m].mean()),
                ensemble_MAE_descriptive=float(ensemble_mae[:,c,m].mean()),
                same_method_full_MSE=float(meanseed_mse[:,0,m].mean()),
                MSE_ratio=float(meanseed_mse[:,c,m].mean()/meanseed_mse[:,0,m].mean()),
                MSE_ratio_upper_one_sided_98_333333=risk_upper,
                MAE_ratio_secondary=float(meanseed_mae[:,c,m].mean()/meanseed_mae[:,0,m].mean()),
                mean_anchor_wall_seconds=float(anchor_wall[:,c].mean()),
                mean_single_seed_inference_wall_seconds=float(inference_wall[:,c,m].mean()),
                mean_single_seed_inference_CPU_seconds=float(inference_cpu[:,c,m].mean()),
                mean_whole_system_wall_seconds=float(whole_wall[:,c,m].mean()),
                mean_ensemble_whole_wall_seconds_descriptive=float(ensemble_wall[:,c,m].mean()),
                whole_wall_saving=float(1-whole_wall[:,c,m].mean()/whole_wall[:,0,m].mean()),
                wall_saving_lower_one_sided_98_333333=saving_lower,
                noninferiority_1pct_gate=bool(risk_gate),wall_reduction_20pct_gate=bool(cost_gate),
                CDEV_pass=bool(config != BASELINE and risk_gate and cost_gate),phase_risks={})
            for p,name in enumerate(PHASES):
                count=int(phase_counts[:,p].sum())
                row['phase_risks'][name]=dict(cells=count,
                    MSE=float(phase_sse[:,c,m,:,p].mean(1).sum()/count) if count else None,
                    MAE=float(phase_sae[:,c,m,:,p].mean(1).sum()/count) if count else None,
                    ensemble_MSE=float(ensemble_phase_sse[:,c,m,p].sum()/count) if count else None,
                    ensemble_MAE=float(ensemble_phase_sae[:,c,m,p].sum()/count) if count else None)
            summaries.append(row)
            if row['CDEV_pass']:passing.append(row)
    passing.sort(key=lambda row:(row['mean_whole_system_wall_seconds'],row['config'],row['method']))
    chosen = ({key:passing[0][key] for key in ('config','method','mean_whole_system_wall_seconds',
                 'MSE_ratio_upper_one_sided_98_333333','wall_saving_lower_one_sided_98_333333')} if passing else None)
    np.savez_compressed(destination/'PARENT_LOSSES.npz',parent=np.array(parents),configs=np.array(CONFIG_NAMES),methods=np.array(METHODS),
        seed=np.array(SEEDS),MSE_seedwise=mse,MAE_seedwise=mae,MSE_meanseed=meanseed_mse,MAE_meanseed=meanseed_mae,
        MSE_ensemble=ensemble_mse,MAE_ensemble=ensemble_mae,anchor_wall=anchor_wall,inference_wall=inference_wall,
        inference_CPU=inference_cpu,whole_wall_meanseed=whole_wall,whole_wall_ensemble=ensemble_wall,
        phases=np.array(PHASES),phase_counts=phase_counts,phase_squared_error_sum_seedwise=phase_sse,
        phase_absolute_error_sum_seedwise=phase_sae,phase_squared_error_sum_ensemble=ensemble_phase_sse,
        phase_absolute_error_sum_ensemble=ensemble_phase_sae)
    with (destination/'PARENT_METRICS.csv').open('w') as stream:
        writer=csv.writer(stream);writer.writerow(('parent','config','method','MSE_meanseed','MAE_meanseed','MSE_ensemble','MAE_ensemble','anchor_wall','whole_wall_single_seed'))
        for n,parent in enumerate(parents):
            for c,config in enumerate(CONFIG_NAMES):
                for m,method in enumerate(METHODS):writer.writerow((parent,config,method,meanseed_mse[n,c,m],meanseed_mae[n,c,m],ensemble_mse[n,c,m],ensemble_mae[n,c,m],anchor_wall[n,c],whole_wall[n,c,m]))
    with (destination/'PARENT_PHASE_METRICS.csv').open('w') as stream:
        writer=csv.writer(stream);writer.writerow(('parent','config','method','phase','cells','MSE_meanseed','MAE_meanseed','MSE_ensemble','MAE_ensemble'))
        for n,parent in enumerate(parents):
            for c,config in enumerate(CONFIG_NAMES):
                for m,method in enumerate(METHODS):
                    for p,name in enumerate(PHASES):
                        count=int(phase_counts[n,p])
                        values=((phase_sse[n,c,m,:,p].mean()/count,phase_sae[n,c,m,:,p].mean()/count,
                                 ensemble_phase_sse[n,c,m,p]/count,ensemble_phase_sae[n,c,m,p]/count)
                                if count else ('','','',''))
                        writer.writerow((parent,config,method,name,count,*values))
    result=dict(status='CDEV_FROZEN_HEAD_PASS' if chosen else 'CDEV_NO_FROZEN_HEAD_PASS',parents=64,
        bootstrap_replicates=20000,bootstrap_unit='independent parent; seed losses averaged before resampling',
        bootstrap_rng_seed=202610270608,one_sided_confidence=1-ALPHA,primary_metric='full-grid complex MSE',
        primary_estimator='mean seed loss, not seed-ensemble loss',selection=chosen,all_arms=summaries,
        selection_rule='minimum measured whole wall among arms with risk ratio upper<=1.01 and wall saving lower>=0.20',
        development_selection_only=True,locked_D_GT_read=False,old_R5B_GT_read=False,new_physics_RHS=0,
        frozen_head_matched_status='NOT_RUN',frozen_learning_cost=frozen_learning_cost(),
        seal=str(destination/'ALL_PREDICTIONS_SEALED.json'),model_load_costs_separate=seal['model_loads'],
        inference_timing_qualification='MEASURED_CURRENT_LOAD; other CPU1 research may run concurrently; clean 30-repeat matched benchmark required before final wall claim',
        scope='Warm loaded single-parent CPU1 inference; anchor wall includes root end-to-end context/report persistence. Loading/training amortization separate.',
        optimizer_paths=evaluate_paths(parents,prior,truth,phases,destination))
    write(destination/'CDEV_RESULTS.json',result)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=OUT/'track_c/evaluation_CDEV')
    parser.add_argument('--predict-only',action='store_true')
    parser.add_argument('--score-only',action='store_true')
    args=parser.parse_args()
    if args.predict_only and args.score_only:raise ValueError('choose at most one phase')
    destination=args.output.resolve()
    if not destination.is_relative_to((OUT/'track_c').resolve()):raise ValueError('Track C output only')
    parents=expected_parents();reports=require_complete(parents)
    if args.score_only:
        seal=json.loads((destination/'ALL_PREDICTIONS_SEALED.json').read_text())
    else:seal=seal_predictions(destination,parents)
    snapshots=validate_seal(destination,seal,parents)
    if args.predict_only:
        print(json.dumps(dict(status='ALL_PREDICTIONS_SEALED',parents=64,GT_read=False)),flush=True);return
    result=score(destination,parents,reports,snapshots,seal)
    print(json.dumps(dict(status=result['status'],selection=result['selection'],output=str(destination/'CDEV_RESULTS.json'))),flush=True)


if __name__=='__main__':main()
