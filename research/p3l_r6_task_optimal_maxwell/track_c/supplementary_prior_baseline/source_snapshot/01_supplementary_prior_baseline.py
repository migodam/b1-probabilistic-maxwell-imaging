"""Closed-by-default prospective strong-baseline supplement,64 per noise condition.

PREPARED ONLY. Root may choose NOT_RUN. A root SUPPLEMENT_FREEZE.json must
select direct-Probit cap/mode from CDEV ONLY before any D GT unlock. The
unchanged primary D family, source and predictions are references, never edited.
No new truth generation. Root alone runs central paid D-stage physical paths.
"""
from __future__ import annotations
import os
for _name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_name]='1'
import argparse,csv,gc,hashlib,io,json,shutil,sys,time
from pathlib import Path
import numpy as np
OUT=Path(__file__).resolve().parents[1];TRACK=OUT/'track_c';PRIMARY=OUT/'track_d'
ROOT=TRACK/'supplementary_prior_baseline';sys.path.insert(0,str(OUT))
GROUPS={'DLOW':.5,'DHIGH':2.};SEEDS=(2026102601,2026102602,2026102603)
LABEL='SUPPLEMENTARY_PRE_GT_EXPLORATORY'


def write(path,value):
    from shared.physics import write as original_write
    original_write(path,value)


def identity(path):
    path=Path(path);content=path.read_bytes()
    return dict(path=str(path.resolve()),bytes=len(content),blake2b=hashlib.blake2b(content).hexdigest())


def matching_identity(value):
    if identity(value['path'])!=value:raise ValueError(f'changed frozen artifact: {value["path"]}')


def source_parents(group):
    rows=json.loads((OUT/f'configs/parents_{group}.json').read_text())['parents']
    parents=[str(row['parent']) for row in rows]
    if len(parents)!=256 or len(set(parents))!=256:raise ValueError('Original256-parent noise manifest required')
    return parents[:64]


def no_GT_yet():
    if (PRIMARY/'GT_UNLOCK.json').exists() or (PRIMARY/'scores').exists():
        raise RuntimeError('CLOSED: supplement freeze/physical/predict cannot start after any main D GT unlock/score')
    if (ROOT/'scores').exists():raise RuntimeError('supplement scoring already started')


def primary_seals(primary_freeze):
    path=PRIMARY/'ALL_512_PREDICTIONS_SEALED.json'
    seal=json.loads(path.read_text())
    if seal.get('status')!='ALL_512_PREDICTIONS_SEALED_BEFORE_D_GT' or seal.get('GT_read') is not False:
        raise RuntimeError('all512 PRIMARY predictions must be sealed before supplement operations')
    if seal['freeze_identity']!=identity(PRIMARY/'D_FREEZE.json'):raise ValueError('primary freeze identity')
    for row in seal['group_seals']:matching_identity(row)
    groups={}
    for group in GROUPS:
        value=json.loads((PRIMARY/'predictions'/group/'PREDICTIONS_SEALED.json').read_text())
        allparents=[str(row['parent']) for row in json.loads((OUT/f'configs/parents_{group}.json').read_text())['parents']]
        if value['parents']!=allparents or value['freeze_identity']!=seal['freeze_identity']:raise ValueError('primary group parent/freeze identity')
        expected={(arm,seed) for arm in primary_freeze['estimator_arms'] for seed in SEEDS}
        if {(r['method'],r['seed']) for r in value['records']}!=expected or len(value['records'])!=len(expected):
            raise ValueError('complete unchanged primary arm/seed family required')
        for row in value['records']:matching_identity(row['identity'])
        matching_identity(value['phases_identity']);groups[group]=value
    return seal,groups


def freeze(operation):
    path=ROOT/'SUPPLEMENT_FREEZE.json'
    if not path.exists():raise RuntimeError('NOT_RUN: root SUPPLEMENT_FREEZE.json absent')
    d=json.loads(path.read_text())
    if (d.get('status')!='ROOT_SUPPLEMENT_FROZEN_PRE_GT' or d.get('execution_authorized') is not True or
            d.get('scope')!=LABEL or d.get('CDEV_selected_baseline_competes') is not True):
        raise RuntimeError('NOT_RUN: explicit root pre-GT strong-baseline selection required')
    if d.get('selection_scope')!='CDEV_ONLY' or d.get('selected_cap') not in (64,128) or d.get('selected_mode') not in ('actual','q'):
        raise ValueError('cap64/128 and actual/q must be frozen from CDEV only')
    if d.get('method')!='deterministic_matched' or d.get('seeds')!=list(SEEDS) or d.get('groups')!=list(GROUPS):raise ValueError('original frozen deterministic3seed,two conditions only')
    if d.get('parents_per_condition')!=64 or d.get('primary_D_unchanged') is not True:raise ValueError('supplement64/condition; primary contract unchanged')
    expected_stats=dict(bootstrap=20000,rng_seed=202610270798,per_comparison_alpha=.05/8,
        NI_fraction=.01,wall_fraction=.8,qualification='separate exploratory supplementary family; not a full512 confirmation or combined-family FWER claim')
    if d.get('statistics')!=expected_stats:raise ValueError('explicit supplementary exploratory statistical contract required')
    if d['primary_D_freeze_identity']!=identity(PRIMARY/'D_FREEZE.json') or d['primary_global_seal_identity']!=identity(PRIMARY/'ALL_512_PREDICTIONS_SEALED.json'):
        raise ValueError('supplement must bind exact unchanged primary D freeze/global seal')
    selection_path=TRACK/'probit_prior_baseline/evaluation_CDEV/RESULTS.json'
    if d['CDEV_selection_results_identity']!=identity(selection_path):raise ValueError('exact completed direct-prior CDEV selection results required')
    for row in d.get('source_identities',[]):matching_identity(row)
    if identity(Path(__file__)) not in d.get('source_identities',[]):raise ValueError('root must freeze this standalone supplement source')
    # Original dependencies were already fixed before the direct CDEV queries.
    dev_freeze=json.loads((TRACK/'probit_prior_baseline/DEVELOPMENT_SOURCE_FREEZE.json').read_text())
    for row in dev_freeze['files']:matching_identity(row)
    for group in GROUPS:
        if d['selected_parents'].get(group)!=source_parents(group):raise ValueError('first64 metadata-index parents only; no GT-based selection')
        matching_identity(d['manifest_identities'][group])
    from track_d.common import freeze as primary_contract
    primary=primary_contract();main_seal,groups=primary_seals(primary)
    if operation in ('anchors','predict'):
        no_GT_yet()
        if d.get('GT_unlock_absent_at_freeze') is not True:raise ValueError('root pre-GT freeze attestation required')
    return d,primary,main_seal,groups


def complete_context(d,group,parent):
    cap=d['selected_cap'];folder=ROOT/'contexts'/group/parent/f'prior_x0_probit{cap}'
    report_path=folder/'REPORT.json'
    if not report_path.exists():return False
    report=json.loads(report_path.read_text());closed=OUT/'raw/jobs'/f'SUPP_prior_{parent}_P{cap}'/'CLOSED.json'
    if (not closed.exists() or json.loads(closed.read_text()).get('status')!='COMPLETE' or report.get('GT_read') is not False or
            report.get('path_complete') is not True or report.get('context_identity')!=identity(folder/'CONTEXT.npz') or
            report.get('supplement_freeze_identity')!=identity(ROOT/'SUPPLEMENT_FREEZE.json')):
        raise RuntimeError('HOLD: existing supplementary context requires explicit interrupted-job/identity audit; never replay')
    return True


def budget(d):
    from shared.physics import Book
    ledger=Book.totals();authorization=json.loads((OUT/'configs/AUTHORIZATION.json').read_text())
    cap=int(d['selected_cap']);worst=40*(cap+1);pending=[]
    for group in GROUPS:
        for parent in source_parents(group):
            if not complete_context(d,group,parent):pending.append((group,parent))
    c_used=ledger['stages'].get('C',{}).get('RHS',0);d_used=ledger['stages'].get('D',{}).get('RHS',0)
    if ledger.get('active_jobs') or ledger.get('active_actions'):raise RuntimeError('HOLD: central physical/receipt job still active')
    c_limit=min(3200000,authorization['stage_caps']['C']);d_limit=min(4100000,authorization['stage_caps']['D']);total_limit=min(8000000,authorization['total_complex_RHS'])
    available=min(d_limit-d_used,total_limit-ledger['RHS'])
    record=dict(scope=LABEL,selected_cap=cap,pending=len(pending),worst_RHS_per_parent=worst,
        worst_remaining_RHS=len(pending)*worst,available_D_and_total_RHS=available,
        C_used=c_used,C_hard_cap=c_limit,D_used=d_used,D_hard_cap=d_limit,total_used=ledger['RHS'],total_hard_cap=total_limit,
        decision='RUN' if c_used<=c_limit and len(pending)*worst<=available else 'HOLD',new_truth_generation_RHS=0)
    write(ROOT/'PRE_RUN_BUDGET.json',record)
    if record['decision']!='RUN':raise RuntimeError('HOLD: authorized C3.2m/D4.1m/total8m hard constraints cannot reserve complete supplement')
    return pending


def snapshot_sources(d):
    destination=ROOT/'source_snapshot'
    if destination.exists():return
    destination.mkdir(parents=True)
    copied=[]
    paths=[ROOT/'SUPPLEMENT_FREEZE.json',Path(__file__),TRACK/'probit_prior_baseline.py',TRACK/'probit_prior_baseline/DEVELOPMENT_PROTOCOL.json']
    for index,path in enumerate(paths):
        target=destination/f'{index:02d}_{path.name}';shutil.copy2(path,target)
        copied.append(dict(original=identity(path),snapshot=identity(target)))
    write(destination/'SOURCE_IDENTITY.json',dict(scope=LABEL,records=copied,primary_D_source_edited=False))


def anchors(d):
    from shared.physics import Book,Target,SIGMA
    from track_c.probit_prior_baseline import run
    pending=budget(d);snapshot_sources(d);cap=int(d['selected_cap']);name=f'prior_x0_probit{cap}'
    for index,(group,parent) in enumerate(pending):
        no_GT_yet()
        with np.load(OUT/'raw/online'/f'{parent}.npz',allow_pickle=False) as z:
            if str(z['parent'])!=parent:raise ValueError('existing observation identity')
            y,sigma=z['observation'].copy(),float(z['sigma'])
        if not np.isclose(sigma/SIGMA,GROUPS[group],rtol=0,atol=1e-13):raise ValueError('existing noise-scale contract')
        destination=ROOT/'contexts'/group/parent/name
        book=Book('D',f'SUPP_prior_{parent}_P{cap}',parent,'SUPPLEMENTARY_PRE_GT_EXPLORATORY',max_wall=120)
        status,error='COMPLETE',None;tick=time.perf_counter()
        try:
            target=Target(book,sigma=sigma);context,report=run(target,y,parent,destination,cap)
            cleanup=time.perf_counter();target.backend()._transfers.clear();del target,context;collected=gc.collect()
            write(destination/'REPORT.json',dict(report,scope=LABEL,group=group,sigma=sigma,GT_read=False,
                supplement_freeze_identity=identity(ROOT/'SUPPLEMENT_FREEZE.json'),context_identity=identity(destination/'CONTEXT.npz'),
                primary_D_freeze_identity=d['primary_D_freeze_identity'],cleanup_seconds=time.perf_counter()-cleanup,GC_collected_objects=collected,
                end_to_end_anchor_wall_seconds=time.perf_counter()-tick,
                wall_scope='standalone prior-state+Probit path+context persistence+cleanup; no cached primary baseline cost subtraction'))
        except BaseException as exc:status,error='FAILED',repr(exc);raise
        finally:book.close(status,error)
        if (index+1)%16==0:print(json.dumps(dict(event='SUPPLEMENT_ANCHORS',completed=index+1,pending_at_start=len(pending),RHS=Book.totals()['RHS'])),flush=True)
    for group in GROUPS:
        for parent in source_parents(group):
            if not complete_context(d,group,parent):raise ValueError('complete supplementary context family')
        write(ROOT/f'ANCHORS_COMPLETE_{group}.json',dict(status='COMPLETE',parents=source_parents(group),n=64,
              scope=LABEL,config=name,freeze_identity=identity(ROOT/'SUPPLEMENT_FREEZE.json'),GT_read=False))


def predict(d):
    from track_c.frozen_heads import load_action,predict_contexts,ONLINE_KEYS
    from track_c.kkt_closure.closure_action import Action as QAction
    cap=d['selected_cap'];mode=d['selected_mode'];name=f'prior_x0_probit{cap}';global_path=ROOT/'ALL_128_SUPPLEMENT_PREDICTIONS_SEALED.json'
    if global_path.exists():raise FileExistsError('supplement predictions already globally sealed; never rewrite')
    for group in GROUPS:
        marker=json.loads((ROOT/f'ANCHORS_COMPLETE_{group}.json').read_text())
        if marker.get('parents')!=source_parents(group) or marker.get('freeze_identity')!=identity(ROOT/'SUPPLEMENT_FREEZE.json'):raise ValueError('both64-condition contexts must complete before prediction')
        for parent in source_parents(group):
            if not complete_context(d,group,parent):raise RuntimeError('complete paid supplementary context receipts required')
    group_seals=[]
    for group in GROUPS:
        no_GT_yet();ids=source_parents(group);destination=ROOT/'predictions'/group
        if destination.exists():raise RuntimeError('partial or existing supplement prediction directory: explicit audit required')
        destination.mkdir(parents=True);contexts=[];costs=[];context_ids=[]
        for parent in ids:
            path=ROOT/'contexts'/group/parent/name/'CONTEXT.npz';report=json.loads(path.with_name('REPORT.json').read_text())
            if report['context_identity']!=identity(path) or report['supplement_freeze_identity']!=identity(ROOT/'SUPPLEMENT_FREEZE.json'):raise ValueError('supplement context identity/freeze')
            with np.load(path,allow_pickle=False) as z:contexts.append({key:z[key].copy() for key in (*ONLINE_KEYS,'parent')})
            costs.append(report['end_to_end_anchor_wall_seconds']);context_ids.append(dict(parent=parent,context=identity(path),report=identity(path.with_name('REPORT.json'))))
        records=[]
        for seed in SEEDS:
            tick=time.perf_counter()
            action=load_action('deterministic_matched',seed) if mode=='actual' else QAction('deterministic_matched',seed,'q',GROUPS[group])
            load=time.perf_counter()-tick;tick,cpu=time.perf_counter(),time.process_time()
            prediction=predict_contexts(contexts,action,ids,batch_size=64);wall,CPU=time.perf_counter()-tick,time.process_time()-cpu
            if prediction.shape!=(64,576) or not np.isfinite(prediction).all():raise ValueError('finite supplementary prediction; no fallback')
            path=destination/f'priorP{cap}_{mode}_det_s{seed}.npz'
            np.savez_compressed(path,parent=np.array(ids),prediction=prediction)
            records.append(dict(seed=seed,identity=identity(path),mode=mode,load_wall_seconds=load,prediction_wall_seconds=wall,prediction_CPU_seconds=CPU,
                per_parent_warm_learning_seconds=wall/64,per_parent_cold_anchor_seconds=costs,GT_read=False,new_RHS=0))
            del action
        seal=dict(status='SUPPLEMENT_GROUP_SEALED_PRE_GT',scope=LABEL,group=group,parents=ids,n=64,records=records,contexts=context_ids,
                  freeze_identity=identity(ROOT/'SUPPLEMENT_FREEZE.json'),GT_read=False)
        write(destination/'PREDICTIONS_SEALED.json',seal);group_seals.append(identity(destination/'PREDICTIONS_SEALED.json'))
    no_GT_yet()
    write(global_path,dict(status='ALL_128_SUPPLEMENT_PREDICTIONS_SEALED_PRE_GT',scope=LABEL,group_seals=group_seals,
        freeze_identity=identity(ROOT/'SUPPLEMENT_FREEZE.json'),primary_global_seal_identity=d['primary_global_seal_identity'],GT_read=False,created_utc=time.time()))


def validate_supplement(d):
    path=ROOT/'ALL_128_SUPPLEMENT_PREDICTIONS_SEALED.json';seal=json.loads(path.read_text())
    if seal.get('status')!='ALL_128_SUPPLEMENT_PREDICTIONS_SEALED_PRE_GT' or seal.get('GT_read') is not False or seal['freeze_identity']!=identity(ROOT/'SUPPLEMENT_FREEZE.json'):
        raise ValueError('exact complete128 supplementary seal required')
    for row in seal['group_seals']:matching_identity(row)
    snapshots={};groups={}
    for group in GROUPS:
        value=json.loads((ROOT/'predictions'/group/'PREDICTIONS_SEALED.json').read_text());ids=source_parents(group)
        if value['parents']!=ids or {r['seed'] for r in value['records']}!=set(SEEDS):raise ValueError('supplement parent/seed family')
        for row in value['contexts']:
            matching_identity(row['context']);matching_identity(row['report'])
        for row in value['records']:
            matching_identity(row['identity'])
            with np.load(row['identity']['path'],allow_pickle=False) as z:
                if z['parent'].astype(str).tolist()!=ids or z['prediction'].shape!=(64,576):raise ValueError('supplement prediction schema')
                snapshots[group,row['seed']]=z['prediction'].copy()
        groups[group]=value
    return seal,groups,snapshots


def bootstrap(values,seed,replicates=20000):
    rng=np.random.default_rng(seed);draws=[]
    for start in range(0,replicates,256):
        indices=rng.integers(len(values),size=(min(256,replicates-start),len(values)))
        draws.append(values[indices].mean(1))
    return np.concatenate(draws)


def score(d,primary,main_groups):
    unlock_path=PRIMARY/'GT_UNLOCK.json';unlock=json.loads(unlock_path.read_text())
    if unlock.get('status')!='ROOT_UNLOCKED_AFTER_SEAL' or unlock.get('global_seal_identity')!=identity(PRIMARY/'ALL_512_PREDICTIONS_SEALED.json'):
        raise RuntimeError('root main GT unlock required; supplement cannot open truth itself')
    seal,groups,supplement=validate_supplement(d)
    if seal['created_utc']>unlock_path.stat().st_mtime:raise RuntimeError('supplement seal must predate main GT unlock')
    if (ROOT/'scores').exists():raise FileExistsError('supplement scores already exist; no implicit rerun')
    # Snapshot only the two unchanged main references for fixed first64 parents;
    # all main family identities were checked by primary_seals() before this.
    references={};reference_cost={}
    for group,value in main_groups.items():
        ids=source_parents(group)
        for method in (primary['primary_selected'],primary['primary_full']):
            for row in value['records']:
                if row['method']!=method:continue
                with np.load(row['identity']['path'],allow_pickle=False) as z:
                    order={str(parent):i for i,parent in enumerate(z['parent'])};indices=np.array([order[p] for p in ids]);references[group,method,row['seed']]=z['prediction'][indices].copy()
                reference_cost[group,method,row['seed']]=np.array(row['per_parent_cold_anchor_seconds'])[indices]+row['per_parent_warm_learning_seconds']
    (ROOT/'scores').mkdir();write(ROOT/'scores/GT_OPEN_STARTED.json',dict(scope=LABEL,supplement_seal=identity(ROOT/'ALL_128_SUPPLEMENT_PREDICTIONS_SEALED.json'),
        primary_seal=d['primary_global_seal_identity'],all128_supplement_and_all512_primary_seals_verified=True,utc=time.time()))
    records=[];alpha=d['statistics']['per_comparison_alpha']
    for gi,group in enumerate(GROUPS):
        ids=source_parents(group);truth=[]
        # FIRST NEW D TRUTH OPEN: root has unlocked both globally sealed phases.
        for parent in ids:
            with np.load(OUT/'raw/offline_truth'/f'{parent}.npz',allow_pickle=False) as z:
                if str(z['parent'])!=parent:raise ValueError('supplement truth parent identity')
                truth.append(z['chi'].copy())
        truth=np.array(truth);names=['supplement',primary['primary_selected'],primary['primary_full']]
        mse=np.zeros((64,3,3));mae=np.zeros_like(mse);seconds=np.zeros_like(mse);ensemble=np.zeros((64,3))
        for j,name in enumerate(names):
            predictions=[]
            for k,seed in enumerate(SEEDS):
                if j==0:
                    pred=supplement[group,seed];row=next(r for r in groups[group]['records'] if r['seed']==seed)
                    wall=np.array(row['per_parent_cold_anchor_seconds'])+row['per_parent_warm_learning_seconds']
                else:pred=references[group,name,seed];wall=reference_cost[group,name,seed]
                error=abs(pred-truth);mse[:,j,k]=(error**2).mean(1);mae[:,j,k]=error.mean(1);seconds[:,j,k]=wall;predictions.append(pred)
            ensemble[:,j]=(abs(np.mean(predictions,0)-truth)**2).mean(1)
        means=mse.mean(2);times=seconds.mean(2);values=np.column_stack((means,times));draws=bootstrap(values,d['statistics']['rng_seed']+gi)
        comparisons={}
        for j,name in enumerate(names[1:],1):
            risk=draws[:,0]-1.01*draws[:,j];wall=draws[:,3]-.8*draws[:,3+j];difference=draws[:,0]-draws[:,j]
            comparisons[name]=dict(MSE_ratio=float(means[:,0].mean()/means[:,j].mean()),
                MSE_difference=float((means[:,0]-means[:,j]).mean()),MSE_difference_interval=np.quantile(difference,[alpha/2,1-alpha/2]).tolist(),
                NI1pct_linear_contrast_upper=float(np.quantile(risk,1-alpha)),
                wall20pct_linear_contrast_upper=float(np.quantile(wall,1-alpha)),
                NI1pct_diagnostic=bool(np.quantile(risk,1-alpha)<0),wall20pct_diagnostic=bool(np.quantile(wall,1-alpha)<0))
        np.savez_compressed(ROOT/'scores'/f'{group}_PARENT_LOSSES.npz',parent=np.array(ids),arms=np.array(names),seeds=np.array(SEEDS),
            MSE_seedwise=mse,MAE_seedwise=mae,MSE_meanseed=means,MAE_meanseed=mae.mean(2),MSE_ensemble=ensemble,whole_wall_meanseed=times)
        record=dict(status=LABEL,group=group,n=64,selection='first64 metadata-index parents; cap/mode chosen on CDEV only before GT',
            selected_cap=d['selected_cap'],selected_mode=d['selected_mode'],mean_MSE=dict(zip(names,map(float,means.mean(0)))),
            mean_MAE=dict(zip(names,map(float,mae.mean((0,2))))),mean_whole_seconds=dict(zip(names,map(float,times.mean(0)))),
            comparisons=comparisons,statistics=d['statistics'],seed_mean_loss_primary=True,ensemble_separate_descriptive=True,
            full512_qualification=False,primary_D_family_or_source_changed=False,no_noise_condition_pooling=True,
            probability_or_risk_certificate_claim=False,new_training=0,new_scoring_RHS=0,
            wall_scope='separately paid standalone supplement path+cleanup vs already paid primary references; batch64 warm head timing; no repeated hardware uncertainty estimated')
        write(ROOT/'scores'/f'{group}_RESULTS.json',record);records.append(record)
    write(ROOT/'scores/SCORING_COMPLETE.json',dict(status=LABEL,conditions=list(GROUPS),n_per_condition=64,
        primary512_qualification=False,global_primary_D_results_unchanged=True,new_truth_generation=0,new_training=0,new_scoring_RHS=0))
    print(json.dumps(records),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('operation',choices=['plan','anchors','predict','score']);args=parser.parse_args()
    if args.operation=='plan':
        print(json.dumps(dict(status='PREPARED_ONLY_CLOSED_UNTIL_ROOT_SUPPLEMENT_FREEZE',scope=LABEL,parents_per_condition=64,conditions=list(GROUPS),
            selected_cap_options=[64,128],selected_feature_options=['actual','q'],worst128_parent_RHS_at_cap128=660480,
            new_truth_generation=0,primary_D_source_changed=False),indent=2));return
    d,primary,_,groups=freeze(args.operation)
    if args.operation=='anchors':anchors(d)
    elif args.operation=='predict':predict(d)
    else:score(d,primary,groups)


if __name__=='__main__':main()
