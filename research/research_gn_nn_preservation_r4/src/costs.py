"""Collect actual paid stages without treating reused GN as free or double work."""
from pathlib import Path
import csv,json

def read(path):return json.loads(Path(path).read_text())

def collect(study):
    study=Path(study);rows=[];seen_setups=set()
    def emit(stage,identifier,wall=None,counts=None,scope='actual execution',additive=True,**extra):
        c=counts or {}
        rows.append(dict(stage=stage,identifier=identifier,wall_seconds=wall,
                         LU=c.get('LU'),complex_RHS=c.get('actual_complex_RHS',c.get('complex_RHS')),
                         additive_execution=additive,scope=scope,all_counters=json.dumps(c),**extra))
    f=study/'data/GENERATION_RECEIPT.json'
    if f.exists():
        d=read(f);emit('full_VIE_data_generation','all690mothers',d.get('wall_seconds'),d.get('counts',{}))
    f=study/'data/reference_jacobian.json'
    if f.exists():
        d=read(f);seen_setups.add(json.dumps(d,sort_keys=True));emit('public_background_J_setup','local_reference',d.get('wall_seconds'),d.get('counts',{}),scope='actual public reference once on this host')
    for f in (study/'anchors').glob('BATCH_*.json'):
        d=read(f);setup=d.get('shared_reference_setup',{})
        identity=json.dumps(setup,sort_keys=True)
        if setup.get('status')=='PREPARED_PUBLIC_REFERENCE'and identity not in seen_setups:
            seen_setups.add(identity)
            emit('public_background_J_setup',f.stem,setup.get('wall_seconds'),setup.get('counts',{}),scope='actual public reference setup in this batch; reused setup records excluded')
    for f in sorted((study/'anchors').glob('obs_*.json')):
        d=read(f);emit('GN_or_traditional_anchor',f.stem,d.get('solve_wall_seconds'),d.get('counts'),
                     status=d.get('status'),stop_reason=d.get('stop_reason'),origin=d.get('execution_origin',d.get('origin')))
        emit('anchor_public_features',f.stem,d.get('feature_wall_seconds'),d.get('feature_counts',{}),
             scope='gb actually computed once here; later method budgets may reference its cost')
    for f in (study/'models').glob('*/summary.json'):
        d=read(f);emit('NN_training',f.parent.name,d.get('training_seconds'),{},
                     scope='training, not online; optimizer and validation are included',epochs=d.get('epochs_run'))
    f=study/'models/shared_setup.json'
    if f.exists():
        d=read(f);emit('TRAIN_normalization_and_operator_setup','shared',d.get('seconds'),{},scope='shared preparation, not multiplied by model count')
    for f in (study/'predictions').glob('*/inference_setup.json'):
        d=read(f);emit('NN_model_load',f.parent.name,d.get('checkpoint_and_operator_setup_seconds'),{},
                     scope='model/checkpoint preparation; operator cache reuse status in source receipt')
    for f in (study/'predictions').glob('*/obs_*.json'):
        d=read(f);c=d.get('preparation_cost',{});wall=sum(c.get(k,0.)for k in ('feature_normalization_seconds','inference_seconds','output_transfer_seconds'))
        # Posthoc projections actually happened too; the original NN prediction is a shared source.
        if 'posthoc_projection_seconds'in c:wall=c['posthoc_projection_seconds']
        emit('NN_prediction',f.parent.name+'/'+f.stem,wall,{},
             scope='actual normalization/inference/transfer only; inherited gb is already counted at anchor',
             method_budget_preparation_seconds=c.get('wall_seconds'))
    for f in (study/'evaluation').glob('*/obs_*.json'):
        d=read(f);c=d.get('verification_costcounts',d.get('costcounts',{}))
        wall=d.get('evaluator_wall_seconds',d.get('verification_wall_seconds',d.get('wall_seconds')))
        emit('physical_verification_correction',f.parent.name+'/'+f.stem,wall,c,
             scope='actual evaluator execution; prior inference is separate',status=d.get('status'),
             method_budget_seconds=d.get('wall_seconds'),budget_scope=d.get('cost_scope','post_GN_increment'),
             wall_overshoot=d.get('wall_overshoot_seconds'))
    for sub in ('e0','test_energy'):
        f=study/'results'/sub/'RECEIPT.json'
        if f.exists():
            d=read(f);emit('offline_energy_diagnostic',sub,d.get('wall_seconds'),d.get('counts',d.get('model_counts',{})),scope='scientific diagnostic, excluded from deployment cost')
    f=study/'audit/ANCHOR_PROTOCOL_CORRECTION.json'
    if f.exists():
        d=read(f)
        # Completed attempts already exist as suffixed canonical anchor receipts.
        emit('preparation_protocol_correction','17completed_reused_receipts',d.get('sum_completed_solve_seconds'),d.get('completed_counts'),
             scope='nonadditive reference: completed attempts counted individually as TRAIN discrepancy anchors',additive=False)
        emit('interrupted_preparation','up_to4_inflight',None,{},scope='unknown partial costs, not zero',status='UNKNOWN_PARTIAL_COST')
    keys=list(dict.fromkeys(k for r in rows for k in r));dest=study/'results/COST_LEDGER.csv';dest.parent.mkdir(exist_ok=True)
    with dest.open('w',newline='')as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
    notes={'rows':len(rows),'wall_warning':'sum of task seconds is not elapsed wall; parallel jobs overlap',
           'actual_vs_deployment':'inherited common GN/features counted once in actual execution; each method must still pay them for independent deployment',
           'unknown_costs':'interrupted preparation partial work is unknown, never zero',
           'shared_setup_warning':'operator/model shared receipt costs are not multiplied by number of predictions',
           'training_hardware':'local MPS; source receipts authoritative; remote GPU not taken from other jobs',
           'other_current_costs':'CIE Fourier/SVD/Green current actions retained in all_counters or original current_phase receipt; not relabeled as LU/RHS',
           'precision':'NN f32 real-pair/complex64; final Maxwell f64/complex128'}
    (dest.parent/'COST_NOTES.json').write_text(json.dumps(notes,indent=2)+'\n')
    return dest
