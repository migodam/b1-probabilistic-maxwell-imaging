"""Collect actual paid stages without treating reused GN as free or double work."""
from pathlib import Path
import csv,json

def read(path):return json.loads(Path(path).read_text())

def collect(study):
    study=Path(study);rows=[];seen_setups=set()
    online=read(study/'data/manifests/online.json')['records']
    case_meta={r['obs_id']:r for r in online}
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
        d=read(f);obs=f.stem.split('_discrepancy')[0].split('_weakreg')[0]
        split=case_meta.get(obs,{}).get('split');primary=obs==f.stem
        purpose='offline_TRAIN_VAL_anchor_preparation'if split in('train','val')else'TEST_target_solver'
        emit('GN_or_traditional_anchor',f.stem,d.get('solve_wall_seconds'),d.get('counts'),
                     split=split,primary_anchor=primary,purpose=purpose,
                     status=d.get('status'),stop_reason=d.get('stop_reason'),origin=d.get('execution_origin',d.get('origin')))
        emit('anchor_public_features',f.stem,d.get('feature_wall_seconds'),d.get('feature_counts',{}),
             split=split,purpose='NN_public_input_features',
             scope='gb actually computed once here; later method budgets may reference its cost')
    for f in (study/'models').rglob('summary.json'):
        d=read(f);emit('NN_training',f.parent.name,d.get('training_seconds',d.get('training_wall_seconds')),
                     d.get('bridge_cost',{}).get('counts',{}),
                     scope='training, not online; optimizer and validation included, bridge wall is part of this training wall',epochs=d.get('epochs_run'))
    setup_files=list((study/'models').rglob('shared_setup_*.json'))
    setup_files += [f for f in (study/'models').rglob('shared_setup.json') if not list(f.parent.glob('shared_setup_*.json'))]
    for f in setup_files:
        d=read(f);emit('TRAIN_normalization_and_operator_setup',str(f.relative_to(study)),d.get('seconds',d.get('wall_seconds')),{},scope='actual preparation per invocation; latest aliases not counted twice')
    inference_files=list((study/'predictions').glob('*/inference_setup_*.json'))
    inference_files += [f for f in (study/'predictions').glob('*/inference_setup.json') if not list(f.parent.glob('inference_setup_*.json'))]
    for f in inference_files:
        d=read(f);emit('NN_model_load',f.parent.name,d.get('checkpoint_and_operator_setup_seconds'),{},
                     scope='model/checkpoint preparation; operator cache reuse status in source receipt')
    for f in (study/'predictions').glob('prediction_runs_*.json'):
        d=read(f);emit('NN_shared_operator_upload',f.stem,d.get('global_operator_setup_seconds'),{},scope='one immutable operator upload for this actual bulk inference invocation')
    for f in (study/'predictions').glob('*/obs_*.json'):
        d=read(f);c=d.get('preparation_cost',{});wall=sum(c.get(k,0.)for k in ('feature_normalization_seconds','inference_seconds','output_transfer_seconds'))
        # Posthoc projections actually happened too; the original NN prediction is a shared source.
        if 'posthoc_projection_seconds'in c:wall=c['posthoc_projection_seconds']
        emit('NN_prediction',f.parent.name+'/'+f.stem,wall,c.get('fixed_operator_cost',{}),
             split=case_meta.get(f.stem,{}).get('split'),
             scope='actual normalization/inference/transfer only; inherited gb is already counted at anchor',
             method_budget_preparation_seconds=c.get('wall_seconds'))
    for f in (study/'evaluation').glob('*/obs_*.json'):
        d=read(f);c=dict(d.get('verification_costcounts',d.get('costcounts',{})))
        for key,value in d.get('current_phase',{}).get('cost',{}).items():
            if isinstance(value,(int,float))and key!='wall_seconds':c['CIE_'+key]=value
        wall=d.get('evaluator_wall_seconds',d.get('verification_wall_seconds',d.get('wall_seconds')))
        emit('physical_verification_correction',f.parent.name+'/'+f.stem,wall,c,
             split=d.get('split',case_meta.get(f.stem,{}).get('split')),
             scope='actual evaluator execution; prior inference is separate',status=d.get('status'),
             method_budget_seconds=d.get('wall_seconds'),budget_scope=d.get('cost_scope','post_GN_increment'),
             wall_overshoot=d.get('wall_overshoot_seconds'))
    for sub in ('e0','test_energy'):
        f=study/'results'/sub/'RECEIPT.json'
        if f.exists():
            d=read(f);emit('offline_energy_diagnostic',sub,d.get('total_recorded_wall_seconds',d.get('wall_seconds')),d.get('physical_costcounts',d.get('model_counts',d.get('counts',{}))),scope='scientific diagnostic, excluded from deployment cost; measured phase sum if whole process was not recorded')
    f=study/'cache/projector_setup_creation.json'
    if f.exists():
        d=read(f);emit('fixed_Green_projector_preparation','creation',d.get('setup_seconds'),{},scope='actual once-only projector creation; later model receipts only load this cache')
    for f in (study/'results/recurrent_iterations').glob('*/RECEIPT.json'):
        d=read(f);emit('offline_recurrent_depth_diagnostic',f.parent.name,d.get('wall_seconds'),d.get('physical_costcounts',{}),scope='VAL public-input depth diagnostic, not deployment or extra training')
    for f in (study/'audit').glob('anchor_grad_*/RECEIPT.json'):
        d=read(f);emit('offline_projected_gradient_diagnostic',f.parent.name,d.get('wall_seconds'),d.get('costcounts',{}),scope='public anchor adjoint diagnostic; separate from original solver and deployment',status=d.get('status'))
    for f in (study/'audit').glob('REMOTE_TRANSFER_*.json'):
        d=read(f);action=d.get('action','unknown');queue=action.startswith('queue_')
        emit('remote_queue_elapsed'if queue else'remote_transfer_or_packaging',f.stem,d.get('wall_seconds'),{},
             scope='overlapping queue elapsed, not additive physical task time'if queue else'actual transport/setup overhead, not Maxwell computation',
             additive=not queue,status=d.get('status'))
    f=study/'audit/PUBLICATION_METADATA_ATTEMPT.json'
    if f.exists():
        d=read(f);emit('stopped_publication_metadata_probe',f.stem,None,{'LU':0,'actual_complex_RHS':0},scope='nonphysical metadata probe; exact elapsed not recorded',status=d['status'])
    f=study/'audit/ANCHOR_PROTOCOL_CORRECTION.json'
    if f.exists():
        d=read(f)
        # Completed attempts already exist as suffixed canonical anchor receipts.
        emit('preparation_protocol_correction','17completed_reused_receipts',d.get('sum_completed_solve_seconds'),d.get('completed_counts'),
             scope='nonadditive reference: completed attempts counted individually as TRAIN discrepancy anchors',additive=False)
        emit('interrupted_preparation','up_to4_inflight',None,{},scope='unknown partial costs, not zero',status='UNKNOWN_PARTIAL_COST')
    for f in (study/'audit').glob('data_recovery_*.json'):
        d=read(f);emit('interrupted_preparation',f.stem,None,{},scope='unknown interrupted local attempts beyond the earlier protocol correction, not zero',status='UNKNOWN_PARTIAL_COST',details=d.get('interrupted_local_inflight'))
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
