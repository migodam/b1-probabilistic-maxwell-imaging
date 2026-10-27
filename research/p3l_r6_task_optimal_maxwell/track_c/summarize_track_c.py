"""Auditable Track C artifact/cost closeout; no truth or physical calls."""
import json
from pathlib import Path
TRACK=Path(__file__).resolve().parent


def main():
    read=lambda path:json.loads(path.read_text())
    original=read(TRACK/'evaluation_CDEV/CDEV_RESULTS.json')
    revision=read(TRACK/'evaluation_CDEV_revision1/CDEV_RESULTS.json')
    native=read(TRACK/'evaluation_matched512_full_native/MATCHED_CDEV_RESULTS.json')
    low=read(TRACK/'evaluation_matched512_full_low/MATCHED_CDEV_RESULTS.json')
    ledger=[json.loads(line) for line in (TRACK/'MATCHED_TRAINING_LEDGER.jsonl').read_text().splitlines()]
    ends=[row for row in ledger if row['event']=='END']
    active={row['job'] for row in ledger if row['event']=='BEGIN'}-{row['job'] for row in ends}
    reports={}
    for family,root,configset in [
        ('CDEV_initial5',TRACK/'contexts/CDEV',{'native96_probit64','native24_probit16','native48_probit32','native72_probit48','adaptive_native72_probit48'}),
        ('CDEV_revision1_3',TRACK/'contexts/CDEV',{'native24_probit64','native48_probit64','native96_probit32'}),
        ('historical_TRAIN512_low48P32',TRACK/'matched_anchors_TRAIN512',{'native48_probit32'})]:
        selected=[read(path) for path in root.glob('*/*/REPORT.json') if path.parent.name in configset]
        reports[family]=dict(parents_config_paths=len(selected),new_RHS=sum(r['path_resource_cost']['RHS'] for r in selected),
            states=sum(r['path_resource_cost']['states'] for r in selected),LU=sum(r['path_resource_cost']['LU'] for r in selected),
            physics_wall_seconds=sum(r['path_resource_cost']['physics_wall_seconds'] for r in selected),
            end_to_end_anchor_wall_seconds=sum(r['end_to_end_anchor_wall_seconds'] for r in selected))
    fits={}
    for arm in ('full','native','low'):
        records=[read(path) for path in (TRACK/'matched_heads_TRAIN512').glob(f'{arm}_s*/RESULT.json')]
        fits[arm]=dict(fits=len(records),capacity=6418,updates_per_fit=4000,train_n=512,provider_updates=0,
            seed_selected_updates={str(r['seed']):r['selected_update'] for r in records},
            recorded_training_wall_seconds=sum(r['training_wall_seconds'] for r in records),
            recorded_training_CPU_seconds=sum(r['training_CPU_seconds'] for r in records))
    benchmark=read(TRACK/'clean_inference_benchmark30/BENCHMARK.json')
    failed=read(TRACK/'evaluation_CDEV_failed_feature_dtype/FAILURE_RECEIPT.json')
    result=dict(status='TRACK_C_DEVELOPMENT_COMPLETE_LOCKED_D_OWNED_BY_ROOT',new_RHS_from_this_worker=0,
        root_paid_anchor_path_costs=reports,CDEV_truth_generation_RHS_separate=64*8,
        original_five_arm_gate=original['status'],revision1_gate=revision['status'],revision1_selection=revision['selection'],
        matched_native_risk_gates=[{k:r[k] for k in ('system','upper98_333_vs_newFull512','upper98_333_vs_originalFull1024','development_risk_candidate')} for r in native['systems']],
        matched_low_risk_gates=[{k:r[k] for k in ('system','upper98_333_vs_newFull512','upper98_333_vs_originalFull1024','development_risk_candidate')} for r in low['systems']],
        matched_training=dict(total_fits=len(ends),active_jobs=list(active),fits_by_arm=fits,
            GPU_reserved_whole_job_seconds=sum(r['GPU_reserved_whole_job_seconds'] for r in ends),
            CPU_seconds=sum(r['CPU_seconds'] for r in ends),assigned_GPU_wall_cap_seconds=1800,new_RHS=0,
            ledger=str(TRACK/'MATCHED_TRAINING_LEDGER.jsonl')),
        teacher_export=read(TRACK/'matched_inputs_TRAIN512/EXPORT_RECEIPT.json'),
        native_cache_reuse={k:v for k,v in read(TRACK/'native96_only_contexts/BUILD_RECEIPT.json').items() if k!='records'},
        clean_benchmark={key:benchmark[key] for key in ('status','repetitions','parents','seeds','arms','job_wall_seconds','job_CPU_seconds','scope')},
        retained_failures=[read(TRACK/'EVALUATOR_LAUNCH_RECEIPT.json'),failed,read(TRACK/'PHASE_CORRECTION_RECEIPT.json')],
        original_bad_phase_outputs_preserved=str(TRACK/'BUG_GT_PHASE_DEFINITION'),
        cost_limits='Original frozen training/provider/R4 data physics amortization is separate; offline scoring/plot orchestration CPU wall not separately recorded; no complete campaign cost total inferred here.',
        scientific_limits=['CDEV64 selected variants and checkpoints; no independent confirmation claim',
            'All new Full512/Native512/Low512 heads fail originalFull1024 1% noninferiority; no threshold relaxation',
            'Native cache reuse has zero incremental RHS, while Original native cold phase cost remains retained',
            'Head microbenchmark excludes anchor computation; standalone D physical timing belongs to root',
            'All phase figures use corrected Original full-anchor q0 masks; first erroneous GT-phase outputs retained'],
        D_GT_read_by_worker=False,figures=str(TRACK/'figures_CDEV_revision2/FIGURE_SOURCES.json'))
    (TRACK/'TRACK_C_CLOSEOUT.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status=result['status'],paid_anchors=reports,training=result['matched_training'])),flush=True)


if __name__=='__main__':main()
