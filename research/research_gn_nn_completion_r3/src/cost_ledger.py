"""Collect saved measured receipts. Rows explicitly distinguish reused costs."""
from pathlib import Path
import argparse,csv,json

def main():
    p=argparse.ArgumentParser();p.add_argument('--study',type=Path,default=Path(__file__).resolve().parents[1]);a=p.parse_args();s=a.study
    rows=[]
    def add(phase,path,obj,wall,counts=None,scope='actual_unique_execution',status=None):
        counts=counts or{}
        rows.append(dict(phase=phase,id=str(path.relative_to(s)),scope=scope,status=status or obj.get('status','RECORDED'),wall_seconds=wall,
            LU=counts.get('LU'),complex_RHS=counts.get('actual_complex_RHS'),FullJ=counts.get('FullJ',counts.get('FULL_ANALYTIC_J')),
            counts_json=json.dumps(counts,separators=(',',':')),time_note='task elapsed; concurrent task times must not be summed as overall elapsed'))
    read=lambda f:json.loads(f.read_text())
    for f in [s/'audit/GENERATION_RECEIPT_INITIAL_STRESS_DESIGN.json',s/'audit/GENERATION_RECEIPT_BALANCED_STRESS_DESIGN.json',*sorted((s/'data').glob('GENERATION_EXPAND*RECEIPT.json'))]:
        if f.exists():
            d=read(f);add('fullwave_data_generation',f,d,d.get('wall_seconds'),d.get('counts'))
    for f in sorted((s/'anchors').glob('obs*.json')):
        d=read(f)
        wall=d.get('anchor_plus_feature_wall_seconds')
        if wall is None and d.get('solve_wall_seconds') is not None:
            wall=d['solve_wall_seconds']+d.get('feature_wall_seconds',0.)
        add('weakreg45_anchor'if '_weakreg'in f.stem else'GN_anchor',f,d,wall,d.get('counts'))
    for f in sorted((s/'anchors').glob('BATCH*.json')):
        d=read(f);add('anchor_batch_elapsed',f,d,d.get('wall_seconds'),d.get('aggregate_model_counts'),scope='aggregate_of_anchor_rows_do_not_double_count')
    for f in sorted((s/'audit').glob('ANCHOR_FAILURE_INITIAL*.json')):
        d=read(f)
        add('failed_anchor_initial_attempt',f,d,d.get('failed_wall_seconds'),scope='additional_failed_execution_operation_counts_unknown',status='FAILED_COUNTS_UNKNOWN')
    for f in sorted((s/'models').glob('**/summary.json')):
        d=read(f);add('NN_training',f,d,d.get('training_seconds'),scope='actual_unique_training_MPS_LU_RHS_not_applicable')
    for f in sorted((s/'predictions').glob('*/inference_setup.json')):
        d=read(f);reused=(f.parent/'REUSED_PROPOSALS.json').exists()or'_filter_'in f.parent.name
        add('NN_inference_cold_setup',f,d,d.get('checkpoint_and_operator_setup_seconds'),scope='reused_source_setup_no_new_execution_do_not_double_count'if reused else'shared_once_per_run_not_in_warm_arm_time')
    for f in sorted((s/'evaluation').glob('*/obs*.json')):
        d=read(f)
        if d.get('schema')!='gn-nn-completion-r3-evaluation-v1':continue
        scope='warm_arm_includes_incremental_feature_and_inference'
        if f.parent.name=='weakreg45_anchor':scope='reuses_weakreg45_anchor_plus_new_verification_do_not_double_count_anchor'
        add('online_evaluation',f,d,d.get('wall_seconds'),d.get('costcounts'),scope=scope)
    for f in sorted((s/'evaluation').glob('RUN_MANIFEST*.json')):
        d=read(f);add('evaluation_batch_elapsed',f,d,d.get('batch_wall_seconds'),scope='aggregate_elapsed_includes_IO_do_not_add_to_case_times')
    f=s/'cache/reference_jacobian.json'
    if f.exists():
        d=read(f);add('local_shared_reference',f,d,d.get('wall_seconds'),d.get('counts'))
    # Oracle costs are offline diagnostics, never an online inference advantage.
    seen_oracles=set()
    for f in [s/'evaluation/DEV_INTEGRATION_SCORES.json',
              s/'evaluation/PILOT256_VAL_SCORES.json',
              s/'evaluation/CONFIRMATION_SCORES.json']:
        if not f.exists():continue
        for d in read(f).get('oracle_rows',[]):
            key=d.get('obs_id')
            if key in seen_oracles:continue
            seen_oracles.add(key)
            add('offline_GT_current_oracle',f,d,d.get('wall_seconds'),d.get('costcounts'),
                scope='offline_diagnostic_not_online_cost')
            rows[-1]['id']+='Object:'+str(key)
    out=s/'results';out.mkdir(exist_ok=True)
    with (out/'COST_LEDGER.csv').open('w',newline='')as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (out/'COST_ACCOUNTING_NOTES.json').write_text(json.dumps(dict(rows=len(rows),aggregate_warning='Anchor batch and evaluation batch rows are elapsed summaries, not additional computations. Weakreg45 evaluation reuses its anchor. Do not sum them twice.',unknown_costs='Initial failed TRAIN SVD attempt has recorded child wall but unknown LU/RHS. Exact all-attempt total operation count is unavailable.',precision='All final physics f64/c128; NN training f32 real-pair complex64',hardware_scope='CPU anchor batches span i9 and M4; post-GN arms on M4. No same-hardware end-to-end acceleration claim.',remote_reference='Remote reference setup is separately retained in remote batch receipt; charge once when cold.'),indent=2))
    print(len(rows))

if __name__=='__main__':main()
