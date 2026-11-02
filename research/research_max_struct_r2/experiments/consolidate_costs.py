"""Aggregate existing measured receipts once per outer process/task; no physics."""
from pathlib import Path
from collections import Counter
import json,csv

ROOT=Path(__file__).resolve().parents[1]
def read(path):return json.loads((ROOT/path).read_text())

def main():
    rows=[];details=[]
    def add(name,path,j,wall='wall_seconds',cpu='cpu_seconds',counts=None,phase='validation',classification='AUTHORIZED_R2'):
        c=counts or {};rows.append(dict(run_id=name,source=path,phase=phase,classification=classification,
            wall_seconds=j.get(wall),cpu_seconds=j.get(cpu),physical_LU=c.get('LU',0),source_RHS=c.get('source_RHS',c.get('forward_RHS',0)),
            adjoint_RHS=c.get('adjoint_RHS',0),receiver_RHS=c.get('receiver_RHS',0),tangent_RHS=c.get('tangent_RHS',0),
            actual_complex_RHS=c.get('actual_complex_RHS',c.get('forward_RHS',0)),current_Cholesky=c.get('current_Cholesky',0),
            current_RHS=c.get('current_RHS',0),scope=j.get('cost_scope',j.get('timer_scope','See source outer-process receipt')),
            missing_cpu_reason='Not recorded in raw receipt'if j.get(cpu)is None else None))
    gen=read('audit/GENERATION_RECEIPT.json');add('BLIND_GENERATION','audit/GENERATION_RECEIPT.json',gen,counts=gen['physical_counts'],phase='generation')
    gradient=read('audit/PROFILE_GRADIENT_WIRING.json');cc=Counter(gradient['physical_model_cost']);cc.update(gradient['result']['cost'])
    add('PROFILE_GRADIENT_WIRING','audit/PROFILE_GRADIENT_WIRING.json',gradient,counts=cc,phase='implementation_validation')
    for path in sorted((ROOT/'results/blind').glob('*/*/RECEIPT.json')):
        j=json.loads(path.read_text());counts=Counter(j['physical_model_cost'])
        if j['arm']=='profile':
            costs=j['auxiliary']['cost'];counts.update({k:costs.get(k,0)for k in ('current_Cholesky','current_RHS')})
        add(j['run_id'],str(path.relative_to(ROOT)),j,counts=counts,phase='inverse',classification='FRESH_BLIND_SEALED')
        details.append(dict(run_id=j['run_id'],arm=j['arm'],auxiliary_stop=j['auxiliary'].get('stop',j['auxiliary'].get('status')),
            corrector_stop=j['corrector_result']['stop'],current_steps=j['auxiliary'].get('cost',{}).get('PR_CG_steps'),
            nested_auxiliary_cost=j['auxiliary'].get('cost'),nested_corrector_cost=j['corrector_result'].get('cost'),
            nested_timers_rule='Already included in task wall/CPU and physical_model_cost; never add again'))
    for stage,record in [('plan','PLAN'),('acquire','ACQUISITION'),('fit','FIT'),('score','SCORE'),('coverage','COVERAGE')]:
        path=f'results/active/{stage}_PROCESS.json';j=read(path);r=read(f'results/active/{record}.json')
        add('ACTIVE_'+stage,path,j,wall='wall_s',cpu='cpu_s',counts={'LU':r['LU'],'forward_RHS':r['forward_RHS']},
            phase='active_'+stage,classification='RESTRICTED_1D_FAMILY_FINITE_DDA')
    for name in ['ENDPOINTS','SOURCE_WIRING']:
        path=f'results/diagnostics/{name}.json'
        if (ROOT/path).exists():
            j=read(path);add(name,path,j,counts=j['physical_model_cost'],phase='post_prediction_no_GT_validation')
    if (ROOT/'results/diagnostics/SUMMARY.json').exists():
        j=read('results/diagnostics/SUMMARY.json').get('summary_assembly')
        if j:
            add('DIAGNOSTIC_SUMMARY_ASSEMBLY','results/diagnostics/SUMMARY.json#summary_assembly',j,
                counts={'LU':0,'forward_RHS':0},phase='scoring_rendering')
    for path,name in [('audit/PREDICTIONS_SEALED.json','PREDICTION_SEAL'),('results/BLIND_SCORES.json','POST_CLOSURE_SCORING'),
                      ('audit/BLIND_PLOT_RECEIPT.json','BLIND_PLOT')]:
        j=read(path);add(name,path,j,phase='scoring_rendering',counts={'LU':0,'forward_RHS':0})
    if (ROOT/'audit/ACTIVE_PLOT.json').exists():
        j=read('audit/ACTIVE_PLOT.json');add('ACTIVE_PLOT','audit/ACTIVE_PLOT.json',j,wall='wall_s',cpu='cpu_s',phase='scoring_rendering',counts={'LU':0,'forward_RHS':0})
    # Post-hoc runs, if later sealed, are a separate development component.
    posthoc_paths=list((ROOT/'results/posthoc').glob('**/RECEIPT.json'))+list((ROOT/'results/capacity_control').glob('**/RECEIPT.json'))
    for path in sorted(posthoc_paths):
        j=json.loads(path.read_text())
        if j.get('status')not in ('SEALED_PREDICTION','SEALED_POST_HOC_PREDICTION','FAILED'):continue
        add(j.get('run_id','POST_HOC_CAPACITY'),str(path.relative_to(ROOT)),j,counts=j.get('physical_model_cost',{}),
            phase='posthoc_development',classification='POST_HOC_DEVELOPMENT_NOT_NEW_BLIND')
    if (ROOT/'results/POSTHOC_CAPACITY_SCORE.json').exists():
        j=read('results/POSTHOC_CAPACITY_SCORE.json');add('POSTHOC_CAPACITY_SCORING','results/POSTHOC_CAPACITY_SCORE.json',j,
            phase='scoring_rendering',classification='POST_HOC_DEVELOPMENT_NOT_NEW_BLIND',counts={'LU':0,'forward_RHS':0})
    count_names=['physical_LU','source_RHS','adjoint_RHS','receiver_RHS','tangent_RHS','actual_complex_RHS','current_Cholesky','current_RHS']
    totals={k:sum(x[k]for x in rows)for k in count_names}
    timed=[x for x in rows if x['wall_seconds']is not None];cpurows=[x for x in rows if x['cpu_seconds']is not None]
    numerics=[x for x in rows if x['phase']not in ('scoring_rendering',)]
    out=dict(status='MEASURED_RECEIPTS_WITH_EXPLICIT_UNMEASURED_OVERHEAD',rows=rows,counts=totals,
       measured_wall_process_sum=sum(x['wall_seconds']for x in timed),measured_cpu_known_component_sum=sum(x['cpu_seconds']for x in cpurows),
       complete_cpu_process_sum=None,complete_cpu_reason='ACTIVE_PLOT CPU and other uninstrumented overhead unavailable; known component sum is not complete total',
       numerical_wall_process_sum=sum(x['wall_seconds']for x in numerics if x['wall_seconds']is not None),
       numerical_cpu_process_sum=sum(x['cpu_seconds']for x in numerics if x['cpu_seconds']is not None),
       nested_inverse_details=details,
       accounting='Outer inverse task, generation, gradient-check, active stage,diagnostic and scoring/rendering timers each added once; nested current/corrector/setup/raw state timers not added again. Process sums are not elapsed research time.',
       physical_dimensions='Blind/generation/endpoints use576; source-wiring/gradient root use16; active coarse DDA192 and fine DDA1536. Aggregate LU/RHS are counted operations, not equal FLOPs or speed proof.',
       source_metric_Cholesky=None,source_metric_Cholesky_reason='Raw global/source-wiring receipts retain these small factorizations, not comprehensively merged into this physical/current cost table',
       small_SVD=None,small_SVD_reason='Raw basis/rank SVD/QR counters retained, but norm/least-squares/legality LAPACK internals not comprehensively instrumented',
       Green_actions=None,Green_actions_reason='Raw CIE/profile operator column counters retained; complete Green/Gram/regularizer actions not exhaustively counted',
       FFT=0,elliptic_solves=0,GPU_seconds=0,monetary_cost=None,monetary_cost_reason='No billed-service meter/local hardware amortization price available',
       unmeasured=['interpreter imports/startup','shell/tool orchestration','scientific reasoning/literature/API work','R1 read-only audit','source preparation/report composition',
                   'some receipt/scoring JSON/CSV serialization','historical operator assembly','packaging/publication unless separately receipted'],
       failed_attempts='No failed numerical inverse dispositions in12 sealed main receipts; successful optimizer MAXITER/early stops are not erased. Any other failed tooling attempts remain outside numerical receipts and unmeasured.',
       no_hash_checks=True)
    (ROOT/'audit/COST_TOTALS.json').write_text(json.dumps(out,indent=2,allow_nan=False))
    with (ROOT/'audit/COST_LEDGER.csv').open('w',newline='')as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    lines=['# R2 paid cost ledger','',
      f"Current numerical receipts: **{totals['physical_LU']:,} physical LU / {totals['actual_complex_RHS']:,} complex RHS**; **{totals['current_Cholesky']:,} current Cholesky / {totals['current_RHS']:,} current RHS** are separate.",'',
      f"Numerical process sums: {out['numerical_wall_process_sum']:.6f}s wall / {out['numerical_cpu_process_sum']:.6f}s CPU. Measured rendering/scoring-inclusive wall sum: {out['measured_wall_process_sum']:.6f}s. Known CPU component sum: {out['measured_cpu_known_component_sum']:.6f}s; full CPU total is unknown because some components were not timed.",'',
      '|Work|Wall/s|CPU/s|Physical LU|Source RHS|Adjoint RHS|All complex RHS|Current Cholesky|Current RHS|','|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    for x in rows:
        ws='null'if x['wall_seconds']is None else f"{x['wall_seconds']:.6f}";cs='null'if x['cpu_seconds']is None else f"{x['cpu_seconds']:.6f}"
        lines.append(f"|{x['run_id']}|{ws}|{cs}|{x['physical_LU']}|{x['source_RHS']}|{x['adjoint_RHS']}|{x['actual_complex_RHS']}|{x['current_Cholesky']}|{x['current_RHS']}|")
    lines+=['','Counts sum measured work across different matrix dimensions and decompositions. They are not a FLOP-equivalence certificate. Profile current Cholesky matrices are576×576; the gradient wiring check uses16×16. They are not free, nor are they additional physical Maxwell forward LU.',
      '',out['accounting'],'','No online post-closure reselection changes the12 main predictions. Post-hoc controls, if included, are labelled development and excluded from blind success rates.',
      '','Unknown overhead, incomplete phase partitions and monetary fees remain null. Existing R1 costs are referenced historical work and never added as new R2 executions. See [structured totals](audit/COST_TOTALS.json) and [CSV ledger](audit/COST_LEDGER.csv).']
    (ROOT/'COST_LEDGER.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps(dict(counts=totals,numerical_wall=out['numerical_wall_process_sum'],wall_measured=out['measured_wall_process_sum'])))

if __name__=='__main__':main()
