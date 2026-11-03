"""Arithmetic-only final R3 receipt consolidation; old snapshot stays unchanged."""
from pathlib import Path
from collections import Counter
from datetime import datetime, timezone
import csv, json, time
ROOT=Path(__file__).resolve().parents[1]

def main():
    tick,cpu=time.perf_counter(),time.process_time()
    rows=[];total=Counter();auxiliary={};seen=set()
    physical=list(sorted((ROOT/'results/development').glob('*/RECEIPT.json')))
    physical += [ROOT/p for p in ['results/diagnostics/FIXED_MATERIAL_H1_H2.json','audit/PASSIVE_GRADIENT_WIRING.json',
        'results/diagnostics/NOISE_FLOOR_DISTRIBUTION.json','results/fiber/RECEIPT.json',
        'results/onecell/RECEIPT.json','results/vector_scope/RECEIPT.json']]
    for p in physical:
        j=json.loads(p.read_text());assert j['status'] in ('SEALED_DEVELOPMENT_PREDICTION','COMPLETED')
        c=j.get('physical_model_cost',j.get('physical_cost',j.get('cost',{})))
        assert 'LU' in c and 'actual_complex_RHS' in c, p
        total.update({k:v for k,v in c.items() if isinstance(v,(int,float))})
        source=str(p.relative_to(ROOT));auxiliary[source]=j.get('auxiliary_cost',{})
        rows.append(dict(run=p.parent.name if p.name=='RECEIPT.json' else p.stem,source=source,
            classification=j.get('classification'),status=j['status'],LU=c['LU'],complex_RHS=c['actual_complex_RHS'],
            wall_seconds=j.get('wall_seconds'),cpu_seconds=j.get('cpu_seconds'),accounting_group='PHYSICAL_OUTER_RECEIPT'))
    reads=[ROOT/'audit'/s for s in ['INHERITED_READ_COVERAGE.json','R2_HISTORY_READBACK.json',
        'FINITE_SCAN_READBACK.json','ALL_PAIR_SCORE_READBACK.json','FINAL_MECHANISM_ARITHMETIC_READBACK.json','FINAL_ARITHMETIC_READBACK.json','RECONSTRUCTION_FIGURE_RECEIPT.json']]
    scores=list(sorted((ROOT/'results/development_scoring_passes').glob('*.json')))+[ROOT/'results/DEVELOPMENT_SCORES.json']
    for p in reads+scores:
        if not p.exists():continue
        j=json.loads(p.read_text());identity=j.get('scored_utc',str(p))
        if identity in seen:continue
        seen.add(identity)
        if j.get('auxiliary_cost'):auxiliary[str(p.relative_to(ROOT))]=j['auxiliary_cost']
        rows.append(dict(run=p.stem,source=str(p.relative_to(ROOT)),classification=j.get('classification','READ_ONLY_AUDIT_SCORING'),
            status=j.get('status','COMPLETED'),LU=0,complex_RHS=0,wall_seconds=j.get('wall_seconds'),
            cpu_seconds=j.get('cpu_seconds'),accounting_group='READ_ONLY_OUTER_RECEIPT'))
    for p in [ROOT/'audit/ALL_PAIR_STARTUP_FAILURE.json',ROOT/'audit/FINAL_SCORE_STARTUP_FAILURE.json',ROOT/'audit/FINAL_ARITHMETIC_INITIAL_FAILURE.json']:
        if not p.exists():continue
        j=json.loads(p.read_text());rows.append(dict(run=p.stem,source=str(p.relative_to(ROOT)),
            classification=j['classification'],status=j.get('status','FAILED_BEFORE_NUMERICS'),LU=0,complex_RHS=0,
            wall_seconds=j.get('wall_seconds'),cpu_seconds=j.get('cpu_seconds'),accounting_group='FAILED_ATTEMPT_UNMEASURED_TIME'))
    sourcegram=ROOT/'results/development/globaldual_restart/SMALL_GRAM_WORK.json'
    out=dict(status='FINAL_DEVELOPMENT_ACCOUNTING_NO_FRESH_CONFIRMATION',created_utc=datetime.now(timezone.utc).isoformat(),
        physical_outer_receipt_count=len(physical),rows=rows,physical_LU=total['LU'],physical_complex_RHS=total['actual_complex_RHS'],
        RHS_components={k:total.get(k,0) for k in ('source_RHS','adjoint_RHS','receiver_transpose_RHS','resolvent_RHS')},
        physical_counters=dict(total),auxiliary_by_receipt=auxiliary,globaldual_small_Gram=json.loads(sourcegram.read_text()),
        wall_process_sum=sum(x['wall_seconds'] for x in rows if x['wall_seconds'] is not None),
        cpu_known_component_sum=sum(x['cpu_seconds'] for x in rows if x['cpu_seconds'] is not None),
        accounting='Each outer physical receipt once. Nested auxiliary model/stage/corrector counters are retained as provenance but never added twice. Archived scoring passes deduplicated by scored_utc. Wall sum is sum of instrumented processes, not elapsed campaign duration.',
        snapshot_preserved='audit/DEVELOPMENT_COST_AUDIT.json is historical six-arm snapshot, not final total.',
        matrix_dimensions='Scalar development/noise/H1H2 576; passive wiring 16; fiber 6; onecell 1; vector scope 24. LU/RHS totals across unequal dimensions are workload counts, not FLOPs or fair speed comparison.',
        fresh_generation='NOT_RUN',fresh_inverse='NOT_RUN',fresh_scoring='NOT_RUN',
        inherited_R2_new_cost=0,inherited_R2_reason='Published reference reused; no historical rerun',
        unmeasured=['startup/imports','pre-instrumentation failures and interrupted arithmetic reviews','tool/reasoning/literature','source/contract/report writing','historical assembly','some serialization','this consolidation output write'],
        monetary_cost=None,monetary_cost_reason='No verified price or billing meter',hash_checks=False,
        consolidation_wall_seconds=time.perf_counter()-tick,consolidation_cpu_seconds=time.process_time()-cpu)
    assert sum(out['RHS_components'].values())==out['physical_complex_RHS'],out['RHS_components']
    (ROOT/'audit/FINAL_DEVELOPMENT_COST_AUDIT.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    with (ROOT/'audit/FINAL_DEVELOPMENT_COST_LEDGER.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    print(json.dumps({k:out[k] for k in ['physical_outer_receipt_count','physical_LU','physical_complex_RHS','RHS_components','wall_process_sum','cpu_known_component_sum']}))
if __name__=='__main__':main()
