"""Read sealed development/diagnostic counters. No physical computation."""
from pathlib import Path
from collections import Counter
import json,csv
ROOT=Path(__file__).resolve().parents[1]

def main():
    rows=[];counts=Counter()
    for p in sorted((ROOT/'results/development').glob('*/RECEIPT.json')):
        j=json.loads(p.read_text());c=j.get('physical_model_cost',j.get('cost',{}));counts.update(c)
        rows.append(dict(run=p.parent.name,source=str(p.relative_to(ROOT)),classification=j.get('classification'),status=j['status'],
          LU=c.get('LU'),complex_RHS=c.get('actual_complex_RHS'),wall_seconds=j.get('wall_seconds'),cpu_seconds=j.get('cpu_seconds')))
    aux=Counter()
    for p in [ROOT/'results/diagnostics/FIXED_MATERIAL_H1_H2.json',ROOT/'audit/PASSIVE_GRADIENT_WIRING.json']:
        if p.exists():
            j=json.loads(p.read_text());c=j['physical_model_cost'];counts.update(c);aux.update(j.get('auxiliary_cost',{}))
            rows.append(dict(run=p.stem,source=str(p.relative_to(ROOT)),classification=j.get('classification'),status=j['status'],
              LU=c.get('LU'),complex_RHS=c.get('actual_complex_RHS'),wall_seconds=j.get('wall_seconds'),cpu_seconds=j.get('cpu_seconds')))
    for p in [ROOT/'audit/INHERITED_READ_COVERAGE.json',ROOT/'results/DEVELOPMENT_SCORES.json']:
        if p.exists():
            j=json.loads(p.read_text());rows.append(dict(run=p.stem,source=str(p.relative_to(ROOT)),classification='READ_ONLY_AUDIT_SCORING',status='COMPLETED',
                LU=0,complex_RHS=0,wall_seconds=j.get('wall_seconds'),cpu_seconds=j.get('cpu_seconds')))
    sourcegram=ROOT/'results/development/globaldual_restart/SMALL_GRAM_WORK.json'
    gram=json.loads(sourcegram.read_text())if sourcegram.exists()else None
    result=dict(status='DEVELOPMENT_ACCOUNTING_SNAPSHOT_NO_FRESH_CONFIRMATION',rows=rows,
      physical_LU=counts.get('LU'),physical_complex_RHS=counts.get('actual_complex_RHS'),
      source_RHS=counts.get('source_RHS'),adjoint_RHS=counts.get('adjoint_RHS'),receiver_transpose_RHS=counts.get('receiver_transpose_RHS'),
      auxiliary_counts=dict(aux),globaldual_small_Gram=gram,
      wall_process_sum=sum(x['wall_seconds']for x in rows if x['wall_seconds']is not None),
      cpu_known_component_sum=sum(x['cpu_seconds']for x in rows if x['cpu_seconds']is not None),
      accounting='Each outer receipt once; nested stage/corrector/raw solver/current Gram timing already included. Published R2 reference not charged.',
      matrix_dimensions='Development and H1/H2 physical systems576; passive-gradient root16. Current/Gram systems have separate dimensions/counters; counts are not equal FLOPs.',
      fresh_generation='NOT_RUN',fresh_inverse='NOT_RUN',fresh_scoring='NOT_RUN',
      unmeasured=['imports/startup','tool/reasoning/literature','source/contract/report writing','some receipt serialization','historical assembly'],
      monetary_cost=None,monetary_cost_reason='No verified price/meter',hashes=None,hashes_reason='User forbids hash checks')
    (ROOT/'audit/DEVELOPMENT_COST_AUDIT.json').write_text(json.dumps(result,indent=2))
    with(ROOT/'audit/DEVELOPMENT_COST_LEDGER.csv').open('w',newline='')as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    print(json.dumps({k:result[k]for k in ['physical_LU','physical_complex_RHS','wall_process_sum','fresh_generation']}))

if __name__=='__main__':main()
