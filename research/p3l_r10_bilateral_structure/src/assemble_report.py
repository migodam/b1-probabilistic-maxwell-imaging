"""Assemble documented results from paid receipts. No physics or model selection."""
from pathlib import Path
import csv, hashlib, json, datetime, subprocess
P=Path(__file__).resolve().parents[1]
R=P.parents[1]
def read(s): return json.loads((P/s).read_text())
def write(s,x):
    (P/s).parent.mkdir(parents=True,exist_ok=True)
    (P/s).write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
cost=read('audit/independent_qa/FINAL_COST_RECONCILIATION.json')
assert cost['actual_RHS256']==306096 and cost['actual_LU256']==9273
fresh=read('results/candidates_fresh10/STATISTICS.json')
transfer=read('results/transfer/FINAL_SUMMARY.json')
reduced=read('results/reduced_response_fresh10/RESULTS.json')
prob=read('results/reduced_response_probability_workload/RESULTS.json')
decisions={
 'FIXED_SHAPE_PHYSICS_REUSE':'SCOPED_GO',
 'BILATERAL_TARGET_COVERAGE':'NO_GO_EXPANSION',
 'SHAPE_TARGET_COVERAGE_ADVANTAGE':'NOT_ESTABLISHED',
 'FIXED_SHAPE_CONDITIONAL_TASK':'NUMERICAL_TASK_PASS',
 'MINIMAL_REALIZATION':'NOT_CERTIFIED',
 'FULL66_POSTERIOR':'UNQUALIFIED',
 'EXACT_MATERIAL_FIBER':'NOT_ESTABLISHED',
 'FULL_IMAGING_SYSTEM_ACCELERATION':'NOT_ESTABLISHED',
 'STRONG_IEEE_TAP_CONTRIBUTION':'HOLD',
 'NEXT_PRIORITY':'shape-update response reuse with residual-controlled refresh; bounded paired path experiment'}
write('EXPERIMENT_RESULTS.json',dict(status='SCIENTIFIC_WORK_COMPLETE',updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),source_base_commit='62fb7f0616849cd0f011e033fcf19acbcc6125f2',decisions=decisions,cost=cost,transfer=transfer,fresh10=dict(summary=fresh['summary'],scope='new material parents; small mechanism study; sameK not samecost across old/new arms'),heldout_shape_reduction=reduced,conditional_workload=prob,statistics='results/candidates_fresh10/CANONICAL_PAIRED_BOOTSTRAP.json',not_run=['nonlinear variable projection','Full66 posterior sampling','SDE','new neural models','global minimal realization','new acquisition/noise','Research W transfer','continuous PDE certification']))
write('RESULTS_SUMMARY.json',dict(status='COMPLETE_SCOPED',decisions=decisions,cost=cost,reference='EXPERIMENT_RESULTS.json'))
# Preserve old DEV table, expose a dataset-labelled combined view.
old=P/'METHOD_COMPARISON.csv'; saved=P/'results/candidates/METHOD_COMPARISON_DEV.csv'
if not saved.exists(): saved.write_bytes(old.read_bytes())
rows=[]
for label,path in [('DEV8_EXPOSED',saved),('FRESH10',P/'results/candidates_fresh10/METHOD_COMPARISON.csv')]:
    with path.open() as f:
        for row in csv.DictReader(f): rows.append(dict(dataset=label,**row))
keys=list(dict.fromkeys(k for row in rows for k in row))
with old.open('w') as f:
    w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
# Raw ledgers remain append-only; the master is a namespaced read-only assembly.
ledgers=['results/transfer/COST_LEDGER.jsonl','results/candidates/RESOURCE_LEDGER.jsonl','results/candidates_fresh10/RESOURCE_LEDGER.jsonl','results/reduced_response/COST_LEDGER.jsonl','results/reduced_response_fresh10/COST_LEDGER.jsonl','results/reduced_response_probability_workload/COST_LEDGER.jsonl']
with (P/'COST_LEDGER.jsonl').open('w') as out:
    out.write(json.dumps(dict(kind='INDEX',interpretation='Source-local cumulative snapshots are not additive. Use final receipt reconciliation; every paid action retained once per source.',canonical_totals=cost))+'\n')
    for source in ledgers:
        path=P/source
        if not path.exists(): raise FileNotFoundError(path)
        for i,line in enumerate(path.open(),1):
            if line.strip():out.write(json.dumps(dict(source=source,line=i,record=json.loads(line)),ensure_ascii=False,allow_nan=False)+'\n')
files=[]
for sub in ('src','config','theory','audit'):
 for p in sorted((P/sub).rglob('*')):
    if p.is_file() and p.suffix in ('.py','.json','.md','.txt','.csv') and '__pycache__' not in p.parts:
      files.append(dict(path=str(p.relative_to(R)),bytes=p.stat().st_size,blake2b_256=hashlib.blake2b(p.read_bytes(),digest_size=32).hexdigest()))
write('SOURCE_MANIFEST.json',dict(id='R10_BILATERAL_SOURCE_IDENTITY',base_commit='62fb7f0616849cd0f011e033fcf19acbcc6125f2',R9_scientific_commit='7a80e7ac44d723eab4cbab54f028d258780e12e7',R9_release='p3lr9-stabilized-volume-20261030',runtime='Python3.13 numpy2.5.3 scipy1.18.1; float64/complex128; singleBLAS',files=files,dependency_manifest='Release R10_dependency_inputs.tar.zst',historical_assets='read-only; final publication commit identifies report/tree; this manifest does not hash itself'))
print(json.dumps(dict(total_rows=len(rows),source_files=len(files),cost=cost)))
