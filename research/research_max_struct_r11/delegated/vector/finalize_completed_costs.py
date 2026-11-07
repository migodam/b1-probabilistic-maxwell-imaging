"""Finalize allvector cost entries after CSV-only completed report production."""
import json,time,sys,subprocess,datetime
from pathlib import Path
ENTRY=time.perf_counter();V=Path(__file__).resolve().parent;ROOT=V.parents[1];O=V/'completed_analysis';O.mkdir(exist_ok=True)
p=subprocess.run([sys.executable,str(ROOT/'experiments/vector_compile_costs.py')],cwd=ROOT.parents[1],capture_output=True,text=True)
(O/'BASE_COST_COMPILER.log').write_text(p.stdout+p.stderr)
if p.returncode:raise RuntimeError('Base metadata costcompiler failed:'+p.stdout+p.stderr)
x=json.loads((V/'COSTS.json').read_text());runs=[json.loads(s) for s in (O/'ANALYSIS_RUNS.jsonl').read_text().splitlines()]
def row(ident,cpu,wall,source,scope):return dict(id=ident,platform='MAC_CPU',CPU_seconds=cpu,CPU_additive=True,CPU_measurement='producer_measured_own_scope_near_end',wall_seconds=wall,wall_additive=False,status='COMPLETE',source=source,scope=scope,physical_counts=[dict(matrix_order=5184,physical_LU=0,primal_complex_RHS=0,adjoint_complex_RHS=0)],notes=['Metadata only,noFW/GTarrays/newplots;not root theory/publication cost'])
x['entries'].append(row('completed_CSV_analysis_and_REPORT_all_runs',sum(r['CPU_seconds'] for r in runs),sum(r['wall_seconds'] for r in runs),'completed_analysis/ANALYSIS_RUNS.jsonl','Own CSVsummary/REPORTmd+HTML/COMPLETED_REPRODUCE generation;allrecordedruns once'))
record=dict(UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),own_CPU_seconds=time.process_time(),own_wall_seconds=time.perf_counter()-ENTRY,scope='Own finalcostwrapper only;child basecompiler already counted byCOST_COMPILER_RUNS,do not addchild again',new_forward_RHS=0,GT_array_reads=0)
with (O/'COST_FINALIZER_RUNS.jsonl').open('a') as f:f.write(json.dumps(record)+'\n')
fins=[json.loads(s) for s in (O/'COST_FINALIZER_RUNS.jsonl').read_text().splitlines()]
x['entries'].append(row('completed_cost_finalizer_own_all_runs',sum(r['own_CPU_seconds'] for r in fins),sum(r['own_wall_seconds'] for r in fins),'completed_analysis/COST_FINALIZER_RUNS.jsonl','Own wrapper only;child compiler already counted separately'))
byhost={}
for r in x['entries']:
 if r['CPU_additive'] and r['CPU_seconds'] is not None:byhost[r['platform']]=byhost.get(r['platform'],0)+r['CPU_seconds']
x.update(status='COMPLETED_VECTOR_ALL50_SCORED_AND_REPORTED_KNOWN_COSTS_WITH_EXPLICIT_NULL_GAPS',generated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),known_CPU_lower_bound_by_platform_seconds=byhost,completed_reporting_producers_included=True,all50_fresh_scoring_status='COMPLETE_ROOT_SEPARATELY_AUTHORIZED_AFTER_PUBLIC_BARRIER',no_further_scientific_writes=True)
for name in ('REPORT.md','REPORT.html'):
 p=V/name;s=p.read_text().replace('@MAC_CPU_COST@',f"{byhost['MAC_CPU']:.6f}").replace('@WINDOWS_CPU_COST@',f"{byhost['WINDOWS_CUDA']:.6f}");p.write_text(s)
(V/'COSTS.json').write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
(O/'COST_FINALIZATION_RECEIPT.json').write_text(json.dumps(dict(status='COMPLETE',platform_entry_sums=byhost,entries=len(x['entries']),physical_totals_by_matrix_order=x['physical_totals_by_matrix_order'],null_gaps=len(x['null_gaps']),scientific_RHS_added_by_reporting=0),indent=2)+'\n')
print(json.dumps(dict(status=x['status'],CPU_by_platform=byhost,entries=len(x['entries']),null_gaps=len(x['null_gaps']))),flush=True)
