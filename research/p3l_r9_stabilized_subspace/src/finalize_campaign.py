"""Synthesize completed, frozen evidence only. Never call physics or modify predictions."""
from pathlib import Path
import json,time
P=Path(__file__).resolve().parents[1];V=P/'results/possibility_volume'
def read(x):return json.loads((P/x).read_text())
def write(x,v):(P/x).write_text(json.dumps(v,ensure_ascii=False,indent=2)+'\n')
assert (V/'offline_statistics/SCORING_COMPLETE.json').exists()
basepath=P/'audit/RESOURCE_SUMMARY_BEFORE_E5.json'
if not basepath.exists():basepath.write_text((P/'RESOURCE_SUMMARY.json').read_text())
r=json.loads(basepath.read_text());e=read('results/possibility_volume/RESOURCE_SNAPSHOT.json')
# Import every original event only once. Source journals remain authoritative.
ledger=P/'COST_LEDGER.jsonl';rel='results/possibility_volume/RESOURCE_LEDGER.jsonl';seen=set()
if ledger.exists():
 with ledger.open() as f:
  for line in f:
   a=json.loads(line)
   if a.get('source')==rel:seen.add(a['source_line'])
peak=0
with ledger.open('a') as out,(P/rel).open() as src:
 for n,line in enumerate(src,1):
  a=json.loads(line);peak=max(peak,a.get('peak_RSS_bytes',0))
  if n not in seen:out.write(json.dumps(dict(source=rel,source_line=n,branch='possibility_volume',entry=a),ensure_ascii=False)+'\n')
r.update(reduced_RHS256=r['reduced_RHS256']+e['reduced_RHS256'],reduced_LU256=r['reduced_LU256']+e['reduced_LU256'],new_physics_after_formula_steering=dict(full_RHS576=0,reduced_RHS256=e['reduced_RHS256']),CPU_seconds_sum=r['CPU_seconds_sum']+e['CPU_seconds'],job_wall_seconds_sum=r['job_wall_seconds_sum']+e['controlled_wall_seconds'],E5_peak_RSS_bytes=peak,E5_cap_RHS256=400000,E5_cap_LU256=12000,E5_cap_wall_seconds=7200,E5=e,scientific_execution_complete=True)
r['branches'].append(dict(branch='possibility_volume',full_RHS=0,full_LU=0,reduced_RHS256=e['reduced_RHS256'],reduced_LU256=e['reduced_LU256'],CPU_seconds=e['CPU_seconds'],controlled_wall_seconds=e['controlled_wall_seconds'],source=rel))
r['scope']='Scientific controlled jobs only; documentation, reading, offline statistics and packaging reported separately if timed; full576/reduced256 counts are not interchangeable. All rejected physical attempts included.'
write('RESOURCE_SUMMARY.json',r)
s=read('POSSIBILITY_SUMMARY_COMPACT.json');s.update(status='COMPLETE_WITH_SCOPED_POSITIVE_AND_NEGATIVE_EVIDENCE',scientific_frozen_source_commit='75184170',published_scientific_base='9cd323cbc88af89c5dcbdef8b39049218f8dfcdf',go=dict(material_valid_candidate_cloud='SCOPED_GO',high_probability_GT_coverage='NOT_ESTABLISHED',fixed_subspace_20pct_advantage='NO_GO',dynamic_necessary='NOT_ESTABLISHED_IN_THIS_RANGE',exact_ROI_reuse='SCOPED_KERNEL_GO_STANDARD_IDENTITY',conditional_probability='NUMERICAL_TASK_PASS_SCOPED',strict_probability_certificate='NOT_ESTABLISHED',Full66_posterior='UNQUALIFIED_NOT_MAIN_GATE',strong_TAP='NOT_YET_ESTABLISHED',W_transfer='NOT_RUN'),resource=r,posthoc_geometry=read('results/possibility_volume/posthoc_geometry/SUMMARY.json'))
write('RESULTS_SUMMARY.json',s)
print(json.dumps({'full576_RHS':r['new_full_RHS'],'reduced256_RHS':r['reduced_RHS256'],'CPU':r['CPU_seconds_sum'],'wall':r['job_wall_seconds_sum'],'E5_peakRSS':peak}))
