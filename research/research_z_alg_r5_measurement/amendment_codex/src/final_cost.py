from complete_quotient import *
from collections import Counter
orig=json.loads((BASE/'results/cost_accounting.json').read_text());stages={'original_inclusive_charged_before_plots':orig['charged_CPU_seconds_before_plots']};counts=Counter();rss=[orig['peak_case_rss_bytes']];count_files=[]
for f in ['plot_receipt.json','plot_revision_receipt.json','FINAL_VERIFICATION.json']:
 o=json.loads((BASE/'results'/f).read_text());stages[f]=o.get('CPU_seconds',0);rss.append(o.get('peak_rss_bytes',0))
for sub in ['results','results/development_v0_before_recovery_rank_audit']:
 for p in sorted((OUT/sub).glob('*.json')):
  o=json.loads(p.read_text())
  if not isinstance(o,dict):continue
  if p.name.startswith('complete_T'):
   # Count operations from individual roots, CPU only from enclosing batch.
   pass
  elif 'CPU_seconds_new_batch' in o:stages[str(p.relative_to(OUT))]=o['CPU_seconds_new_batch']
  elif 'CPU_seconds' in o:stages[str(p.relative_to(OUT))]=o['CPU_seconds']
  elif 'CPU_charge_upper_bound' in o:stages[str(p.relative_to(OUT))]=o['CPU_charge_upper_bound']
  for key in ['receipt','physics_receipt']:
   if key in o:
    counts.update(o[key].get('counts',{}));rss.append(o[key].get('peak_rss_bytes',0));count_files.append(str(p.relative_to(OUT))+':'+key)
# Conservative explicit provisions for earlier review diagnostics, amendment starts/brokers,
# one prototype before the original v0 batch, document/metadata work and packaging.
reserves={'prior_review_numerical_diagnostics':10.,'amendment_process_startups_and_broker_children':30.,'initial_single_case_prototype_before_v0_batch':3.,'document_metadata_final_packaging_allowance':60.}
charged=sum(stages.values())+sum(reserves.values());rootchecks=[]
for p in (OUT/'results').glob('complete_T*.json'):
 o=json.loads(p.read_text());rootchecks.append(dict(file=p.name,complex_roots=len(o['certificates']),root_radius_trials=sum(len(c['trials']) for c in o['certificates']),leading_dimension=o['leading']['dimension'],precision_bits=o['bits']))
report=dict(charged_CPU_seconds=charged,CPU_budget_seconds=7200,stage_inclusive_CPU_seconds=stages,additional_conservative_allowances_seconds=reserves,largest_recorded_process_rss_bytes=max(rss),memory_limit_bytes=8*1024**3,original_physics_counts=orig['physical_and_reduced_operations'],original_certification_counts=orig['certification_operations_separate'],amendment_counted_physics_operations=dict(counts),amendment_count_provenance=count_files,final_quotient_certificate_work=rootchecks,unmetered='LLM remote inference and tool-server/network CPU are not measured as local scientific CPU. Small local startup/metadata/packaging and earlier diagnostics have explicit conservative allowances; these are provisions, not exact timings or proven upper bounds.',accounting_rule='Inclusive stage times only; never sum nested ledger event times. Original tiny retries are in original charged total; amendment failure and early v0 included. New source policies are counterfactual alternatives, not a single deployment with all sources.',limits_satisfied_numerically=bool(charged<7200 and max(rss)<8*1024**3))
save(BASE/'results/FINAL_COST.json',report);print({k:report[k] for k in ['charged_CPU_seconds','largest_recorded_process_rss_bytes','limits_satisfied_numerically']})
