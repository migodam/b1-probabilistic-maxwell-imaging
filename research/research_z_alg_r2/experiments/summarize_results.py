"""Consolidate existing receipts only. This script performs no physical solves."""
import json
from pathlib import Path
P=Path(__file__).resolve().parents[1]
def read(name): return json.loads((P/name).read_text())
def write(name,data): (P/name).write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n')
h=read('results/quotient_holdout_summary.json'); scores=read('results/quotient_holdout_scores.json'); resource=read('results/RESOURCE_SUMMARY.json'); imaging=read('results/highfreq48_scores.json')
setup_kinds={'GREEN_SETUP','RECEIVER_SVD','FIXED_ACTIVE_SETUP'}
methods=[]
for r in h['methods']:
 setup=sum(e['wall_seconds'] for e in r['accounting']['events'] if e['kind'] in setup_kinds)
 methods.append(dict(method=r['method'],targets=r['target_count'],setup_seconds=setup,
  target_average_seconds=(r['total_wall_including_target_logging']-setup)/r['target_count'],
  total_wall_seconds=r['total_wall_including_target_logging'],full_complex_RHS=r['accounting']['total_complex_rhs'],
  qualified_material_roots=r['qualified_roots'],method_exceptions=r['failures'],fallbacks=r['fallbacks']))
q=methods[0]; comparisons=[]
import math
for r in methods[1:]:
 denominator=r['target_average_seconds']-q['target_average_seconds']
 n=(q['setup_seconds']-r['setup_seconds'])/denominator if denominator>0 else None
 comparisons.append(dict(control=r['method'],continuous_intersection_targets=n,
  smallest_positive_integer_advantage=(max(1,math.floor(n)+1) if n is not None else None),
  observed_four_target_wall_ratio=r['total_wall_seconds']/q['total_wall_seconds']))
write('results/CROSSOVER.json',dict(status='EMPIRICAL_AFFINE_PROJECTION_FROM_FOUR_TARGET_BATCH',
 formula='T(N)=T_fixed_setup + N*T_target_average',setup_event_kinds=sorted(setup_kinds),
 methods=methods,comparisons=comparisons,
 scope='One measured four-target batch in fixed method order. Observation-dependent polynomial/SVD/recovery/validation/logging costs remain in target average. This is not a separately measured one-target crossover; Python process/import time excluded equally from all method timings.'))
nearest=max(min(x['material_RMS'] for x in row['roots']) for row in scores['rows'])
summary=dict(schema='Z-ALG-R2-RESULTS-SUMMARY-v1',scientific_status='QUALIFIED_METHOD_RESULT_WITH_LIMITS',
 model='Finite scalar-TM disk-pulse nonlinear Hankel-Green VIE; not Born, continuous PDE, or 3D vector certification',
 theorem=dict(status='PROVED_UNDER_EXPLICIT_ASSUMPTIONS',unknown_real_materials=4,receivers=1,independent_sources=2,
  quotient_length=6,requirements=['fixed Green/geometry','rank(q)=2','invertible cubic leading block T',
  'physical conjugate-swap coefficient structure','simple roots and separating eigenform for numerical extraction'],
  excluded=['nodes','resolvent poles','non-real material','out-of-box candidates'],
  numerical_completeness='NOT_RIGOROUSLY_CERTIFIED',reference='MATHEMATICAL_THEORY.md'),
 fresh_holdout=dict(status=scores['status'],targets=4,methods=methods,cases=h['cases'],
  returned_sets_contain_GT=sum(row['GT_in_returned_set_RMS_1e5'] for row in scores['rows']),
  total_method_target_sets=len(scores['rows']),maximum_nearest_GT_RMS=nearest,
  source_frozen_before_generation=True,all16_predictions_frozen_before_scoring=True,
  truth_used_for_initialization_or_selection=False,roots_identical_for_all_four_methods=True,
  unique_GT_selection='NOT_ESTABLISHED: case03 has another legitimate material about1.968RMS away',
  one_batch_timing=True,crossover='results/CROSSOVER.json',
  receipts=['results/quotient_holdout_summary.json','results/quotient_holdout_FREEZE_RECEIPT.json','results/quotient_holdout_scores.json']),
 coverage=dict(exposed_quotient_tasks=101,nonempty_tasks=24,qualified_endpoints=27,
  all_root_sets_match_six_path=True,complex_candidates=606,
  E1_same_slice_distinct_material_RMS=0.6280685,
  connected_component_count='UNKNOWN',strong_64_start_control_recovers_previously_omitted_roots=True,
  all_material_domain_search='NOT_ESTABLISHED'),
 spatial_imaging=dict(unknown_materials=48,sources=8,receivers=96,frequency_k='20*pi',
  geometric_model='48 independent disk cells; all material contrasts unknown',
  method_totals=imaging['method_totals'],setup_components=imaging['setup_components'],
  clean_original_gates='3/3 per arm',noisy_raw_SVD='FAIL: large material error and box violations',
  noisy_refined_methods='Noise-scale fits with low material error; original exact-response gate remains FALSE',
  novelty='STANDARD_SVD_PLUS_FULL_J_NOT_CLAIMED_ORIGINAL',quotient_scaling_to48='NOT_RUN'),
 costs=dict(full_complex_RHS=resource['actual_distinct_full_complex_RHS'],
  breakdown=resource['full_RHS_breakdown'],small_reduced_RHS=resource['small_algebraic_counts_saved_receipts_only']['reduced_complex_RHS'],
  canonical_receipt_rows=resource['canonical_receipt_rows'],paid_optimizer_nonconvergences=len(resource['failures']['optimizer_nonconverged_paid_attempts']),
  nonconverged_holdout_starts=22,full_physics_operation_errors=0,
  wall_scope=resource['timing'],hardware=resource['hardware'],agent_token_and_service_cost='NOT_AVAILABLE'),
 failures=['Six-path root solving slower than FullJ8 despite fewer full RHS',
  '24-path endpoint collisions in8 slices; no numerical completeness certificate',
  'Most coarse E1 slices empty withinbox',
  'Low-frequency48-cell receiver matrix rank/conditioning fails direct inverse',
  'Raw receiver SVD fails at fixed relative complex noise1e-3',
  '28development RHS lack raw endpoint receipts',
  'ScholarQA retrieval/verify returned429; finite primary-source audit cannot establish global novelty'],
 contribution_decision=dict(method_paper_candidate='KEEP_WITH_EXPLICIT_SCOPE',
  global_novelty='HOLD',globally_certified_enumeration='NOT_ESTABLISHED',
  superior_coverage_to_strongest_control='NOT_SUPPORTED',Full66='NOT_RUN'),
 reproducibility=['run_solver.py','replay_holdout.py','REPRODUCE.md'],
 formula_compendium='FORMULA_COMPENDIUM.md',independent_review='audit/FINAL_SCIENTIFIC_REVIEW.md')
g=read('results/general6_summary.json'); gs=read('results/general6_scores.json')
gm=[]
for r in g['methods']:
 setup=sum(e['wall_seconds'] for e in r['accounting']['events'] if e['kind'] in setup_kinds)
 gm.append(dict(method=r['method'],total_wall_seconds=r['total_wall_including_logging'],
  fixed_setup_seconds=setup,target_average_seconds=(r['total_wall_including_logging']-setup)/3,
  full_complex_RHS=r['accounting']['total_complex_rhs'],qualified_material_roots=r['qualified_roots']))
summary['theorem'].update(family='2s independent real materials,s independent sources,one complex receiver',
  general_quotient_length='binom(2s,s)',general_basis='all monomials of total degree<=s',
  implemented_instances=[dict(materials=4,sources=2,quotient_dimension=6),dict(materials=6,sources=3,quotient_dimension=20)],
  general_reference='theory/general_quotient.md')
summary['general6_fresh_holdout']=dict(status=gs['status'],targets=3,all6_materials_unknown=True,
 methods=gm,cases=g['cases'],score_rows=gs['rows'],
 new_full_RHS=48,generation_RHS=9,all60_complex_candidates_distinct=True,polish_steps=0,
 additional_qualified_root_vs_frozen64_starts=True,extra_material_pair_RMS=1.1587903026954673,
 solver_frozen_before_GT_generation=True,predictions_frozen_before_scoring=True,
 comparison_scope='One fixed64-start budget and one three-target timing batch; specialized algebraic solvers not compared',
 wall_ratio=gm[1]['total_wall_seconds']/gm[0]['total_wall_seconds'],
 source='results/general6_summary.json')
summary['coverage']['general6_frozen64_start_missed_one_qualified_root']=True
summary['contribution_decision']['superior_coverage_to_strongest_control']='OBSERVED_FOR_GENERAL6_FROZEN64_STARTS_ONLY;NOT_UNIVERSAL'
summary['reproducibility'].append('replay_general6.py')
summary['costs']['prior_stage_RHS']=53903
summary['costs']['general6_new_RHS']=48
write('RESULTS_SUMMARY.json',summary)

print('Consolidated receipts, costs, and empirical crossover;0newRHS.')
