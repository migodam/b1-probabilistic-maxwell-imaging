"""Read sealed result/receipt summaries; never opens truth or runs physics."""
from pathlib import Path
import json,csv,time
P=Path(__file__).resolve().parent
def read(path):return json.loads((P/path).read_text())
def write(path,value):(P/path).write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
def main():
 ledger=read('RESOURCE_LEDGER_SNAPSHOT.json');closure=read('CLOSURE.json')
 assert not ledger['active_jobs'] and not ledger['active_actions'] and ledger['RHS']==closure['actual_RHS']
 b=read('track_b/TRAINING_JOB.json');c=read('track_c/TRACK_C_CLOSEOUT.json')['matched_training'];mps=b['wall_seconds']+c['GPU_reserved_whole_job_seconds']
 resources=dict(status='CLOSED_WITH_RECEIPTS',new_complex_RHS=ledger['RHS'],hard_complex_RHS=8000000,unused_complex_RHS=8000000-ledger['RHS'],
  stages=ledger['stages'],new_full_states=ledger['states'],LU=ledger['LU'],forward_RHS=ledger['forward_RHS'],derivative_RHS=ledger['derivative_RHS'],
  controlled_physics_job_wall_seconds=ledger['controlled_wall_seconds'],solver_kernel_wall_seconds=ledger['physics_wall_seconds'],physical_CPU_seconds=ledger['CPU_seconds'],
  conservative_MPS_training_job_seconds=mps,conservative_MPS_device_hours=mps/3600,MPS_cap_hours=48,
  new_provider_fits=len(b['records']),new_matched_B_head_fits=len(b['records']),new_C_head_fits=c['total_fits'],
  B_training_wall_seconds=b['wall_seconds'],C_training_wall_seconds=c['GPU_reserved_whole_job_seconds'],C_training_CPU_seconds=c['CPU_seconds'],
  physical_action_failures=ledger['failures'],administrative_failed_jobs=closure['administrative_failed_jobs'],historical_R5B_RHS_counted_as_new=0,
  hardware_memory='24GiB RAM; bounded target/LU lifetime; per-job sampled RSS and MPS memory in receipts, not a proof of maximum whole-machine memory',
  storage='lossless per-query q/prediction, final chi/current/field/J; no per-query LU archive',
  receipt_coverage='Central JSONL covers every attempted physical action. B and C training, B6, inference/analysis and delivery have separate receipts. Metadata/report editing and literature I/O are not individually CPU-profiled; physical job wall is not mislabeled as the complete campaign CPU total.',
  nonphysical_receipts=['track_b/TRAINING_JOB.json','track_c/MATCHED_TRAINING_LEDGER.jsonl','track_b/b6_moments/JOB.json','track_c/clean_inference_benchmark30','track_c/kkt_closure/','track_d/predictions/','track_d/figures/','track_a/COPULA_FIBRE_REFINEMENT_QA.json'])
 write('RESOURCE_SUMMARY.json',resources)
 result=dict(status='SCIENCE_COMPLETE',full66_posterior_qualified=False,root_scientific_judgments=dict(
  GO_A='REDIRECT: exact material fibre/compensation diagnostics and local Q1/Q2 moments; no global reference or task-cost advantage',
  GO_B='STOP expanding the tested N1024 objective family; preserve the historical trained functional estimator',
  GO_C='SCOPED_GO versus original frozen Full-actual baseline; stronger Full-q low-noise NI and cross-window20pct lower bound remain HOLD; ordinary direct Probit is competitive',
  GO_D='GO for the two tested noise scales; acquisition generalization NOT_RUN',
  GO_TAP='REDIRECT: useful computations and counterexamples, but no demonstrated uniquely Maxwell-specific algorithmic advantage'),
  stage_A=dict(RHS=ledger['stages']['A']['RHS'],model='active 2D/3D truncated-Gaussian submodels, not Full66 marginals',reference='local GH15/GH21 Q1/Q2; global tail not certified',cost_separation_established=False,dimension4='NOT_RUN',copula=read('track_a/COPULA_COUNTEREXAMPLE_RESULTS.json')),
  stage_B=dict(new_provider_fits=12,new_head_fits=12,train_parent_count=1024,validation_parent_count=256,validation_previously_exposed=True,objective_comparisons='feature, GT, feature+action, GT+action; same capacity/schedule/data',new_family_beats_original_N32768=False,Full66_posterior_claim=False),
  stage_C=dict(development_parents=64,new_matched_head_fits=9,original_full_baseline_preserved=True,direct_prior_development=read('track_c/probit_prior_baseline/evaluation_CDEV/RESULTS.json'),common_window_replay=read('track_c/common_window_timing/RESULTS.json')),
  conditions={},supplement={},resources=resources)
 table=[]
 for group in ('DLOW','DHIGH'):
  d=read(f'track_d/scores/{group}_RESULTS.json');r=d['mean_raw_full_grid_complex_MSE'];w=d['mean_anchor_plus_learning_seconds'];a=d['primary_selected'];f=d['primary_full']
  phases={name:d['secondary_seed_mean_metrics'][name]['raw']['anchor_phase'] for name in (f,a,d['actualgradient_reference'])}
  result['conditions'][group]=dict(n=d['n'],sigma_ratio=d['sigma_scale'],MSE=r,whole_wall_seconds=w,relative_risk_change_selected_vs_full=r[a]/r[f]-1,
   whole_wall_saving=1-w[a]/w[f],comparisons=d['comparisons'],cost_GO=d['cost_GO'],mechanism_GO=d['mechanism_GO'],stronger_full_q_posthoc=read('track_d/FULL_Q_REFERENCE_POSTHOC.json')[group],
   physical_point_baselines=d['physical_point_baselines'],optimizer=d['optimizer_and_fallback'],phase_metrics=phases,
   phase_ROI_counts=d['phase_ROI_counts_host_transition_pure'],selected_phase_risk_minus_full={k:phases[a][k]['full_grid_mse_contribution']-phases[f][k]['full_grid_mse_contribution'] for k in phases[f]},
   seed_disagreement={k:r[k]-d['secondary_ensemble_metrics'][k]['raw']['mse_full'] for k in r},
   qualification='separate256-parent zero-shot noise condition; mean three seed losses; global8 Bonferroni; no condition pooling or posterior claim')
  result['supplement'][group]=read(f'track_c/supplementary_prior_baseline/scores/{group}_RESULTS.json')
  for name in r:table.append(dict(group=group,sigma_ratio=d['sigma_scale'],method=name,MSE_meanseed=r[name],whole_wall_seconds=w[name]))
 with (P/'MAIN_RISK_COST_TABLE.csv').open('w',newline='') as f:o=csv.DictWriter(f,fieldnames=table[0].keys());o.writeheader();o.writerows(table)
 low=result['conditions']['DLOW'];old=low['comparisons']['whole_system_wall_20pct']['one_sided_upper'];delta=read('track_d/ADMINISTRATIVE_TIMING_SENSITIVITY_PRE_GT.json')['known_difference_upper_seconds']
 sensitivity=dict(status='DESCRIPTIVE_PRE_GT_SPECIFIED_SENSITIVITY',one_recovered_parent_Full_wall_range_width=delta,primary_wall_contrast_upper=old,
  conservative_upper_after_any_bootstrap_resampling=old+.8*delta,still_below_zero=old+.8*delta<0,
  proof='Each bootstrap mean can add at most0.8*delta when all draws select the one affected parent; actual average effect is smaller. No primary result rewritten.',new_RHS=0)
 write('track_d/ADMINISTRATIVE_TIMING_SENSITIVITY_RESULT.json',sensitivity)
 if (P/'RESULTS_SUMMARY.json').exists():
  previous=read('RESULTS_SUMMARY.json')
  if 'site' in previous:result['site']=previous['site']
 write('RESULTS_SUMMARY.json',result)
 print(json.dumps(dict(RHS=ledger['RHS'],MPS_hours=mps/3600,conditions={g:{k:v for k,v in x.items() if k in ('relative_risk_change_selected_vs_full','whole_wall_saving','cost_GO','mechanism_GO')} for g,x in result['conditions'].items()}),indent=2))
if __name__=='__main__':main()
