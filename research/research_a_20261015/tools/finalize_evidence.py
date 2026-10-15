from pathlib import Path
import json,subprocess,platform,sys
import numpy as np,scipy
from b1_joint.accounting import dump
from research_a.budget import Budget
from research_a.target import OUT,ROOT
ledger=Budget(OUT)
with ledger.job('final-numerical-evidence-postprocess-repair1'):
 e1=ledger.summary();e2=json.loads((OUT/'results/E2/RESULTS.json').read_text());post=json.loads((OUT/'results/E2/E2_POSTPROCESS_JOB.json').read_text())
 summary=json.loads((OUT/'results/SUMMARY_E1.json').read_text());timing=json.loads((OUT/'results/COLD_WARM_SEARCH_TIMING.json').read_text())
 times={}
 for arm in summary['methods']:
  r=[t for t in timing['records'] if t['arm']==arm]
  times[arm]=dict(cold_initial_first4_median=float(np.median([t['cold_initial_search_seconds'] for t in r])),warm_replayed_first4_five_rotations_median=float(np.median([t['warm_wall_seconds'] for t in r])),new_full_warm_attempts=0,exact_endpoint_cache_requests=len(r)*12,scope='entireproposal/correction/selection warmreplaywithvalidendpointcache; base/modelconditioning separate; coldincludesendpointfullcalls')
 dump(OUT/'results/TIMING_SUMMARY.json',times)
 pairs=json.loads((OUT/'results/PAIRWISE_AMBIGUITY.json').read_text());export=json.loads((OUT/'results/ACTUAL_PACKAGE_EXPORT.json').read_text())
 gate=dict(research_stage='RESEARCH-A-20261015',overall='KEEP_RADIAL_CANDIDATE_BENCHMARK_AND_B_HANDOFF_STOP_CURRENT_FEEDBACK_ADAPTIVE_AND_E3',E0='PASS',E1_physical_worlds='PASS_IN_RADIAL_SUPPORT',E1_Best4='POSITIVE_LIMITED_RESULT',E1_strong_conjunctive_algorithm_gate='FAIL',E2='NO_GO_CURRENT_MECHANISM',E3='NOT_RUN_NOT_AUTHORIZED',posterior_weights=None,IEEE_TAP='INSUFFICIENT_EVIDENCE',archive_worlds=193,observation_count=32,witness=pairs['first_fixed_order_witness'],method_statistics=summary,all_historical_decisions_unchanged=True)
 dump(OUT/'GATE_DECISION.json',gate)
 paths=dict(B1=ROOT,C1=ROOT.parents[1]/'C1/c1-joint-p2',A23=ROOT.parents[1]/'A23/A23_research',A24=ROOT.parents[1]/'A24/a24-source-anchored-om-feedback',A30=ROOT.parents[1]/'A30/a30-maxwell-material-manifold')
 repos={}
 for name,path in paths.items():
  repos[name]=dict(path=str(path),commit=subprocess.check_output(['git','-C',str(path),'rev-parse','HEAD'],text=True).strip(),status=subprocess.check_output(['git','-C',str(path),'status','--short'],text=True).strip())
 dump(OUT/'UPSTREAM_MANIFEST.json',dict(repos=repos,A17='closure_final source snapshot; no independentgit; see delegated/history/EVIDENCE.json',assets=dict(package_received='Research_A_Maxwell_Compensation_Package.zip extractedpro_original',bank='B1 JointP2 COMMON_RESPONSE_BANK1250paidstates',accepted_reference='P2accepted continuousreference retained; notonlinefeature',checkpoint='C1 P2 full_distill-26101321.pt original4Dproxy target; no newtraining',A30='inputs/a29_release saved24x24TM backgrounds andN3; oldledgerreadonly'),digest_policy='NoSHA256computed/checkedorinvented; Git revisionsareprovenance identifiers',original_P2_count=3459,original_B1_count=384,new_authorization='E12304/1800s,E2conditional384/9216RHS/3600s,nontransferable',numerical_versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,platform=platform.platform()),source_access_gaps='SOURCE_ACCESS_GAPS.md'))
 dump(OUT/'TRUTH_ISOLATION_AUDIT.json',dict(status='PASS_AT_DECLARED_PROCESS_BOUNDARY',online_files=['DEV_OBSERVATIONS.json','CONFIRM_OBSERVATIONS.json'],offline_files=['offline_truth/DEV_TRUTH.json','offline_truth/CONFIRM_TRUTH.json','A30 OFFLINE_TARGETS truth openedafterallsearches'],search_signature='observation/sigma/base/support/history/seed/budget only',forbidden_import_check='E0 AST search.py hasnooffline import',selection='baseplusfarthestfield; minimumresidual singleimage; nottruth',packageexport='readsCONFIRM_CAMPAIGN only, nooffline truth',fresh_protocol_file='configs/CONFIRM_FREEZE.json writtenbeforetruthjobBEGIN',seedstreams='independentparent261015102/noise261015103/proposal261015101',scope='codepaths/files andsavedreceipt audit; notOSsandboxsecretprotection',old_locked_truth_reads=0))
 # Indexed hierarchy, not duplicate reservations; E2 retainedonly inits ownledger.
 if not any(r.get('event')=='EXTERNAL_SUBLEDGER' for r in ledger.rows()):ledger.append(dict(event='EXTERNAL_SUBLEDGER',external_stage='E2',path='results/E2/COST_LEDGER.jsonl',state_attempts=e2['resources']['new_full_state_attempts'],complex_RHS=e2['resources']['new_complex_LU_RHS'],counting='referenceonly; do notsumasadditionalBEGIN'))
 print(json.dumps(times,indent=2))
# Closeafterthepostprocessjobreceipt exists; writingreportisnotanother numericaljob.
e1=ledger.summary()
resource=dict(E1=e1,E2=e2['resources'],E2_postprocess=post,zero_physics_E0=dict(packagealgebra_seconds=json.loads((OUT/'results/E0_PACKAGE_ALGEBRA.json').read_text())['runtime_seconds'],interface='0 newphysics'),historical_shared_setup=dict(P2_ledger_closed=3459,originalB1_closed=384,bank1250='already paidandshared allmethods',C1_P2_training='15paidCPUmodels historical; checkpointproposalextra supervision disclosed',A30='paidknownbackgroundgeometry+N3 stateLU reused, freshLUaudit+32adjointRHS charged'),GPU=0,deepseek_calls=0,old_locked_truth_reads=0,publication_IO_included_in_numericwall=False,unmetered_small_initial_postprocessing=dict(scoredev1_wall_seconds=.671970167,scoredev2_wall_seconds=.371983583,scoredev3_wall_seconds=.338385834,scoreconfirm_wall_seconds=.521929792,note='observedtoolwall accounted separately; not physics; existingreadonly preparation also reported E0'),all_failure_and_rejection_costs_retained=True)
resource['E1']['controlled_execution_wall_including_initial_scoretools']=e1['controlled_execution_wall']+sum(resource['unmetered_small_initial_postprocessing'][k] for k in resource['unmetered_small_initial_postprocessing'] if k.endswith('_seconds'))
dump(OUT/'results/RESOURCE_SUMMARY_FINAL.json',resource)
dump(OUT/'RUN_CLOSED.json',dict(scientific_run_closed=True,E1_new_forward_requires_new_explicit_authorization=True,E2_completed_no_rerun=True,reason='intendedstudycomplete; E2RHSbudgeteffectivelyexhausted; noextraphysicsforsupplementalpolish',final_E1_attempts=e1['attempts'],final_E2_attempts=e2['resources']['new_full_state_attempts']))
