"""After50publicseals: explicit selections/attempt paths,never GT or new physics."""
import time,json,datetime,shutil
from pathlib import Path
ENTRY=time.perf_counter();ROOT=Path(__file__).resolve().parents[1];V=ROOT/'delegated/vector';OUT=V/'confirmation_host_union'
def read(p):return json.loads(Path(p).read_text())
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
barrier=read(OUT/'ALL50_PUBLIC_PREDICTION_AND_PREFIX_SEALS.json');assert barrier['global_tasks']==50 and not barrier['fresh_GT_access']
queue=read(ROOT/'configs/vector_confirmation_queue.json')['tasks'];assert len({x['task_id'] for x in queue})==50
registry=read(V/'ATTEMPT_PATH_REGISTRY.json');registry['relative_path_base']=str(ROOT.resolve());registry['path_resolution_rule']='Absolute paths remain absolute;olderdelegated/vector paths resolve against relative_path_base;originalreceipt strings are not rewritten';mapfile=V/'gpu_recovery_attempt2_transport/CANONICAL_COLLECTION_MAP.json';mapping=read(mapfile)['mappings'];assert len(mapping)==9
if not registry.get('recovery2_mappings_finalized'):
 for m in mapping:
  registry['entries'].append(dict(task_id=m['task_id'],scientific_attempt=0,execution_generation='gpu_recovery_attempt2',status='FIRST_SCIENTIFIC_ATTEMPT_NATURALLY_SEALED',historical_local_root=str((ROOT.parents[1]/m['original_archive_path']).resolve()),historical_literal_archive_path=m['original_archive_path'],current_local_root=str((ROOT.parents[1]/m['canonical_collection_path']).resolve()),prediction_sealed=True,old_receipt_paths_modified=False))
 registry.update(status='ALL50_SEALED_FINAL_ATTEMPT_PATH_REGISTRY',recovery2_mappings_finalized=True,pending_mapping=None,updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),scientific_confirmation_attempts_including_failure=51,controlled_no_FW_directories_are_not_scientific_attempts=True)
 dump(V/'ATTEMPT_PATH_REGISTRY.json',registry)
# Running publicmonitor loaded an earlier source snapshot; preserve originals and
# correct only flattening metadata,not any state/Phi/cost value or scientific seal.
flatpath=OUT/'ALL_PUBLIC_COST_PREFIX_CHOICES.json';costpath=OUT/'NESTED_ARM_COST_BY_HOST.json'
flat=read(flatpath);cost=read(costpath);selections=[]
for task in queue:
 ident=task['task_id'];d=V/'inversions'/ident;r=read(d/'RECEIPT.json');p=read(d/'COST_PREFIX_SEAL.json');retry=bool(r.get('infrastructure_attempt_metadata'))
 for row in flat:
  if row['task_id']==ident:row.update(cost_scope=p.get('cost_scope','NORMAL_FIRST_ATTEMPT_END_TO_END_ARM'),end_to_end_equal_cost_qualified=not retry)
 for row in cost['rows']:
  if row['task_id']==ident:row.update(infrastructure_attempt_metadata=r.get('infrastructure_attempt_metadata'),cost_scope='RETRY_ATTEMPT_ONLY_NOT_END_TO_END' if retry else 'NORMAL_FIRST_ATTEMPT_ARM')
 selections.append(dict(task_id=ident,public_id=task['public_id'],host=r['host_dispatch_role'],method=task['method_label'],
  selected_state_id=r['physical_corrector']['selected_state_id'],public_Phi=r['physical_Phi'],selection=r['physical_corrector']['selection'],
  prediction_path=str((d/'prediction.npz').resolve()),prediction_seal_path=str((d/'SEALED.json').resolve()),
  public_paid_prefix_seal_path=str((d/'COST_PREFIX_SEAL.json').resolve()),public_paid_prefix_material_path=str((d/'COST_PREFIX_MATERIALS.npz').resolve()),
  corrector_accepted=r['physical_corrector']['accepted'],corrector_requests=r['physical_corrector']['requests'],corrector_stop=r['physical_corrector']['stop'],
  current_update='ordinary3Dreal-isotropic CIE; unchanged inheritedcurrent block',material_readout=r['material_metric'],
  Green_basis_rank=120 if r['material_metric']!='ordinary' else None,Green_basis='Publicbackgroundsource-driven orth([Einc,G Einc]),notfullGapproximation' if r['material_metric']!='ordinary' else None,
  physical_model='True3DdyadicCM+RR;5184complexcurrent/1728realmaterial;all60sources',not_2D_G19=True,not_BothBlocks=True,
  infrastructure_attempt_metadata=r.get('infrastructure_attempt_metadata'),end_to_end_equal_cost_qualified=not retry,
  original_failed_attempt_path=str((V/'infrastructure_failed_attempts'/ident/'attempt0').resolve()) if retry else None,GT_access=False))
registry['canonical_predictions']=[dict(task_id=x['task_id'],prediction_path=x['prediction_path'],prediction_seal_path=x['prediction_seal_path'],host=x['host'],infrastructure_retry=bool(x['infrastructure_attempt_metadata'])) for x in selections]
dump(V/'ATTEMPT_PATH_REGISTRY.json',registry)
for path,value in ((flatpath,flat),(costpath,cost)):
 prior=path.with_name(path.stem+'_before_retry_scope_metadata_correction.json')
 if not prior.exists():shutil.copyfile(path,prior)
 dump(path,value)
dump(OUT/'SELECTIONS.json',dict(status='ALL50_PUBLIC_SELECTIONS_SEALED_BEFORE_FRESH_GT',tasks=selections,unique_task_ids=50,matched_tasks=40,crossgrid_tasks=10,independent_materials=8,scientific_inverse_attempts_including_failed_attempt0=51,
 primary_selector='Lowest public Phi among paidcorrectorstates,including rejectedtrials;noGT',fixed_prefix_selector='Lowest public Phi among allpaidphysicalstates at completedCPU/wall/LUbudgets;retryaffectedone isattempt-only,notendtoend',
 retry_affected_object='c17',all8_objects_retained=True,seven_object_sensitivity_excludes_whole_c17=True,original_failed_prefix_history_complete=False,
 attempt_path_registry=str((V/'ATTEMPT_PATH_REGISTRY.json').resolve()),fresh_GT_access=False,fresh_GT_authorized=False,new_forward_RHS=0))
dump(OUT/'PUBLIC_FINALIZATION_COST.json',dict(status='COMPLETE',CPU_seconds=time.process_time(),wall_seconds=time.perf_counter()-ENTRY,scope='Onlypublicselections/attemptregistry/flattenedmetadata;originalvalues+receipts retained,zeroGT/FW',new_forward_RHS=0,GT_access=False))
print(json.dumps(dict(status='ALL50_PUBLIC_DELIVERY_READY_GT_NOT_READ',selections=str(OUT/'SELECTIONS.json'))),flush=True)
