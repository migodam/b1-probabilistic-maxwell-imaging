"""Public-metadata-only authoritative vector cost ledger; never loads GT/NPZ."""
import time,json,datetime,statistics
from pathlib import Path
ENTRY=time.perf_counter();ROOT=Path(__file__).resolve().parents[1];V=ROOT/'delegated/vector'
rows=[];gaps=[];used=set()
def read(rel,required=False):
 p=V/rel
 if not p.exists():
  if required:raise FileNotFoundError(p)
  return None
 used.add(rel);return json.loads(p.read_text())
def physical(x,dimension=5184):
 return dict(matrix_order=dimension,physical_LU=x.get('physical_LU',0),primal_complex_RHS=x.get('primal_RHS',x.get('physical_source_RHS',0)),adjoint_complex_RHS=x.get('adjoint_RHS',x.get('physical_adjoint_RHS',0)))
def add(ident,platform,cpu,wall,source,*,status='COMPLETE',scope='',counts=None,additive=True,notes=None,exact=False):
 rows.append(dict(id=ident,platform=platform,CPU_seconds=cpu,CPU_additive=additive,CPU_measurement='whole_OS_external' if exact else 'producer_measured_scope_near_end_or_function',wall_seconds=wall,wall_additive=False,status=status,source=source,scope=scope,physical_counts=counts or [],notes=notes or []))
def gap(ident,platform,scope,reason,source):gaps.append(dict(id=ident,platform=platform,CPU_seconds=None,scope=scope,reason=reason,source=source))
def cpus(x,keys):return sum(x[k] for k in keys if x.get(k) is not None)
def no_physics():return [dict(matrix_order=5184,physical_LU=0,primal_complex_RHS=0,adjoint_complex_RHS=0)]
# Independent public setup, one-off derivative wiring and fixed-state readout trials.
for name in ('top32_audit','source_basis','source_basis_audit','derivative_check','z_derivative_check','fixed_state_cuda','source_gn_cuda'):
 p=name+'/RECEIPT.json';x=read(p,True);cuda=name in ('fixed_state_cuda','source_gn_cuda');host='WINDOWS_CUDA' if cuda else 'MAC_CPU'
 add(name,host,x.get('CPU_seconds'),x.get('wall_seconds'),p,status=x['status'],scope='Independent measured diagnostic/build scope; zero new physical Maxwell solves',counts=no_physics(),notes=['Auxiliary/Green actions remain in producer receipt; physicalLU0 does not mean zero algebra cost.'])
 rows[-1]['reported_operator_counts']=x.get('counts',{})
 if cuda:
  add(name+'_query_helpers',host,x.get('gate_helper_CPU_seconds'),None,p,scope='Distinct helper CPU,not in Windows parent CPU',counts=no_physics())
  gap(name+'_import_exit',host,'Whole scientific child CPU beyond recorded post-import timer','WholeOS external child CPU was never measured; historical template corrected',p)
# Local transports; SSH/SCP CPU is local,never remote scientific CPU.
transports=('fixed_state_transport','source_gn_transport','dev_ordinary_transport','remaining_transport','c18_append_transport','gpu_recovery_attempt1_transport','gpu_recovery_attempt2_transport')
for name in transports:
 p=name+'/TRANSPORT_COST.json';x=read(p)
 if not x:gap(name,'MAC_CPU','Final transport cost','Not yet final',p);continue
 status=x['status']
 if name=='c18_append_transport':status='CONTROLLED_CANCEL_WAITING_ZERO_FW_RAW_FAILURE_RETAINED'
 if name=='gpu_recovery_attempt1_transport':status='CONTROLLED_NEXT_ARM_MKDIR_DRAIN_RAW_FAILURE_RETAINED'
 add(name,'MAC_CPU',cpus(x,['local_controller_CPU_seconds','local_transport_child_CPU_seconds']),x.get('wall_seconds'),p,status=status,scope='Local own process+local SSH/SCP children only;remote CPU excluded;wall overlaps remote work',counts=no_physics())
 logs=sorted((V/name).glob('pack*.log'));found=False
 for log in logs:
  for line in log.read_text(errors='replace').splitlines():
   try:r=json.loads(line)
   except ValueError:continue
   if 'CPU_seconds' in r:
    found=True;add(name+'_remote_pack','WINDOWS_CUDA',r['CPU_seconds'],r.get('wall_seconds'),str(log.relative_to(V)),scope='Distinct remote archive producer,zero physical solve;not local SSH CPU',counts=no_physics())
 if not found:gap(name+'_remote_pack','WINDOWS_CUDA','Remote zip archive producer','Remote pack CPU not separately instrumented;local SSH wall is not CPU',name+'/pack*.log')
 gap(name+'_remote_extract','WINDOWS_CUDA','Remote source extraction/import finalization','Not separately instrumented;local transport wall must not be substituted',name+'/extract.log')
# Every explicit future-source/host metadata overlay was a separate paid action.
for p in sorted((V/'remaining_transport').glob('*COST.json')):
 if p.name=='TRANSPORT_COST.json':continue
 rel=str(p.relative_to(V));x=read(rel);cpu=cpus(x,['local_controller_CPU_seconds','local_child_CPU_seconds'])
 add('overlay_'+p.stem,'MAC_CPU',cpu,x.get('wall_seconds'),rel,scope=x.get('scope','Local engineering overlay transport'),counts=no_physics())
 if 'local_child_CPU_seconds' not in x:gap('overlay_'+p.stem+'_children','MAC_CPU','Overlay child transport CPU','Explicitly unmeasured in producer',rel)
 gap('overlay_'+p.stem+'_remote_extract','WINDOWS_CUDA','Remote overlay extraction CPU','No remote external CPU receipt',rel)
# Ordinary DEV6: parent process timer and all arm/paid clocks are nested.
x=read('dev_ordinary_process/PROCESS_COST.json',True)
add('ordinary_DEV6','WINDOWS_CUDA',x['CPU_seconds'],x['wall_seconds'],'dev_ordinary_process/PROCESS_COST.json',status=x['status'],scope='Whole ordinary6 shared process timer;replaces all6nestedarmCPU',counts=[physical(x['counts'])])
rows[-1]['reported_operator_counts']=x['counts']
add('ordinary_DEV6_query_helpers','WINDOWS_CUDA',x['gate_helper_CPU_seconds'],None,'dev_ordinary_process/PROCESS_COST.json',scope='Distinct nvidia-smi helpers once',counts=no_physics())
# DEV6 and failed original confirmation share an external measured controller.
ph=read('remaining_controller/PHASE_COSTS.json',True)['phases']
for name in ('dev_denominator','confirmation'):
 phase=next(p for p in ph if p['phase']==name);proc=read(name+'_process/PROCESS_COST.json',True)
 ident='new_DEV6_including_one6LU_library_probe' if name=='dev_denominator' else 'original_GPU25_failed_after15seals_plus315accepted'
 add(ident,'WINDOWS_CUDA',phase['CPU_seconds'],phase['wall_seconds'],'remaining_controller/PHASE_COSTS.json#'+name,status=proc['status'],scope='GetProcessTimes whole scientific child replaces process/arm CPU;originalfailedphase includes15sealed+failedarm,not15armplusphase',counts=[physical(proc['counts'])],exact=True,notes=['6libraryLU physically included in newDEV6modelcounts;probeCPU also nested,never add again.'] if name=='dev_denominator' else ['320extra5184LU from failedc17attempt0 are included here exactly once.'])
 rows[-1]['reported_operator_counts']=proc['counts']
 add(ident+'_query_helpers','WINDOWS_CUDA',proc['gate_helper_CPU_seconds'],None,name+'_process/PROCESS_COST.json',scope='Windows GetProcessTimes parent excludes nvidia-smi child CPU;add once',counts=no_physics())
add('remaining_controller_own','WINDOWS_CUDA',read('remaining_controller/PROCESS_COST.json',True)['CPU_seconds'],None,'remaining_controller/PROCESS_COST.json',scope='Controller own only;phase children above separately counted',counts=no_physics())
gap('original_failed_timeout_helper','WINDOWS_CUDA','The one inherited nvidia-smi helper throwing TimeoutExpired','R10 exception occurred before helperCPU accounting;unknown is null,not zero','infrastructure_query_timeout_failure/FAILED_ATTEMPT_COST_AND_RECOVERY_SCOPE.json')
# Infrastructure retry1 and cached-query continuation2: separate paid scopes.
for name in ('gpu_recovery_attempt1','gpu_recovery_attempt2'):
 phase=read(name+'_controller/SCIENTIFIC_PROCESS_COST.json');proc=read(name+'_process/PROCESS_COST.json')
 if not phase or not proc:gap(name,'WINDOWS_CUDA','Recovery scientific child','Not yet terminal',name+'_controller/SCIENTIFIC_PROCESS_COST.json');continue
 add(name+'_scientific','WINDOWS_CUDA',phase['CPU_seconds'],phase['wall_seconds'],name+'_controller/SCIENTIFIC_PROCESS_COST.json',status='CONTROLLED_ZERO_FW_NEXT_ARM_DRAIN' if name.endswith('1') else proc['status'],scope='Additional whole external child CPU,not replacement of originalfailedphase;replaces ownnestedarmCPU',counts=[physical(proc['counts'])],exact=True)
 rows[-1]['reported_operator_counts']=proc['counts']
 add(name+'_query_helpers','WINDOWS_CUDA',proc['gate_helper_CPU_seconds'],None,name+'_process/PROCESS_COST.json',scope='All successful and retried query-helper CPU returned by wrapper;raw timeout logs overlap and are not added again',counts=no_physics())
 ctrl=read(name+'_controller/CONTROLLER_COST.json',True)
 add(name+'_controller_own','WINDOWS_CUDA',ctrl['CPU_seconds'],ctrl['wall_seconds'],name+'_controller/CONTROLLER_COST.json',scope='Own finite controller only;scientific child separate',counts=no_physics())
 add(name+'_preflight_helper','WINDOWS_CUDA',ctrl.get('preflight_query_helper_CPU_seconds'),None,name+'_controller/CONTROLLER_COST.json',scope='One standalone public resource snapshot helper,not scientific gatehelper',counts=no_physics())
# Generation: 8coarse LU and 2fine iterative forward objects,not600equal-sizeRHS.
x=read('confirmation_generation/PROCESS_COST.json',True)
fine=[r for r in x['actual_solve_records'] if r['label'].endswith('_grid18')]
fine_labels={r['label'] for r in fine}
if len(fine_labels)!=2 or len(fine)!=120 or any({r['source'] for r in fine if r['label']==label}!=set(range(60)) for label in fine_labels):raise ValueError('Finegeneration2all60RHSbatch provenance differs')
add('confirmation_generation','MAC_CPU',x['CPU_seconds'],x['wall_seconds'],'confirmation_generation/PROCESS_COST.json',scope='Whole measured generation process;private writer only,not scoring or inverse',counts=[dict(matrix_order=5184,physical_LU=8,primal_complex_RHS=480,adjoint_complex_RHS=0),dict(matrix_order=17496,physical_LU=0,primal_complex_RHS=len(fine),adjoint_complex_RHS=0,iterative_forward_batches=len(fine_labels),GMRES_iterations=x['counts']['gmres_iterations'])])
rows[-1]['reported_operator_counts']=x.get('counts',{})
# Mac wait4 includes the entire scientific process tree: do not add RSS helpers.
for name in ('mac_confirmation_controller','mac_confirmation_spooled_controller','mac_c18_return_controller'):
 jobs=read(name+'/PROCESS_RECEIPTS.json')
 if jobs:
  for j in jobs:
   proc=read('mac_confirmation_processes/'+j['public_id']+'/PROCESS_COST.json',True)
   add('mac_object_'+j['public_id'],'MAC_CPU',j['whole_child_CPU_seconds'],j['whole_child_wall_seconds'],name+'/PROCESS_RECEIPTS.json#'+j['public_id'],status='COMPLETE' if j['exit']==0 else 'FAILED',scope='Macwait4 complete scientific process and ALL children;replacesnestedprocess,all5armCPU,RSShelperCPU',counts=[physical(proc['counts'])],exact=True)
   rows[-1]['reported_operator_counts']=proc['counts']
 ctrl=read(name+'/CONTROLLER_COST.json')
 if ctrl:add(name+'_own','MAC_CPU',cpus(ctrl,['CPU_seconds','monitor_helper_CPU_seconds']),ctrl['wall_seconds'],name+'/CONTROLLER_COST.json',status='CONTROLLED_IDLE_PARENT_DRAIN' if name=='mac_confirmation_controller' else ctrl['status'],scope='Own controller+its ps/vm_stat helpers only;scientificwait4 children separate',counts=no_physics())
 else:gap(name+'_own','MAC_CPU','Controller own+monitorhelper CPU','Not yet terminal',name+'/CONTROLLER_COST.json')
x=read('mac_spool_controller_handoff/PROCESS_COST.json',True)
add('mac_handoff_helper_own','MAC_CPU',cpus(x,['CPU_seconds','query_helper_CPU_seconds']),x['wall_seconds'],'mac_spool_controller_handoff/PROCESS_COST.json',scope='Handoff helper own+querychildren only;PHASE2_INCLUSIVE_WAIT4_VIEW overlaps allspooledcontroller andobjectcosts,neveradd',counts=no_physics())
# c18 out-and-back cancellation: repeated physical task count is exactly0.
x=read('c18_return_mac_engineering/CANCELLATION_ENGINEERING_RECEIPT.json',True)
add('cancelled_waiting_c18_GPU_controller','WINDOWS_CUDA',x['waiting_controller_CPU_seconds'],x['waiting_controller_wall_seconds'],'c18_return_mac_engineering/CANCELLATION_ENGINEERING_RECEIPT.json',status=x['status'],scope='Whole terminated OWN waiting controller,0scientific child,measured GetProcessTimes;duplicatecopy in c18_gpu_append_controller excluded',counts=no_physics(),exact=True)
x=read('c18_return_mac_engineering/CANCELLATION_ACTION_COST.json',True)
add('c18_cancel_action','WINDOWS_CUDA',cpus(x,['own_CPU_seconds','powershell_inventory_child_CPU_seconds']),x['wall_seconds'],'c18_return_mac_engineering/CANCELLATION_ACTION_COST.json',scope='Separate cancellation action own+its powershell child,not targetcontroller',counts=no_physics())
for rel,keys in [('c18_return_mac_engineering/CANCELLATION_TRANSPORT_COST.json',['local_own_CPU_seconds','local_transport_children_CPU_seconds']),('gpu_recovery_controlled_drain/PLACEMENT_TRANSPORT_COST.json',['own_CPU_seconds','local_transport_children_CPU_seconds'])]:
 x=read(rel,True);add(rel.split('/')[0]+'_local_action','MAC_CPU',cpus(x,keys),x['wall_seconds'],rel,scope='Own local action+SSH/SCPchildren,not remote actionCPU',counts=no_physics())
sentinel='gpu_recovery_attempt1_inversions/c17_matched12_globalDual_CONF/CONTROLLED_DRAIN_SENTINEL.json';x=read(sentinel,True)
add('remote_controlled_drain_sentinel_action','WINDOWS_CUDA',x['action_CPU_seconds'],x['action_wall_seconds'],sentinel,scope='Distinct own metadataaction;not scientificworker;nextarm0FW provenfinallycounts',counts=no_physics())
# DEV scorer produces its plots in the same process: count external scorer once.
x=read('dev_early_postseal/PROCESS_COST.json',True)
add('DEV_early_postseal_controller_own','MAC_CPU',x['CPU_seconds'],x['wall_seconds'],'dev_early_postseal/PROCESS_COST.json',scope='Own finite local controller only',counts=no_physics())
steps=x['steps'];add('DEV_early_postseal_transport_children','MAC_CPU',sum(s.get('local_child_CPU_seconds',0) for s in steps if s['operation'] in ('remote_pack','download_extract')),None,'dev_early_postseal/PROCESS_COST.json',scope='Local SSH/SCP childCPU,not remote pack or scorer',counts=no_physics())
s=next(s for s in steps if s['operation']=='postseal_scorer');add('DEV_score_and_all3Dplots','MAC_CPU',s['whole_child_CPU_seconds'],s['wall_seconds'],'dev_early_postseal/PROCESS_COST.json#postseal_scorer',scope='Externalwait4 scorer includes allplots/exports;replacesnesteddev_postseal_scores/SCORING_COST',counts=no_physics(),exact=True)
x=read('dev_early_postseal/REMOTE_PACK_COST.json',True);add('DEV_early_remote_pack','WINDOWS_CUDA',x['CPU_seconds'],x['wall_seconds'],'dev_early_postseal/REMOTE_PACK_COST.json',scope='Distinct remote archive CPU',counts=no_physics())
x=read('dev_postseal_scores/METADATA_CORRECTION.json',True);add('DEV_scope_metadata_correction','MAC_CPU',x['CPU_seconds'],x['wall_seconds'],'dev_postseal_scores/METADATA_CORRECTION.json',scope='Only postscoremetadatafield correction,no new GT/FW',counts=no_physics())
x=read('confirmation_host_union/MERGE_COST.json')
if x:add('fresh50_public_merge_monitor','MAC_CPU',x['CPU_seconds'],x['wall_seconds'],'confirmation_host_union/MERGE_COST.json',scope='Publiconly collectionwait/shape+box/copy validation,noGT,walloverlaps science',counts=no_physics())
else:gap('fresh50_public_merge_monitor','MAC_CPU','Own publiccollectionmonitor CPU','Not yet all50terminal','confirmation_host_union/MERGE_COST.json')
x=read('confirmation_host_union/CONTROLLED_35_15_OBSERVER_RESTART.json',True)
observed=sum(int(z['observed_self_CPU_time_ps'].split(':')[0])*60+float(z['observed_self_CPU_time_ps'].split(':')[1]) for z in x['stopped_owned_public_monitors'])
add('two_stopped_public_observers','MAC_CPU',observed,None,'confirmation_host_union/CONTROLLED_35_15_OBSERVER_RESTART.json',scope='pscentisecond ownCPU observations before controlledstop,lowerbound;onlymetadata monitors,0FW',counts=no_physics());gap('stopped_observer_exit_tail','MAC_CPU','Stopped publicobserver CPU beyondpssample','Not externally measured','confirmation_host_union/CONTROLLED_35_15_OBSERVER_RESTART.json')
x=read('confirmation_host_union/PUBLIC_FINALIZATION_COST.json')
if x:add('public50_selection_and_registry_finalization','MAC_CPU',x['CPU_seconds'],x['wall_seconds'],'confirmation_host_union/PUBLIC_FINALIZATION_COST.json',scope=x['scope'],counts=no_physics())
# Fresh scorer is intentionally absent before root separate GT authorization.
x=read('confirmation_postseal_scores/SCORING_COST.json');external=read('fresh_postseal_external/PROCESS_COST.json')
if external:
 add('fresh50_score_and_all3Dplots','MAC_CPU',external['whole_scorer_child_CPU_seconds'],external['whole_scorer_child_wall_seconds'],'fresh_postseal_external/PROCESS_COST.json',scope='Externalwait4 scorer includingallplots;replacesnested SCORING_COST',counts=no_physics(),exact=True)
 add('fresh50_score_controller_own','MAC_CPU',external['controller_own_CPU_seconds'],external['controller_wall_seconds'],'fresh_postseal_external/PROCESS_COST.json',scope='Own finite postseal controller only,not childscorer CPU',counts=no_physics())
elif x:add('fresh50_score_and_all3Dplots','MAC_CPU',x['CPU_seconds'],x['wall_seconds'],'confirmation_postseal_scores/SCORING_COST.json',scope='Measured scorer includingallplots;near-endproducerCPU,not extraFW',counts=no_physics())
else:gap('fresh50_score_and_plots','MAC_CPU','Future root-authorized postseal scorer/plots','Not yet run/authorized;not a zero-cost completed task','confirmation_postseal_scores/SCORING_COST.json')
# Nested evidence for costs/CPUcomparisons,never additional charges.
nested=[]
for path in sorted((V/'inversions').glob('*/RECEIPT.json')):
 x=json.loads(path.read_text());nested.append(dict(task_id=path.parent.name,source=str(path.relative_to(V)),CPU_seconds=x['CPU_seconds'],wall_seconds=x['wall_seconds'],CPU_additive=False,wall_additive=False,host=x.get('host_dispatch_role'),CUDA_event_stream_interval_ms=x.get('CUDA_event_stream_interval_ms'),CUDA_event_scope=x.get('CUDA_event_scope'),counts=x['counts'],infrastructure_attempt_metadata=x.get('infrastructure_attempt_metadata'),replaced_by='Whole appropriateWindowsphase or Macwhole-objectwait4 above'))
for path in sorted(V.glob('*_process/SHARED_BACKGROUND_COST.json')):
 x=json.loads(path.read_text());add('nested_setup_'+path.parent.name,'MAC_CPU' if path.parent.name.startswith('mac') else 'WINDOWS_CUDA',x.get('CPU_seconds'),x.get('wall_seconds'),str(path.relative_to(V)),scope='Cold sharedsetup/basis/libraryprobe view already nested in scientificparent,NOT additive',additive=False)
for path in sorted((V/'mac_confirmation_processes').glob('*/SHARED_BACKGROUND_COST.json')):
 x=json.loads(path.read_text());add('nested_setup_'+path.parent.name,'MAC_CPU',x.get('CPU_seconds'),x.get('wall_seconds'),str(path.relative_to(V)),scope='Mac cold sharedsetup/basis is inside objectwait4,NOT additive',additive=False)
spool_views=[];batch_profiles=[];memory_peaks=[]
for path in sorted((V/'inversions').glob('*/MATERIAL_HISTORY_SPOOL.json')):
 x=json.loads(path.read_text());spool_views.append(dict(task_id=path.parent.name,source=str(path.relative_to(V)),CPU_additive=False,wall_additive=False,write_includes_flush_subset=True,reported=x))
for path in sorted((V/'inversions').glob('*/PHYSICAL_BATCH_PROFILE.json')):
 x=json.loads(path.read_text());grouped={}
 for r in x['records']:
  t=grouped.setdefault(r['operation'],dict(CPU_seconds=0.,wall_seconds=0.,records=0));t['records']+=1
  for k in ('CPU_seconds','wall_seconds'):t[k]+=r.get(k,0.) or 0.
 batch_profiles.append(dict(task_id=path.parent.name,source=str(path.relative_to(V)),groups=grouped,CPU_additive=False,wall_additive=False,scope='Producerbatchtiming only;operatornames preserved because CPU LU+primal and CUDA separateLU/solve have different batching'))
for name in ('dev_ordinary_process','dev_denominator_process','confirmation_process','gpu_recovery_attempt1_process','gpu_recovery_attempt2_process'):
 x=read(name+'/PROCESS_COST.json')
 if x:memory_peaks.append(dict(source=name+'/PROCESS_COST.json',platform='WINDOWS_CUDA',GPU_peak_allocated_bytes=x.get('GPU_peak_allocated_bytes'),GPU_peak_reserved_bytes=x.get('GPU_peak_reserved_bytes'),peak_bytes_not_additive=True))
for name in ('mac_confirmation_controller','mac_confirmation_spooled_controller','mac_c18_return_controller'):
 x=read(name+'/CONTROLLER_COST.json')
 if x:memory_peaks.append(dict(source=name+'/CONTROLLER_COST.json',platform='MAC_CPU',peak_sampled_concurrent_sum_RSS_bytes=x.get('peak_sum_RSS_observed_bytes'),sampling_seconds=20,peak_bytes_not_additive=True))
# Common data provenance and accounting policy. No heterogeneous speedup inferred.
barrier=read('confirmation_host_union/ALL50_PUBLIC_PREDICTION_AND_PREFIX_SEALS.json');status='ALL50_SCIENCE_SEALED_KNOWN_COSTS_WITH_EXPLICIT_NULL_GAPS' if barrier else 'RUNNING_COMPLETED_COST_LOWER_BOUND'
byhost={}
for x in rows:
 if x['CPU_additive'] and x['CPU_seconds'] is not None:byhost[x['platform']]=byhost.get(x['platform'],0)+x['CPU_seconds']
small_totals={}
for x in rows:
 if not x['CPU_additive']:continue
 raw=x.get('reported_operator_counts',{});n=raw.get('denominator_small_LU',0)
 if n:
  rank=raw['denominator_small_LU_dimension_sum']//n
  if rank*n!=raw['denominator_small_LU_dimension_sum']:raise ValueError('Mixed small-LU orders need explicit split')
  t=small_totals.setdefault(str(rank),dict(small_LU=0,small_triangular_complex_RHS=0));t['small_LU']+=n;t['small_triangular_complex_RHS']+=raw.get('denominator_small_triangular_complex_RHS',0)
physical_totals={}
for x in rows:
 if not x['CPU_additive']:continue
 for c in x['physical_counts']:
  key=str(c['matrix_order']);t=physical_totals.setdefault(key,dict(physical_LU=0,primal_complex_RHS=0,adjoint_complex_RHS=0,iterative_forward_batches=0,GMRES_iterations=0))
  for k in t:t[k]+=c.get(k,0)
# This compiler's own repeatedruns are separately retained rather than overwritten.
compiler=dict(UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),CPU_seconds=time.process_time(),wall_seconds=time.perf_counter()-ENTRY,scope='Own metadata-only costcompiler,zeroFW/GT,near-write timer')
with (V/'COST_COMPILER_RUNS.jsonl').open('a') as f:f.write(json.dumps(compiler)+'\n')
compiler_runs=[json.loads(s) for s in (V/'COST_COMPILER_RUNS.jsonl').read_text().splitlines()]
add('cost_compiler_all_recorded_runs','MAC_CPU',sum(x['CPU_seconds'] for x in compiler_runs),sum(x['wall_seconds'] for x in compiler_runs),'COST_COMPILER_RUNS.jsonl',scope='Each metadata-only compiler run ownCPU exactly once;history costs also appear here so entries sum matches platform totals',counts=no_physics())
byhost['MAC_CPU']=byhost.get('MAC_CPU',0)+sum(x['CPU_seconds'] for x in compiler_runs)
gap('near_end_timer_exit_tails','MAC_CPU_AND_WINDOWS_CUDA','Producerfunction/near-endtimerswithout externalwholeOSmeters','Unmeasured interpreter/finalwrite/exit outside each statedproducer scope;platformsums are knownlowerbounds','entries.CPU_measurement')
gap('interactive_orchestration_metadata','MAC_CPU','Uninstrumented interactive orchestration/read/edit actions','No wholeOS CPU meter for individual agent shell/fileedit actions;not counted as numericalzero','COSTS.json scope')
result=dict(status=status,generated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),scope='VectorR11 only: all instrumented generation/diagnostic/DEV/freshscientific/failed/recovery/helper/controller/transport/scoring/plot scopes; explicit unknown gaps retained',compiler_GT_arrays_read=0,compiler_new_forward_RHS=0,GT_scope='OldexposedDEVpostseal scoring included;freshGTscorer only when separately rootauthorized receipt exists',monetary_cost=None,monetary_cost_null_reason='No tariff or monetary receipt; measured CPU/wall/LU/RHS are reported',excluded_producers=[dict(source='audit/MAC_POSTRUN_PACKAGE_METADATA.json',known_CPU_lower_bound_seconds=.107974,reason='Already counted separately by root total ledger'),dict(source='Root sharedWindowsenvironment inventory,main48 ledger entry',reason='Shared environment inventory already counted in root main48,not duplicated here'),dict(source='Root scientific theory/review/publication and all2D producers',reason='Outside vector worker producer scope;root ledger owns these')],entries=rows,nested_arm_receipts=nested,nonadditive_spool_views=spool_views,nonadditive_batch_profiles=batch_profiles,memory_peaks=memory_peaks,compiler_runs=compiler_runs,known_CPU_lower_bound_by_platform_seconds=byhost,heterogeneous_CPU_seconds_not_an_algorithm_speed_comparison=True,wall_durations_not_additive=True,physical_totals_by_matrix_order=physical_totals,denominator_counts_by_small_matrix_order=small_totals,null_gaps=gaps,physical_LU_dtype='complex128',real_material_dtype='float64',inverse_material_coordinates=1728,inverse_current_shape=[5184,60],receiver_shape=[180,60],small_LU_dtype='complex128',scientific_confirmation_unique_ids=50,scientific_confirmation_sealed_local=sum((V/'inversions'/t['task_id']/'SEALED.json').exists() for t in json.loads((ROOT/'configs/vector_confirmation_queue.json').read_text())['tasks']),original_failed_attempt_extra_LU=320,original_failed_phase_replaces15sealed_arm_CPU_plus_failedarm=True,public_library_probe6LU_nested_in_DEV6=True,Mac_wait4_includes_all_children_RSShelpers=True,source_and_primal_RHS_are_aliases=True,denominator_small_LU_order=120,top32_trial_small_LU_order=32,attempt_path_registry='ATTEMPT_PATH_REGISTRY.json',accounting_rules=['Never add nestedarm/setup/paid/profile/spool timers to wholephase or Macwait4','Never addMacPHASE2inclusivewait4view to objectwait4+controllerown','OriginalfailedGPUphase AND independentretryphase are both paid;320failedLUincludedexactlyonce','ReturnedwrapperhelperCPUalreadycontains timeouthelperCPU;rawtimeoutrows are nonadditive','Cancelledc18GPUwaitingcontroller0FW;its duplicatecostreceipt paths are aliases','Mac/WindowsCPUwork totals are separate;concurrentwall andwaitingwall never summed as researchmakespan','Affectedc17sourceDualretry prefixes are ATTEMPT_ONLY_NOT_END_TO_END;original315cost retained;all8main plus7object sensitivity','No new scientificexperiment/benchmark/GTaccess by thiscompiler'])
(V/'COSTS.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(dict(status=status,entries=len(rows),null_gaps=len(gaps),known_CPU_by_platform=byhost,sealed=result['scientific_confirmation_sealed_local'])),flush=True)
