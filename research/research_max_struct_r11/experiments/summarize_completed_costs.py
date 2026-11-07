"""Read-only R11 cost consolidation. Stdlib only; no NPZ/GT, science or remote calls."""
from pathlib import Path
from datetime import datetime,timezone
from collections import defaultdict
import argparse,json,time


def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--study-root',type=Path,default=Path(__file__).absolute().parents[1]);ap.add_argument('--output-dir',type=Path);args=ap.parse_args();root=args.study_root.absolute();out=(args.output_dir or root).absolute();out.mkdir(parents=True,exist_ok=True);entries=[];warnings=[]
 def read(rel):
  p=root/rel
  if not p.exists():return None
  try:return json.loads(p.read_text())
  except (ValueError,OSError) as e:warnings.append(dict(source=rel,error=str(e),action='Not counted; rerun after producer completes atomic update'));return None
 def source(rel):
  p=root/rel
  if not p.exists():return dict(path=rel,exists=False)
  st=p.stat();return dict(path=rel,exists=True,bytes=st.st_size,mtime_ns=st.st_mtime_ns)
 def physical(dimension,raw=None,zero=False,**overrides):
  raw=raw or {};r=dict(dimension=dimension,LU=0 if zero else raw.get('LU',raw.get('physical_LU')),complex_RHS=0 if zero else raw.get('actual_complex_RHS',raw.get('complex_RHS')),source_RHS=0 if zero else raw.get('source_RHS',raw.get('physical_source_RHS',raw.get('primal_RHS'))),adjoint_RHS=0 if zero else raw.get('adjoint_RHS',raw.get('physical_adjoint_RHS')),iterative_forward_batches=0 if zero else None)
  r.update(overrides)
  if r['complex_RHS'] is None and r['source_RHS'] is not None and r['adjoint_RHS'] is not None:r['complex_RHS']=r['source_RHS']+r['adjoint_RHS']
  r['counts_known_null_reason']='Known counters only; missing RHS splits/batch counts were not recorded and are not inferred. Small denominator LU and RHS are excluded from this physical view.' if any(v is None for k,v in r.items() if k!='dimension') else 'Explicit reported counters or a producer explicitly doing no physical solve.'
  return r
 def add(id,scope,src,status,cpu=None,basis='UNMEASURED',phys=None,aux=None,notes=None,additive=True,coverage='COMPLETE_KNOWN_SCOPE',nested=None):
  reasons=[];timing=[]
  for rel in src:
   if not rel.endswith('.json'):continue
   raw=read(rel)
   if isinstance(raw,dict):
    fields={k:raw[k] for k in ('wall_seconds','campaign_shared_wall_seconds','CUDA_event_stream_interval_ms','CUDA_event_scope','GPU_peak_allocated_bytes','GPU_peak_reserved_bytes','scope','CPU_scope') if k in raw}
    if fields:timing.append(dict(source=rel,reported=fields,additive=False))
  if cpu is None:reasons.append('CPU not measured, still running, or only nested into another counted producer; not imputed as zero.')
  if phys is None:reasons.append('Physical counts unavailable for this scope; not imputed as zero.')
  entries.append(dict(id=id,scope=scope,status=status,sources=[source(s) for s in src],CPU_known_lower_bound_seconds=cpu,CPU_measurement_scope=basis,CPU_additive=additive,coverage=coverage,physical_counts=phys,reported_timing_context=timing,reported_auxiliary_counts=aux or {},nested_or_excluded_sources=nested or [],notes=notes or [],counts_known_null_reasons=reasons or ['All fields stated below retain their producer scope; no missing cost is synthesized.']))
 def zero():return [physical(None,zero=True)]
 # Authoritative upstream ledgers already deduplicate their own transfer/controller scopes.
 for label,path,kind in [('twod_remote_DEV','delegated/twod/TWOD_DEV_COSTS.json','inversions'),('twod_independent48','delegated/twod/independent48/COSTS.json','jobs')]:
  v=read(path)
  if not v:continue
  jobs=v.get(kind,[]);raw=dict(LU=sum(x.get('LU',0) or 0 for x in jobs),complex_RHS=sum(x.get('complex_RHS',0) or 0 for x in jobs));ph=[physical(576,raw)]
  notes=[v.get('cost_rule','Producer ledger disjoint component sums are reused verbatim.'),f'Physical completed-job counts cover {len(jobs)} recorded tasks; not pending or live CPU.',*v.get('unknown_costs',[])]
  if kind=='jobs':
   gen=read('delegated/twod/independent48/GENERATION_RECEIPT.json');
   if gen:ph.append(physical(576,gen.get('physical_cost',{})));notes.append('Generation CPU is already in producer known_CPU_lower_bound; its separate physical counters are included once.')
   progress=read('delegated/twod/dev_transport/INDEPENDENT48_STATUS.json');notes.append('Latest progress snapshot may be ahead of synchronized counters; do not add its full CPU to this ledger. See live_status_views.')
  add(label,'Remote 2D DEV attempts/setup/transport/failed metadata' if kind=='inversions' else 'Remote 48x6 confirmation completed lower bound, generation and transport', [path],v['status'],v.get('known_CPU_lower_bound_seconds'),'Upstream disjoint known-CPU lower bound; full child GetProcessTimes plus distinct controller/transport components',ph,notes=notes,coverage='RUNNING_COMPLETED_LOWER_BOUND' if v['status'] not in ('COMPLETE','COMPLETED_KNOWN_CPU_LOWER_BOUND') else 'COMPLETED_WITH_UNMEASURED_GAPS')
 # Local DEV campaigns. Scalar-seed/readout checks and scoring are separate processes.
 for name,path,qpath in [('contracted_DEV3','delegated/diagnostics/CONTRACTED_DEV_COST_LEDGER.json',None),('BothBlocks_DEV3','delegated/diagnostics/BOTHBLOCKS_DEV_COST_LEDGER.json','delegated/diagnostics/BOTHBLOCKS_QUEUE.json'),('current_only_DEV3','delegated/diagnostics/current_only_dev/CURRENT_ONLY_DEV_COST_LEDGER.json','delegated/diagnostics/current_only_dev/CURRENT_ONLY_QUEUE.json')]:
  v=read(path)
  if not v:continue
  cpu=v.get('external_whole_child_CPU_seconds_sum',v.get('whole_child_CPU_seconds_sum'));add(name,'Three actual local exposed-DEV inversions',[path],v['status'],cpu,'Exact external wait4 child CPU' if v.get('external_whole_child_CPU_seconds_sum') is not None else 'Near-exit process_time lower bound; final external exit CPU unavailable',[physical(576,v.get('physical_cost',{}))],notes=[v.get('CPU_precision_note',''), 'Setup/imports/auxiliary/corrector already included in child CPU; do not add nested readout/solver timers.'])
  if qpath and (q:=read(qpath)):add(name+'_controller','Distinct local queue controller',[qpath],q['status'],q.get('controller_CPU_seconds'),'Own queue lifetime process_time; excludes children',zero())
  add(name+'_scoring','Standalone cached postseal scoring',[path],v['status'],v.get('scoring_CPU_seconds'),'Own scoring function timer, imports/exit not sampled; lower bound',zero())
 # Completed16 has distinct exact children, generator, controller and postprocess producers.
 path='delegated/diagnostics/bothblocks_confirmation/COST_LEDGER.json';v=read(path)
 if v:
  add('BothBlocks16_inversions','80 actual local confirmation children',[path],v['status'],v['exact_inversion_child_CPU_seconds'],'Exact external wait4 complete child CPU',[physical(576,v['physical_inversion_cost'])],v.get('auxiliary_inversion_cost_count_sums'),notes=['All5 methods same Mac; rawCPU is not pooled as a performance comparison with Windows48.','Portfolio overlapping CIE attempts are already in these80 children; no portfolio add-on charge.'])
  add('BothBlocks16_generation','16 generation children/exports',[path],v['status'],v['exact_generation_child_CPU_seconds'],'Exact external wait4 generator CPU',[physical(576,v['generation_physical_cost'])])
  add('BothBlocks16_controller','Distinct controller, no science',[path],v['status'],v['controller_CPU_seconds'],'Own controller process_time; shared campaign wall never added',zero())
  post='delegated/diagnostics/bothblocks_confirmation/POSTSEAL_EXTERNAL_PROCESS_COSTS.json';ps=read(post)
  if ps:add('BothBlocks16_postprocess','Main scorer/frontier/two plot renders incl layout refinement',[post],'COMPLETE',sum(x['external_complete_child_CPU_seconds'] for x in ps),'Exact external wait4 per standalone child, each once',zero(),notes=['Includes main scoring, frontier and both plot invocations. FIGURES plot receipt and embedded ledger postprocess are nested mirrors, not extra CPU.'])
  for tag,path in [('prediction_seal','delegated/diagnostics/bothblocks_confirmation/PREDICTION_SEAL_RECEIPT.json'),('common_LU_selection','delegated/diagnostics/bothblocks_confirmation/COMMON_PHYSICS_CAP_SELECTION_BEFORE_GT.json')]:
   x=read(path)
   if x:add('BothBlocks16_'+tag,'Public-only selection/seal process',[path],x['status'],x.get('CPU_seconds'),'Own main-function CPU timer, excludes import/exit; lower bound',zero())
 # Diagnostics with a corrected attempt: physical counters are already combined, CPU timers are disjoint.
 path='delegated/diagnostics/DIAGNOSTIC_RECEIPT.json';v=read(path);repairpath='delegated/diagnostics/ORDINARY_REPAIR_RECEIPT.json';repair=read(repairpath)
 if v:
  cpu=v.get('cpu_seconds');cpu=cpu+repair.get('cpu_seconds',0) if cpu is not None and repair else cpu
  add('fixed_state10_counterfactuals','10 exposed endpoints; final30 candidates plus10 archived seed-control attempts',[path,repairpath],v['status'],cpu,'Original main diagnostic timer plus distinct repair-process timer; imports/exit unmeasured',[physical(576,v['physical_cost'])],v.get('auxiliary_cost'),notes=['Combined receipt already contains50 physical FW/LU and720 RHS, including the repaired ordinary10 FW. Do NOT add repair or SEED_ONLY_ATTEMPT physical counts again.'],nested=['delegated/diagnostics/SEED_ONLY_ATTEMPT_RECEIPT.json',repairpath+'#physical_cost already merged'])
 for id,path,cpu_key,aux_key in [('metric_split','delegated/diagnostics/METRIC_SPLIT_RECEIPT.json','cpu_seconds','operation_cost_total_including_corrected_attempt'),('BothBlocks_wiring','delegated/diagnostics/BOTHBLOCKS_WIRING_CHECK.json','CPU_seconds','operation_cost'),('current_curvature','RESULTS/current_curvature_exposed/RECEIPT.json','CPU_seconds','auxiliary_counts'),('current_remainder','RESULTS/current_remainder_exposed/RECEIPT.json','CPU_seconds','auxiliary_counts'),('exposed_DEV_scoring','RESULTS/exposed_dev_scores/lagged/RECEIPT.json','CPU_seconds',''),('stage_compare','delegated/diagnostics/BOTHBLOCKS_STAGE_COMPARISON_RECEIPT.json','CPU_seconds','')]:
  x=read(path)
  if x:add(id,'Independent zero-FW diagnostic / cached scoring',[path],x['status'],x.get(cpu_key),'Own diagnostic CPU timer, import/exit coverage limited',zero(),x.get(aux_key,{}),notes=['Corrected metric split includes20 low-rank actions/240 auxiliary RHS but only final attempt CPU was timed; earlier attempt CPU unknown.'] if id=='metric_split' else [])
 path='delegated/diagnostics/CONTRACTED_WIRING_CHECK.json';x=read(path)
 if x:add('contracted_wiring','Standalone contracted material closure check',[path],x['status'],None,'Whole helper CPU not recorded; nested readout11,121microseconds is not called whole process CPU',zero(),x.get('auxiliary_cost'))
 paths=['FIGURES/exposed_development/PLOT_PRIOR_RENDER_RECEIPT.json','FIGURES/exposed_development/PLOT_RECEIPT.json'];plots=[x for p in paths if (x:=read(p))]
 if plots:add('exposed_figures','Two distinct exposed figure render processes',[p for p in paths if read(p)],'COMPLETE',sum(x.get('CPU_seconds_near_exit',0) for x in plots),'Absolute near-exit render child process_time; exit tail unavailable',zero(),notes=['These are separate render invocations, not the independent16 plots included in postprocess.'])
 path='delegated/diagnostics/bothblocks_confirmation/finite_noise_prior_audit/EXTERNAL_PROCESS_COST.json';x=read(path)
 if x:add('finite_noise_prior_audit','All16x5 cached finite-noise/prior audit',[path],'COMPLETE',x['external_complete_child_CPU_seconds'],'Exact external wait4 core audit; editorial formatting CPU unmeasured',zero(),notes=[x.get('measurement_note','')],nested=['finite_noise_prior_audit/AUDIT_RECEIPT.json own CPU is inside this child'])
 # Vector independent audits / transport. No cross-size LU/RHS mixing.
 vbase='delegated/vector/'
 for name in ('top32_audit','source_basis','source_basis_audit','derivative_check','z_derivative_check','fixed_state_cuda','source_gn_cuda'):
  path=vbase+name+'/RECEIPT.json';x=read(path)
  if not x:continue
  cpu=x.get('CPU_seconds');notes=['Public saved-field/background audit; zero new physical Maxwell solve.']
  if name in ('fixed_state_cuda','source_gn_cuda'):notes+=['Recorded post-import CPU is only a lower bound. Despite an early template sentence, no whole-OS external CUDA-child CPU was actually measured; imports/exit remain null.']
  add('vector_'+name,'Standalone vector diagnostic/setup',[path],x['status'],cpu,'Recorded process or function CPU lower bound, producer scope retained',[physical(5184,zero=True)],x.get('counts'),notes=notes)
  if x.get('gate_helper_CPU_seconds') is not None:add('vector_'+name+'_gate_helper','Separate cooperative-gate helper process',[path],x['status'],x['gate_helper_CPU_seconds'],'Helper process CPU separate from scientific parent, counted once',zero())
 for name in ('fixed_state_transport','source_gn_transport','dev_ordinary_transport','remaining_transport'):
  path=vbase+name+'/TRANSPORT_COST.json';x=read(path)
  add('vector_'+name,'Local transport/controller only; remote scientific CPU excluded',[path],x.get('status','COMPLETE') if x else 'RUNNING_NO_FINAL_RECEIPT',sum(x.get(k,0) for k in ('local_controller_CPU_seconds','local_transport_child_CPU_seconds')) if x else None,'Distinct Mac controller+SSH/SCP child CPU, not remote science',zero(),notes=['Transport wall overlaps remote science; do not add wall durations.'],coverage='RUNNING_COMPLETED_LOWER_BOUND' if not x else 'COMPLETE_KNOWN_SCOPE')
 for name in ('MEASUREMENT_OVERLAY_COST','LINALG_PROTOCOL_OVERLAY_COST','LINALG_TIMING_AMENDMENT_COST'):
  path=vbase+'remaining_transport/'+name+'.json';x=read(path)
  if x:add('vector_'+name,'Separate local source/measurement overlay transport',[path],'COMPLETE_KNOWN_COMPONENTS',x.get('local_controller_CPU_seconds',0)+x.get('local_child_CPU_seconds',0),'Recorded local components only; absent child/extraction CPU is unknown',zero(),notes=[x.get('scope','')])
 path=vbase+'dev_ordinary_process/PROCESS_COST.json';x=read(path)
 if x:
  add('vector_ordinary_DEV6','Six ordinary DEV arms in one shared scientific process',[path],x['status'],x.get('CPU_seconds'),'Single whole process_time producer; six arm timers are nested',[physical(5184,x.get('counts',{}))],x.get('counts'),notes=['Arm receipts, corrector timers, paid events and arm CPU sum are not added.','primal_RHS and physical_source_RHS are aliases; likewise adjoint aliases.'],nested=['delegated/vector/inversions/v04_*DEV/RECEIPT.json','delegated/vector/inversions/v09_*DEV/RECEIPT.json'])
  add('vector_ordinary_DEV6_gate_helper','Distinct gate helper of ordinary shared process',[path],x['status'],x.get('gate_helper_CPU_seconds'),'Separate helper CPU once, not the six arm count',zero())
 path=vbase+'confirmation_generation/PROCESS_COST.json';x=read(path)
 if x:add('vector_confirmation_generation','10 independent vector observations generated, no inverse score',[path],x['status'],x.get('CPU_seconds'),'Whole generation process_time including imports/setup/exports',[physical(5184,LU=8,source_RHS=480,adjoint_RHS=0,complex_RHS=480),physical(17496,LU=0,source_RHS=120,adjoint_RHS=0,complex_RHS=120,iterative_forward_batches=2)],x.get('counts'),notes=['Eight coarse5184 LU/all60source =480 RHS; two fine17496 GMRES/all60 =120 RHS. Raw primal_RHS600 spans sizes and is NOT added to physical_source_RHS480.'])
 # Future completed phase whole-OS measurement replaces nested same-scope receipt CPU.
 phasefile=vbase+'remaining_controller/PHASE_COSTS.json';phases=read(phasefile) or {};phaseix={p['phase']:p for p in phases.get('phases',[])}
 for name,count in [('dev_denominator',6),('confirmation',50)]:
  procpath=vbase+name+'_process/PROCESS_COST.json';proc=read(procpath);phase=phaseix.get(name);engpath=vbase+'engineering_library_selection/RECEIPT.json';eng=read(engpath) if name=='dev_denominator' else None
  if phase or proc:
   cpu=phase['CPU_seconds'] if phase else proc.get('CPU_seconds');raw=(proc or {}).get('counts',{})
   add('vector_'+name,'Shared vector scientific phase including setup/probe and task arms',[phasefile,procpath],(proc or {}).get('status','PHASE_TERMINAL'),cpu,'Whole OS GetProcessTimes replaces same-scope process/arm CPU' if phase else 'Recorded phase process_time; future PHASE_COSTS replaces this CPU',[physical(5184,raw)] if raw else None,raw,notes=['Engineering6LU probe is nested in dev_denominator; never add its CPU or physical counts on top.'] if eng else [])
   if proc and proc.get('gate_helper_CPU_seconds') is not None:add('vector_'+name+'_gate_helper','Distinct phase gate helper',[procpath],proc['status'],proc['gate_helper_CPU_seconds'],'Separate helper once',zero())
  else:
   events={};logpath=vbase+'remaining_transport/run.log';lp=root/logpath
   if lp.exists():
    for line in lp.read_text(errors='replace').splitlines():
     try:event=json.loads(line)
     except ValueError:continue
     if event.get('event')=='PHASE_TASK_SEALED' and event.get('phase')==name:events[event['task_id']]=event
   known=(eng.get('CPU_seconds',0) if eng else 0)+sum(e.get('CPU_seconds',0) for e in events.values())
   raw=dict((eng or {}).get('counts_delta',{}));raw['physical_LU']=raw.get('physical_LU',0)+sum(e.get('physical_LU',0) or 0 for e in events.values())
   add('vector_'+name,'Running/queued vector phase; completed nested producer lower bound only',[procpath,phasefile,logpath]+([engpath] if eng else []),'RUNNING_OR_QUEUED_NOT_SCORED',known if eng or events else None,'Engineering probe plus disjoint sealed-arm CPU timers only; whole phase imports/setup/helper tail unknown',[physical(5184,raw)] if eng or events else None,raw,coverage='RUNNING_COMPLETED_LOWER_BOUND',notes=[f'{len(events)} sealed arm timers seen in public controller log; RHS for those arm events unavailable, only reported LU added.','These partial counters/CPU are replaced, not added, when whole phase receipt arrives.','Scientific child continues; no remote query or interference by this compiler.'])
 path=vbase+'remaining_controller/PROCESS_COST.json';x=read(path)
 if x:add('vector_remaining_controller','Controller own lifetime, no science',[path],x['status'],x.get('CPU_seconds'),'Controller only; PHASE_COSTS are not repeated',zero())
 else:add('vector_remaining_controller','Still-running controller, own final CPU not yet received',[path],'RUNNING',None,'Unknown until final own-process receipt',zero(),coverage='RUNNING_COMPLETED_LOWER_BOUND')
 # Early exposed-DEV postseal producer: independent metadata/pack/scorer costs, not freshGT.
 ep=vbase+'dev_early_postseal/PROCESS_COST.json';early=read(ep);scpath=vbase+'dev_postseal_scores/SCORING_COST.json';score=read(scpath)
 if early:
  add('vector_DEV_early_postseal_controller','Independent local waiting/export controller',[ep],early.get('status','COMPLETE'),early.get('CPU_seconds'),'Own controller CPU only, child/scorer CPU separate',zero())
  steps=early.get('steps',[]);local=sum(x.get('local_child_CPU_seconds',0) for x in steps if x.get('operation') in ('remote_pack','download_extract'))
  add('vector_DEV_early_postseal_transport_children','Distinct local pack/download transport children',[ep],early.get('status','COMPLETE'),local,'Local child OS CPU only, remote pack/scorer CPU excluded',zero())
  sc=[x for x in steps if x.get('operation')=='postseal_scorer'];scpu=sum(x.get('whole_child_CPU_seconds',0) for x in sc) if sc else (score or {}).get('CPU_seconds')
  if sc or score:add('vector_DEV_early_scorer','Old exposed DEV postseal scorer only',[ep,scpath],(score or {}).get('status','COMPLETE'),scpu,'Whole external scorer child CPU replaces nested SCORING_COST timer' if sc else 'Nested scorer timer lower bound; replaced when external child receipt arrives',zero())
 else:
  add('vector_DEV_early_postseal_controller','Independent early-DEV waiting controller, final costs not yet available',[ep],'RUNNING_WAITING_DEV6_TERMINAL',None,'Await final controller/transport/scorer receipts',zero(),coverage='RUNNING_COMPLETED_LOWER_BOUND')
  if score:add('vector_DEV_early_scorer','Completed exposed DEV scorer, external final measurement not yet synced',[scpath],score.get('status','COMPLETE'),score.get('CPU_seconds'),'Scorer timer lower bound; future external cost replaces this',zero())
 rp=vbase+'dev_early_postseal/REMOTE_PACK_COST.json';remote=read(rp)
 if remote:add('vector_DEV_early_remote_pack','Independent remote archive producer, no physics',[rp],remote.get('status','COMPLETE'),remote.get('CPU_seconds'),'Own remote Python archive CPU',zero())
 fh='delegated/diagnostics/bothblocks_confirmation/first_cn99_hit/EXTERNAL_PROCESS_COSTS.json';firsthit=read(fh)
 if firsthit:add('first_CN99_hit_posthoc','Two standalone zero-FW public selection then GT-scoring processes',[fh],'COMPLETE_POSTHOC_OFFLINE_DIAGNOSTIC',sum(x['external_complete_child_CPU_seconds'] for x in firsthit),'Exact external wait4 each once; source selections and original confirmation unchanged',zero())
 bc='delegated/diagnostics/bothblocks_confirmation/posthoc_budget_control/COST_LEDGER.json';bcontrol=read(bc)
 if bcontrol:
  add('posthoc_budget_control_children','Only triggered restarted original-Maxwell corrector children; original80 costs excluded',[bc],bcontrol['status'],bcontrol['new_exact_child_CPU_seconds'],'Exact external wait4 complete continuation children',[physical(576,bcontrol.get('new_physical_cost',{}))],notes=['Original CIE/Dual/G19 costs are deployment reference views already included in original80, not paid again.','Only eligible public-Phi states within per-object remaining allowance selected; atomic/exit overrun charged; actual spend not equal targets.'])
  add('posthoc_budget_control_controller','Distinct finite continuation controller',[bc],bcontrol['status'],bcontrol['controller_CPU_seconds'],'Own controller CPU separate from child sums',zero())
  bp='delegated/diagnostics/bothblocks_confirmation/posthoc_budget_control/POSTPROCESS_EXTERNAL_COSTS.json';post=read(bp)
  if post:add('posthoc_budget_control_postprocess','Public32pair seal then separate postseal scorer',[bp],'COMPLETE',sum(x['external_complete_child_CPU_seconds'] for x in post),'Exact external wait4 replaces nested scoring timer',zero())
 oracle='delegated/diagnostics/bothblocks_confirmation/oracle_trace_boundary/RECEIPT.json';o=read(oracle)
 if o:add('oracle_trace_boundary','Root saved-scored158626node CSV audit, no new GT arrays',[oracle],o['status'],o['CPU_seconds'],'Independent diagnostic own CPU timer; no new optimizer/FW',zero(),notes=[o.get('warning','')])
 # Coverage gaps and excluded aliases remain visible; absent costs are not invented.
 add('unmeasured_setup_and_failures','Recorded failed environment starts, source preparation and uninstrumented remote extraction',['delegated/vector/IMPLEMENTATION_NOTES.md','delegated/vector/COST_METADATA_CORRECTIONS.md','delegated/twod/TWOD_DEV_COSTS.json'],'KNOWN_GAPS',None,'No OS CPU receipt exists',None,notes=['Default python missing numpy and an incorrect venv path failed before science; CPU unmeasured.','Some source packaging, remote extraction/vendor-copy and metadata probe CPU were not instrumented.','Human/editorial changes and unreceipted control calls are outside any exact total; this ledger is a known lower bound.'],additive=False)
 add('historical_R10_reference','Read-only reused R10 states/operators/scores, not rerun R10 campaign',[],'REFERENCE_ONLY_NOT_CHARGED',0,'R11 only charges its actual new calls; old R10 campaign cost excluded',zero(),additive=False)
 add('portfolio_views','Overlapping deployment/secondary views, no new solver run',['delegated/diagnostics/bothblocks_confirmation/PORTFOLIO_RULE_BEFORE_DATA_SCORE.json','delegated/twod/independent48/secondary/REGISTRATION_RECEIPT.json'],'NONADDITIVE_VIEW',0,'Shared underlying attempts counted exactly once; deployment sums are separate scenarios',zero(),additive=False)
 # Normalize duplicate identity mistakes before output; counts aggregate by physical matrix dimension.
 assert len({e['id'] for e in entries})==len(entries)
 totals=defaultdict(float);dim=defaultdict(lambda:dict(LU_known=0,complex_RHS_known=0,source_RHS_known=0,adjoint_RHS_known=0,iterative_forward_batches_known=0,entries_with_unknown_counters=[]))
 for e in entries:
  if not e['CPU_additive']:continue
  if e['CPU_known_lower_bound_seconds'] is not None:totals['known_nonoverlapping_CPU_lower_bound_seconds']+=e['CPU_known_lower_bound_seconds']
  for p in e['physical_counts'] or []:
   if p['dimension'] is None:continue
   d=dim[str(p['dimension'])]
   for key in ('LU','complex_RHS','source_RHS','adjoint_RHS','iterative_forward_batches'):
    if p.get(key) is None:
     if e['id'] not in d['entries_with_unknown_counters']:d['entries_with_unknown_counters'].append(e['id'])
    else:d[key+'_known']+=p[key]
 status48=read('delegated/twod/dev_transport/INDEPENDENT48_STATUS.json');vectorstatus=read('delegated/vector/CURRENT_RUN.json');live=dict(independent48=status48,vector=vectorstatus)
 payload=dict(status='REPEATABLE_PARTIAL_KNOWN_COST_LEDGER',generated_utc=datetime.now(timezone.utc).isoformat(),scope='R11 actual paid work only; engineering bookkeeping, no scientific conclusion',monetary_cost=None,monetary_cost_null_reason='No billing or hardware tariff record; CPU/wall/counts are not converted to currency',sources_are_metadata_only=True,compiler_CPU_seconds_near_write=time.process_time(),compiler_cost_scope='This invocation own CPU shown separately, not cumulatively added across unlogged recompiles',new_forward=0,GT_files_read=0,totals=dict(totals),physical_totals_by_dimension=dict(dim),live_status_views=live,entries=entries,warnings=warnings,accounting_rules=['Never add nested arm/corrector/aux timers to their whole process.','576,5184 LU and17496 iterative solves remain separate.','Count source/primal and adjoint aliases only once.','Fixed-state merged50FW already includes ordinary repair10FW.','Portfolio scenarios overlap underlying attempts; never add them to campaign total.','GPU stream interval, wall time, active kernel time and CPU are not interchangeable.','Missing/running/unmeasured values stay null; displayed totals are lower bounds, not complete exact project cost.'])
 (out/'COST_LEDGER.json').write_text(json.dumps(payload,indent=2,allow_nan=False)+'\n')
 text=['# R11 cost ledger','',f"Generated {payload['generated_utc']}. Metadata-only compilation; no science, forward, optimizer, GT read or remote command.",'',f"Known nonoverlapping recorded CPU lower bound: **{totals['known_nonoverlapping_CPU_lower_bound_seconds']:.6f} seconds**. This mixes documented exact child measurements and partial timers as a lower bound; it is not an exact total or a cross-hardware speed comparison. Ongoing48/3D work and unmeasured setup/import/exit costs prevent a complete total.",'','|Physical matrix dimension|Known LU|Known complex RHS|Known source RHS|Known adjoint RHS|','|---:|---:|---:|---:|---:|']
 for size,d in sorted(dim.items(),key=lambda z:int(z[0])):text.append(f"|{size}|{d['LU_known']}|{d['complex_RHS_known']}|{d['source_RHS_known']}|{d['adjoint_RHS_known']}|")
 text+=['','Counters are lower bounds when producer splits are missing. 17496 uses two fine-grid iterative forward batches, not LU. Small denominator LU/RHS remain auxiliary counters in JSON.','', '|Entry|Status|Recorded additive CPU s|Scope/source|','|---|---|---:|---|']
 for e in entries:
  c=e['CPU_known_lower_bound_seconds'];cs='null' if c is None else f'{c:.6f}';cs+=' (view only)' if not e['CPU_additive'] else '';src=e['sources'][0]['path'] if e['sources'] else 'no new execution source';text.append(f"|{e['id']}|{e['status']}|{cs}|{e['scope']}; `{src}`|")
 text+=['','## Running and incomplete coverage','',f"48 progress snapshot: {status48.get('status')} / {status48.get('completed')} terminal / {status48.get('pending')} pending." if status48 else '48 status unavailable.',f"3D: {(vectorstatus or {}).get('status','UNKNOWN')}. Whole running-phase/helper/transport CPU is null until final receipt; completed engineering probe is only a replaceable lower bound.",'','## Deduplication and gaps','']+[f'- {r}' for r in payload['accounting_rules']]+['','- Two early CUDA audit receipts contained a template suggestion of external OS CPU measurement. The producer confirms none was captured: post-import CPU is retained as a lower bound, imports/exit stay unknown.','- DEV/48 transport scopes are divided by their authoritative upstream producer ledgers; this compiler never adds their raw transport logs again.','- Engineering6LU probe is included in its future denominator shared process, and six ordinary vector arm receipts are inside one ordinary process.','- CPU for ordinary repair is a separate process timer, while its10FW are already merged in the fixed-state50FW count.','- Plot layout-refinement attempts are included; mathematical/scientific experiments are not rerun to fill accounting gaps.','', 'Recompile after refreshing local producer ledgers:', '', '```', 'python research/research_max_struct_r11/experiments/twod_compile_costs.py', 'python research/research_max_struct_r11/experiments/twod_independent_costs.py', 'python research/research_max_struct_r11/experiments/summarize_completed_costs.py', '```', '', 'The compiler does not fetch remote results or stop/restart running queues. For exact scope, null reasons, raw auxiliary counts and source byte/mtime provenance, use COST_LEDGER.json.']
 (out/'COST_LEDGER.md').write_text('\n'.join(text)+'\n');print(json.dumps(dict(entries=len(entries),known_CPU_lower_bound_seconds=totals['known_nonoverlapping_CPU_lower_bound_seconds'],physical_dimensions=sorted(dim),warnings=len(warnings))))
if __name__=='__main__':main()
