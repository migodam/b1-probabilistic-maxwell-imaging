#!/usr/bin/env python3
"""Approved Heun mechanism repair; small-h gate precedes task-scale physics."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[name]='1'
from pathlib import Path
import sys,json,argparse,time,traceback,gc,resource
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from p3r4.core import OUT,Book,Target,BudgetStop,save,write,plain
from p3r4.dynamic_geometry import (Chart,GeometryInvalid,material_frame,projector_field,material_real,
    distance_match,phase_diagnostic,finite_material_normal,heun_corrector)

JOB='dynamic_geometry_heun_r1'
RESULT=OUT/'results/III'/JOB
FIRST=OUT/'results/III/dynamic_geometry_first'
PROTOCOL=OUT/'configs/GEOMETRY_HEUN_PROTOCOL.json'
PARENTS=tuple(f'R4DEV_{i:05d}' for i in range(8))
HS=(.001,.002,.004)
RADII=(.02,.04,.08)
CAP=168192


def protocol():
 return dict(campaign='P3LR4-20261025',stage='III',job=JOB,revision=1,
   authority='Parent explicitly approved Heunrevision1 after verified Euler O(h^2) self-defect; first job must be CLOSED COMPLETE',
   parents=PARENTS,selection='first8fixed newDEV, not GT-selected',dimensions=[64,66],ranks=[2,4,8],
   fields='same frozen seeds and Euclidean generalized-material projector fields from first pilot; normalized by exact decoder metric',
   integrator='explicit Heun: predictor=x+tF(x); corrected=x+t/2*(F(x)+F(predictor)); each field refresh uses actual booked Full+J',
   regularity='at C1 only second-order Taylor consistency with o(h^2); no global O(h^3) at seams; empirical convergence required',
   orientations=dict(forward=[[0,1],[1,1],[0,-1],[1,-1]],reverse=[[1,1],[0,1],[1,-1],[0,-1]]),
   small_h=HS,reference_bracket_h=.001,reference_acquisition='already paid firstpilot 3h centralFD;480RHS/config=23040RHS for selected48; retained in first ledger',
   gate=dict(both_fine_h_relative_material_error_max=.2,
     fine_mid_change_max_factor_times_mid_coarse=1.25,change_absolute_floor=1e-8,
     require_no_material_rank_change=True,require_no_relative_gap_le_001=True,
     seams='flag separately; seam crossing does not imply high regularity or automatically pass gate'),
   gate_failure='no task-scale physics, no additional h, no new seeds or amplitude adjustment',
   task_radii=RADII,stress_radius=.08,
   task_step='sqrt(requested_RMS/||D reference_bracket||), both orientations; report actual endpoint distance, do not force match',
   direct_control='both signed reference bracket directions, exact decoder-distance match requested RMS, then freshly paid Full endpoint',
   latent_displacement_cap=8.,per_curve_worst_cost=dict(Full_states=8,J_builds=7,complex_RHS=288),
   all48_small_h_RHS=82944,all48_task_and_direct_RHS=85248,
   max_new_complex_RHS=CAP,parent_additional_cap=250000,max_controlled_wall_seconds=1800,
   teacher_eligible=dict(radius_relative_error_max=.01,signed_material_progress_min_fraction=.5,
     Full_DeltaPhi_max=2.,prior_energy_delta_max=2.,backward_residual_max=1e-10,prior_legal=True),
   teacher_semantics='truth-free local computational teacher only; no posterior sample/GT-risk qualification',
   no_GT=True,continuous_path_certificate=False,Heun_map_exactly_reversible=False,
   source_before_pay=True)


def append(path,row):
 with Path(path).open('a') as f:f.write(json.dumps(plain(row),allow_nan=False)+'\n');f.flush();os.fsync(f.fileno())


def rss():
 value=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
 return value if sys.platform=='darwin' else value*1024


def run():
 if json.loads(PROTOCOL.read_text())!=plain(protocol()):raise ValueError('Heun frozen protocol mismatch')
 first_receipt=json.loads((OUT/'III/jobs/dynamic_geometry_first/CLOSED.json').read_text())
 if first_receipt['status']!='COMPLETE':raise ValueError('firstpilot not sealed')
 if RESULT.exists():raise FileExistsError('no implicit Heun rerun')
 RESULT.mkdir(parents=True)
 write(RESULT/'FROZEN_PROTOCOL.json',protocol())
 for source in ('tools/run_r4_geometry_heun.py','src/p3r4/dynamic_geometry.py','src/p3r4/core.py'):
  dest=RESULT/'source'/source;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes((ROOT/source).read_bytes())
 book=Book('III',JOB,max_wall=1800);target=Target(book,keep_fields=True)
 costs=dict(RHS=0,states=0,LU=0);owners={};owner=['common'];original=book._update
 def bounded(event,**kw):
  if costs['RHS']+kw.get('rhs',0)>CAP:raise BudgetStop('Heun hard168192 RHS cap')
  if rss()>3*1024**3:raise BudgetStop('Heun3GB RSS cap')
  key='|'.join(map(str,owner));arm=owners.setdefault(key,dict(RHS=0,states=0,LU=0))
  result=original(event,**kw)
  for name,arg in [('RHS','rhs'),('states','states'),('LU','lu')]:costs[name]+=kw.get(arg,0);arm[name]+=kw.get(arg,0)
  return result
 book._update=bounded
 status='COMPLETE';error=None;started=time.perf_counter();summaries=[];teacher_rows=[]
 try:
  for parent in PARENTS:
   owner[:]=[parent,'common_anchor'];anchor=target.load_cache(OUT/'results/contexts/cache'/f'{parent}.npz','R4_StageI_common_context')
   y=target.observation(parent);phi0=target.phi(anchor,y)
   for dimension in (64,66):
    chart=Chart(anchor.q,dimension);x0=chart.x0
    for rank in (2,4,8):
     name=f'd{dimension}_r{rank}';dest=RESULT/parent/name;dest.mkdir(parents=True)
     before=costs.copy();case_started=time.perf_counter();firstfolder=FIRST/parent/name
     with np.load(firstfolder/'BASE_FRAME.npz',allow_pickle=False) as a:base={key:a[key].copy() for key in a.files}
     if not np.array_equal(base['q'],anchor.q) or not np.allclose(base['x'],x0,rtol=0,atol=1e-14):raise ValueError('firstpilot base mismatch')
     with np.load(firstfolder/'BRACKET_h0.001.npz',allow_pickle=False) as a:B=a['bracket'].copy()
     D0=base['Dchi'];M0=material_real(D0);ref=M0@B;amplitude=float(np.linalg.norm(ref));seeds=base['seeds'];X0,Y0=base['X'],base['Y']
     def fields(frame,D):return (projector_field(frame,D,seeds[:,0]),projector_field(frame,D,seeds[:,1]))
     def evaluate(x,phase,label,need_frame=True):
      owner[:]=[parent,name,phase]
      if np.linalg.norm(x-x0)>8*(1+1e-12):raise GeometryInvalid('LATENT_DISPLACEMENT_CAP_BEFORE_PHYSICS')
      point=target.evaluate(chart.native(x),purpose=f'{parent}_{name}_{phase}_{label}')
      phi=target.phi(point,y)
      row=dict(parent=parent,dimension=dimension,rank=rank,phase=phase,label=label,endpoint=point.key,
       Phi=phi,DeltaPhi=phi-phi0,prior_delta=chart.prior_energy(x)-chart.prior_energy(x0),
       actual_RMS=target.prior.material_rms(point.state.chi,anchor.state.chi),
       prior_legal=bool(target.prior.legal(point.q)),backward_residual=float(point.state.residual),
       **phase_diagnostic(target.prior,point.q,anchor.q))
      save(dest/'coordinates'/f'{point.key}.npz',x=x,q=point.q)
      if need_frame:
       Jq=target.jacobian(point);D=chart.derivative(target.prior,x);J=chart.jacobian(Jq,x)
       save(dest/'derivatives'/f'{point.key}.npz',x=x,q=point.q,Jq=Jq,Jcoordinate=J)
       frame=material_frame(J,D,rank)
       row.update(material_rank=frame['material_rank'],relative_gap=frame['relative_gap'],near_crossing=frame['near_crossing'])
       append(dest/'ENDPOINT_ROWS.jsonl',row)
       return point,row,fields(frame,D)
      append(dest/'ENDPOINT_ROWS.jsonl',row);return point,row,None
     def curve(step,orientation,phase,label):
      sequence=((0,1),(1,1),(0,-1),(1,-1)) if orientation=='forward' else ((1,1),(0,1),(1,-1),(0,-1))
      x=x0.copy();current_fields=(X0,Y0);path=[];curve_before=costs.copy();tick=time.perf_counter()
      for segment,(index,sign) in enumerate(sequence):
       v=current_fields[index];predictor=x+sign*step*v
       p,row,pfields=evaluate(predictor,phase,f'{label}_s{segment}_predictor',True)
       row['stage_role']='Heun_predictor_field_evaluation';path.append(row);del p
       corrected=heun_corrector(x,v,pfields[index],sign*step)
       p,row,cfields=evaluate(corrected,phase,f'{label}_s{segment}_corrected',segment<3)
       row['stage_role']='Heun_corrected_endpoint';path.append(row)
       x=corrected
       if segment<3:current_fields=cfields;del p
      delta=p.state.chi-anchor.state.chi;material=np.r_[delta.real,delta.imag]/24
      sign_ref=1 if orientation=='forward' else -1
      scaled=material/(step*step)
      expected=sign_ref*ref
      out=dict(orientation=orientation,step=step,endpoint=p.key,x=x.copy(),
       actual_RMS=float(np.linalg.norm(material)),actual_endpoint_normal=finite_material_normal(delta,D0,base),
       signed_material_progress=float(material@(expected/max(amplitude,1e-30))),
       Phi=row['Phi'],DeltaPhi=row['DeltaPhi'],prior_delta=row['prior_delta'],
       prior_legal=row['prior_legal'],backward_residual=row['backward_residual'],
       relative_material_bracket_error=float(np.linalg.norm(scaled-expected)/max(amplitude,1e-30)),
       scaled_material_vector=scaled,
       relative_latent_bracket_error=float(np.linalg.norm((x-x0)/(step*step)-sign_ref*B)/max(np.linalg.norm(B),1e-30)),
       any_rank_change=any(z.get('material_rank',base['material_rank'])!=int(base['material_rank']) for z in path),
       any_near_crossing=bool(base['near_crossing']) or any(z.get('near_crossing',False) for z in path),
       any_seam_crossing=any(z['seam_crossed_cells']>0 for z in path),
       stage_evaluations=len(path),stage_endpoint_keys=[z['endpoint'] for z in path],
       maximum_stage_Phi=max(z['Phi'] for z in path),maximum_stage_prior_delta=max(z['prior_delta'] for z in path),
       intermediate_stage_points_are_continuous_path_certificate=False,
       actual_cost={key:costs[key]-curve_before[key] for key in costs},wall_seconds=time.perf_counter()-tick)
      save(dest/'curve_vectors'/f'{label}_{orientation}.npz',x_end=x,scaled_material_vector=scaled,reference_material_vector=expected)
      del p
      return out
     small={};failures=[]
     for orientation in ('forward','reverse'):
      for h in HS:
       try:
        row=curve(h,orientation,'small_h_'+orientation,f'h{h}_{orientation}')
        small[(orientation,h)]=row
        append(dest/'SMALL_H_ROWS.jsonl',{key:value for key,value in row.items() if key not in ('x','scaled_material_vector')})
       except GeometryInvalid as exc:
        fail=dict(phase='small_h',orientation=orientation,h=h,reason=str(exc));failures.append(fail);append(dest/'SKIPS.jsonl',fail)
     gate_reasons=[];convergence=[]
     if len(small)!=6:gate_reasons.append('INCOMPLETE_PREDECLARED_SMALL_H_CURVES')
     for orientation in ('forward','reverse'):
      if all((orientation,h) in small for h in HS):
       fine,mid,coarse=[small[(orientation,h)] for h in HS]
       first=float(np.linalg.norm(fine['scaled_material_vector']-mid['scaled_material_vector']))
       second=float(np.linalg.norm(mid['scaled_material_vector']-coarse['scaled_material_vector']))
       contraction=first<=max(1e-8,1.25*second)
       convergence.append(dict(orientation=orientation,fine_h_error=fine['relative_material_bracket_error'],fine_mid_change=first,mid_coarse_change=second,change_contraction_pass=contraction))
       if fine['relative_material_bracket_error']>.2:gate_reasons.append(orientation+'_FINE_MATERIAL_ERROR_GT_0.2')
       if not contraction:gate_reasons.append(orientation+'_NO_EMPIRICAL_H_CONTRACTION')
     if any(z['any_rank_change'] for z in small.values()):gate_reasons.append('MATERIAL_RANK_CHANGE')
     if any(z['any_near_crossing'] for z in small.values()):gate_reasons.append('NEAR_CLUSTER_CROSSING')
     gate_pass=not gate_reasons
     gate=dict(parent=parent,dimension=dimension,rank=rank,passed=gate_pass,reasons=gate_reasons,convergence=convergence,
       any_small_h_seam=any(z['any_seam_crossing'] for z in small.values()),
       no_global_smooth_order_claim=True,task_scale_will_run=gate_pass)
     write(dest/'SMALL_H_GATE.json',gate)
     tasks=[]
     if gate_pass:
      for radius in RADII:
       for orientation in ('forward','reverse'):
        sign=1 if orientation=='forward' else -1
        try:
         x,match=distance_match(target.prior,chart,anchor.state.chi,x0,sign*B,radius)
         direct_before=costs.copy();tick=time.perf_counter()
         p,row,_=evaluate(x,'direct_'+orientation,f'r{radius}_{orientation}',False)
         delta=p.state.chi-anchor.state.chi;material=np.r_[delta.real,delta.imag]/24
         direct=dict(arm='direct_'+orientation,radius=radius,endpoint=p.key,actual_RMS=row['actual_RMS'],
          radius_matched=match['radius_matched'],signed_material_progress=float(material@(sign*ref/max(amplitude,1e-30))),
          Phi=row['Phi'],DeltaPhi=row['DeltaPhi'],prior_delta=row['prior_delta'],prior_legal=row['prior_legal'],backward_residual=row['backward_residual'],
          actual_endpoint_normal=finite_material_normal(delta,D0,base),
          actual_cost={key:costs[key]-direct_before[key] for key in costs},wall_seconds=time.perf_counter()-tick)
         tasks.append(direct);del p
        except GeometryInvalid as exc:append(dest/'SKIPS.jsonl',dict(phase='task_direct',orientation=orientation,radius=radius,reason=str(exc)))
        try:
         if amplitude<=1e-12:raise GeometryInvalid('NO_OP_REFERENCE_BRACKET')
         row=curve(float(np.sqrt(radius/amplitude)),orientation,'task_heun_'+orientation,f'r{radius}_{orientation}')
         row={key:value for key,value in row.items() if key not in ('x','scaled_material_vector')}
         row.update(arm='Heun_'+orientation,radius=radius,radius_matched=bool(abs(row['actual_RMS']/radius-1)<=.01),actual_RMS_relative_target=row['actual_RMS']/radius)
         tasks.append(row)
        except GeometryInvalid as exc:append(dest/'SKIPS.jsonl',dict(phase='task_Heun',orientation=orientation,radius=radius,reason=str(exc)))
      for row in tasks:
       row['teacher_eligible']=bool(abs(row['actual_RMS']/row['radius']-1)<=.01 and row['signed_material_progress']>=.5*row['radius'] and
           row['DeltaPhi']<=2 and row['prior_delta']<=2 and row['prior_legal'] and row['backward_residual']<=1e-10)
       append(dest/'TASK_ROWS.jsonl',row)
       if row['teacher_eligible']:teacher_rows.append(dict(parent=parent,dimension=dimension,rank=rank,**row))
     else:
      append(dest/'SKIPS.jsonl',dict(phase='all_task_scales',reason='PREDECLARED_SMALL_H_GATE_FAILED',radii=RADII,gate_reasons=gate_reasons))
     summary=dict(parent=parent,dimension=dimension,rank=rank,gate=gate,completed_small_h_curves=len(small),completed_task_endpoints=len(tasks),
        eligible_teachers=sum(z['teacher_eligible'] for z in tasks),reference_material_bracket_amplitude=amplitude,
        previously_paid_shared_FD_RHS=480,actual_new_cost={key:costs[key]-before[key] for key in costs},wall_seconds=time.perf_counter()-case_started)
     write(dest/'SUMMARY.json',summary);summaries.append(summary)
     write(RESULT/'PROGRESS.json',dict(completed_configurations=len(summaries),passed_gates=sum(z['gate']['passed'] for z in summaries),last=summary,costs=costs))
     print(json.dumps(plain(dict(parent=parent,dimension=dimension,rank=rank,gate=gate_pass,reasons=gate_reasons,task_endpoints=len(tasks),total_RHS=costs['RHS']))),flush=True)
     gc.collect()
   del anchor;target.backend()._transfers.clear();gc.collect()
  write(RESULT/'ELIGIBLE_LOCAL_TEACHERS.json',dict(qualification='LOCAL_COMPUTATIONAL_TEACHERS_NOT_POSTERIOR_SAMPLES',frozen_criteria=protocol()['teacher_eligible'],rows=teacher_rows))
  write(RESULT/'RESULTS.json',dict(status='COMPLETE',configurations=summaries,passed_gates=sum(z['gate']['passed'] for z in summaries),failed_gates=sum(not z['gate']['passed'] for z in summaries),
     new_costs=costs,arm_costs=owners,previously_paid_selected_shared_FD_RHS=23040,
     eligible_local_teacher_count=len(teacher_rows),controlled_elapsed_at_summary=time.perf_counter()-started,RSS_peak_bytes=rss(),
     no_GT=True,posterior_samples=False,continuous_path_certificate=False,Heun_map_exactly_reversible=False,
     no_global_O_h3_claim_at_C1_seams=True,scientific_review='PARENT_OWNER'))
 except BaseException as exc:
  status='STOPPED' if isinstance(exc,BudgetStop) else 'ERROR';error=repr(exc)
  write(RESULT/'ERROR.json',dict(status=status,error=error,traceback=traceback.format_exc(),new_costs=costs))
  raise
 finally:
  book.close(status,error)
  write(RESULT/'FINAL_JOB_COST.json',dict(status=status,error=error,new_costs=costs,total_wall_seconds=time.perf_counter()-started,RSS_peak_bytes=rss(),receipt=str(book.jobdir/'CLOSED.json')))


def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--prepare',action='store_true');parser.add_argument('--run',action='store_true');args=parser.parse_args()
 if args.prepare:
  if PROTOCOL.exists():raise FileExistsError('frozen Heun protocol exists')
  write(PROTOCOL,protocol());print(str(PROTOCOL));return
 if args.run:run();return
 parser.error('choose --prepare or --run')

if __name__=='__main__':main()
