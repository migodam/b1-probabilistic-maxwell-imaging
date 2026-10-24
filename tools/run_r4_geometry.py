#!/usr/bin/env python3
"""Frozen R4 StageIII local geometry experiment; preparation imports no backend."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[name]='1'
from pathlib import Path
import sys,argparse,json,time,traceback,gc,resource
from types import SimpleNamespace
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from p3r4.core import OUT,Book,Target,BudgetStop,save,write,plain
from p3r4.dynamic_geometry import (Chart,GeometryInvalid,material_frame,projector_field,
    central_bracket,material_normal,rotation,distance_match,phase_diagnostic,finite_material_normal)
from p3r2.geometry import strong_correct

JOB='dynamic_geometry_first'
PROTOCOL=OUT/'configs/GEOMETRY_PROTOCOL.json'
RESULT=OUT/'results/III'/JOB
PARENTS=tuple(f'R4DEV_{i:05d}' for i in range(32))
SEED=20261025030
RANKS=(2,4,8)
HS=(.001,.002,.004)
RADII=(.02,.04,.08)
CAP=800000
STATE_CAP=96
ARMS=('finite_difference','circuit_forward','circuit_reverse','direct_bracket_forward','direct_bracket_reverse',
    'frozen_weak','dynamic_weak','ordinary_field_curvature','ordinary_field_secant','strong_compensated')


def protocol():
    return dict(campaign='P3LR4-20261025',stage='III',job=JOB,parents=PARENTS,seed=SEED,
        release_marker='results/DEV_CONTEXT_COMPLETE.json',coordinate_dimensions=[64,66],ranks=RANKS,
        frame='realified decoder D/sqrt576 SVD cutoff max(smax*1e-8,1e-14); material-whiten; eig(T.T(J.TJ+I)T); Euclidean projector of selected span',
        fields='fixed independent Gaussian ci per parent/dimension/rank from SeedSequence; P(x)ci normalized by actual decoder material RMS; one pair only',
        central_FD_h=HS,FD_states_per_h=4,reference_FD_h=.002,
        bracket='DY X-DX Y including material normalization; retain all h, no best-h selection',
        radii=RADII,stress_radius=.08,arms=ARMS,
        circuit='four explicit Euler steps +X,+Y,-X,-Y; reverse +Y,+X,-Y,-X; not exact inverses or MH maps',
        circuit_step='s=sqrt(requested_RMS/||D bracket_h=.002||); actual endpoint RMS may differ; no adaptive radius rescue',
        direct='decoder-only first dyadic distance root, latent displacement cap8, no redraw',
        dynamic_weak='two segments: matched radius/2 frozen X then refreshed X to matched radius',
        curvature_control='central self derivative D_X X at h=.002',
        secant_control='one-sided (Y(x+hX)-Y0)/h at h=.002',
        strong='existing R2 strong_correct, only fixed-u64; material radius/progress preserved, latent cap2, maxiter20;66 NOT_RUN',
        states_cap_per_parent_per_arm=STATE_CAP,arm_caps_shared_across_all_dimensions_ranks_radii=True,
        max_new_complex_RHS=CAP,max_controlled_wall_seconds=3600,RSS_cap_bytes=3*1024**3,
        no_GT=True,GT_filter=False,all_endpoints_persist_current_field=True,
        frame_gap_threshold_diagnostic=1e-3,material_field_noop_cutoff=1e-10,
        literal_noiseball_white_radius=23.785449004057032,
        noiseball_role='historical diagnostic threshold, not posterior credible region',
        continuous_path_certificate=False,posterior_samples=False,full_Lie_rank_claim=False,
        cost='charge all Full/J/failed paths; shared FD cost explicitly retained for every-arm cold comparisons; report unequal marginal costs',
        prior_context_owner='R4 StageI paid native/probit 96+64 contexts; no stationarity claim',
        expected_RHS_upper_bound_before_skips=293632,
        stop='invalid material rank, no-op, missing finite probit coordinate, cannot match radius within cap; retain failures; no physics retuning')


def append(path,row):
    with Path(path).open('a') as f:f.write(json.dumps(plain(row),allow_nan=False)+'\n');f.flush();os.fsync(f.fileno())


def rss():
    value=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return value if sys.platform=='darwin' else value*1024


def run():
    frozen=json.loads(PROTOCOL.read_text())
    if frozen!=plain(protocol()):raise ValueError('frozen geometry protocol mismatch')
    release=json.loads((OUT/'results/DEV_CONTEXT_COMPLETE.json').read_text())
    if release.get('status')!='COMPLETE' or release.get('parents')!=128:raise ValueError('context worker release marker missing')
    for parent in PARENTS:
        for path in (OUT/'results/contexts/cache'/f'{parent}.npz',OUT/'inputs/online'/f'{parent}.npz'):
            if not path.is_file():raise FileNotFoundError(path)
    if RESULT.exists():raise FileExistsError('existing geometry run; never reset/replay')
    RESULT.mkdir(parents=True)
    write(RESULT/'FROZEN_PROTOCOL.json',frozen)
    book=Book('III',JOB,max_wall=3600);target=Target(book,keep_fields=True)
    counters=dict(RHS=0,states=0,LU=0);arm_counts={};owner=[None,None];original=book._update
    def bounded(event,**kw):
        if counters['RHS']+kw.get('rhs',0)>CAP:raise BudgetStop('geometry job RHS cap')
        if rss()>3*1024**3:raise BudgetStop('geometry RSS cap')
        key='|'.join(str(v) for v in owner)
        arm=arm_counts.setdefault(key,dict(RHS=0,states=0,LU=0))
        if arm['states']+kw.get('states',0)>STATE_CAP:raise BudgetStop('per-parent-arm geometry state cap')
        result=original(event,**kw)
        for name,arg in [('RHS','rhs'),('states','states'),('LU','lu')]:
            counters[name]+=kw.get(arg,0);arm[name]+=kw.get(arg,0)
        return result
    book._update=bounded
    before_job=time.perf_counter();status='COMPLETE';error=None;config_summaries=[];all_endpoint_rows=[];fd_rows=[]
    try:
      for parent_index,parent in enumerate(PARENTS):
        folder=RESULT/parent;folder.mkdir()
        owner[:]=[parent,'common_anchor']
        anchor=target.load_cache(OUT/'results/contexts/cache'/f'{parent}.npz','R4_StageI_context_common_96native_64probit')
        y=target.observation(parent)
        Jq0=target.jacobian(anchor);phi0=target.phi(anchor,y)
        for dimension in (64,66):
          try:chart=Chart(anchor.q,dimension)
          except GeometryInvalid as exc:
            append(folder/'SKIPS.jsonl',dict(dimension=dimension,reason=str(exc),paid=False));continue
          x0=chart.x0;D0=chart.derivative(target.prior,x0);J0=chart.jacobian(Jq0,x0)
          for rank in RANKS:
            name=f'd{dimension}_r{rank}';dest=folder/name;dest.mkdir()
            cost_start=counters.copy();case_started=time.perf_counter()
            config=dict(parent=parent,dimension=dimension,rank=rank,status='COMPLETE')
            try:
              frame0=material_frame(J0,D0,rank)
              rng=np.random.default_rng(np.random.SeedSequence([SEED,parent_index,dimension,rank]))
              seeds=rng.normal(size=(dimension,2));X0=projector_field(frame0,D0,seeds[:,0]);Y0=projector_field(frame0,D0,seeds[:,1])
              save(dest/'BASE_FRAME.npz',x=x0,q=anchor.q,J=J0,Dchi=D0,seeds=seeds,X=X0,Y=Y0,**frame0)
              def fields(frame,D):return projector_field(frame,D,seeds[:,0]),projector_field(frame,D,seeds[:,1])
              def evaluate(x,arm,label,need_frame=False):
                owner[:]=[parent,arm]
                if np.linalg.norm(x-x0)>8*(1+1e-12):raise GeometryInvalid('LATENT_DISPLACEMENT_CAP_BEFORE_PHYSICS')
                point=target.evaluate(chart.native(x),purpose=f'{parent}_{name}_{arm}_{label}')
                q=point.q;chi=point.state.chi;phi=target.phi(point,y)
                delta=chi-anchor.state.chi
                row=dict(parent=parent,dimension=dimension,rank=rank,arm=arm,label=label,endpoint=point.key,
                    Phi=phi,DeltaPhi=phi-phi0,white_residual=float(np.sqrt(2*phi)),
                    prior_delta=chart.prior_energy(x)-chart.prior_energy(x0),
                    material_RMS=target.prior.material_rms(chi,anchor.state.chi),latent_displacement=float(np.linalg.norm(x-x0)),
                    prior_legal=bool(target.prior.legal(q)),backward_residual=float(point.state.residual),
                    **phase_diagnostic(target.prior,q,anchor.q))
                append(dest/'ENDPOINT_ROWS.jsonl',row);all_endpoint_rows.append(row)
                save(dest/'coordinates'/f'{point.key}.npz',x=x,q=q)
                if need_frame:
                    Jq=target.jacobian(point);D=chart.derivative(target.prior,x);J=chart.jacobian(Jq,x)
                    save(dest/'derivatives'/f'{point.key}.npz',x=x,q=q,Jq=Jq,Jcoordinate=J)
                    frame=material_frame(J,D,rank)
                    return point,row,frame,D
                return point,row
              probes={};brackets={}
              for h in HS:
                around={};meta=[]
                for axis,v in (('X',X0),('Y',Y0)):
                  for sign in (-1,1):
                    x=x0+sign*h*v
                    point,row,frame,D=evaluate(x,'finite_difference',f'h{h}_{axis}{sign:+}',True)
                    Fx,Fy=fields(frame,D)
                    around[(axis,sign)]=dict(X=Fx,Y=Fy)
                    meta.append(dict(axis=axis,sign=sign,endpoint=point.key,rotation=rotation(frame,frame0),phase=phase_diagnostic(target.prior,point.q,anchor.q)))
                    del point,frame,D
                bracket=central_bracket(around[('Y',1)]['X'],around[('Y',-1)]['X'],around[('X',1)]['Y'],around[('X',-1)]['Y'],h)
                curvature=(around[('X',1)]['X']-around[('X',-1)]['X'])/(2*h)
                secant=(around[('X',1)]['Y']-Y0)/h
                brackets[h]=bracket;probes[h]=dict(curvature=curvature,secant=secant)
                diag=material_normal(bracket,D0,frame0)
                fd=dict(parent=parent,dimension=dimension,rank=rank,h=h,**diag,probe_metadata=meta,
                    base_material_rank=frame0['material_rank'],base_relative_gap=frame0['relative_gap'],
                    max_projector_rotation=max(m['rotation']['projector_operator'] for m in meta),
                    any_rank_change=any(m['rotation']['material_rank']!=frame0['material_rank'] for m in meta),
                    any_near_crossing=frame0['near_crossing'] or any(m['rotation']['near_crossing'] for m in meta),
                    any_seam_crossing=any(m['phase']['seam_crossed_cells']>0 for m in meta),
                    full_Lie_rank_proven=False,continuous_path_verified=False)
                append(dest/'BRACKET_ROWS.jsonl',fd);fd_rows.append(fd)
                save(dest/f'BRACKET_h{h}.npz',bracket=bracket,curvature=curvature,secant=secant)
              write(dest/'FD_CONVERGENCE.json',[dict(h1=a,h2=b,latent_difference=float(np.linalg.norm(brackets[a]-brackets[b])),
                  relative_difference=float(np.linalg.norm(brackets[a]-brackets[b])/max(np.linalg.norm(brackets[b]),1e-30)),
                  material_difference=float(np.linalg.norm(np.r_[D0.real,D0.imag]@(brackets[a]-brackets[b]))/24))
                  for a,b in zip(HS[:-1],HS[1:])])
              B=brackets[.002];amplitude=material_normal(B,D0,frame0)['material_amplitude']
              arm_rows=[]
              def endpoint_summary(arm,radius,point,row,match,extra=None):
                data=dict(arm=arm,radius=radius,stress=radius==.08,endpoint=point.key,
                    actual_RMS=row['material_RMS'],radius_matched=bool(abs(row['material_RMS']-radius)<=max(1e-8,radius*.01)) if match else False,
                    Phi=row['Phi'],DeltaPhi=row['DeltaPhi'],prior_delta=row['prior_delta'],
                    seam_crossed_cells=row['seam_crossed_cells'],continuous_path_verified=False,
                    actual_endpoint_normal=finite_material_normal(point.state.chi-anchor.state.chi,D0,frame0))
                if extra:data.update(extra)
                append(dest/'ARM_ROWS.jsonl',data);arm_rows.append(data)
              for radius in RADII:
                predictor=None
                direct=[('direct_bracket_forward',B),('direct_bracket_reverse',-B),('frozen_weak',X0),
                    ('ordinary_field_curvature',probes[.002]['curvature']),('ordinary_field_secant',probes[.002]['secant'])]
                for arm,vector in direct:
                  try:
                    before=counters.copy();tick=time.perf_counter()
                    x,match=distance_match(target.prior,chart,anchor.state.chi,x0,vector,radius)
                    point,row=evaluate(x,arm,f'radius{radius}')
                    if arm=='frozen_weak':predictor=SimpleNamespace(q=point.q.copy())
                    endpoint_summary(arm,radius,point,row,True,dict(decoder_match=match,
                        marginal_cost={k:counters[k]-before[k] for k in counters},wall_seconds=time.perf_counter()-tick))
                    del point
                  except GeometryInvalid as exc:append(dest/'SKIPS.jsonl',dict(arm=arm,radius=radius,reason=str(exc)))
                # Two-segment refreshed weak field; midpoint and endpoint are paid.
                try:
                  arm='dynamic_weak';before=counters.copy();tick=time.perf_counter()
                  xmid,_=distance_match(target.prior,chart,anchor.state.chi,x0,X0,radius/2)
                  mid,mrow,fmid,Dmid=evaluate(xmid,arm,f'radius{radius}_mid',True)
                  Xmid,_=fields(fmid,Dmid)
                  x,match=distance_match(target.prior,chart,anchor.state.chi,xmid,Xmid,radius)
                  point,row=evaluate(x,arm,f'radius{radius}_end')
                  endpoint_summary(arm,radius,point,row,True,dict(decoder_match=match,
                      midpoint=mrow,marginal_cost={k:counters[k]-before[k] for k in counters},wall_seconds=time.perf_counter()-tick))
                  del mid,point,fmid,Dmid
                except GeometryInvalid as exc:append(dest/'SKIPS.jsonl',dict(arm='dynamic_weak',radius=radius,reason=str(exc)))
                # Explicit Euler circuits, each direction owns a separate <=96 state arm.
                for arm,sequence in (('circuit_forward',((0,1),(1,1),(0,-1),(1,-1))),
                                     ('circuit_reverse',((1,1),(0,1),(1,-1),(0,-1)))):
                  before=counters.copy();tick=time.perf_counter();path=[]
                  try:
                    if amplitude<=1e-12:raise GeometryInvalid('NO_OP_BRACKET_MATERIAL_AMPLITUDE')
                    step=float(np.sqrt(radius/amplitude));x=x0.copy();Fx,Fy=X0,Y0
                    for segment,(field_index,sign) in enumerate(sequence):
                      x=x+sign*step*(Fx if field_index==0 else Fy)
                      if segment<3:
                        point,row,frame,D=evaluate(x,arm,f'radius{radius}_segment{segment}',True)
                        Fx,Fy=fields(frame,D);del frame,D
                      else:point,row=evaluate(x,arm,f'radius{radius}_segment{segment}')
                      path.append(row)
                      if segment<3:del point
                    endpoint_summary(arm,radius,point,row,False,dict(euler_step=step,
                        theoretical_bracket_RMS=radius,actual_RMS_relative_target=row['material_RMS']/radius,
                        numerical_map_is_exact_reverse=False,ideal_flow_Jacobian_used=False,
                        path_sample_count=4,max_sampled_Phi=max(z['Phi'] for z in path),
                        max_sampled_prior_delta=max(z['prior_delta'] for z in path),
                        sampled_outside_reference_ball=sum(z['white_residual']>23.785449004057032 for z in path),
                        sampled_seam_crossing=any(z['seam_crossed_cells']>0 for z in path),
                        intermediate_endpoints=[z['endpoint'] for z in path],
                        marginal_cost={k:counters[k]-before[k] for k in counters},wall_seconds=time.perf_counter()-tick))
                    del point
                  except GeometryInvalid as exc:append(dest/'SKIPS.jsonl',dict(arm=arm,radius=radius,reason=str(exc),
                      partial_paid_endpoints=[z['endpoint'] for z in path],marginal_cost={k:counters[k]-before[k] for k in counters}))
                if dimension==64 and predictor is not None:
                  arm='strong_compensated';before=counters.copy();tick=time.perf_counter()
                  try:
                    velocity=target.prior.derivative(anchor.q)[:,:64]@X0
                    q,fit=strong_correct(target,anchor,predictor,anchor,frame0,velocity,radius,y,maxiter=20)
                    if not fit['radius_and_progress_retained']:
                        append(dest/'SKIPS.jsonl',dict(arm=arm,radius=radius,reason='STRONG_CORRECTOR_CONSTRAINT_FAILURE',fit=fit));continue
                    point,row=evaluate(q[:64],arm,f'radius{radius}')
                    endpoint_summary(arm,radius,point,row,True,dict(corrector=fit,
                      marginal_cost={k:counters[k]-before[k] for k in counters},wall_seconds=time.perf_counter()-tick))
                    del point
                  except (ValueError,GeometryInvalid) as exc:append(dest/'SKIPS.jsonl',dict(arm=arm,radius=radius,reason=str(exc)))
                elif dimension==66:append(dest/'SKIPS.jsonl',dict(arm='strong_compensated',radius=radius,reason='NOT_RUN_EXISTING_HELPER_FIXED_U64_ONLY'))
              config.update(bracket_reference=material_normal(B,D0,frame0),completed_arm_endpoints=len(arm_rows),
                  FD_shared_setup_RHS=480,FD_shared_setup_states=12,
                  arm_rows=arm_rows,material_rank=frame0['material_rank'],relative_gap=frame0['relative_gap'],near_crossing=frame0['near_crossing'])
            except GeometryInvalid as exc:
              config.update(status='SKIPPED_OR_PARTIAL_INVALID_GEOMETRY',reason=str(exc))
              append(dest/'SKIPS.jsonl',dict(phase='configuration',reason=str(exc)))
            config.update(actual_new_costs={k:counters[k]-cost_start[k] for k in counters},wall_seconds=time.perf_counter()-case_started)
            write(dest/'SUMMARY.json',config);config_summaries.append(config)
            write(RESULT/'PROGRESS.json',dict(completed_configurations=len(config_summaries),last={k:v for k,v in config.items() if k!='arm_rows'},costs=counters))
            print(json.dumps(plain(dict(parent=parent,dimension=dimension,rank=rank,status=config['status'],cost=config['actual_new_costs'],total_RHS=counters['RHS']))),flush=True)
            gc.collect()
        del anchor,Jq0
        target.backend()._transfers.clear();gc.collect()
      write(RESULT/'RESULTS.json',dict(status='COMPLETE',configurations=config_summaries,endpoint_rows=len(all_endpoint_rows),bracket_rows=len(fd_rows),
          new_costs=counters,arm_costs=arm_counts,controlled_elapsed_at_summary=time.perf_counter()-before_job,RSS_peak_bytes=rss(),
          no_GT=True,posterior_samples=False,continuous_path_verified=False,full_Lie_rank_claim=False,
          scientific_review='PARENT_OWNER; compare absolute material amplitudes, finite path/data/prior costs and actual arm RHS/wall'))
    except BaseException as exc:
      status='STOPPED' if isinstance(exc,BudgetStop) else 'ERROR';error=repr(exc)
      write(RESULT/'ERROR.json',dict(status=status,error=error,traceback=traceback.format_exc(),new_costs=counters,arm_costs=arm_counts))
      raise
    finally:
      book.close(status,error)
      write(RESULT/'FINAL_JOB_COST.json',dict(status=status,error=error,new_costs=counters,arm_costs=arm_counts,
        total_wall_seconds=time.perf_counter()-before_job,RSS_peak_bytes=rss(),receipt=str(book.jobdir/'CLOSED.json')))


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--prepare',action='store_true');parser.add_argument('--run',action='store_true');args=parser.parse_args()
    if args.prepare:
        if PROTOCOL.exists():raise FileExistsError('frozen protocol exists; inspect rather than overwrite')
        write(PROTOCOL,protocol());print(str(PROTOCOL));return
    if args.run:run();return
    parser.error('choose --prepare or --run; neither default starts physics')

if __name__=='__main__':main()
