"""P0-C fixed-pair full-Maxwell phase table and distorted-wave Born comparison.

Known-reference physical existence only. Import performs no physics. The main
entry requires an independently qualified *joint* witness before opening C's
budget; no B2 allowance or online-data subtraction is used.
"""
from __future__ import annotations

import argparse
import csv
import fcntl
import json
import time
import uuid
from pathlib import Path

import numpy as np
from p0.geometry import probes

from .accounting import ROOT, Book, BudgetStop, write, encode
from .native_world import LinearWorld
from .physical_witness import mesh_for
from .physics import Forward, coefficient_jacobian


IDENTITY = 'P0_C_FIXED_VALIDATED_PAIR_KNOWN_REFERENCE_NOT_ONLINE'
FINE = dict(h=.13,tol=2.5e-5,res=4,n=512)
TX_SUBSET = [1,2]
RX_SUBSET = list(range(3,9))


def _rows(book):
    if not book.path.exists():
        return []
    with book.path.open() as stream:
        fcntl.flock(stream,fcntl.LOCK_SH)
        return [json.loads(line) for line in stream if line.strip()]


def _csv(path, rows):
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _distance(a,b):
    difference=np.asarray(a)-np.asarray(b)
    full=float(np.linalg.norm(difference)/.01)
    central=float(np.linalg.norm(difference[:,TX_SUBSET,:][:,:,RX_SUBSET])/.01)
    if central>full+1e-10*max(1.,full):
        raise RuntimeError('strict nested-aperture norm exceeded full-aperture norm')
    return dict(d_mean_all=full,d_mean_central_nested=central,
                complex_channels_all=int(difference.size),
                complex_channels_central=int(difference.shape[0]*len(TX_SUBSET)*len(RX_SUBSET)))


class _WallLossWorld(LinearWorld):
    """Scale only wall conductivity; inherited probe remains exactly sigma=.02."""
    def __init__(self,seed,probe=None,loss_scale=1.):
        if loss_scale not in (.5,1.,1.5):
            raise ValueError('only fixed common wall-loss factors .5, 1, 1.5 are allowed')
        super().__init__(seed,probe=probe)
        self.loss_scale=float(loss_scale)

    def components(self,x,pts):
        re,loss,b,c0,sigma0,a,c=super().components(x,pts)
        return re,loss*self.loss_scale,b,c0,sigma0*self.loss_scale,a,c

    def material_metadata(self,x):
        row=super().material_metadata(x)
        row['modulated_wall_bounds']['conductivity_max']=.063*self.loss_scale
        row.update(common_wall_loss_scale=self.loss_scale,
                   wall_conductivity_bound=.063*self.loss_scale,
                   probe_conductivity_unchanged=.02,
                   purpose='fixed-pair sensitivity; no parameter reoptimization')
        return row


class _ProbeInsertion:
    """Two affine target parameters about a wall-only distorted-wave background."""
    n=2

    def __init__(self,background,q,probe):
        self.background,self.q,self.probe=background,np.asarray(q).copy(),probe

    def _active(self,pts):
        center,radius,_=self.probe
        return np.linalg.norm(np.asarray(pts)-np.asarray(center),axis=1)<radius

    def coefficient(self,x,f):
        x=np.asarray(x,dtype=float)
        if x.shape!=(2,) or not np.all(np.isfinite(x)) or not 0<=x[0]<=8 or not 0<=x[1]<=.15:
            raise ValueError('target parameters outside original passive material bounds')
        background=self.background.coefficient(self.q,f)
        def eps(pts):
            values=background(pts)
            # Both fixed targets are in strict air for the wall reference.
            active=self._active(pts)
            if np.any(abs(values[active]-1.)>1e-14):
                raise RuntimeError('insertion circle is not air in the fixed wall background')
            values[active]+=x[0]+1j*x[1]/f
            return values
        return eps

    def derivative(self,x,pts,f):
        derivative=np.zeros((len(pts),2),complex)
        active=self._active(pts)
        derivative[active,0]=1.
        derivative[active,1]=1j/f
        return derivative

    def jacobian(self,p,s,x):
        return coefficient_jacobian(p,s,lambda pts,f:self.derivative(x,pts,f))

    def material_metadata(self,x):
        return dict(model='two_parameter_affine_target_insertion_about_fixed_wall',
                    identity=IDENTITY,wall_parameters=self.q,
                    target_center=self.probe[0],target_radius=self.probe[1],
                    x=np.asarray(x),parameters=['target_contrast','target_conductivity'],
                    background=self.background.material_metadata(self.q),
                    approximation='full wall background; first target derivative only',
                    not_online_data_subtraction=True)


class _LiveBook:
    """Mirror existing Book C events/costs to local files, without extra charges."""
    def __init__(self,book,out,receipt):
        self.book,self.out,self.receipt=book,out,receipt

    def __getattr__(self,name):
        return getattr(self.book,name)

    def _persist(self,row):
        with (self.out/'BACKEND_PROGRESS.jsonl').open('a') as stream:
            stream.write(json.dumps(row,default=encode,allow_nan=False)+'\n')
        write(self.out/'LIVE_RHS_CPU_WALL.json',self.receipt())

    def reserve(self,phase,label,rhs=0,factorizations=0,**metadata):
        try:
            reservation=self.book.reserve(phase,label,rhs,factorizations,**metadata)
        except Exception as exc:
            self._persist(dict(event='reservation_failed',phase=phase,label=label,
                               error_type=type(exc).__name__,error=str(exc)))
            raise
        self._persist(dict(event='reservation',id=reservation,phase=phase,label=label,
                           rhs=rhs,factorizations=factorizations))
        return reservation

    def event(self,**row):
        self.book.event(**row)
        self._persist(row)


def _validated_pair(manifest):
    path=Path(manifest)
    if not path.is_absolute():
        path=ROOT/path
    cfg=json.loads(path.read_text())
    if cfg['kind']!='joint' or cfg.get('chart','linear')!='linear' or cfg.get('cavities',[]):
        raise RuntimeError('P0-C requires a validated joint pair in the unchanged LinearWorld wall support')
    identifier=cfg['id']
    if not isinstance(identifier,str) or Path(identifier).name!=identifier or identifier in ('.','..'):
        raise ValueError('invalid validation identifier')
    validation_dir=ROOT/'results/B2/validation'/identifier
    gate=json.loads((validation_dir/'RESULTS.json').read_text())
    if gate.get('numerical_pass') is not True or gate.get('joint_physics_threshold') is not True:
        raise RuntimeError('P0-C blocked: independent numerical_pass AND joint_physics_threshold required')
    validated_manifest=json.loads((validation_dir/'MANIFEST.json').read_text())
    if cfg!=validated_manifest:
        raise RuntimeError('manifest differs from the independently validated manifest')
    with np.load(ROOT/cfg['x_file']) as archive:
        q=archive['x'].copy()
    model=LinearWorld(int(cfg['seed']))
    model._parameters(q)
    # Bind the current file contents to the actual endpoint that validation ran.
    endpoint=validation_dir/'fields/joint8_world1_fine/endpoint_0000/PARAMETERS.json'
    recorded=json.loads(endpoint.read_text())
    if not np.array_equal(np.asarray(recorded['x']),q):
        raise RuntimeError('manifest latent file differs from the validated fine endpoint')
    return cfg,q,path,gate,validation_dir


def main(manifest,include_loss=False):
    """Fixed pair only; include_loss adds common wall factors .5 and 1.5.

    No noise/gain change, fitting, witness reoptimization, or continuum claim.
    Narrow/wide frequency lists are not nested. Central aperture is strictly
    nested within each complete response, with identical .01 whitening.
    """
    cfg,q1,manifest_path,gate,validation_dir=_validated_pair(manifest)
    original=json.loads((ROOT/'configs/ORIGINAL_P0.json').read_text())
    if original['noise']['each_real_quadrature_std']!=.01:
        raise RuntimeError('original noise differs from fixed primary sigma=.01')
    book=Book('C')
    out=ROOT/'results/C/phase_diagram'
    out.mkdir(parents=True,exist_ok=False)
    (out/'responses').mkdir()
    (out/'born').mkdir()
    seed=int(cfg['seed'])
    qs=(np.zeros_like(q1),q1.copy())
    p0=probes()[8]
    p1=((p0[0][0]+.1,p0[0][1]),p0[1],p0[2])
    positions=(p0[0],p1[0])
    run_id=uuid.uuid4().hex
    prefix=f'C_phase_{run_id}_'
    cpu_start,wall_start=time.process_time(),time.perf_counter()
    ledger_start=len(_rows(book))
    phase_rows,born_rows,joint_born_rows=[],[],[]
    full_cache,mesh_cache={},{}
    failure=None
    response_counter=0
    born_results={}

    def receipt():
        rows=[r for r in _rows(book)[ledger_start:]
              if r.get('stage')=='C' and str(r.get('phase','')).startswith(prefix)]
        reserves=[r for r in rows if r.get('event')=='reserve']
        return dict(stage='C',new_charged_complex_rhs=sum(r.get('rhs',0) for r in reserves),
                    new_charged_factorizations=sum(r.get('factorizations',0) for r in reserves),
                    confirmed_backend_rhs=sum(r.get('rhs',0) for r in rows if r.get('event')=='backend_solve'),
                    reserve_ids=[r['id'] for r in reserves],
                    cpu_seconds=time.process_time()-cpu_start,
                    wall_seconds=time.perf_counter()-wall_start,
                    cpu_billing='external serial launcher charges process CPU once; measured here only',
                    ledger=str(book.path),B2_budget_borrowed=False)

    live=_LiveBook(book,out,receipt)
    scales=[1.,.5,1.5] if include_loss else [1.]
    write(out/'MANIFEST.json',cfg)
    write(out/'GATE_PROVENANCE.json',dict(manifest=str(manifest_path),validation=str(validation_dir),
          independently_validated_results=gate,gate_required=['numerical_pass','joint_physics_threshold'],
          fixed_validated_q_bound_to_endpoint=True))
    np.savez_compressed(out/'FIXED_PAIR.npz',q0=qs[0],q1=qs[1],target_centers=positions)
    write(out/'PLAN.json',dict(identity=IDENTITY,seed=seed,stage='C',gamma=1.,sigma=.01,
          frequencies=original['frequencies'],frequency_lists_nested=False,
          distance_definition='Euclidean norm of fixed complex mean-response difference / .01; no channel-count renormalization',
          fixed_pair_not_reoptimized=True,fine=FINE,max_dofs=110000,max_factor_bytes=3500000000,
          all_forward_streaming=True,common_wall_loss_scales=scales,
          wall_conductivity_bound_max=.063*max(scales),target_conductivity_unchanged=.02,
          target0=p0,target1=p1,joint_target_mass_relative_change=0.,
          nested_aperture=dict(tx=TX_SUBSET,rx=RX_SUBSET),
          phase_cases=['cal','joint_contrast3','fixed_target0_contrast3',
                       'fixed_target0_contrast7','joint_contrast7'],
          born=dict(band='narrow',wall_loss_scale=1.,worlds=[0,1],target_positions=[0,1],
                    contrasts=[3.,7.],conductivity=.02,background_parameters=[0.,0.],
                    expected_rhs_per_five_frequency_full_FJ=80,
                    formula='paid_background_data + paid_complex_J @ [contrast,.02]',
                    no_extra_full_solve_for_linear_prediction=True,
                    no_online_data_subtraction=True),
          min_feature='UNRESOLVED_ORIGINAL_DEFINITION',conditional_fixed_support_geometry=True,
          continuous_PDE_certification=False,probabilities=None))

    def mesh(probe):
        key=None if probe is None else (tuple(probe[0]),float(probe[1]))
        if key not in mesh_cache:
            mesh_cache[key]=mesh_for(seed,probe,**FINE)
        return mesh_cache[key]

    def full_response(band,world,probe,scale):
        nonlocal response_counter
        probe_key=None if probe is None else (tuple(probe[0]),float(probe[1]),float(probe[2]))
        key=(band,world,probe_key,scale)
        if key in full_cache:
            return full_cache[key]
        label=f'full_{response_counter:03d}'
        response_counter+=1
        model=_WallLossWorld(seed,probe=probe,loss_scale=scale)
        forward=Forward(model,mesh(probe),original['frequencies'][band],live,
                        f'C/phase_diagram/fields/{label}',phase=prefix+label,identity=IDENTITY,
                        h=.13,pml_strength=6,max_dofs=110000,max_factor_bytes=3500000000,
                        retain_current=False)
        data,_=forward.evaluate(qs[world],jac=False)
        path=out/'responses'/f'{label}.npz'
        np.savez_compressed(path,data=data,q=qs[world],frequencies=original['frequencies'][band],
                            wall_loss_scale=scale)
        write(path.with_suffix('.json'),dict(world=world,band=band,probe=probe,
              wall_loss_scale=scale,full_Maxwell=True,source_folder=str(forward.folder),
              actual_cumulative_cost=receipt(),identity=IDENTITY))
        full_cache[key]=(data.copy(),str(path))
        return full_cache[key]

    def phase_case(band,scale,case,contrast=None,joint=False):
        if contrast is None:
            a=b=None
        else:
            a=(positions[0],p0[1],contrast)
            b=(positions[1] if joint else positions[0],p0[1],contrast)
        d0,path0=full_response(band,0,a,scale)
        d1,path1=full_response(band,1,b,scale)
        row=dict(case=case,band=band,common_wall_loss_scale=scale,
                 target_contrast=contrast,target_shift=.1 if joint else 0.,
                 target_mass_relative_change=0. if joint else None,
                 **_distance(d1,d0),response0=path0,response1=path1,
                 fixed_pair_not_reoptimized=True,identity=IDENTITY)
        phase_rows.append(row)
        write(out/'PHASE_TABLE.json',phase_rows)
        _csv(out/'PHASE_TABLE.csv',phase_rows)
        write(out/'LIVE_RHS_CPU_WALL.json',receipt())

    try:
        # Primary calibration + complete joint pair first, in both bands.
        for band in ('narrow','wide'):
            phase_case(band,1.,'cal')
            phase_case(band,1.,'joint_contrast3',3.,True)
        for band in ('narrow','wide'):
            phase_case(band,1.,'fixed_target0_contrast3',3.)
            phase_case(band,1.,'fixed_target0_contrast7',7.)
            phase_case(band,1.,'joint_contrast7',7.,True)
        for scale in scales[1:]:
            for band in ('narrow','wide'):
                phase_case(band,scale,'cal')
                phase_case(band,scale,'joint_contrast3',3.,True)
                phase_case(band,scale,'fixed_target0_contrast3',3.)
                phase_case(band,scale,'fixed_target0_contrast7',7.)
                phase_case(band,scale,'joint_contrast7',7.,True)

        # Distorted-wave target Born: four independently paid full FJ backgrounds
        # (two fixed walls x two target circles). No historical live LU is loaded.
        band='narrow'
        frequencies=original['frequencies'][band]
        for world in (0,1):
            for target in (0,1):
                probe=(positions[target],p0[1],3.)
                label=f'born_w{world}_target{target}'
                background=_WallLossWorld(seed,loss_scale=1.)
                insertion=_ProbeInsertion(background,qs[world],probe)
                forward=Forward(insertion,mesh(probe),frequencies,live,
                                f'C/phase_diagram/fields/{label}',phase=prefix+label,
                                identity=IDENTITY,h=.13,pml_strength=6,max_dofs=110000,
                                max_factor_bytes=3500000000,retain_current=False)
                baseline,jacobian=forward.evaluate(np.zeros(2),jac=True)
                np.savez_compressed(out/'born'/f'{label}_BASE.npz',data=baseline,J=jacobian,
                                    x=np.zeros(2),wall_q=qs[world],frequencies=frequencies)
                cal_no_circle,_=full_response(band,world,None,1.)
                for contrast in (3.,7.):
                    parameters=np.array([contrast,.02])
                    born=baseline+(jacobian@parameters).reshape(baseline.shape)
                    true,source=full_response(band,world,(positions[target],p0[1],contrast),1.)
                    born_results[(world,target,contrast)]=(born.copy(),true.copy())
                    row=dict(world=world,target_position=target,contrast=contrast,band=band,
                             wall_loss_scale=1.,**_distance(born,true),true_response_path=source,
                             background_vs_no_circle_cal_white=float(np.linalg.norm(baseline-cal_no_circle)/.01),
                             interpretation='Born truncation error against complete Maxwell target response',
                             extra_full_solve_for_Born_prediction=False,identity=IDENTITY)
                    np.savez_compressed(out/'born'/f'{label}_contrast{int(contrast)}.npz',
                                        born=born,true=true,error=born-true,background=baseline,
                                        parameters=parameters,wall_q=qs[world],frequencies=frequencies)
                    born_rows.append(row)
                    write(out/'BORN_ENDPOINTS.json',born_rows)
                    _csv(out/'BORN_ENDPOINTS.csv',born_rows)
        for contrast in (3.,7.):
            b0,t0=born_results[(0,0,contrast)]
            b1,t1=born_results[(1,1,contrast)]
            row=dict(contrast=contrast,band='narrow',
                     born_joint_d_mean=float(np.linalg.norm(b1-b0)/.01),
                     full_joint_d_mean=float(np.linalg.norm(t1-t0)/.01),
                     joint_difference_Born_error_white=float(np.linalg.norm((b1-b0)-(t1-t0))/.01),
                     interpretation='same fixed pair and shifted targets; no q reoptimization',identity=IDENTITY)
            np.savez_compressed(out/'born'/f'JOINT_contrast{int(contrast)}.npz',
                                born0=b0,born1=b1,true0=t0,true1=t1)
            joint_born_rows.append(row)
            write(out/'BORN_JOINT_COMPARISON.json',joint_born_rows)
            _csv(out/'BORN_JOINT_COMPARISON.csv',joint_born_rows)
    except Exception as exc:
        failure=dict(error_type=type(exc).__name__,error=str(exc),
                     status='INCOMPLETE_C_BUDGET_STOP' if isinstance(exc,BudgetStop)
                     else 'INCOMPLETE_NUMERICAL_OR_INTERFACE_FAILURE',
                     no_recovery_physics_attempted=True)
        write(out/'FAILURE.json',failure)
    finally:
        summary=dict(identity=IDENTITY,status='FIXED_PAIR_PHASE_TABLE_COMPLETE' if failure is None
                     else 'INCOMPLETE_PRESERVED_PARTIAL_RESULTS',failure=failure,
                     phase_rows_completed=len(phase_rows),born_endpoint_rows_completed=len(born_rows),
                     born_joint_rows_completed=len(joint_born_rows),actual_cost=receipt(),
                     C_totals=book.totals()['C'],B2_budget_borrowed=False,
                     frequency_lists_non_nested_no_monotonicity_claim=True,
                     aperture_subset_strictly_nested=True,gamma=1.,sigma=.01,
                     fixed_pair_not_reoptimized=True,
                     min_feature_status='UNRESOLVED_ORIGINAL_DEFINITION',
                     conditional_fixed_support_geometry=True,continuous_certification=False,
                     scientific_status='fixed validated pair sensitivity; empirical numerical evidence only',
                     probabilities=None)
        write(out/'RESULTS.json',summary)
        write(out/'RHS_CPU_WALL_RECEIPT.json',receipt())
        write(out/('COMPLETE.json' if failure is None else 'INCOMPLETE.json'),summary)
    print(json.dumps(summary),flush=True)
    return summary


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('manifest')
    parser.add_argument('--include-loss',action='store_true')
    args=parser.parse_args()
    result=main(args.manifest,include_loss=args.include_loss)
    if result['failure'] is not None:
        raise SystemExit(2)
