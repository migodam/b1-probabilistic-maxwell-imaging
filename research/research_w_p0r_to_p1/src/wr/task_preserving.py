"""Known-reference physical-existence corrector preserving finite task response.

Not online inference or a posterior. Import constructs no solver or budget book.
The fixed projection direction is selected once from a bounded linear proposal.
"""
from __future__ import annotations

import argparse
import json
import time
import uuid

import numpy as np
from p0.geometry import probes

from .accounting import ROOT, ScopedBook, BudgetStop, write
from .native_cal import bounded_task
from .native_edit import _cost, _ledger
from .native_world import LinearWorld
from .physical_witness import mesh_for
from .physics import Forward, pack, packj
from .optimization import trust_region


IDENTITY = 'ORACLE_KNOWN_REFERENCE_TASK_PRESERVING_PHYSICAL_EXISTENCE_NOT_ONLINE'


def _metrics(cal, task, cal0, task0, u):
    cr, tr = pack(cal-cal0), pack(task-task0)
    projection = float(u@tr)
    return dict(d_cal=float(np.linalg.norm(cr)), d_T=float(np.linalg.norm(tr)),
                fixed_task_projection=projection, task_projection_target=3.2,
                hinge_residual=20.*min(0., projection-3.2))


def main(seed, band='narrow', use_cached_baseline=True):
    """Use only the existing cal_<seed>_<band> scope; no allowance reset.

    Proposal: metric RMS <=.12, box [-1,1], linear cal norm <=.8. Corrector:
    fixed task-direction hinge plus cal residual, box [-1,1], 8 iter / 20 calls.
    The .12 proposal radius is reported separately from final achieved RMS.
    """
    seed = int(seed)
    if band not in ('narrow', 'wide'):
        raise ValueError('band must be narrow or wide')
    cfg = json.loads((ROOT/'configs/ORIGINAL_P0.json').read_text())
    if cfg['noise']['each_real_quadrature_std'] != .01:
        raise RuntimeError('original noise differs from shared unchanged whitening')
    frequencies = cfg['frequencies'][band]
    book = ScopedBook('B2', f'cal_{seed}_{band}')
    folder = f'B2/task_preserving_{seed}_{band}'
    out = ROOT/'results'/folder
    out.mkdir(parents=True, exist_ok=False)
    run_id = uuid.uuid4().hex
    model = LinearWorld(seed)
    task_model = LinearWorld(seed, probes()[8])
    zero = np.zeros(model.n)
    grid = np.stack(np.meshgrid(np.linspace(-2.,2.,201),
                                np.linspace(0.,3.,151)), axis=-1).reshape(-1,2)
    material0 = model.coefficient(zero,1.)(grid)
    task_material0 = task_model.coefficient(zero,1.)(grid)
    ledger_start = len(_ledger(book))
    cpu_start, wall_start = time.process_time(), time.perf_counter()
    provenance = None
    checkpoint = latest = None
    current_x = zero.copy()
    failure = None
    fit = None
    proposal = u = None
    history = []
    cal0 = task0 = None
    write(out/'PLAN.json', dict(seed=seed, band=band, identity=IDENTITY,
          scope=book.scope, scope_rhs_cap=book.scope_caps['rhs'], run_id=run_id,
          proposal=dict(metric_RMS_max=.12, box=[-1.,1.], linear_cal_norm_max=.8,
                        solver='existing bounded_task; local sequential linear search, not global certificate'),
          corrector=dict(max_iter=8,max_calls=20,prior_lambda=0.,box=[-1.,1.],
                         task_projection_target=3.2,hinge_weight=20.,
                         objective='half squared cal whitened residual plus half squared negative task-projection hinge',
                         metric_RMS_constraint='proposal only; actual final RMS reported'),
          task_probe=probes()[8], gain=1., original_support_noise_material_prior_unchanged=True,
          baseline_parameter_vector='same zero vector for cal and task',
          probability_claim=None, independent_validation_required=True))

    def cost():
        return _cost(book, ledger_start, cpu_start, wall_start)

    def save_point(name, point, gradient=None):
        x = point['x']
        metrics = _metrics(point['cal'], point['task'], cal0, task0, u)
        stationarity = (float(np.linalg.norm(x-np.clip(x-gradient,-1.,1.), ord=np.inf))
                        if gradient is not None else None)
        np.savez_compressed(out/f'{name}.npz', x=x, cal=point['cal'], task=point['task'],
                            cal_base=cal0, task_base=task0, proposal=proposal, u=u,
                            material_grid=grid, material=model.coefficient(x,1.)(grid),
                            task_material=task_model.coefficient(x,1.)(grid),
                            reference_material=material0, reference_task_material=task_material0,
                            frequencies=frequencies,
                            **({'gradient':gradient} if gradient is not None else {}))
        row = dict(**metrics, projected_gradient=stationarity,
                   material_RMS=float(np.sqrt(max(0.,x@metric@x))),
                   actual_cost=cost(), identity=IDENTITY,
                   provisional_actual_thresholds_met=bool(metrics['d_cal']<=1 and metrics['d_T']>=3),
                   scientific_status='AWAITING_INDEPENDENT_VALIDATION', probabilities=None)
        write(out/f'{name}.json', row)
        return row

    try:
        # One common public reference geometry; cal has no probe material, task
        # has fixed probe8. LinearWorld(0) exactly equals the original world.
        mesh = mesh_for(seed, probes()[8])
        cal_fw = Forward(model, mesh, frequencies, book, folder+'/cal',
                         phase=f'B2_task_preserving_cal_{run_id}', identity=IDENTITY,
                         h=.145,pml_strength=6,retain_current=True)
        task_fw = Forward(task_model, mesh, frequencies, book, folder+'/task',
                          phase=f'B2_task_preserving_task_{run_id}', identity=IDENTITY,
                          h=.145,pml_strength=6,retain_current=False)
        source = ROOT/f'results/B2/native_cal_{seed}_{band}/BASE.npz'
        if use_cached_baseline and source.exists():
            with np.load(source) as cached:
                j0, k0 = cached['J'].copy(), cached['K'].copy()
                cal0, task0 = cached['cal'].copy(), cached['task'].copy()
                metric = cached['material_metric'].copy()
                cached_grid, cached_material = cached['grid'].copy(), cached['material0'].copy()
            if (not np.array_equal(cached_grid, grid)
                    or not np.allclose(cached_material, material0, rtol=0.,atol=1e-14)):
                raise RuntimeError('cached baseline native grid/reference material mismatch')
            source_rows = [r for r in _ledger(book) if r.get('scope')==book.scope
                           and r.get('phase') in ('B2_native_cal','B2_native_task')
                           and r.get('event')=='reserve']
            provenance = dict(source=str(source), cached_arrays_only=True,
                              original_setup_not_free=True, live_LU_or_fields_loaded=False,
                              new_FJ_always_passes_through_Forward_and_Book=True,
                              source_scope=book.scope,
                              prior_source_phase_reserved_rhs=sum(r.get('rhs',0) for r in source_rows),
                              source_phase_total_includes_later_nonsetup_native_cal_work=True,
                              source_reserve_ids=[r['id'] for r in source_rows],
                              scope_usage_already_counts_original_setup=True)
        else:
            cal0,j0 = cal_fw.evaluate(zero,jac=True)
            task0,k0 = task_fw.evaluate(zero,jac=True)
            dm = model.derivative(zero,grid,1.)
            metric = np.real(dm.conj().T@dm)/len(grid)
            provenance = dict(source='new paid full-FJ baseline at zero for both states',
                              cached_arrays_only=False, original_setup_not_free=True,
                              actual_setup_cost=cost(), live_LU_or_fields_loaded=False)
        shape=(len(frequencies),len(cfg['tx']),len(cfg['rx']))
        if (cal0.shape!=shape or task0.shape!=shape or j0.shape!=(np.prod(shape),model.n)
                or k0.shape!=j0.shape or metric.shape!=(model.n,model.n)):
            raise RuntimeError('baseline arrays have incompatible shapes')
        if not all(np.all(np.isfinite(a)) for a in (cal0,task0,j0,k0,metric)):
            raise RuntimeError('nonfinite baseline arrays')
        if np.max(abs(metric-metric.T))>1e-12:
            raise RuntimeError('material metric is not symmetric')
        write(out/'SETUP_PROVENANCE.json',provenance)
        np.savez_compressed(out/'BASE.npz',J=j0,K=k0,cal=cal0,task=task0,
                            material_metric=metric,grid=grid,material0=material0)
        proposal,linear = bounded_task(packj(j0),packj(k0),metric,.12)
        write(out/'LINEAR_SEARCH.json',linear)
        predicted_cal=float(np.linalg.norm(packj(j0)@proposal))
        proposal_rms=float(np.sqrt(max(0.,proposal@metric@proposal)))
        if np.max(abs(proposal))>1.+1e-12 or predicted_cal>.8+1e-8 or proposal_rms>.12+1e-9:
            raise RuntimeError('bounded_task proposal failed frozen feasibility checks')
        direction=pack(k0@proposal)
        norm=float(np.linalg.norm(direction))
        if not np.isfinite(norm) or norm<=0:
            raise RuntimeError('bounded proposal has no nonzero finite predicted task direction')
        u=direction/norm
        write(out/'PROPOSAL_PROTOCOL.json',dict(predicted_cal=predicted_cal,
              predicted_task_norm=norm,proposal_material_RMS=proposal_rms,
              fixed_u_norm=float(np.linalg.norm(u)),direction_frozen_before_corrector=True,
              no_gain_fitting=True,no_regularizer=True,no_posterior=True))

        def evaluate(x,jac=False):
            nonlocal current_x,latest
            current_x=np.asarray(x).copy()
            np.savez_compressed(out/'CURRENT_LATENT.npz',current_trial_x=current_x,
                                accepted_x=checkpoint['x'] if checkpoint is not None else zero,
                                accepted_available=checkpoint is not None)
            cal,jc=cal_fw.evaluate(x,jac=jac)
            # Retain this completed half even if the task state's budget stops.
            np.savez_compressed(out/'CURRENT_CAL_PARTIAL.npz',x=x,cal=cal)
            task,jt=task_fw.evaluate(x,jac=jac)
            cr,tr=pack(cal-cal0),pack(task-task0)
            projection=float(u@tr)
            residual=np.r_[cr,20.*min(0.,projection-3.2)]
            if jac:
                if jc is None or jt is None:
                    raise RuntimeError('full Jacobian missing from one of the two paid states')
                task_row=20.*(u@packj(jt)) if projection<3.2 else np.zeros(model.n)
                jacobian=np.vstack((packj(jc),task_row))
                gradient=jacobian.T@residual
            else:
                jacobian=gradient=None
            latest=dict(x=np.asarray(x).copy(),cal=cal.copy(),task=task.copy(),gradient=gradient)
            write(out/'CURRENT_TRIAL.json',dict(**_metrics(cal,task,cal0,task0,u),
                  jac_requested=jac,actual_cost=cost()))
            return residual,jacobian

        def callback(event):
            nonlocal checkpoint
            history.append(event)
            write(out/'ITERATIONS.json',history)
            if event.get('event')=='anchor':
                anchor=np.asarray(event['x'])
                if latest is None or latest['gradient'] is None or not np.array_equal(anchor,latest['x']):
                    raise RuntimeError('accepted step lacks exact paired full-FJ checkpoint')
                checkpoint={k:v.copy() for k,v in latest.items()}
                index=sum(r.get('event')=='anchor' for r in history)-1
                save_point(f'accepted_{index:03d}',checkpoint,checkpoint['gradient'])
                save_point('RETAINED_ACCEPTED',checkpoint,checkpoint['gradient'])

        # Explicit proposal physics is separately recorded and charged setup;
        # no retained fields/J are fabricated from baseline or prediction arrays.
        current_x=proposal.copy()
        evaluate(proposal,jac=False)
        save_point('PROPOSAL',latest)
        write(out/'ITERATIONS.json',history)
        fit=trust_region(evaluate,proposal,model.lower,model.upper,max_iter=8,max_calls=20,
                         radius=.25,max_radius=.75,l2_radius=False,prior_lambda=0.,
                         callback=callback,fatal_exceptions=(Exception,))
        if checkpoint is None or not np.array_equal(fit.x,checkpoint['x']):
            raise RuntimeError('optimizer returned without its paired accepted checkpoint')
        final=save_point('FINAL',checkpoint,checkpoint['gradient'])
    except Exception as exc:
        failure=dict(error_type=type(exc).__name__,error=str(exc),
                     status='INCOMPLETE_BUDGET_STOP' if isinstance(exc,BudgetStop)
                     else 'INCOMPLETE_NUMERICAL_OR_INTERFACE_FAILURE')
        if checkpoint is not None:
            final=save_point('RETAINED_ACCEPTED',checkpoint,checkpoint['gradient'])
        else:
            final=None
        write(out/'FAILURE.json',failure)
        np.savez_compressed(out/'STOPPED_LATENT.npz',current_trial_x=current_x,
                            accepted_x=checkpoint['x'] if checkpoint is not None else zero,
                            accepted_available=checkpoint is not None)
    finally:
        summary=dict(seed=seed,band=band,identity=IDENTITY,
                     status='NUMERICAL_CORRECTOR_COMPLETED' if failure is None
                     else 'INCOMPLETE_NO_CONVERGENCE_CLAIM',
                     failure=failure,setup_provenance=provenance,
                     numerically_converged=bool(fit.success) if fit is not None else False,
                     stopping=str(fit.message) if fit is not None else None,
                     iterations=int(fit.nit) if fit is not None else None,
                     nfev=int(fit.nfev) if fit is not None else None,
                     njev=int(fit.njev) if fit is not None else None,
                     retained_result=final if 'final' in locals() else None,
                     actual_total_cost=cost(),B2_totals=book.totals()['B2'],
                     scientific_status='PROVISIONAL_ACTUAL_THRESHOLDS_ONLY_AWAIT_INDEPENDENT_VALIDATION',
                     probabilities=None,no_new_allowance=True)
        write(out/'RESULTS.json',summary)
        write(out/'RHS_CPU_WALL_RECEIPT.json',cost())
        write(out/('COMPLETE.json' if failure is None else 'INCOMPLETE.json'),summary)
    print(json.dumps(summary),flush=True)
    return summary


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('seed',type=int)
    parser.add_argument('band',nargs='?',default='narrow',choices=['narrow','wide'])
    parser.add_argument('--fresh-baseline',action='store_true')
    args=parser.parse_args()
    result=main(args.seed,args.band,use_cached_baseline=not args.fresh_baseline)
    if result['failure'] is not None:
        raise SystemExit(2)
