"""Conditional Q6/Q20 proposals and shared exact-kernel, full-material controls.

run_multibranch(model,Y,m0,config) is the public, label-free prediction API.
The quotient routines are imported unchanged from the historical Z-ALG R2.
"""
from collections import Counter
from pathlib import Path
import importlib
import sys
import time
import numpy as np
from scipy.optimize import least_squares
from .finite_response import (ExactFiniteResponse, WorkLedger, spatial_supports,
                              pack_complex, full_data_jacobian)


def _frozen_solvers():
    repository = Path(__file__).resolve().parents[3]
    if str(repository) not in sys.path:
        sys.path.insert(0, str(repository))
    q6 = importlib.import_module('research.research_z_alg_r2.src.quotient_solver')
    q20 = importlib.import_module('research.research_z_alg_r2.src.general_quotient_solver')
    return q6.enumerate_quotient, q20.enumerate_general_quotient


def enumerate_material_candidates(kernel, Y, channels=1, box=(.02, 20.), tolerance=1e-7):
    """Keep every rejection, chart failure and all-data residual separately."""
    rows, charts, candidates = [], [], []
    try:
        sources, receivers = kernel.select_channels(channels)
    except ValueError as exc:
        return dict(candidates=[], roots=[], charts=[dict(status='UNSUPPORTED', reason=str(exc))])
    solvers = _frozen_solvers()
    for receiver in receivers:
        view = kernel.quotient_view(sources, int(receiver))
        target = Y[receiver:receiver+1, sources]
        chart = dict(receiver=int(receiver), sources=sources.tolist(), k=kernel.k)
        try:
            result = solvers[kernel.k == 6](view, target)
            system = result['system']
            chart.update(status='NUMERICAL_ENUMERATION', resolved_paths=result['resolved_paths'],
                         unresolved_paths=result['unresolved_paths'],
                         all_paths_numerically_resolved=result['all_paths_numerically_resolved'],
                         certificate='NOT_RIGOROUSLY_CERTIFIED', counters=result['counters'])
        except (ValueError, np.linalg.LinAlgError) as exc:
            chart.update(status='UNRESOLVED_CHART', reason=str(exc))
            charts.append(chart)
            continue
        for path in result['paths']:
            row = dict(receiver=int(receiver), sources=sources.tolist(), path=path['path'],
                       endpoint_status=path['status'], endpoint_residual=path['endpoint_residual'],
                       support=kernel.support.tolist(), status='REJECTED')
            if path['status'] != 'RESOLVED_FINITE_NUMERICAL':
                row['reason'] = 'UNRESOLVED_ENDPOINT'
                rows.append(row)
                continue
            x, s = path['x'], kernel.k//2
            conjugacy = np.linalg.norm(x[s:]-x[:s].conj())/(1+np.linalg.norm(x[:s]))
            row['conjugacy_relative'] = float(conjugacy)
            if conjugacy > 1e-6:
                row['reason'] = 'NONREAL_CONJUGATE_LIFT'
            else:
                w = system.w0+x[:s]@system.N
                v = view.L.reshape(kernel.k)+w@kernel.C
                row['minimum_node_abs'] = float(np.min(np.abs(v)))
                if np.min(np.abs(v)) <= 1e-10:
                    row['reason'] = 'UNRESOLVED_DUAL_NODE'
                else:
                    ratio = w/v
                    delta = np.real(w*v.conj())/np.abs(v)**2
                    constitutive = np.linalg.norm(w-delta*v)/max(np.linalg.norm(w), 1e-300)
                    row.update(maximum_material_imaginary=float(np.max(np.abs(ratio.imag))),
                               dual_constitutive_relative=float(constitutive))
                    if np.max(np.abs(ratio.imag)) > tolerance or constitutive > 1e-8:
                        row['reason'] = 'NONREAL_MATERIAL'
                    else:
                        m = kernel.material(delta)
                        row['material'] = m.tolist()
                        with kernel.ledger.operation('SMALL_POLE_FILTER', small_condition_SVD=1):
                            singular = np.linalg.svd(kernel.coefficients(delta)[0], compute_uv=False)
                        gap = float(singular[-1]/max(singular[0], 1e-300))
                        row['resolvent_relative_singular_gap'] = gap
                        if np.any(m < box[0]-1e-8) or np.any(m > box[1]+1e-8):
                            row['reason'] = 'OUTSIDE_REAL_BOX'
                        elif gap <= 1e-10:
                            row['reason'] = 'PHYSICAL_POLE_OR_UNRESOLVED'
                        elif any(np.linalg.norm(m-c['m'])/np.sqrt(len(m)) < 1e-7 for c in candidates):
                            row['reason'] = 'DUPLICATE_MATERIAL_ROOT'
                        else:
                            prediction = kernel.evaluate(delta)
                            subset = prediction[receiver:receiver+1, sources]
                            row.update(status='LEGAL_SUBSET_ROOT', reason=None,
                                       subset_relative_residual=float(np.linalg.norm(subset-target)/max(np.linalg.norm(target), 1e-300)),
                                       all_data_relative_residual=float(np.linalg.norm(prediction-Y)/max(np.linalg.norm(Y), 1e-300)))
                            candidates.append(dict(m=m, prediction=prediction, provenance=row.copy()))
            rows.append(row)
        charts.append(chart)
    return dict(candidates=candidates, roots=rows, charts=charts)


def _prior(model,m):
    return float(model.reg.evaluate(m)[0]) if hasattr(model,'reg') else 0.


def _prior_residual(model,m):
    return model.regularizer_residual(m) if hasattr(model,'regularizer_residual') else (np.zeros(0),np.zeros((0,len(m))))


def _select(candidates, anchor, quota, scale, model):
    """One best original score, then material-distance clusters, without GT."""
    if quota <= 0:
        return []
    for c in candidates:
        c['score'] = .5*c['relative_residual']**2+_prior(model,c['m'])
    pool = sorted(candidates, key=lambda c:c['score'])
    chosen = []
    for c in pool:
        if not chosen or all(np.linalg.norm(c['m']-p['m'])/np.sqrt(len(anchor))/scale > 1e-5 for p in chosen):
            chosen.append(c)
        if len(chosen) >= quota:
            break
    # Fixed diversity reservation replaces the worst selected seed, when useful.
    if quota >= 2 and len(pool) > 1:
        farthest = max(pool, key=lambda c:np.linalg.norm(c['m']-chosen[0]['m']))
        if all(farthest is not p for p in chosen):
            if len(chosen) < quota:
                chosen.append(farthest)
            else:
                chosen[-1] = farthest
    return chosen


class _BudgetStop(RuntimeError):
    pass


def _full_correct(model, Y, seed, anchor, config, ledger, deadline):
    from .full_maxwell_corrector import corrector
    with ledger.operation('SHARED_FULL_MATERIAL_TRF',corrector_requests=1):
        out=corrector(model,Y,np.clip(seed,config['box'][0]+1e-10,config['box'][1]-1e-10),
                      dict(method='trf',box=config['box'],max_nfev=config['corrector_nfev'],
                           wall_seconds=max(deadline-time.perf_counter(),0.)))
        # Any endpoint re-evaluation is a counted model call, not free scoring.
        state=model.state(out['m'])
        residual=float(np.linalg.norm(state['Y']-Y)/max(np.linalg.norm(Y),1e-300))
    return dict(m=out['m'],prediction=state['Y'].copy(),relative_residual=residual,
                score=.5*residual**2+_prior(model,out['m']),status=out['stop'],
                shared_corrector_cost=out['cost'],history=out['history'])


def run_multibranch(model, Y, m0, config=None):
    """Label-free bounded predictor; returns arrays, diagnostics and full ledger.

    `method`: algebraic, block_trf, block_multistart, random, full_multistart.
    No proposal ever uses an unknown support or material label.  All controls
    use the same all-data finite kernel and full-material corrector.  The
    anchor is retained; returned top1 minimizes the original physical score.
    Wall limits are common hard stop points between paid numerical operations.
    """
    cfg = dict(method='algebraic',box=(.02,20.),sizes=(4,6),max_supports=4,
               channels=1,max_candidates=3,corrector_nfev=40,block_nfev=30,
               block_starts=6,proposal_evals=180,seed=20261101,
               wall_seconds=float('inf'),support_order='cycle')
    cfg.update(config or {})
    if cfg.get('ridge',0.):
        raise ValueError('Use the common model regularizer; method-specific ridge is forbidden')
    begin = time.perf_counter()
    deadline = begin+cfg['wall_seconds']
    ledger = WorkLedger()
    model_initial_cost=Counter(getattr(model,'cost',{}))
    Y, anchor = np.asarray(Y,complex), np.asarray(m0,float).copy()
    if Y.shape != (model.S.shape[0],model.E.shape[1]) or anchor.shape != (model.n,):
        raise ValueError('All-data observation or full-material anchor shape mismatch')
    box, normY = cfg['box'],max(np.linalg.norm(Y),1e-300)
    rng = np.random.default_rng(cfg['seed'])
    with ledger.operation('ANCHOR_FULL_STATE',full_state_requests=1):
        state = model.state(anchor)
    anchor_rel = float(np.linalg.norm(state['Y']-Y)/normY)
    baseline = dict(m=anchor.copy(),prediction=state['Y'].copy(),relative_residual=anchor_rel,
                    score=.5*anchor_rel**2+_prior(model,anchor),status='PRESERVED_ANCHOR')
    supports = cfg.get('supports',spatial_supports(model.shape,cfg['sizes']))
    if cfg['support_order'] == 'gradient':
        J = full_data_jacobian(model,state,ledger)/normY
        gradient = J.T@(pack_complex(state['Y']-Y)/normY)
        scores = [-np.linalg.norm(gradient[np.asarray(a,int)])/np.sqrt(len(a)) for a in supports]
        # Leave one cyclic exploration slot even at a stationary anchor.
        ordering = list(np.argsort(scores,kind='stable'))
        if len(ordering)>cfg['max_supports'] and cfg['max_supports']>1:
            ordering[cfg['max_supports']-1] = next(i for i in range(len(supports)) if i not in ordering[:cfg['max_supports']-1])
        supports = [supports[i] for i in ordering]
    elif cfg['support_order'] == 'random':
        supports = [supports[i] for i in rng.permutation(len(supports))]
    proposal_begin=time.perf_counter()
    proposal_initial_cost=Counter(getattr(model,'cost',{}))
    candidates, root_rows, charts, block_runs = [],[],[],[]
    method = cfg['method']
    if method == 'full_multistart':
        for i in range(cfg['max_candidates']):
            m = anchor.copy() if i == 0 else rng.uniform(box[0],box[1],model.n)
            candidates.append(dict(m=m,relative_residual=float('inf'),provenance=dict(kind='FULL_MATERIAL_MULTISTART',start=i)))
    else:
        for A in supports[:cfg['max_supports']]:
            if time.perf_counter() >= deadline:
                break
            kernel = ExactFiniteResponse(model,anchor,A,state=state,ledger=ledger)
            if method == 'algebraic':
                out = enumerate_material_candidates(kernel,Y,cfg['channels'],box)
                root_rows.extend(out['roots'])
                charts.extend(dict(support=list(A),**c) for c in out['charts'])
                local = out['candidates']
            elif method in ('block_trf','block_multistart','random'):
                local = []
                lower,upper = box[0]-anchor[kernel.support],box[1]-anchor[kernel.support]
                starts = 1 if method=='block_trf' else cfg['block_starts']
                if method == 'random':
                    starts = cfg['proposal_evals']
                for start in range(starts):
                    if time.perf_counter() >= deadline:
                        break
                    delta = np.zeros(kernel.k) if start == 0 else rng.uniform(lower,upper)
                    if method != 'random':
                        cache = dict(d=None)
                        def both(d):
                            if cache['d'] is None or not np.array_equal(cache['d'],d):
                                if time.perf_counter() >= deadline:
                                    raise _BudgetStop()
                                pred,J = kernel.evaluate(d,jacobian=True)
                                residual = pack_complex(pred-Y)/normY
                                prior,prior_J=_prior_residual(model,kernel.material(d))
                                residual=np.r_[residual,prior]
                                J=np.vstack([J/normY,prior_J[:,kernel.support]])
                                cache.update(d=d.copy(),residual=residual,J=J)
                            return cache
                        try:
                            fit = least_squares(lambda d:both(d)['residual'],delta,jac=lambda d:both(d)['J'],
                                                bounds=(lower,upper),max_nfev=cfg['block_nfev'],
                                                ftol=1e-10,xtol=1e-10,gtol=1e-10)
                            delta=fit.x
                            block_runs.append(dict(support=list(A),start=start,nfev=fit.nfev,njev=fit.njev,status=int(fit.status)))
                        except _BudgetStop:
                            break
                    prediction=kernel.evaluate(delta)
                    local.append(dict(m=kernel.material(delta),prediction=prediction,
                                      provenance=dict(kind=method,support=list(A),start=start)))
            else:
                raise ValueError('Unknown method '+method)
            for c in local:
                c['relative_residual']=float(np.linalg.norm(c['prediction']-Y)/normY)
                candidates.append(c)
    proposal_wall=time.perf_counter()-proposal_begin
    proposal_cost=dict(Counter(getattr(model,'cost',{}))-proposal_initial_cost)
    if method == 'full_multistart':
        selected=candidates
    else:
        selected=_select(candidates,anchor,cfg['max_candidates'],box[1]-box[0],model) if candidates else []
    # Every method pays one common anchor correction within the same seed quota.
    if not selected and cfg['max_candidates'] > 0:
        selected=[dict(m=anchor.copy(),provenance=dict(kind='NO_PROPOSAL_ANCHOR_FALLBACK'))]
    corrected=[]
    correction_begin=time.perf_counter()
    correction_initial_cost=Counter(getattr(model,'cost',{}))
    for c in selected:
        if time.perf_counter() >= deadline:
            break
        endpoint=_full_correct(model,Y,c['m'],anchor,cfg,ledger,deadline)
        endpoint['seed_material']=c['m'].copy()
        endpoint['provenance']=c['provenance']
        corrected.append(endpoint)
    pool=[baseline]+[c for c in corrected if np.isfinite(c['score'])]
    top=min(pool,key=lambda c:c['score'])
    reason_counts=Counter(row.get('reason') or row['status'] for row in root_rows)
    selected_ids={id(c) for c in selected}
    proposals=[dict(material=c['m'].tolist(),relative_residual=c.get('relative_residual'),
                    selected=id(c) in selected_ids,provenance=c['provenance']) for c in candidates]
    return dict(method=method,m=top['m'],prediction=top['prediction'],
                relative_residual=top['relative_residual'],score=top['score'],
                initial_relative_residual=anchor_rel,corrected=corrected,
                proposals=proposals,roots=root_rows,charts=charts,block_runs=block_runs,
                filter_counts=dict(reason_counts),proposal_count=len(candidates),
                corrected_count=len(corrected),ledger=ledger.as_dict(),
                model_cost=dict(Counter(getattr(model,'cost',{}))-model_initial_cost),
                phase_wall=dict(reference=proposal_begin-begin,proposal=proposal_wall,
                                corrector=time.perf_counter()-correction_begin),
                phase_model_cost=dict(reference=dict(proposal_initial_cost-model_initial_cost),
                                      proposal=proposal_cost,
                                      corrector=dict(Counter(getattr(model,'cost',{}))-correction_initial_cost)),
                wall_seconds=time.perf_counter()-begin,config=cfg,
                status='WALL_LIMIT' if time.perf_counter() >= deadline else 'COMPLETED_FROZEN_QUOTAS',
                scope='CONDITIONAL_SUBSET_ROOT_PROPOSALS_NOT_FULL_MATERIAL_ALLROOTS',
                certificate='NOT_RIGOROUSLY_CERTIFIED')
