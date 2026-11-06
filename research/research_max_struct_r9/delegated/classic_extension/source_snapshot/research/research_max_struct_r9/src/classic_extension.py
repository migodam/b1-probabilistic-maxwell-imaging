"""R9 original-acquisition, full-pixel classic controls; no truth access.

The frozen R8 disk-pulse FFT-TSOM/current/material implementations are reused.
Only dense fallback assembly is cached by this adapter, without changing the
discrete operator, solver tolerances, optimization objective, or material box.
"""
import json
import time
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.linalg import lu_factor, lu_solve, svd
from scipy.optimize import minimize
from scipy.special import hankel1, jv

from research.research_max_struct_r8.src.paper_physics import (
    MatrixFreeModel, EqualAreaGreen, acquisition,
)
from research.research_max_struct_r8.src.paper_cie_fft import run_cie_fft
from research.research_max_struct_r6.src.fft_tsom_ops import NativeFourierAPCS
from research.research_max_struct_r7.src.aplus_readout import mobius_material_readout


def dump(path, obj):
    def convert(x):
        if isinstance(x, np.ndarray): return x.tolist()
        if isinstance(x, np.generic): return x.item()
        raise TypeError(type(x).__name__)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(json.dumps(obj, indent=2, default=convert, allow_nan=False))
    temporary.replace(path)


class ClassicModel(MatrixFreeModel):
    """Same R8 operator, one retained immutable dense Green after first fallback.

    GMRES remains <=1000 iterations/RHS and uses the inherited warm starts and
    residual checks. Dense LU is retained for all source and adjoint RHS at m.
    """
    def __init__(self, public, prior=True):
        super().__init__(public, prior=prior)
        self._immutable_dense_green = None

    def dense_fallback(self, m, B, adjoint=False):
        if self._dense_m is None or not np.array_equal(m, self._dense_m):
            started, cpu = time.perf_counter(), time.process_time()
            if self._immutable_dense_green is None:
                p = self.G.public
                k, radius, coords = p['k'], p['radius'], p['centers']
                G = np.empty((self.n, self.n), complex, order='F')
                for i in range(0, self.n, 128):
                    stop = min(i + 128, self.n)
                    d = np.linalg.norm(coords[i:stop, None] - coords[None, :], axis=2)
                    d[np.arange(stop-i), np.arange(i, stop)] = 1.
                    G[i:stop] = .5j*np.pi*k*radius*jv(1, k*radius)*hankel1(0, k*d)
                np.fill_diagonal(G, self.green.diagonal)
                self._immutable_dense_green = G
                self.cost['immutable_dense_Green_assemblies'] += 1
                self.cost['immutable_dense_Green_bytes'] = G.nbytes
            A = np.empty_like(self._immutable_dense_green, order='F')
            np.multiply(-m[:, None], self._immutable_dense_green, out=A)
            A[np.diag_indices(self.n)] += 1.
            self._dense_lu = lu_factor(A, overwrite_a=True, check_finite=False)
            self._dense_m = m.copy()
            self.cost['LU'] += 1
            self.cost['dense_fallback_operator_assemblies'] += 1
            self.cost['dense_fallback_setup_wall_seconds'] += time.perf_counter()-started
            self.cost['dense_fallback_setup_CPU_seconds'] += time.process_time()-cpu
        self.cost['dense_fallback_calls'] += 1
        self.cost['actual_complex_RHS'] += B.shape[1]
        self.cost['dense_fallback_complex_RHS'] += B.shape[1]
        Z = lu_solve(self._dense_lu, B, trans=2 if adjoint else 0, check_finite=False)
        residual = (Z-self.green.adjoint(m[:, None]*Z)-B if adjoint
                    else Z-m[:, None]*self.green.apply(Z)-B)
        rel = np.linalg.norm(residual, axis=0)/np.maximum(np.linalg.norm(B, axis=0), 1e-300)
        if not np.all(np.isfinite(Z)) or np.max(rel) > 2e-9:
            raise RuntimeError('Unsolved physical dense-fallback state: '+str(float(np.max(rel))))
        self.cost['maximum_dense_fallback_relative_residual'] = max(
            self.cost.get('maximum_dense_fallback_relative_residual', 0), float(np.max(rel)))
        self._solve_initial[adjoint] = Z.copy()
        return Z


def paper_reference(model, Y, m0, stages, lsie, iteration_callback, stage_callback):
    """Independent NIE2016 Eq3/Eq14 real-box adaptation, with paid telemetry.

    Same equations as R8 paper_nie.py: sourcewise PR+, projected native Fourier
    coefficients, closed-form real modified-contrast LS, no added prior.
    This is an independent implementation, not author software.
    """
    started, cpu = time.perf_counter(), time.process_time()
    counts = Counter()
    E, S = model.E, model.S
    ns = E.shape[1]
    geom = model.paper_geometry
    green = EqualAreaGreen(geom['coords'], model.shape, geom['radius'], geom['k'], counts)
    U, singular, Vh = svd(S, full_matrices=False, check_finite=False)
    counts['receiver_SVD'] += 1
    V = Vh[:15].conj().T
    JD = V@((U[:, :15].conj().T@Y)/singular[:15, None])
    GJD, SJD = green.apply(JD), S@JD
    wd = 1/np.maximum(np.sum(abs(Y)**2, axis=0), 1e-300)
    ws = 1/np.sum(abs(E)**2, axis=0)
    old, m, history, states = {}, np.array(m0).copy(), [], []
    for stage_index, stage in enumerate(stages):
        beta = None if lsie else float(stage['beta_material'])
        basis = NativeFourierAPCS(model.shape, stage['mf'], V, counts)
        A = np.array([old.get(tuple(f), np.zeros(ns, complex)) for f in basis.frequencies])
        previous_gradient = previous_direction = None
        stop = 'MAXITER'
        for it in range(stage['maxiter']):
            factor = np.ones(len(m)) if beta is None else beta/(1+beta*m)
            def action(a):
                f = basis.apply(a)
                return factor[:, None]*(f-m[:, None]*green.apply(f))
            def adjoint(r):
                f = factor[:, None]*r
                return basis.adjoint(f-green.adjoint(m[:, None]*f))
            f = basis.apply(A)
            J, Q = JD+f, E+GJD+green.apply(f)
            D, R = S@J-Y, factor[:, None]*(J-m[:, None]*Q)
            g = basis.adjoint(S.conj().T@(D*wd))+adjoint(R*ws)
            pr = (0 if previous_gradient is None else np.maximum(
                0, np.sum(np.real(g.conj()*(g-previous_gradient)), axis=0)/
                np.maximum(np.sum(abs(previous_gradient)**2, axis=0), 1e-300)))
            d = -g if previous_direction is None else -g+previous_direction*pr
            bad = np.sum(np.real(g.conj()*d), axis=0) >= 0
            d[:, bad] = -g[:, bad]
            dJ = basis.apply(d)
            dD, dR = S@dJ, action(d)
            den = wd*np.sum(abs(dD)**2, axis=0)+ws*np.sum(abs(dR)**2, axis=0)
            alpha = -np.sum(np.real(g.conj()*d), axis=0)/np.maximum(den, 1e-300)
            change = d*alpha
            delta = float(np.sqrt(np.mean(np.sum(abs(change)**2, axis=0)/
                         np.maximum(np.sum(abs(A)**2, axis=0), 1e-300))))
            A += change
            f = basis.apply(A)
            J, Q = JD+f, E+GJD+green.apply(f)
            T, B = (J, Q) if beta is None else (beta*J, Q+beta*J)
            r = np.sum(ws*np.real(B.conj()*T), axis=1)/np.maximum(
                np.sum(ws*abs(B)**2, axis=1), 1e-300)
            if beta is None:
                m = np.clip(r, 0, 20)
            else:
                r = np.clip(r, 0, beta*20/(1+beta*20))
                m = r/(beta*(1-r))
            previous_gradient, previous_direction = g, d
            counts['PR_steps'] += 1
            history.append(dict(stage=stage_index, iteration=it+1, delta2D=delta,
                                cpu_seconds=time.process_time()-cpu,
                                wall_seconds=time.perf_counter()-started))
            iteration_callback(dict(stage=stage_index, iteration=it+1, m=m,
                                    operation_counts=dict(counts)))
            if it >= 1 and delta < stage['paper_delta']:
                stop = 'COEFFICIENT_TOL'
                break
        state = dict(stage=stage_index, m=m.copy(), J=J.copy(), Q=Q.copy(),
                     coefficients=A.copy(), status=stop, beta_material=beta,
                     beta_current=beta, iterations=it+1)
        states.append(state)
        stage_callback(state)
        history.append(dict(stage=stage_index, event='stage_end', status=stop,
                            iterations=it+1, cpu_seconds=time.process_time()-cpu))
        old = {tuple(f): a.copy() for f, a in zip(basis.frequencies, A)}
    return dict(stage_states=states, history=history,
                execution_receipt=dict(status='COMPLETED', cpu_seconds=time.process_time()-cpu,
                    wall_seconds=time.perf_counter()-started, operation_counts=dict(counts),
                    implementation='independent real-box NIE/LSIE; no added prior; sourcewise PR'))


def solve(manifest, task, dest, process_started, process_cpu_started):
    """Save all paid stages and a continuous 300-step physical correction."""
    dest = Path(dest)
    dest.mkdir(parents=True, exist_ok=True)
    dump(dest/'FROZEN_TASK.json', dict(task=task, contract=manifest['contract'],
                                     GT_access='NONE', status='FROZEN_BEFORE_SOLVE'))
    case = next(c for c in manifest['cases'] if c['case_id'] == task['case_id'])
    started, cpu = time.perf_counter(), time.process_time()
    model = None
    auxiliary_receipt = {}
    history, curve, stage_rows, evaluations = [], [], [], []
    try:
        observation_file=Path(case['public_observation_file'])
        if not observation_file.is_absolute():
            observation_file=Path(__file__).resolve().parents[3]/observation_file
        data = np.load(observation_file)
        Y = np.array(data['Y'], complex)
        if Y.shape != (40, 20) or int(data['forward_side']) != 150:
            raise ValueError('Unexpected original public observation schema')
        if float(data['frequency']) != 400e6:
            raise ValueError('Unexpected public frequency')
        setup_t, setup_c = time.perf_counter(), time.process_time()
        public = acquisition(64, 400e6, dense=False)
        prior = not task['arm'].startswith('paper_')
        model = ClassicModel(public, prior=prior)
        setup = dict(cpu_seconds=time.process_time()-setup_c,
                     wall_seconds=time.perf_counter()-setup_t)
        m0 = np.zeros(model.n)
        def tick(phase, **kw):
            return dict(phase=phase, cpu_seconds=time.process_time()-cpu,
                        process_cpu_seconds=time.process_time()-process_cpu_started,
                        wall_seconds=time.perf_counter()-started,
                        physical_cost=dict(model.cost), **kw)
        def progress(phase, **kw): dump(dest/'PROGRESS.json', tick(phase, **kw))
        progress('INITIALIZED', setup=setup)
        def iteration_callback(st):
            # No physical forward per joint iteration: that would change the
            # cost substantially. Stage Phi and every corrector Phi are paid.
            if st['iteration'] == 1 or st['iteration'] % 25 == 0:
                progress('AUXILIARY', stage=st['stage'], iteration=st['iteration'],
                         auxiliary_operation_counts=st['operation_counts'])
        def stage_callback(st):
            material = st['m']
            phi = model.physical_score(material, Y)
            state = model.state(material)
            np.savez_compressed(dest/f"stage_{st['stage']:02d}.npz",
                                m=material, J=st['J'], Q=st['Q'],
                                coefficients=st['coefficients'], Y_prediction=state['Y'])
            row = tick('AUXILIARY_STAGE', stage=st['stage'], Phi=phi,
                       relative_data_residual=float(np.linalg.norm(state['Y']-Y)/np.linalg.norm(Y)),
                       beta_current=st.get('beta_current'), beta_material=st.get('beta_material'))
            stage_rows.append(row)
            curve.append(row)
            dump(dest/'STAGES.json', stage_rows)
            progress('AUXILIARY_STAGE_COMPLETE', stage=st['stage'], Phi=phi)
        stages = [dict(mf=mf, beta_material=b,
                       beta_current=b if task['arm']=='cie_single' else case['beta'][0],
                       maxiter=500, paper_delta=.01 if j==0 else .001,
                       coefficient_tol=(.01 if j==0 else .001)**2)
                  for j, (mf, b) in enumerate(zip(case['mf'], case['beta']))]
        initial_phi = model.physical_score(m0, Y)
        curve.append(tick('INITIAL', Phi=initial_phi))
        if task['arm'] == 'physical_lbfgs':
            result = dict(stage_states=[], history=[], execution_receipt={})
            candidates = [m0]
        elif task['arm'].startswith('paper_'):
            result = paper_reference(model, Y, m0, stages, task['arm']=='paper_lsie',
                                     iteration_callback, stage_callback)
            candidates = [m0]+[s['m'] for s in result['stage_states']]
        else:
            cfg = dict(arm=task['arm'], stages=stages, dpcs_rank=15,
                       material_box=(0, 20), miniter=2,
                       sourcewise_cg=task['arm'] in ('cie_single','cie_dual'),
                       regularizer=lambda z: model.reg.evaluate(z)[:2],
                       fft_geometry=model.paper_geometry, physical_sigma=float(np.linalg.norm(Y)),
                       material_inner_maxiter=60, cache_material_value_gradient=True,
                       aplus_regularizer=model.reg, aplus_sweeps=2, aplus_arc_candidates=True,
                       iteration_callback=iteration_callback, stage_callback=stage_callback)
            if task['arm'] == 'a_plus': cfg['material_readout'] = mobius_material_readout
            result = run_cie_fft(model, Y, m0, cfg)
            candidates = [m0]+[s['m'] for s in result['stage_states']]
        history = result['history']
        auxiliary_receipt = result['execution_receipt']
        dump(dest/'AUXILIARY_HISTORY.json', history)
        dump(dest/'AUXILIARY_RECEIPT.json', auxiliary_receipt)
        values = [model.physical_score(m, Y) for m in candidates]
        selected = int(np.argmin(values))
        auxm = candidates[selected].copy()
        np.savez_compressed(dest/'AUXILIARY_PREDICTION.npz', m=auxm,
                            Y_prediction=model.state(auxm)['Y'])
        dump(dest/'AUXILIARY_SEAL.json', dict(status='SEALED_NO_GT', selected_candidate=selected,
                                             candidate_Phi=values, GT_loaded=False))
        best = [auxm.copy(), float(values[selected])]
        accepted = 0
        curve.append(tick('CORRECTOR_START', Phi=best[1]))
        def fg(z):
            # A failed physical state is a failed run, never a finite dummy
            # objective/zero gradient labelled a converged physical optimum.
            v, g = model.data_value_gradient(z, Y)
            if v < best[1]: best[:] = [z.copy(), float(v)]
            evaluations.append(tick('CORRECTOR_EVALUATION', evaluation=len(evaluations)+1,
                                    Phi=float(v), projected_gradient_inf=float(np.max(abs(
                                        np.where((z<=0)&(g>0)|(z>=20)&(g<0), 0, g))))))
            if len(evaluations) % 10 == 0:
                progress('CORRECTOR', evaluation=len(evaluations), accepted=accepted, Phi=float(v))
                dump(dest/'CORRECTOR_EVALUATIONS.json', evaluations)
            return 1e6*v, 1e6*g
        def callback(z):
            nonlocal accepted
            accepted += 1
            phi = model.physical_score(z, Y)
            if phi < best[1]: best[:] = [z.copy(), float(phi)]
            curve.append(tick('CORRECTOR_ACCEPTED', accepted_iteration=accepted, Phi=phi,
                              incumbent_Phi=best[1]))
            if accepted == 80:
                np.savez_compressed(dest/'CHECKPOINT_080.npz', m=best[0],
                                    Y_prediction=model.state(best[0])['Y'])
                dump(dest/'CHECKPOINT_080.json', tick('CORRECTOR_080', Phi=best[1],
                    accepted_iteration=80, history='continuous L-BFGS-B; no restart', GT_loaded=False))
            if accepted % 10 == 0: dump(dest/'PHI_CPU_CURVE.json', curve)
        opt = minimize(fg, auxm, jac=True, method='L-BFGS-B', bounds=[(0,20)]*model.n,
                       callback=callback, options=dict(maxiter=300, maxfun=1200,
                           ftol=1e-12, gtol=.01, maxls=30))
        final = best[0]
        state = model.state(final)
        np.savez_compressed(dest/'PREDICTION.npz', m=final, m_auxiliary=auxm,
                            Y_prediction=state['Y'], m_optimizer_final=opt.x)
        curve.append(tick('FINAL', Phi=best[1]))
        dump(dest/'PHI_CPU_CURVE.json', curve)
        dump(dest/'CORRECTOR_EVALUATIONS.json', evaluations)
        statuses = [h for h in history if h.get('event')=='stage_end']
        receipt = dict(status='COMPLETED', task_id=task['task_id'], case_id=case['case_id'],
            arm=task['arm'], dataset_status='PUBLIC_PAPER_DEVELOPMENT', GT_loaded=False,
            setup=setup, cpu_seconds=time.process_time()-cpu,
            wall_seconds=time.perf_counter()-started,
            process_cpu_seconds=time.process_time()-process_cpu_started,
            process_wall_seconds=time.perf_counter()-process_started,
            physical_cost=dict(model.cost), auxiliary=auxiliary_receipt,
            stages=statuses, selected_auxiliary_candidate=selected, candidate_Phi=values,
            original_Phi=best[1],
            relative_data_residual=float(np.linalg.norm(state['Y']-Y)/np.linalg.norm(Y)),
            corrector=dict(iterations=int(opt.nit), evaluations=int(opt.nfev),
                           accepted_callbacks=accepted, stop=str(opt.message),
                           optimizer_success=bool(opt.success), cap=300,
                           status='CONVERGED_OPTIMIZER' if opt.success else 'CAP' if opt.nit>=300 else 'STOPPED',
                           checkpoint80='SAVED' if accepted>=80 else 'NOT_REACHED_EARLY_STOP'),
            convergence_qualified=bool(opt.success), all_material_pixels_free=True,
            material_prior='smooth=1e-4,TV=1e-4' if prior else 'NONE',
            source_observation_file=case['public_observation_file'],
            adaptation='real chi[0,20]; disk-cell quadrature; finite auxiliary cap; common corrector',
            timer_scope='process_* includes Python imports since entry; queue wait4 CPU captures full child',
            software_sources=dict(adapter=__file__,
                physics='research/research_max_struct_r8/src/paper_physics.py',
                cie='research/research_max_struct_r8/src/paper_cie_fft.py',
                Aplus='research/research_max_struct_r7/src/aplus_readout.py'))
        dump(dest/'RECEIPT.json', receipt)
        dump(dest/'PREDICTION_SEAL.json', dict(status='SEALED_BEFORE_OFFLINE_SCORING',
            task_id=task['task_id'], prediction_file='PREDICTION.npz', GT_loaded=False,
            UTC=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())))
        return receipt
    except BaseException as exc:
        dump(dest/'FAILURE.json', dict(status='FAILED', task_id=task['task_id'],
            error=repr(exc), cpu_seconds=time.process_time()-cpu,
            wall_seconds=time.perf_counter()-started,
            process_cpu_seconds=time.process_time()-process_cpu_started,
            process_wall_seconds=time.perf_counter()-process_started,
            physical_cost=dict(model.cost) if model is not None else {},
            auxiliary=auxiliary_receipt or getattr(exc, 'cie_execution_receipt', {}),
            GT_loaded=False))
        raise
