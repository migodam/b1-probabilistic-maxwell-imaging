"""Truth-free finite-candidate evaluation using frozen R4 VIE and R2 TRF.

Only public operators, observations, GN anchors, and sealed NN outputs enter
this module. The anchor's original construction is reported separately; every
incremental operation from proposal preparation through selection is charged.
"""
from contextlib import contextmanager
import json
from pathlib import Path
import time
import numpy as np
from scipy.optimize import minimize
from threadpoolctl import threadpool_limits
from research.research_gn_search_space_r2.src.search import (
    Search, Stop, PhaseStop, Model, Regularizer, safe,
)

LOWER, UPPER = .02, 20.
BUDGET = dict(wall_seconds=10., max_LU=100, max_complex_RHS=6000)


def normalize_public(public):
    out = dict(public)
    if 'E0' not in out:
        out['E0'] = out.get('E', out.get('U0'))
    if 'S' not in out:
        out['S'] = out.get('H')
    for key, shape in [('G', (576, 576)), ('E0', (576, 12)), ('S', (48, 576)), ('Y', (48, 12))]:
        if np.shape(out.get(key)) != shape:
            raise ValueError(f'{key} must have shape {shape}; got {np.shape(out.get(key))}')
    return out


def _preparation(cost):
    if cost is None:
        return dict(wall_seconds=0., LU=0, actual_complex_RHS=0, complete=False,
                    missing_reason='No incremental preparation receipt supplied')
    counts = cost.get('costcounts', cost)
    result = dict(cost)
    result.update(wall_seconds=float(cost.get('wall_seconds', 0.)),
                  LU=int(counts.get('LU', 0)),
                  actual_complex_RHS=int(counts.get('actual_complex_RHS', 0)),
                  complete=bool(cost.get('complete', True)))
    if any(result[k] < 0 for k in ('wall_seconds', 'LU', 'actual_complex_RHS')):
        raise ValueError('Preparation costs cannot be negative')
    return result


class CountedModel(Model):
    """Add prospective budget guards without modifying the frozen Model."""
    def __init__(self, public, owner):
        self.owner = owner
        super().__init__(public)

    @contextmanager
    def operation(self, kind, rhs=0, **cost):
        self.owner.guard()
        if self.counts['LU'] + int(cost.get('LU', 0)) > BUDGET['max_LU']:
            raise Stop('MAX_LU_PROSPECTIVE')
        if self.counts['actual_complex_RHS'] + int(rhs) > BUDGET['max_complex_RHS']:
            raise Stop('MAX_RHS_PROSPECTIVE')
        with super().operation(kind, rhs=rhs, **cost):
            yield


class Evaluation(Search):
    """Reuse R2's mature bounded analytic-TRF correction, with 10-second gates."""
    def __init__(self, public, anchor, *, noise, preparation_cost, regularization=0.):
        self.started, self.cpu = time.perf_counter(), time.process_time()
        self.preparation = _preparation(preparation_cost)
        self.noise = str(noise)
        self.data_gate = .001 if self.noise == 'clean' else .003
        self.task = {'noise': self.noise, 'stopping_profile': 'tight'}
        self.stopping_profile = 'tight'
        self.phase_end = np.inf
        self._gate_enabled = False
        self._gate_pending = False
        self.model = CountedModel(normalize_public(public), self)
        self.model.counts['LU'] += self.preparation['LU']
        self.model.counts['actual_complex_RHS'] += self.preparation['actual_complex_RHS']
        self.model.guard = self.guard
        self.target = np.asarray(public['Y'])
        self.scale2 = float(np.vdot(self.target, self.target).real)
        if self.scale2 <= 0 or not np.isfinite(self.scale2):
            raise ValueError('Nonzero finite observed Y required')
        self.anchor = np.asarray(anchor['m0'] if isinstance(anchor, dict) else anchor, float).reshape(576)
        if not np.all(np.isfinite(self.anchor)):
            raise ValueError('Nonfinite anchor')
        self.reg = Regularizer((24, 24), smooth=regularization, tv=regularization)
        self.old_reg = Regularizer((24, 24), smooth=1e-4, tv=1e-4)
        self.best, self.trace, self.attempts, self.arrays = None, [], [], {}
        self.candidates = []

    def elapsed(self):
        return time.perf_counter() - self.started + self.preparation['wall_seconds']

    def guard(self):
        if self.elapsed() >= BUDGET['wall_seconds']:
            raise Stop('WALL_SECONDS')
        if self.model.counts['LU'] > BUDGET['max_LU']:
            raise Stop('MAX_LU')
        if self.model.counts['actual_complex_RHS'] > BUDGET['max_complex_RHS']:
            raise Stop('MAX_RHS')
        if self.elapsed() >= self.phase_end:
            raise PhaseStop('PHASE_SECONDS')
        if self._gate_enabled and self._gate_pending:
            raise PhaseStop('PHYSICAL_DATA_GATE')

    def evaluate(self, m, label):
        state, row = super().evaluate(m, label)
        row['old_data_TV_score'] = .5 * row['relative_residual'] ** 2 + self.old_reg.evaluate(m)[0]
        row['physical_data_feasible'] = bool(row['relative_residual'] <= self.data_gate)
        row['full_VIE_equation_defect'] = state['state_relative_residual']
        # R2 keeps local_best before its next guard; this stops at the first
        # feasible evaluated point and preserves its material and receipt.
        if self._gate_enabled and row['physical_data_feasible']:
            self._gate_pending = True
            self._first_feasible = dict(row)
        return state, row

    def correct(self, m, label, seconds=3., max_nfev=8):
        self._gate_pending = False
        self._first_feasible = None
        self._gate_enabled = True
        try:
            best = super().correct(m, label, seconds=seconds, max_nfev=max_nfev)
            return self._first_feasible if self._first_feasible is not None else best
        finally:
            self._gate_enabled = self._gate_pending = False

    def add_candidate(self, m, label, *, kind, state=None, row=None):
        if row is None:
            state, row = self.evaluate(m, label)
        item = {k:v for k,v in row.items() if k != 'm'}
        item.update(index=len(self.candidates), kind=kind, label=label,
                    distance_to_full_proposal=float(np.linalg.norm(m - self.arrays['m_proposal_unclipped'])))
        self.candidates.append(item)
        self.arrays[f'candidate_{item["index"]}_m'] = np.asarray(m).copy()
        return item

    def auxiliary(self, prediction):
        if 'j_aux' not in prediction:
            return None
        self.guard()
        j = np.asarray(prediction['j_aux'], complex)
        if j.shape != (576, 12) or not np.all(np.isfinite(j)):
            raise ValueError('Finite j_aux shape (576,12) required')
        with self.model.operation('AUXILIARY_CURRENT_DIAGNOSTICS', current_diagnostics=1):
            q = self.model.E0 + self.model.G @ j
            m = self.arrays['m_proposal_clipped']
            defect = j - m[:, None] * q
            out = dict(auxiliary_data_relative_residual=float(np.linalg.norm(self.model.S @ j - self.target) / np.sqrt(self.scale2)),
                       constitutive_relative_defect=float(np.linalg.norm(defect) / max(np.linalg.norm(j), 1e-300)),
                       field_definition='q = E0 + G @ j_aux; no auxiliary data fit is called full VIE')
            if 'q_aux' in prediction:
                qa = np.asarray(prediction['q_aux'], complex)
                if qa.shape != q.shape or not np.all(np.isfinite(qa)):
                    raise ValueError('Finite q_aux shape (576,12) required')
                out['q_aux_field_relative_defect'] = float(np.linalg.norm(qa - q) / max(np.linalg.norm(q), 1e-300))
            return out

    def softfilter_box_qp(self, p):
        """Full fixed-anchor J ridge box-QP; no weak subspace truncation."""
        self.guard()
        pnorm = float(np.linalg.norm(p))
        if pnorm <= 1e-30:
            return self.anchor.copy(), dict(lambda_value=None, status=0, message='Zero displacement: identity proposal, ridge ratio undefined and not evaluated', fixed_reference='GN anchor', no_weak_truncation=True)
        state, _ = self.evaluate(np.clip(self.anchor, LOWER, UPPER), 'softfilter_reference')
        J = self.model.jacobian(state) / np.sqrt(self.scale2)
        epsilon = 1e-4 if self.noise == 'clean' else 1e-3
        lam = (epsilon / pnorm) ** 2
        with self.model.operation('FULL_SOFTFILTER_BOX_QP', box_QP=1):
            gram = J.T @ J
            def objective(h):
                self.guard()
                d = h - p
                return .5 * float(h @ gram @ h) + .5 * lam * float(d @ d), gram @ h + lam * d
            fit = minimize(objective, np.clip(p, LOWER - self.anchor, UPPER - self.anchor),
                           jac=True, method='L-BFGS-B', bounds=list(zip(LOWER-self.anchor, UPPER-self.anchor)),
                           options={'maxiter': 1000, 'ftol': 1e-15, 'gtol': 1e-12})
        return self.anchor + fit.x, dict(lambda_value=lam, epsilon=epsilon, status=int(fit.status),
                                       message=str(fit.message), iterations=int(fit.nit),
                                       fixed_reference='full analytic J at GN anchor', no_weak_truncation=True)

    def run(self, prediction, *, arm, secondary_mode=None, control=None):
        status, reason, auxiliary, filter_receipt = 'COMPLETED', None, None, None
        original = np.asarray(prediction.get('m_proposal', self.anchor), float).reshape(576)
        if not np.all(np.isfinite(original)):
            raise ValueError('Finite m_proposal shape (576,) required')
        self.arrays.update(anchor=self.anchor.copy(), m_proposal_unclipped=original.copy(),
                           m_proposal_clipped=np.clip(original, LOWER, UPPER))
        clipping = dict(bounds=[LOWER, UPPER], proposal_clipped_cells=int(np.count_nonzero(original != np.clip(original, LOWER, UPPER))),
                        proposal_clip_L2=float(np.linalg.norm(original-np.clip(original, LOWER, UPPER))),
                        anchor_clipped_cells=int(np.count_nonzero(self.anchor != np.clip(self.anchor, LOWER, UPPER))))
        try:
            anchor = np.clip(self.anchor, LOWER, UPPER)
            self.add_candidate(anchor, 'anchor', kind='anchor')
            if control:
                best = self.correct(anchor, control, seconds=10., max_nfev=1000)
                if best is not None:
                    self.add_candidate(best['m'], control, kind='corrected_control', row=best)
            else:
                auxiliary = self.auxiliary(prediction)
                proposal = original.copy()
                if secondary_mode == 'full_softfilter_box_qp':
                    proposal, filter_receipt = self.softfilter_box_qp(original-self.anchor)
                    self.arrays['secondary_m_proposal'] = proposal
                # Both raw seeds are verified before correction, so all-best-K
                # has a complete fixed raw candidate set even if TRF stops.
                seeds = [(alpha, np.clip(anchor + alpha*(proposal-anchor), LOWER, UPPER)) for alpha in (.5, 1.)]
                raw = []
                for alpha, seed in seeds:
                    raw.append(self.add_candidate(seed, f'raw_{alpha:g}', kind='raw'))
                for (alpha, seed), before in zip(seeds, raw):
                    if before['physical_data_feasible']:
                        self.attempts.append(dict(label=f'correct_{alpha:g}', status='NOT_NEEDED_DATA_GATE', reason='PHYSICAL_DATA_GATE', nfev=0))
                        best = dict(before, m=seed)
                    else:
                        best = self.correct(seed, f'correct_{alpha:g}', seconds=3., max_nfev=8)
                    if best is not None:
                        self.add_candidate(best['m'], f'corrected_{alpha:g}', kind='corrected', row=best)
        except Stop as exc:
            status, reason = 'BUDGET_STOP', str(exc)
        except Exception as exc:
            status, reason = 'FAILED', f'{type(exc).__name__}: {exc}'
        failures = [a for a in self.attempts if a['status'] in ('FAILED', 'BUDGET_STOP')]
        if any(a['status']=='BUDGET_STOP' for a in failures):
            status, reason = 'BUDGET_STOP', next(a['reason'] for a in reversed(failures) if a['status']=='BUDGET_STOP')
        feasible = [c for c in self.candidates if c['physical_data_feasible'] and c['kind'] != 'anchor']
        anchor_rows = [c for c in self.candidates if c['kind']=='anchor' and c['physical_data_feasible']]
        if control:
            eligible = feasible + anchor_rows
            selected = min(eligible, key=lambda c: c['objective']) if eligible else None
        else:
            selected = min(feasible, key=lambda c: (c['distance_to_full_proposal'], c['index'])) if feasible else (anchor_rows[0] if anchor_rows else None)
        old = min(self.candidates, key=lambda c:c['old_data_TV_score']) if self.candidates else None
        for key, item in [('top1', selected), ('old_score_top1', old)]:
            if item is not None:
                self.arrays[key] = self.arrays[f'candidate_{item["index"]}_m'].copy()
        self.arrays['candidate_materials'] = np.asarray([self.arrays[f'candidate_{c["index"]}_m'] for c in self.candidates]).reshape(-1, 576)
        wall = self.elapsed()
        receipt = dict(schema='gn-nn-completion-r3-evaluation-v1', arm=arm, noise=self.noise,
                       status=status, stop_reason=reason, selection_status='PASS' if selected else 'FAIL_NO_FEASIBLE_CANDIDATE',
                       accounting_status='COMPLETE' if self.preparation['complete'] else 'HOLD_INCOMPLETE_PREPARATION',
                       truth_access='NONE_PUBLIC_INPUT_ONLY', data_gate=self.data_gate,
                       wall_seconds=wall, evaluator_wall_seconds=time.perf_counter()-self.started,
                       cpu_seconds=time.process_time()-self.cpu, preparation_cost=self.preparation,
                       budget=BUDGET, wall_overshoot_seconds=max(0.,wall-BUDGET['wall_seconds']),
                       costcounts=dict(self.model.counts), operation_seconds=dict(self.model.timings),
                       clipping=clipping, auxiliary=auxiliary, secondary_filter=filter_receipt,
                       candidates=self.candidates, attempts=self.attempts, trace=self.trace,
                       top1=selected, old_score_top1=old,
                       online_selection='physical-feasible material closest in L2 to full unclipped NN proposal; otherwise feasible anchor; otherwise FAIL',
                       old_selection='minimum original 0.5 physical relative residual squared + smooth1e-4 + TV1e-4',
                       anchor_construction_accounting='separate common sunk cost; not charged twice to incremental 10s',
                       control=control)
        return safe(receipt), self.arrays


def evaluate_arm(public, anchor, prediction, *, arm, noise='clean', preparation_cost=None, secondary_mode=None):
    with threadpool_limits(limits=1):
        return Evaluation(public, anchor, noise=noise, preparation_cost=preparation_cost).run(prediction, arm=arm, secondary_mode=secondary_mode)


def evaluate_control(public, anchor, *, control='continue_trf', noise='clean'):
    if control not in ('continue_trf', 'weakreg_trf'):
        raise ValueError('Unknown continuation control')
    with threadpool_limits(limits=1):
        evaluator = Evaluation(public, anchor, noise=noise,
                               preparation_cost=dict(wall_seconds=0., LU=0, actual_complex_RHS=0, complete=True),
                               regularization=1e-5 if control=='weakreg_trf' else 0.)
        return evaluator.run({}, arm=control, control=control)


def save_result(directory, obs_id, receipt, arrays):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    output = directory / f'{obs_id}.npz'
    np.savez_compressed(output, **arrays)
    (directory / f'{obs_id}.json').write_text(json.dumps(safe(receipt), indent=2, allow_nan=False)+'\n')
    return output
