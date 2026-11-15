"""Paid post-hoc full-state metric interventions; never a recovery algorithm.

Restore one completed R14 state, switch only its two metric families and
optionally reset PR memory, and compare fixed numbers of native updates.
Central perturbations alter only the real material; A/J/Q and PR are retained.
No truth, corrector, optimization tuning, or new physics backend is used.
"""
import os
for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
             'VECLIB_MAXIMUM_THREADS', 'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_key] = '1'

import argparse
from collections import Counter
from copy import deepcopy
import json
from pathlib import Path
import sys
import time
import traceback

import numpy as np


_PUBLIC_KEYS = {
    'G', 'E', 'E0', 'U0', 'S', 'H', 'Y', 'sigma', 'delta', 'm_init',
    'lower', 'upper', 'grid_shape', 'coords', 'cell_radius', 'k',
    'operators_path', 'centers', 'sources', 'receivers', 'radius',
}
_RESPONSE_FIELDS = ('m', 'J', 'Q', 'A', 'old_grad', 'old_direction')


def _json_value(value):
    if isinstance(value, np.ndarray):
        return _json_value(value.tolist())
    if isinstance(value, np.generic):
        return _json_value(value.item())
    if isinstance(value, dict):
        return {str(k): _json_value(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(v) for v in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, complex):
        return dict(real=_json_value(value.real), imag=_json_value(value.imag))
    if isinstance(value, float) and not np.isfinite(value):
        return dict(nonfinite=repr(value))
    return value


def _dump(path, value):
    Path(path).write_text(json.dumps(_json_value(value), indent=2, allow_nan=False)+'\n', encoding='utf8')


def _delta(after, before):
    return {key: after.get(key, 0)-before.get(key, 0) for key in sorted(set(after)|set(before))
            if isinstance(after.get(key, 0), (int, float, np.number))
            and isinstance(before.get(key, 0), (int, float, np.number))
            and after.get(key, 0) != before.get(key, 0)}


def _norm(value):
    return float(np.linalg.norm(value))


def _relative(a, b):
    return _norm(a-b)/max(_norm(a), _norm(b), 1e-300)


def _steps(value):
    values = sorted(set(int(v) for v in value.split(',')))
    if not values or values[0] < 1:
        raise argparse.ArgumentTypeError('--steps must contain positive comma-separated update counts')
    return values


class InterventionDiagnostics:
    def __init__(self, args):
        self.args = args
        self.root = args.root.resolve()
        self.out = args.out.resolve()
        self.out.mkdir(parents=True, exist_ok=False)
        self.cpu_started, self.wall_started = time.process_time(), time.perf_counter()
        self.model = None
        self.events, self.failures, self.branches = [], [], []
        self.perturbations = []
        self.repeat_control = None
        self.input_keys_read = []
        self.status = 'RUNNING'

    def charge(self, label, function, *, solver=None, metadata=None):
        before_model = dict(self.model.cost) if self.model is not None else {}
        before_aux = dict(solver.cost) if solver is not None else {}
        cpu, wall = time.process_time(), time.perf_counter()
        event = dict(event=len(self.events), label=label, metadata=metadata)
        try:
            result = function()
            event['status'] = 'COMPLETE'
            return result
        except BaseException as exc:
            event.update(status='FAILED', error=repr(exc))
            raise
        finally:
            after_model = dict(self.model.cost) if self.model is not None else {}
            after_aux = dict(solver.cost) if solver is not None else {}
            event.update(CPU_seconds=time.process_time()-cpu, wall_seconds=time.perf_counter()-wall,
                process_CPU_seconds=time.process_time(),
                physical_operation_delta=_delta(after_model, before_model),
                auxiliary_operation_delta=_delta(after_aux, before_aux))
            self.events.append(event)

    def failure(self, label, exc):
        row = dict(label=label, error=repr(exc), traceback=traceback.format_exc(),
                   process_CPU_seconds=time.process_time())
        self.failures.append(row)
        _dump(self.out/'FAILURES.json', self.failures)

    def setup(self):
        sys.path[:0] = [str(self.root/'src'), str(self.args.legacy_repo.resolve())]
        from research.research_fbg_r3_experiments.src.physics import MaxwellModel
        from dualmetric_core import StatefulDualMetric, _load_state_npz
        from run_experiment import build_config
        self.Core = StatefulDualMetric
        def load_public():
            with np.load(self.args.input, allow_pickle=False) as archive:
                keys = [key for key in archive.files if key in _PUBLIC_KEYS]
                public = {key: archive[key].copy() for key in keys}
            self.input_keys_read.extend(keys)
            if 'operators_path' in public:
                operator_path = self.args.input.parent/str(public['operators_path'].item())
                with np.load(operator_path, allow_pickle=False) as archive:
                    keys = [key for key in archive.files if key in _PUBLIC_KEYS]
                    operators = {key: archive[key].copy() for key in keys}
                self.input_keys_read.extend('operators:'+key for key in keys)
                public = {**operators, **public}
            if 'E0' not in public:
                public['E0'] = public.get('E', public.get('U0'))
            return public
        self.public = self.charge('load public observations/operators only', load_public)
        def build_model():
            self.model = MaxwellModel(self.public)
        self.charge('original MaxwellModel setup', build_model)
        self.Y = np.asarray(self.public['Y'], complex)
        self.normY = _norm(self.Y)
        if not np.isfinite(self.normY) or self.normY <= 0:
            raise ValueError('Physical data normalization requires nonzero finite Y')
        self.config, self.rank_decision = self.charge('original build_config including DPCS selection',
            lambda: build_config(self.model, self.public, self.out))
        self.base_state = self.charge('read full source checkpoint NPZ',
            lambda: _load_state_npz(self.args.checkpoint))
        base = self.restore(self.base_state, label='base full-state restore')
        self.base_state = self.charge('detach complete base state', base.state_dict, solver=base)
        self.lo, self.hi = base.lo.copy(), base.hi.copy()
        self.base_iteration = int(base.iterations)
        self.base_pair = base.pair
        self.charge('save detached base checkpoint', lambda: base.checkpoint(self.out/'BASE_STATE.npz'), solver=base)
        _dump(self.out/'CONFIGURATION.json', dict(truth_access='NONE',
            configuration_source='original run_experiment.build_config plus saved complete checkpoint numerical flags',
            checkpoint=str(self.args.checkpoint.resolve()), input=str(self.args.input.resolve()),
            checkpoint_pair=base.pair, checkpoint_global_iteration=base.iterations,
            checkpoint_stage=base.stage, checkpoint_iteration=base.iteration,
            checkpoint_stage_pending=base.stage_pending, checkpoint_done=base.done,
            rank_decision=self.rank_decision, numerical_flags=self.base_state['frozen_config'],
            initialization='restore actual complete state; never restart from m alone',
            fixed_update_counts=self.args.steps, diagnostic_only=True))
        self.base_snapshot = self.snapshot(base, self.out/'BASE_PHYSICAL.npz', row=None,
            label='base paid physical state')
        # Diagnostics must not silently alter the state from which interventions
        # start.  The saved base contains exactly the pre-diagnostic counters.
        return base

    def restore(self, state, *, pair=None, reset=False, label='full-state restore'):
        return self.charge(label, lambda: self.Core.from_state_dict(self.model, self.Y,
            self.config, state, pair=pair, reset_memory=reset),
            metadata=dict(pair=state['pair'] if pair is None else pair, reset_PR=reset,
                          source_global_iteration=int(state['iterations']),
                          actions_recomputed='NONE; restore cached DPCS/FFT/block/state arrays'))

    def active_masks(self, m):
        return np.asarray(m <= self.lo+1e-8), np.asarray(m >= self.hi-1e-8)

    def branch_signature(self, row, solver):
        signature = {key: row.get(key) for key in (
            'direction_restarted', 'material_fallback', 'material_inner_success',
            'material_inner_nit', 'aplus_sweeps', 'aplus_inherited_gtol',
            'stage', 'iteration', 'stage_end', 'stage_stop')}
        signature['PR_plus_zero'] = np.asarray(row.get('PR_plus', 0.)) == 0
        old = row.get('material_old_block_loss')
        seed = row.get('material_scalar_seed_block_loss')
        signature['material_scalar_seed_preferred'] = None if old is None or seed is None else bool(seed < old)
        lo, hi = self.active_masks(solver.m)
        signature['lower_active_mask'] = lo
        signature['upper_active_mask'] = hi
        return signature

    @staticmethod
    def signature_difference(a, b):
        return [key for key in sorted(set(a)|set(b))
                if not np.array_equal(np.asarray(a.get(key)), np.asarray(b.get(key)))]

    def snapshot(self, solver, path, *, row, label):
        """Original material Phi/gradient, physical fields, and cached defects."""
        path = Path(path)
        def action():
            arrays = {}
            before_metrics = None
            angle = cosine = None
            angle_reason = 'initial snapshot has no material step'
            if row is not None:
                m_before = np.asarray(row['m_before'], float)
                phi_before, gradient_before = self.model.data_value_gradient(m_before, self.Y)
                physical_before = self.model.state(m_before)
                arrays.update(m_before=m_before.copy(), physical_gradient_before=gradient_before.copy(),
                    physical_j_before=physical_before['j'].copy(), physical_Q_before=physical_before['Q'].copy(),
                    physical_Y_before=physical_before['Y'].copy())
                residual_before = physical_before['Y']-self.Y
                R_before = float(self.model.reg.evaluate(m_before)[0])
                before_metrics = dict(Phi=float(phi_before), D=.5*_norm(residual_before)**2/self.normY**2, R=R_before)
                step = solver.m-m_before
                denominator = _norm(step)*_norm(gradient_before)
                if denominator > 0:
                    cosine = float(np.clip(np.dot(step, -gradient_before)/denominator, -1., 1.))
                    angle = float(np.degrees(np.arccos(cosine)))
                    angle_reason = None
                else:
                    angle_reason = 'zero material step or zero original-Phi gradient'
            Phi, gradient = self.model.data_value_gradient(solver.m, self.Y)
            physical = self.model.state(solver.m)
            residual = physical['Y']-self.Y
            D = .5*_norm(residual)**2/self.normY**2
            R = float(self.model.reg.evaluate(solver.m)[0])
            port = solver.S@solver.J-self.Y
            solver.cost['S_columns'] += solver.ns
            defect = solver.J-solver.m[:, None]*solver.Q
            transformed_current = solver._P(solver.m, defect, solver.pair[0])
            transformed_material = solver._P(solver.m, defect, solver.pair[1])
            state_current = .5*solver.rho*float(np.vdot(transformed_current, transformed_current@solver.Os).real)
            state_material = .5*solver.rho*float(np.vdot(transformed_material, transformed_material@solver.Os).real)
            aux_data = .5*float(np.vdot(port, port@solver.Od).real)
            solver.cost['source_metric_products'] += 3
            lower, upper = self.active_masks(solver.m)
            arrays.update(m=solver.m.copy(), physical_j=physical['j'].copy(), physical_Q=physical['Q'].copy(),
                physical_Y=physical['Y'].copy(), physical_residual=residual.copy(), physical_Phi_gradient=gradient.copy(),
                auxiliary_J=solver.J.copy(), auxiliary_Q=solver.Q.copy(), auxiliary_port_residual=port,
                auxiliary_constitutive_defect=defect, auxiliary_current_transformed_defect=transformed_current,
                auxiliary_material_transformed_defect=transformed_material, lower_active_mask=lower, upper_active_mask=upper)
            np.savez_compressed(path, **arrays)
            return dict(Phi=float(Phi), D=D, R=R, physical_relative_residual=_norm(residual)/self.normY,
                Phi_decomposition_error=float(Phi-D-R), physical_Phi_gradient_norm=_norm(gradient),
                auxiliary_port_relative_residual=_norm(port)/self.normY,
                auxiliary_port_original_normalized_D=.5*_norm(port)**2/self.normY**2,
                auxiliary_data_source_metric_loss=aux_data,
                auxiliary_constitutive_defect_norm=_norm(defect),
                auxiliary_constitutive_relative_to_J=_norm(defect)/max(_norm(solver.J), 1e-300),
                auxiliary_constitutive_relative_to_E=_norm(defect)/max(_norm(solver.E), 1e-300),
                auxiliary_state_current_metric_loss=state_current, auxiliary_state_material_metric_loss=state_material,
                lower_active=int(lower.sum()), upper_active=int(upper.sum()),
                material_step_angle_to_negative_physical_Phi_gradient_degrees=angle,
                material_step_cosine_to_negative_physical_Phi_gradient=cosine,
                material_step_angle_null_reason=angle_reason, physical_before=before_metrics,
                legal=bool(np.isfinite(solver.m).all() and np.all(solver.m >= self.lo) and np.all(solver.m <= self.hi)),
                global_iteration=int(solver.iterations), stage=int(solver.stage),
                stage_iteration=int(solver.iteration), pair=solver.pair, arrays=str(path))
        metrics = self.charge(label, action, solver=solver)
        _dump(path.with_suffix('.json'), metrics)
        return metrics

    def run_branch(self, pair, reset):
        name = 'pair'+pair+('_resetPR' if reset else '_retainPR')
        folder = self.out/'branches'/name
        folder.mkdir(parents=True)
        summary = dict(name=name, pair=pair, reset_PR=reset, diagnostic_only=True,
            update_budget=max(self.args.steps), completed_updates=0, snapshots={}, status='RUNNING', steps=[])
        solver = None
        last_row = None
        try:
            solver = self.restore(self.base_state, pair=pair, reset=reset, label=name+' restore')
            self.charge(name+' initial state export', lambda: solver.checkpoint(folder/'step_00000_state.npz'), solver=solver)
            summary['snapshots']['0'] = self.snapshot(solver, folder/'step_00000_physical.npz',
                row=None, label=name+' initial paid physical audit')
            for count in range(1, max(self.args.steps)+1):
                try:
                    row = self.charge(name+' native update', solver.step, solver=solver,
                                      metadata=dict(relative_update=count))
                except StopIteration as exc:
                    summary.update(status='SCHEDULE_EXHAUSTED', schedule_stop=str(exc))
                    break
                last_row = row
                summary['completed_updates'] = count
                signature = self.branch_signature(row, solver)
                summary['steps'].append(dict(relative_update=count, global_iteration=int(solver.iterations),
                    stage=int(row['stage']), stage_iteration=int(row['iteration']),
                    native_current_objective_change=row['current_step_objective_change'],
                    material_step_norm=float(row['readout_material_change']),
                    branch_flags={k:v for k,v in signature.items() if not k.endswith('_mask')}))
                if count in self.args.steps:
                    stem = folder/f'step_{count:05d}'
                    self.charge(name+' full completed state export',
                        lambda: solver.checkpoint(str(stem)+'_state.npz'), solver=solver,
                        metadata=dict(relative_update=count))
                    summary['snapshots'][str(count)] = self.snapshot(solver, str(stem)+'_physical.npz',
                        row=row, label=name+f' paid physical audit after {count} updates')
            if summary['status'] == 'RUNNING':
                summary['status'] = 'COMPLETE_FIXED_UPDATES'
            count = summary['completed_updates']
            if str(count) not in summary['snapshots']:
                summary['snapshots'][str(count)] = self.snapshot(solver, folder/'final_physical.npz',
                    row=last_row, label=name+' final paid physical audit')
            summary['final_snapshot'] = summary['snapshots'][str(count)]
            summary['initial_snapshot'] = summary['snapshots']['0']
            summary['physical_Phi_change'] = summary['final_snapshot']['Phi']-summary['initial_snapshot']['Phi']
        except Exception as exc:
            summary.update(status='FAILED', error=repr(exc))
            self.failure(name, exc)
        finally:
            if solver is not None:
                try:
                    self.charge(name+' final full-state export', lambda: solver.checkpoint(folder/'FINAL_STATE.npz'), solver=solver)
                    summary['final_state'] = str(folder/'FINAL_STATE.npz')
                except Exception as exc:
                    self.failure(name+' final checkpoint', exc)
            _dump(folder/'BRANCH.json', summary)
            self.branches.append(summary)

    def repeated_identical_control(self):
        folder = self.out/'repeat_identical_step'
        folder.mkdir()
        result = dict(pair='01', reset_PR=False, status='RUNNING')
        try:
            twins = []
            for index in (1, 2):
                solver = self.restore(self.base_state, pair='01', label=f'repeat control {index} restore')
                row = self.charge(f'repeat control {index} identical native update', solver.step, solver=solver)
                self.charge(f'repeat control {index} full-state export',
                    lambda: solver.checkpoint(folder/f'twin_{index}_state.npz'), solver=solver)
                twins.append((solver, row))
            a, b = twins[0][0], twins[1][0]
            exact = {key: bool(np.array_equal(getattr(a, key), getattr(b, key))) for key in _RESPONSE_FIELDS}
            sig_a, sig_b = [self.branch_signature(row, solver) for solver, row in twins]
            flags = self.signature_difference(sig_a, sig_b)
            result.update(status='COMPLETE', arrays_bitwise_equal=exact, branch_differences=flags,
                all_equal=bool(all(exact.values()) and not flags),
                relative_errors={key:_relative(getattr(a,key),getattr(b,key)) for key in _RESPONSE_FIELDS})
            if not result['all_equal']:
                raise RuntimeError('Repeated identical native step control did not reproduce exactly')
        except Exception as exc:
            result.update(status='FAILED', error=repr(exc))
            self.failure('repeated identical step', exc)
        self.repeat_control = result
        _dump(folder/'CONTROL.json', result)

    def one_step_response(self, state, pair, folder, label):
        folder.mkdir(parents=True)
        solver = self.restore(state, pair=pair, label=label+' restore exact perturbed state')
        self.charge(label+' initial full-state export', lambda: solver.checkpoint(folder/'INITIAL_STATE.npz'), solver=solver)
        row = self.charge(label+' one native response update', solver.step, solver=solver)
        # Preserve next-step state before observational diagnostic counters.
        response = {key:getattr(solver, key).copy() for key in _RESPONSE_FIELDS}
        self.charge(label+' completed full-state export', lambda: solver.checkpoint(folder/'FINAL_STATE.npz'), solver=solver)
        signature = self.branch_signature(row, solver)
        metrics = self.snapshot(solver, folder/'PHYSICAL.npz', row=row,
                                label=label+' paid physical/gradient audit')
        result = dict(label=label, pair=pair, status='COMPLETE', signature=signature,
            physical=metrics, initial_state=str(folder/'INITIAL_STATE.npz'), final_state=str(folder/'FINAL_STATE.npz'))
        _dump(folder/'RESPONSE.json', result)
        return response, result

    def perturbation_responses(self):
        folder = self.out/'perturbations'
        folder.mkdir()
        base_m = self.base_state['m']
        free = np.minimum(base_m-self.lo, self.hi-base_m) > 1e-5
        def gradient_action():
            value, gradient = self.model.data_value_gradient(base_m, self.Y)
            return float(value), gradient.copy()
        Phi, gradient = self.charge('perturbation original-Phi gradient at non-root state', gradient_action)
        h1 = np.where(free, -gradient, 0.)
        rng = np.random.default_rng(14014)
        h2 = np.where(free, rng.normal(size=base_m.shape), 0.)
        norms = (_norm(h1), _norm(h2))
        reason = None
        if not np.any(free):
            reason = 'all coordinates are active or within 1e-5 of a box boundary'
        elif norms[0] == 0:
            reason = 'original-Phi gradient restricted to strictly interior coordinates is zero'
        elif norms[1] == 0:
            reason = 'seeded random interior direction has zero norm'
        direction_info = dict(Phi=Phi, full_physical_gradient_norm=_norm(gradient),
            restricted_physical_gradient_norm=norms[0], strictly_interior_coordinates=int(free.sum()),
            random_seed=14014, direction1='negative original-Phi physical gradient restricted to strictly interior coordinates, unit Euclidean norm',
            direction2='seeded standard-normal material direction on the same strictly interior coordinates, unit Euclidean norm',
            direction_sum='h1+h2, deliberately not renormalized', null_reason=reason)
        if reason is not None:
            np.savez_compressed(folder/'DIRECTIONS.npz', m=base_m, physical_Phi_gradient=gradient, free=free)
            _dump(folder/'DIRECTIONS.json', direction_info)
            self.perturbation_summary = dict(status='NULL_DIRECTIONS', reason=reason, direction_info=direction_info)
            return
        h1 /= norms[0]
        h2 /= norms[1]
        directions = dict(h1=h1, h2=h2, hsum=h1+h2)
        np.savez_compressed(folder/'DIRECTIONS.npz', m=base_m, physical_Phi_gradient=gradient,
                            free=free, h1=h1, h2=h2, hsum=h1+h2)
        _dump(folder/'DIRECTIONS.json', direction_info)
        derivatives, jvp_arrays, outcomes = {}, {}, []
        epsilons = (1e-3, 3e-4)
        for pair in ('01', '11'):
            baseline_response = baseline_info = None
            try:
                baseline_response, baseline_info = self.one_step_response(self.base_state, pair,
                    folder/('pair'+pair)/'unperturbed', 'pair'+pair+' unperturbed response')
            except Exception as exc:
                self.failure('pair'+pair+' unperturbed response', exc)
            for direction_name, h in directions.items():
                maxh = float(np.max(abs(h)))
                nonzero = abs(h) > 0
                if maxh == 0 or not np.any(nonzero):
                    outcomes.append(dict(pair=pair, direction=direction_name, status='NULL_DIRECTION'))
                    continue
                boundmargin = float(np.min(np.minimum(base_m-self.lo, self.hi-base_m)[nonzero]))
                for declared in epsilons:
                    actual = min(declared, .25*boundmargin/maxh)
                    entry = dict(pair=pair, direction=direction_name, declared_epsilon=declared,
                        actual_epsilon=actual, bound_margin=boundmargin, max_abs_direction=maxh,
                        epsilon_rule='min(declared, .25*minimum participating box margin/maxabs(h)); no clipping',
                        same_A_J_Q_PR=True, status='RUNNING')
                    responses, infos = {}, {}
                    for sign, multiplier in (('plus', 1.), ('minus', -1.)):
                        label = f'pair{pair} {direction_name} eps{declared:g} {sign}'
                        destination = folder/('pair'+pair)/direction_name/f'eps_{declared:g}'/sign
                        try:
                            def perturb_state():
                                state = deepcopy(self.base_state)
                                state['m'] = base_m+multiplier*actual*h
                                if np.any(state['m'] < self.lo) or np.any(state['m'] > self.hi):
                                    raise ValueError('Declared feasible perturbation left original material box')
                                for key in ('A', 'J', 'Q', 'old_grad', 'old_direction'):
                                    if not np.array_equal(state[key], self.base_state[key]):
                                        raise RuntimeError('Material-only perturbation altered '+key)
                                return state
                            state = self.charge(label+' material-only full-state copy', perturb_state,
                                                metadata=dict(actual_epsilon=actual, sign=multiplier))
                            responses[sign], infos[sign] = self.one_step_response(state, pair, destination, label)
                        except Exception as exc:
                            self.failure(label, exc)
                            entry.setdefault('failed_signs', []).append(dict(sign=sign, error=repr(exc)))
                    if len(responses) == 2:
                        derivative = {key:(responses['plus'][key]-responses['minus'][key])/(2*actual)
                                      for key in _RESPONSE_FIELDS}
                        derivatives[(pair, direction_name, declared)] = derivative
                        differences = self.signature_difference(infos['plus']['signature'], infos['minus']['signature'])
                        baseline_differences = None if baseline_info is None else {
                            sign:self.signature_difference(infos[sign]['signature'], baseline_info['signature'])
                            for sign in ('plus','minus')}
                        nonlinear_even = None if baseline_response is None else {
                            key:_norm(.5*(responses['plus'][key]+responses['minus'][key])-baseline_response[key])
                            for key in _RESPONSE_FIELDS}
                        nonsmooth = bool(differences or (baseline_differences and any(baseline_differences.values())))
                        entry.update(status='COMPLETE', central_JVP_norms={key:_norm(value) for key,value in derivative.items()},
                            plus_minus_branch_differences=differences, baseline_branch_differences=baseline_differences,
                            nonsmooth_branch_or_active_change=nonsmooth, central_even_response_norms=nonlinear_even,
                            branch_observability='available PR/readout/stage flags and full material active masks; internal unrecorded branches cannot be certified absent')
                        prefix=f'pair{pair}__{direction_name}__eps_{declared:g}'.replace('.', 'p').replace('-', 'm')
                        for key, value in derivative.items():
                            jvp_arrays[prefix+'__'+key] = value
                    else:
                        entry['status'] = 'FAILED_RESPONSE'
                    outcomes.append(entry)
                    self.perturbations = outcomes
                    _dump(folder/'RESPONSES.json', outcomes)
        comparisons = []
        outcome_lookup = {(v['pair'],v['direction'],v['declared_epsilon']):v for v in outcomes
                          if 'declared_epsilon' in v}
        for pair in ('01', '11'):
            for direction_name in directions:
                coarse = derivatives.get((pair, direction_name, epsilons[0]))
                fine = derivatives.get((pair, direction_name, epsilons[1]))
                if coarse is not None and fine is not None:
                    metadata = [outcome_lookup[(pair,direction_name,eps)] for eps in epsilons]
                    actuals = [entry['actual_epsilon'] for entry in metadata]
                    comparisons.append(dict(kind='central_JVP_epsilon_convergence', pair=pair,
                        direction=direction_name, declared_epsilons=list(epsilons),
                        actual_epsilons=actuals, actual_epsilon_refined=bool(actuals[1] < actuals[0]),
                        nonsmooth_branch_or_active_change=any(entry['nonsmooth_branch_or_active_change'] for entry in metadata),
                        interpretation='epsilon comparison cannot establish convergence if the feasible caps make the actual epsilons identical',
                        relative_difference={key:_relative(coarse[key],fine[key]) for key in _RESPONSE_FIELDS},
                        diagnostic_only=True))
            for declared in epsilons:
                values = [derivatives.get((pair,name,declared)) for name in ('h1','h2','hsum')]
                if all(value is not None for value in values):
                    first, second, combined = values
                    metadata = [outcome_lookup[(pair,name,declared)] for name in ('h1','h2','hsum')]
                    result = dict(kind='directional_additivity', pair=pair, declared_epsilon=declared,
                        direction='h1+h2 without renormalization',
                        actual_epsilons={name:entry['actual_epsilon'] for name,entry in zip(('h1','h2','hsum'),metadata)},
                        nonsmooth_branch_or_active_change=any(entry['nonsmooth_branch_or_active_change'] for entry in metadata),
                        relative_defect={key:_norm(combined[key]-first[key]-second[key])/
                            max(_norm(combined[key]),_norm(first[key])+_norm(second[key]),1e-300)
                            for key in _RESPONSE_FIELDS}, diagnostic_only=True)
                    comparisons.append(result)
                    prefix=f'pair{pair}__additivity__eps_{declared:g}'.replace('.', 'p').replace('-', 'm')
                    for key in _RESPONSE_FIELDS:
                        jvp_arrays[prefix+'__'+key] = combined[key]-first[key]-second[key]
        self.charge('save central-JVP and additivity arrays',
                    lambda: np.savez_compressed(folder/'JVP_ARRAYS.npz', **jvp_arrays))
        self.perturbation_summary = dict(status='COMPLETE' if not any(v['status'].startswith('FAILED') for v in outcomes) else 'PARTIAL_WITH_FAILURES',
            direction_info=direction_info, responses=outcomes, comparisons=comparisons,
            arrays=str(folder/'JVP_ARRAYS.npz'),
            interpretation='post-hoc finite-update response only; branch changes flag possible nonsmoothness, small defects do not prove a theorem')
        _dump(folder/'SUMMARY.json', self.perturbation_summary)

    def finish(self):
        auxiliary_total = Counter()
        for event in self.events:
            auxiliary_total.update(event['auxiliary_operation_delta'])
        completed = [b for b in self.branches if b['status'] == 'COMPLETE_FIXED_UPDATES']
        qa = dict(
            what_was_changed='Only current/material metric pair; retain/reset control additionally clears old PR gradient/direction.',
            were_states_identical_before_switch='All branches restore the same saved full m,A,J,Q,DPCS/source metric/stage/PR state; reset branches only clear PR memory.',
            is_this_a_new_recovery_algorithm=False,
            are_comparisons_equal_CPU=False,
            comparison_scope='fixed native update counts from one exposed post-hoc checkpoint; all actual costs charged separately',
            GT_access='NONE; input reads use an explicit physical-public-field allowlist',
            completed_fixed_count_branches=len(completed), requested_branches=8,
            repeated_identical_step_control=self.repeat_control,
            original_physical_objective='Phi=.5||Y(m)-Yobs||F^2/||Yobs||F^2+original material regularizer',
            material_step_angle='latest material displacement versus negative original-Phi gradient at m_before, using paid original Maxwell adjoint',
            physical_Phi_changes={b['name']:b.get('physical_Phi_change') for b in self.branches},
            perturbation_status=getattr(self,'perturbation_summary',{}).get('status','NOT_COMPLETED'),
            root_scope='finite non-root one-update response, subject to recorded branch/active-set changes; no basin or GT certificate',
            failures=len(self.failures))
        _dump(self.out/'Q_AND_A.json', qa)
        _dump(self.out/'BRANCHES.json', self.branches)
        _dump(self.out/'FAILURES.json', self.failures)
        _dump(self.out/'COST_LEDGER.json', self.events)
        receipt = dict(status=self.status, diagnostic_only=True, truth_access='NONE',
            checkpoint=str(self.args.checkpoint.resolve()), input=str(self.args.input.resolve()),
            root=str(self.root), steps=self.args.steps, process_CPU_seconds=time.process_time(),
            diagnostic_CPU_seconds=time.process_time()-self.cpu_started,
            diagnostic_wall_seconds=time.perf_counter()-self.wall_started,
            physical_operation_counts=dict(self.model.cost) if self.model is not None else {},
            diagnostic_auxiliary_operation_counts=dict(auxiliary_total),
            setup_operation_counts=dict(receiver_SVD=int(getattr(self,'rank_decision',{}).get('new_receiver_SVD',0))),
            historical_checkpoint_counts_not_recharged=True, input_keys_read=self.input_keys_read,
            CPU_scope='whole process imports plus setup/restore/copies/actions/physical audits/failed attempts/output; end-of-receipt serialization follows this timestamp',
            fixed_update_branches=len(self.branches), failures=len(self.failures),
            perturbation_summary='perturbations/SUMMARY.json', cost_ledger='COST_LEDGER.json',
            Q_and_A='Q_AND_A.json', currency='no metered external service', GPU_seconds=0)
        _dump(self.out/'RECEIPT.json', receipt)
        print(json.dumps(dict(status=self.status, out=str(self.out), branches=len(self.branches),
                              failures=len(self.failures), CPU_seconds=receipt['diagnostic_CPU_seconds'])))

    def run(self):
        try:
            self.setup()
            for pair in ('00','01','10','11'):
                for reset in (False,True):
                    self.run_branch(pair, reset)
            self.repeated_identical_control()
            self.perturbation_responses()
            self.status = 'COMPLETE' if not self.failures else 'PARTIAL_WITH_FAILURES'
        except BaseException as exc:
            self.status = 'FAILED' if isinstance(exc, Exception) else 'INTERRUPTED'
            self.failure('top-level diagnostic', exc)
        finally:
            self.finish()
        return 0 if self.status == 'COMPLETE' else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True, help='R14 experiment directory')
    parser.add_argument('--legacy-repo', type=Path, required=True)
    parser.add_argument('--input', type=Path, required=True, help='public observation/operator NPZ')
    parser.add_argument('--checkpoint', type=Path, required=True, help='actual completed full-state NPZ')
    parser.add_argument('--out', type=Path, required=True, help='new output directory; existing paths are rejected')
    parser.add_argument('--steps', type=_steps, default=_steps('1,8,32'))
    args = parser.parse_args()
    return InterventionDiagnostics(args).run()


if __name__ == '__main__':
    sys.exit(main())
