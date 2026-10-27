"""Instrument original R4 optimizers without changing their objectives or steps.

No Book is constructed here. The caller supplies the sole paid R6 Target. The
temporary SciPy wrapper only observes accepted iterations and optionally raises
the documented StopIteration callback signal. It never evaluates the forward
model, Jacobian, gradient, or a counterfactual future endpoint for a decision.
The original modules retain their own bounded Point caches. Weak references here
cannot retain historical LU factors. Physical calls remain on the supplied Target.
"""
from __future__ import annotations
import os
for _name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
              'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_name] = '1'
import json
import time
import weakref
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from pathlib import Path
import numpy as np


@dataclass(frozen=True)
class AnchorConfig:
    name: str
    native_cap: int
    probit_cap: int
    adaptive: bool = False
    minimum_native_states: int = 24
    minimum_probit_states: int = 16
    stable_iterations: int = 3
    residual_relative_change: float = 1e-3
    material_relative_change: float = 5e-3


CONFIGS = {c.name: c for c in (
    AnchorConfig('native96_probit64', 96, 64),
    AnchorConfig('native24_probit16', 24, 16),
    AnchorConfig('native48_probit32', 48, 32),
    AnchorConfig('native72_probit48', 72, 48),
    AnchorConfig('adaptive_native72_probit48', 72, 48, adaptive=True),
)}


def _plain(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(k): _plain(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [_plain(v) for v in value]
    return value


def _write(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(_plain(data), indent=2, allow_nan=False) + '\n')


class _BookView:
    def __init__(self, original, jobdir):
        self.original, self.jobdir = original, Path(jobdir)

    def __getattr__(self, name):
        return getattr(self.original, name)


class _ObservedTarget:
    def __init__(self, original, destination):
        self.original = original
        self.book = _BookView(original.book, Path(destination) / 'original_optimizer')
        self.points = weakref.WeakValueDictionary()
        self.attempts = 0
        self.jacobian_invocations = 0

    def __getattr__(self, name):
        return getattr(self.original, name)

    def evaluate(self, q, purpose='target'):
        self.attempts += 1
        point = self.original.evaluate(q, purpose=purpose)
        self.points[np.asarray(point.q, float).tobytes()] = point
        return point

    def jacobian(self, point):
        self.jacobian_invocations += 1
        self.points[np.asarray(point.q, float).tobytes()] = point
        return self.original.jacobian(point)


class _Trace:
    def __init__(self, target, y, parent, folder, config, decision_predictor=None):
        self.target, self.y, self.parent = target, y, str(parent)
        self.folder, self.config = Path(folder), config
        self.folder.mkdir(parents=True, exist_ok=True)
        self.stream = (self.folder / 'PATH.jsonl').open('x')
        self.start = time.perf_counter()
        self.ledger_start = self.target.book.totals()
        self.phase = None
        self.phase_start = 0
        self.iteration = 0
        self.previous = None
        self.stable = 0
        self.stops = []
        self.decision_predictor = decision_predictor

    def begin(self, phase):
        self.phase, self.phase_start = phase, self.target.attempts
        self.iteration, self.previous, self.stable = 0, None, 0

    def observe(self, intermediate_result):
        # SciPy supplies accepted x/fun/cost/nfev, but not generally its current
        # gradient/Jacobian. Only already cached Point.J is used below.
        x = np.asarray(intermediate_result.x, float)
        if self.phase == 'probit':
            from shared.physics import from_probit
            q = from_probit(x)
        else:
            q = x.copy()
        point = self.target.points.get(q.tobytes())
        self.iteration += 1
        row = dict(event='ACCEPTED_ITERATION', parent=self.parent, phase=self.phase,
                   iteration=self.iteration, fresh_state_attempts=self.target.attempts,
                   phase_fresh_state_attempts=self.target.attempts-self.phase_start,
                   wall_seconds=time.perf_counter()-self.start,
                   scipy_nfev=int(intermediate_result.nfev),
                   scipy_cost=float(intermediate_result.cost),
                   scipy_gradient_available=hasattr(intermediate_result, 'grad'),
                   cached_point_available=point is not None,
                   residual_relative_change=None, material_relative_change=None,
                   scaled_gradient_inf=None, gradient_source='MISSING_NO_EXTRA_SOLVE',
                   current_output_source=None, stable_iterations=0,
                   online_stop=False, future_endpoint_used=False, GT_used=False)
        ledger = self.target.book.totals()
        row['ledger_sequence'] = ledger['sequence']
        row['path_resource_cost'] = {key: ledger[key]-self.ledger_start[key]
                                    for key in ('RHS', 'forward_RHS', 'derivative_RHS',
                                                'states', 'LU', 'physics_wall_seconds')}
        current = None
        if point is not None:
            residual = np.asarray(self.target.white(point.state.data-self.y))
            chi = np.asarray(point.state.chi)
            output = chi
            source = 'decoded_material_update_proxy_NOT_RISK_CALIBRATED'
            if self.decision_predictor is not None and point.J is not None:
                output = np.asarray(self.decision_predictor(self.y, point, self.parent))
                source = 'frozen_material_action_update_proxy_NOT_RISK_CALIBRATED'
            current = (residual.copy(), output.copy(), q.copy())
            snapshot = f'{self.phase}_{self.iteration:04d}.npz'
            np.savez_compressed(self.folder / snapshot, q=q, x=x, chi=chi,
                                prediction=point.state.data, residual=self.y-point.state.data,
                                current_output=output)
            row.update(endpoint=point.key, output_file=snapshot,
                       whitened_residual_squared=float(residual@residual),
                       phi=float(self.target.phi(point, self.y)),
                       native_objective=float(self.target.objective(point, self.y)),
                       current_output_source=source, cached_J_available=point.J is not None)
            if point.J is not None:
                if self.phase == 'probit':
                    from shared.physics import probit_derivative
                    J = point.J * probit_derivative(x)[None, :]
                    prior_residual, prior_J = x, np.eye(66)
                else:
                    J = point.J
                    prior_residual = q[:64]
                    prior_J = np.column_stack((np.eye(64), np.zeros((64, 2))))
                augmented_J = np.vstack((J, prior_J))
                gradient = augmented_J.T @ np.r_[residual, prior_residual]
                if self.phase == 'native':
                    # Exact same native property-box projected gradient as R4.
                    for j in (64, 65):
                        if ((q[j] <= 1e-8 and gradient[j] > 0) or
                                (q[j] >= 1-1e-8 and gradient[j] < 0)):
                            gradient[j] = 0
                scale = np.sqrt(np.maximum((augmented_J**2).sum(0), 1.))
                row.update(scaled_gradient_inf=float(np.max(np.abs(gradient)/scale)),
                           gradient_source='CURRENT_CACHED_J_ARITHMETIC_ONLY')
            if self.previous is not None:
                r0, o0, q0 = self.previous
                row['residual_relative_change'] = float(abs(residual@residual-r0@r0) /
                                                       max(1., float(r0@r0)))
                row['material_relative_change'] = float(np.linalg.norm(output-o0) /
                                                       max(1e-12, np.linalg.norm(o0)))
                row['native_step_l2'] = float(np.linalg.norm(q-q0))
                stable = (row['residual_relative_change'] <= self.config.residual_relative_change
                          and row['material_relative_change'] <= self.config.material_relative_change)
                self.stable = self.stable+1 if stable else 0
            minimum = (self.config.minimum_native_states if self.phase == 'native'
                       else self.config.minimum_probit_states)
            stop = (self.config.adaptive and self.target.attempts-self.phase_start >= minimum
                    and self.stable >= self.config.stable_iterations)
            row.update(stable_iterations=self.stable, online_stop=stop)
            self.previous = current
        self.stream.write(json.dumps(_plain(row), allow_nan=False)+'\n')
        self.stream.flush()
        if row['online_stop']:
            self.stops.append(dict(phase=self.phase, iteration=self.iteration,
                                   fresh_state_attempts=self.target.attempts))
            raise StopIteration('R6 frozen prefix-only residual/material-update rule')

    def close(self):
        self.stream.close()


@contextmanager
def _accepted_iteration_hook(function, trace, phase):
    """Patch only the original function's SciPy entry while called serially."""
    namespace = function.__globals__
    original = namespace['least_squares']
    trace.begin(phase)

    def instrumented(*args, **kwargs):
        callback = kwargs.pop('callback', None)

        def combined(intermediate_result):
            if callback is not None:
                callback(intermediate_result)
            trace.observe(intermediate_result)

        return original(*args, callback=combined, **kwargs)

    namespace['least_squares'] = instrumented
    try:
        yield
    finally:
        namespace['least_squares'] = original


def build_online_context(target, y, parent, point, x=None, H=None):
    """One final J is required by the original anchor/head contract and billed."""
    J = target.jacobian(point)
    context = dict(observation=np.asarray(y).copy(), q0=point.q.copy(),
                   chi0=point.state.chi.copy(), J=J.copy(),
                   prediction=point.state.data.copy(), residual=np.asarray(y)-point.state.data,
                   phi=np.array(target.phi(point, y)), parent=np.array(str(parent)))
    if x is not None:
        context['x'] = np.asarray(x).copy()
    if H is not None:
        context['H'] = np.asarray(H).copy()
    return context


def run(target, y, parent, destination, config='native96_probit64', decision_predictor=None):
    """Run exactly one paid anchor path on a caller-owned R6 Target.

    Returns (context, report). Persistence is exclusively under destination.
    All physics/attempt costs must be read from the caller's singleton ledger;
    the local counters are path diagnostics, never a substitute resource ledger.
    """
    from shared.physics import reconstruct, probit_optimize
    config = CONFIGS[config] if isinstance(config, str) else config
    if not isinstance(config, AnchorConfig):
        raise TypeError('AnchorConfig or frozen config name required')
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    if (destination/'ANCHOR_REPORT.json').exists() or (destination/'path/PATH.jsonl').exists():
        raise FileExistsError('existing Track C path: no implicit physical replay')
    y = np.asarray(y)
    if y.shape != (32, 8) or not np.iscomplexobj(y) or not np.isfinite(y).all():
        raise ValueError('original finite complex32x8 observation required')
    observed = _ObservedTarget(target, destination)
    trace = _Trace(observed, y, parent, destination/'path', config, decision_predictor)
    tick = time.perf_counter()
    ledger_start = target.book.totals()
    try:
        with _accepted_iteration_hook(reconstruct, trace, 'native'):
            native, native_report = reconstruct(observed, y, str(parent), method='full_TRF',
                                                 max_evaluations=config.native_cap)
        # Do not retain an extra old LU once the probit optimizer has started.
        native_chi, native_q = native.state.chi.copy(), native.q.copy()
        with _accepted_iteration_hook(probit_optimize, trace, 'probit'):
            point, x, H, probit_report = probit_optimize(observed, y, str(parent),
                destination/'original_optimizer/probit', native, cap=config.probit_cap)
        context = build_online_context(observed, y, parent, point, x, H)
        context.update(native_chi=native_chi, native_q=native_q)
        np.savez_compressed(destination/'CONTEXT.npz', **context)
        report = dict(parent=str(parent), config=asdict(config), GT_read=False,
                      exposed_R5B_GT_used=False, original_optimizer_math=True,
                      original_native_report=native_report, original_probit_report=probit_report,
                      attempted_forward_states=observed.attempts,
                      jacobian_invocations_including_hits=observed.jacobian_invocations,
                      resource_cost_source='caller singleton shared.physics.Book ledger',
                      wall_seconds=time.perf_counter()-tick, adaptive_stops=trace.stops,
                      stopping_rule_label='PREFIX_RESIDUAL_AND_MATERIAL_UPDATE_PROXY_NOT_RISK_CALIBRATED',
                      current_decision_uses_no_new_physics=True, final_context_J_required=True,
                      path_complete=True, output_context=str(destination/'CONTEXT.npz'))
        ledger = target.book.totals()
        report['path_resource_cost'] = {key: ledger[key]-ledger_start[key]
                                       for key in ('RHS', 'forward_RHS', 'derivative_RHS',
                                                   'states', 'LU', 'physics_wall_seconds')}
        _write(destination/'ANCHOR_REPORT.json', report)
        return context, report
    except BaseException as exc:
        _write(destination/'FAILURE.json', dict(parent=str(parent), config=asdict(config),
                    error=repr(exc), attempted_forward_states=observed.attempts, GT_read=False,
                    path_complete=False, no_implicit_retry=True))
        raise
    finally:
        trace.close()
