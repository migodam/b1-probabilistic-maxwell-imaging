#!/usr/bin/env python3
"""Fit two saved observations with the fixed, nonnegative ridge-GN protocol."""
from collections import Counter
from pathlib import Path
import sys
import json
from time import perf_counter, process_time

HERE = Path(__file__).resolve().parent
STUDY = HERE.parents[1]
HIGH = STUDY.parent / 'z_alg_secant_tasks' / 'highdim'
sys.path.insert(0, str(HIGH))
from replay import ScalarTM, load_operators  # Sets a single BLAS thread before NumPy.

import numpy as np
from scipy.fft import fft2, ifft2
from scipy.linalg import eigvalsh
from scipy.sparse.linalg import LinearOperator, gmres

RING = 'r10_twod_0001'
NATIVE = 'optfresh_01__pw180_10x20__n02'
CASES = (RING, NATIVE)
APERTURE = 'PW180_10x20'
RTOL = 1.e-10


def read_json(path):
    return json.loads(path.read_text())


def write_json(path, data):
    pending = path.with_suffix('.partial')
    pending.write_text(json.dumps(data, indent=2) + '\n')
    pending.replace(path)


def clocks():
    return perf_counter(), process_time()


def elapsed(start):
    return dict(zip(('wall_seconds', 'cpu_seconds'),
                    (current - old for current, old in zip(clocks(), start, strict=True))))


class FourierTM:
    """Native crop/pad scattering model for both state and transpose solves."""

    def __init__(self, geometry):
        self.S = geometry[f'{APERTURE}__S']
        self.B = geometry[f'{APERTURE}__B']
        self.n, self.ns = self.B.shape
        self.nr = self.S.shape[0]
        self.side = int(np.sqrt(self.n))
        self.spectrum = geometry['spectrum']
        self.pad = self.spectrum.shape[0]
        opposite = -np.arange(self.pad) % self.pad
        self.transpose_spectrum = self.spectrum[np.ix_(opposite, opposite)]
        self.counts = Counter(cached_operator_imports=1)

    def green(self, value, transpose=False):
        self.counts['green_fft_applications'] += 1
        spectrum = self.transpose_spectrum if transpose else self.spectrum
        shape = (self.pad, self.pad)
        padded = fft2(value.reshape(self.side, self.side), s=shape)
        return ifft2(spectrum * padded)[:self.side, :self.side].ravel()

    def solve(self, m, rhs, transpose=False):
        if not np.any(m):
            return rhs.copy()
        if transpose:
            def action(x):
                return x - self.green(m * x, transpose=True)
        else:
            def action(x):
                return x - m * self.green(x)
        A = LinearOperator((self.n, self.n), matvec=action, dtype=complex)
        result = []
        for vector in rhs.T:
            counts = self.counts
            counts['actual_complex_RHS'] += 1
            counts['receiver_transpose_RHS' if transpose else 'source_RHS'] += 1

            def completed_step(_):
                counts['gmres_iterations'] += 1

            sol, info = gmres(A, vector, rtol=RTOL, atol=0., restart=64,
                              maxiter=12, callback=completed_step, callback_type='pr_norm')
            if info:
                raise ArithmeticError(f'Native GMRES did not converge: info={info}.')
            result.append(sol)
        return np.column_stack(result)

    def forward(self, m):
        self.counts['FULL_FORWARD'] += 1
        rhs = m[:, None] * self.B
        current = self.solve(m, rhs)
        total = self.B + np.column_stack([self.green(j) for j in current.T])
        defect = current - m[:, None] * total
        return dict(m=m, j=current, Q=total, y=self.S @ current,
                    state_relative_residual=float(np.linalg.norm(defect) /
                                                  max(np.linalg.norm(rhs), 1.e-300)))

    def jacobian(self, state):
        self.counts['FULL_ANALYTIC_J'] += 1
        receiver = self.solve(state['m'], self.S.T, transpose=True).T
        jac = np.einsum('ri,is->sri', receiver, state['Q']).reshape(-1, self.n)
        return np.stack((jac.real, jac.imag), axis=1).reshape(-1, self.n)


def load_case(case):
    if case == RING:
        with np.load(HIGH / 'results' / (case + '.npz'), allow_pickle=False) as raw:
            measurement = raw['Y_archive']
        return ScalarTM(load_operators(case)), measurement
    with np.load(STUDY / 'synthetic/public/fine_geometry.npz', allow_pickle=False) as raw:
        model = FourierTM(raw)
    with np.load(STUDY / 'synthetic/public' / (case + '.npz'), allow_pickle=False) as raw:
        return model, raw['Y']


def objective(state, y, strength):
    return .5 * (np.linalg.norm(state['y'] - y)**2 + strength * np.dot(state['m'], state['m']))


def iterate(model, y):
    m = np.zeros(model.n)
    state = model.forward(m)
    jac = model.jacobian(state)
    sigma2 = eigvalsh(jac @ jac.T, subset_by_index=(jac.shape[0]-1, jac.shape[0]-1))[0]
    strength = 0.1 * sigma2
    packed = ScalarTM.pack(y)
    count_updates = 0
    accepted = [0]
    trials = []

    def record(state, step_fraction, outer, keep):
        trials.append(dict(outer_iteration=outer, fraction=step_fraction, accepted=keep,
                           residual_frobenius=float(np.linalg.norm(state['y'] - y)),
                           peak=float(np.max(state['m'])), objective=float(objective(state, y, strength)),
                           state_relative_residual=state['state_relative_residual']))
        snapshots.append((state['m'], state['y'], state['j']))

    snapshots = []
    record(state, 0., 0, True)
    stop = '12 outer updates'
    for outer in range(1, 13):
        old_objective = objective(state, y, strength)
        linear_target = packed - ScalarTM.pack(state['y']) + jac @ m
        x = z = m.copy()
        momentum = 1.
        for _ in range(128):
            gradient = jac.T @ (jac @ z - linear_target) + strength * z
            x_next = np.maximum(z - gradient / (sigma2 + strength), 0.)
            momentum_next = (1. + np.sqrt(1. + 4. * momentum**2)) / 2.
            z = x_next + ((momentum - 1.) / momentum_next) * (x_next - x)
            x, momentum = x_next, momentum_next
        for backtrack in range(10):
            fraction = 2.**-backtrack
            trial_m = m + fraction * (x - m)
            trial = model.forward(trial_m)
            new_objective = objective(trial, y, strength)
            keep = bool(new_objective < old_objective)
            record(trial, fraction, outer, keep)
            if keep:
                accepted.append(len(trials)-1)
                break
        else:
            stop = 'full-wave regularized objective failed to decrease on ten trials'
            break
        count_updates += 1
        m, state = trial_m, trial
        if old_objective - new_objective < 1.e-6 * old_objective:
            stop = 'relative full-wave regularized objective decrease < 1e-6'
            break
        if outer < 12:
            jac = model.jacobian(state)
            sigma2 = eigvalsh(jac @ jac.T, subset_by_index=(jac.shape[0]-1, jac.shape[0]-1))[0]
    row = dict(lambda_ridge=float(strength), stopping_reason=stop,
               accepted_updates=count_updates, accepted_indices=accepted, trials=trials,
               operations=dict(model.counts), projected_gradient_iterations=128*(outer))
    return row, snapshots


def main():
    for case in CASES:
        start = clocks()
        model, y = load_case(case)
        receipt, states = iterate(model, y)
        materials, predictions, currents = map(np.stack, zip(*states, strict=True))
        np.savez_compressed(HERE / (case + '_path.npz'), Y=y, materials=materials,
                            predictions=predictions, currents=currents)
        receipt.update(case=case, state_cells=model.n, source_receiver_counts=[model.ns, model.nr],
                       load_solve_store_seconds=elapsed(start))
        write_json(HERE / (case + '_online.json'), receipt)
        endpoint = receipt['trials'][receipt['accepted_indices'][-1]]
        print(case, receipt['stopping_reason'], endpoint, elapsed(start), flush=True)


if __name__ == '__main__':
    main()
