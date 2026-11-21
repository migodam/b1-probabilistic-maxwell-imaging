"""Replay every saved high-dimensional endpoint using only the public package."""

import argparse
from collections import Counter
from contextlib import contextmanager
import json
import os
from pathlib import Path
import time

for _variable in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_variable] = '1'

import numpy as np
from scipy.linalg import lu_factor, lu_solve

HERE = Path(__file__).resolve().parent
SOURCE = HERE / 'results'
OPERATORS = HERE / 'public_vie_inputs.npz'
OBSERVATION_RADII = HERE / 'observation_radii.json'


class ScalarTM:
    """Discrete VIE on stored binary64 operators, with complete analytic Jacobian."""

    def __init__(self, public):
        self.G, self.E0, self.S = (public[name] for name in ('G', 'E0', 'S'))
        self.n, self.ns = self.E0.shape
        self.nr = self.S.shape[0]
        self.identity = np.eye(self.n, dtype=complex)
        self.factor = self.state = None
        self.counts = Counter({'cached_operator_imports': 1})
        self.timings, self.cpu_timings = Counter(), Counter()

    @contextmanager
    def operation(self, name, rhs=0, **counts):
        self.counts.update({name: 1, 'actual_complex_RHS': int(rhs),
                            **{key: int(value) for key, value in counts.items()}})
        wall, cpu = time.perf_counter(), time.process_time()
        try:
            yield
        finally:
            self.timings[name] += time.perf_counter() - wall
            self.cpu_timings[name] += time.process_time() - cpu

    @staticmethod
    def pack(y):
        port_order = np.asarray(y).T
        return np.stack((port_order.real, port_order.imag), axis=-1).ravel()

    def forward(self, m, label='replay'):
        with self.operation('FULL_FORWARD', rhs=self.ns, LU=1, source_RHS=self.ns,
                            full_factorization_dimension=self.n):
            A = self.identity - m[:, None] * self.G
            b = m[:, None] * self.E0
            self.factor = lu_factor(A, check_finite=False)
            j = lu_solve(self.factor, b, check_finite=False)
            q = self.E0 + self.G @ j
            self.state = {'m': m, 'j': j, 'Q': q, 'y': self.S @ j, 'label': label,
                          'state_relative_residual': float(np.linalg.norm(A @ j - b) /
                                                           max(np.linalg.norm(b), 1e-300))}
        return self.state

    def jacobian(self, state):
        if self.state is not state:
            raise ValueError('Jacobian requires the matching retained VIE system')
        with self.operation('FULL_ANALYTIC_J', rhs=self.nr, receiver_transpose_RHS=self.nr,
                            FullJ=1, material_J_builds=1):
            t = lu_solve(self.factor, self.S.T, trans=1, check_finite=False).T
            j = np.einsum('ri,is->sri', t, state['Q']).reshape(self.ns * self.nr, self.n)
            return np.stack((j.real, j.imag), axis=1).reshape(2 * self.ns * self.nr, self.n)


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def load_operators(case, path=OPERATORS):
    with np.load(path, allow_pickle=False) as data:
        green_key = str(data[f'{case}__G_key'].item())
        return {'G': data[green_key], 'E0': data[f'{case}__E0'], 'S': data[f'{case}__S']}


def residual(y, reference):
    return float(np.linalg.norm(y - reference))


def interval_status(row, delta):
    lo, hi = row
    return 'inside' if hi <= delta else 'outside' if lo > delta else 'uncertain'


def archive_radii(case, document, original_noise):
    result = {'frozen_observation_ball': document['archive_radius']}
    original = original_noise[case].get('original_r8_noise_absolute_radius_reconstructed_offline')
    if original is not None:
        result['original_R8_noise_ball'] = original
    return result


def summarize_fields(y_minus, y_plus, target, center, radii):
    radii = {name: dict(radius=float(value), both_inside=bool(max(residual(y_minus, target),
                                                                 residual(y_plus, target)) <= value))
             for name, value in radii.items()}
    return {
        'field_separation': residual(y_plus, y_minus),
        'archived_residual_norm_minus_plus': [residual(y_minus, target), residual(y_plus, target)],
        'relative_archived_residual_minus_plus': [residual(y_minus, target) / np.linalg.norm(target),
                                                 residual(y_plus, target) / np.linalg.norm(target)],
        'center_residual_norm_minus_plus': [residual(y_minus, center), residual(y_plus, center)],
        'balls': radii,
    }


class ReplayCheck:
    """Track the difference between fresh forwards and the sealed stored fields."""

    def __init__(self, model):
        self.model = model
        self.maximum = 0.

    def evaluate(self, material, stored):
        y = self.model.forward(material)['y']
        relative = residual(y, stored) / max(np.linalg.norm(stored), 1e-300)
        self.maximum = max(relative, self.maximum)
        if relative > 2e-10:
            raise ArithmeticError(f'Stored full-wave prediction fails replay: {relative:.3g}')
        return y


def replay_case(case, checker, source, original_noise):
    doc = read_json(source / f'{case}.json')
    radii = archive_radii(case, doc, original_noise)
    rows = []
    with np.load(source / f'{case}.npz', allow_pickle=False) as data:
        target, task = data['Y_archive'], data['task_vector']
        centers = {}
        for name in doc['centers']:
            key = f'history__{name}'
            centers[name] = checker.evaluate(data[f'{key}__m'], data[f'{key}__Y'])
        for row in doc['pairs']:
            key = row['array_prefix']
            minus, plus = (data[f'{key}__m_{end}'] for end in ('minus', 'plus'))
            fields = [checker.evaluate(m, data[f'{key}__Y_{end}'])
                      for m, end in zip((minus, plus), ('minus', 'plus'), strict=True)]
            gap = float(abs(task @ (plus - minus)))
            if not np.isclose(gap, row['task_gap'], rtol=1e-10, atol=1e-12):
                raise ArithmeticError('Stored task gap fails replay')
            fresh = summarize_fields(*fields, target, centers[row['center']], radii)
            expected = row['archive_ball']['both_inside']
            if 'uncertain' not in (row['archive_ball']['plus'], row['archive_ball']['minus']):
                if fresh['balls']['frozen_observation_ball']['both_inside'] != expected:
                    raise ArithmeticError('Archived feasibility fails replay')
            rows.append({'candidate': key, 'task_gap': gap, **fresh})
    return rows


def replay_comparison(case, checker, source, comparison):
    if case not in comparison['result']:
        return None
    doc = comparison['result'][case]
    methods = {}
    with np.load(source / doc['arrays'], allow_pickle=False) as data:
        y, task = data['Y_archive'], data['task_vector']
        for method, arm in doc['methods'].items():
            candidates = []
            for item in arm['candidates']:
                key = item['array_prefix']
                m = data[f'{key}__material']
                field = checker.evaluate(m, data[f'{key}__prediction'])
                candidates.append({'candidate': key, 'task': float(task @ m), 'archived_residual': residual(field, y),
                                   'relative_archived_residual': residual(field, y) / np.linalg.norm(y)})
            spans = {}
            for ball, contract in doc['spans'].items():
                delta = contract['radius']
                inside = [row['task'] for row in candidates if row['archived_residual'] <= delta]
                width = max(inside) - min(inside) if inside else None
                expected = contract['spans_by_method'][method]['feasible_span']
                ambiguous = any(interval_status(item['archived_distance_interval'], delta) == 'uncertain'
                                for item in arm['candidates'])
                if not ambiguous and ((width is None) != (expected is None) or
                                      (width is not None and not np.isclose(width, expected, atol=1e-12))):
                    raise ArithmeticError('Archived comparison span fails replay')
                spans[ball] = {'radius': delta, 'endpoints_inside': len(inside), 'task_span': width}
            methods[method] = {'candidates': candidates, 'balls': spans}
    return methods


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, help='write replay figures as JSON to this path')
    args = parser.parse_args()
    seal = read_json(SOURCE / 'ENDPOINT_SEAL.json')
    comparison = read_json(SOURCE / 'GN_BUDGET_COMPARISON.json')
    original_noise = read_json(SOURCE / 'OFFLINE_SCORE.json')
    result = {}
    for case in seal['cases']:
        checker = ReplayCheck(ScalarTM(load_operators(case)))
        pairs = replay_case(case, checker, SOURCE, original_noise)
        matched_comparison = replay_comparison(case, checker, SOURCE, comparison)
        result[case] = {'pairs': pairs, 'matched_comparison': matched_comparison,
                        'max_relative_discrepancy_from_sealed_predictions': checker.maximum}
    body = json.dumps(result, indent=2, allow_nan=False) + '\n'
    if args.output:
        args.output.write_text(body, encoding='utf-8')
        print(f'Wrote {sum(len(case["pairs"]) for case in result.values())} fresh endpoint-pair replays and '
              'the corresponding budget comparisons.')
    else:
        print(body, end='')


if __name__ == '__main__':
    main()
