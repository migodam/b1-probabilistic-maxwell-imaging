"""Find rigorous material-diameter lower bounds using only the stored old reads."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations
import json
from time import perf_counter, process_time

import numpy as np
from flint import acb_mat, arb
from scipy.optimize import minimize
from scipy.stats import qmc

from data14 import ACTIONS, FrequencyTM, Fit, ROOT, load_cases, ri
from pivot_bounds import ball, lower_float, mat, row, sqnorm, upper_float

ARMS = ('repeat', 'one_minus_0.02', 'one_plus_0.02', 'split_0.02', 'split_0.08')
CASES = ('spectral-00', 'spectral-01', 'spectral-12', 'spectral-19', 'fresh6-00', 'fresh6-06',
         'fresh6-09', 'confirm6-00', 'confirm6-04', 'confirm6-05')
MARGIN = .997


def complex_pair(v):
    return [float(v.real), float(v.imag)]


def rms(v):
    return float(np.linalg.norm(v) / np.sqrt(v.size))


class PairSearch:
    def __init__(self, case):
        self.case = case
        self.model = FrequencyTM(case.geometry)
        self.eps, self.gamma = case.additive_radius, case.gain_radius
        self.lo, self.hi = case.box
        self.hard = case.hard_gain_radius or self.gamma
        self.verified_solves = 0
        self.old_cache = {}

    def anchor(self, material):
        f, _ = self.model.response(material)
        z = np.vdot(f, self.case.old - f) / (np.vdot(f, f).real + (self.eps / self.gamma)**2)
        if abs(z) > self.hard:
            z *= (.999 * self.hard / abs(z))
        x = np.r_[material, z.real / self.gamma, z.imag / self.gamma]
        score = np.linalg.norm(((1 + z) * f - self.case.old) / self.eps)**2 + abs(z / self.gamma)**2
        return x, float(score)

    def anchors(self, choices):
        points = []
        scores = []

        def add(material):
            x, score = self.anchor(np.clip(material, self.lo, self.hi))
            scores.append(score)
            if score < MARGIN and all(rms(x[:6] - p[:6]) > 1e-6 for p in points):
                points.append(x)

        for start in choices['old_starts']:
            add(np.asarray(start))
        if not points or len(choices['exact_central_roots']) < 2:
            fit = Fit(self.model, [0.], [self.case.old], gain_calibration=((self.eps / self.gamma)**2, 1.))
            samples = qmc.Sobol(6, seed=646427, scramble=True).random_base2(2)
            starts = [np.full(6, 12.), *[np.asarray(x) for x in choices['old_starts']],
                      *list(self.lo + (self.hi - self.lo) * samples)]
            for start in starts:
                add(fit.solve(start, 140)[0])
        return points, scores

    def response(self, p, delta):
        m = p[:6]
        z = self.gamma * complex(*p[6:])
        f, j = self.model.response(m, delta)
        F = (1 + z) * f
        J = np.column_stack(((1 + z) * j, self.gamma * f, 1j * self.gamma * f))
        return F, J

    def fun_parts(self, p, action, r):
        f, j = self.response(p, 0.)
        old = np.r_[ri((f - self.case.old) / self.eps), p[6:]]
        old_j = np.vstack((ri(j / self.eps), np.eye(8)[6:]))
        blocks = []
        derivatives = []
        for delta, energy in zip(action.deltas, action.energies):
            y, Jy = (f, j) if delta == 0 else self.response(p, delta)
            weight = np.sqrt(energy) / r
            blocks.extend(ri(y * weight))
            derivatives.extend(ri(Jy * weight))
        return old, old_j, np.array(blocks), np.array(derivatives)

    def constraints(self, a, b, action, r):
        oa, ja, fa, Ja = self.fun_parts(a, action, r)
        ob, jb, fb, Jb = self.fun_parts(b, action, r)
        half_diff = (fa - fb) / 2
        cap = (self.hard / self.gamma)**2
        vals = np.asarray([MARGIN - oa @ oa, MARGIN - ob @ ob,
                           MARGIN - half_diff @ half_diff,
                           cap * MARGIN - a[6:] @ a[6:], cap * MARGIN - b[6:] @ b[6:]])
        jac = np.zeros((len(vals), 16))
        jac[0, :8] = -2 * (oa @ ja)
        jac[1, 8:] = -2 * (ob @ jb)
        jac[2, :8] = -half_diff @ Ja
        jac[2, 8:] = half_diff @ Jb
        jac[3, 6:8] = -2 * a[6:]
        jac[4, 14:] = -2 * b[6:]
        return vals, jac

    def inside(self, a, b, action, r):
        if min(np.min(a[:6]), np.min(b[:6])) < self.lo or max(np.max(a[:6]), np.max(b[:6])) > self.hi:
            return False
        return bool(np.min(self.constraints(a, b, action, r)[0]) >= -.0001)

    def transform(self, x, action, r):
        _, Dj, _, J = self.fun_parts(x, action, r)
        _, s, vh = np.linalg.svd(np.vstack((Dj, J)), full_matrices=False)
        return vh.T / np.maximum(s, .002), float(s[-1])

    def line(self, center, direction, action, r):
        max_step = 2.
        directions = (-direction.copy(), direction.copy())
        for vv in directions:
            for j, (m, d) in enumerate(zip(center[:6], vv[:6])):
                if d == 0:
                    continue
                distance = self.hi - m if d > 0 else m - self.lo
                if distance < 1e-8:
                    vv[j] = 0.
                else:
                    max_step = min(max_step, distance / abs(d))
        if max_step <= 0.:
            return center, center
        a, b = [center + v * max_step for v in directions]
        if self.inside(a, b, action, r):
            return a, b
        lo, hi = 0., max_step
        for _ in range(17):
            step = (lo + hi) / 2
            a, b = [center + v * step for v in directions]
            if self.inside(a, b, action, r):
                lo = step
            else:
                hi = step
        return tuple(center + v * lo for v in directions)

    def optimize(self, a, b, action, r, maxiter):
        tic, calls0 = process_time(), self.model.factorizations
        p0 = np.r_[a, b]
        T = np.zeros((16, 16))
        T[:8, :8], weak1 = self.transform(a, action, r)
        T[8:, 8:], weak2 = self.transform(b, action, r)
        origin = None
        con = None
        Jcon = None
        checked = []
        scale = max(.03, rms(a[:6] - b[:6]))**2
        I = np.r_[np.arange(6), 8+np.arange(6)]
        M = T[I, :]

        def objective(y):
            p = p0 + T @ y
            diff = p[:6] - p[8:14]
            g = np.r_[diff, 0., 0., -diff, 0., 0.]
            return -np.dot(diff, diff) / (6 * scale), -(T.T @ g) / (3 * scale)

        def constraint(y):
            nonlocal origin, con, Jcon
            if origin is None or not np.array_equal(y, origin):
                p = p0 + T @ y
                box = p[I]
                q, B = self.constraints(p[:8], p[8:], action, r)
                con = np.r_[q, box - self.lo, self.hi - box]
                Jcon = np.vstack((B @ T, M, -M))
                origin = y.copy()
                if np.min(con) >= -1e-8:
                    checked.append((rms(p[:6] - p[8:14]), p[:8].copy(), p[8:].copy()))
            return con

        res = minimize(objective, np.zeros(16), jac=True, method='SLSQP',
                       constraints=[dict(type='ineq', fun=constraint, jac=lambda y: (constraint(y), Jcon)[1])],
                       options=dict(maxiter=maxiter, ftol=1e-10))
        constraint(res.x)
        return checked, dict(status=int(res.status), success=bool(res.success), iterations=int(res.nit),
                            full_forward_with_derivative_calls=self.model.factorizations - calls0,
                            cpu_s=process_time() - tic, minimum_joint_singular_values=[weak1, weak2],
                            final_min_constraint=float(np.min(con)),
                            best_float_RMS=max((x[0] for x in checked), default=None))

    def verify(self, a, b, action, radius):
        """Use a fixed explicit gain for each material, and prove compatibility jointly."""
        K, S = len(a) - 2, len(self.case.old)
        zeta = [self.gamma * complex(*p[6:]) for p in (a, b)]
        if any(np.min(p[:K]) < self.lo or np.max(p[:K]) > self.hi for p in (a, b)):
            return None
        eps, gamma, r = [arb(float(q)) for q in (self.eps, self.gamma, radius)]
        hard = arb(float(self.hard))
        old_scores = []
        calibrated = []
        initial = self.verified_solves
        for p, z in zip((a, b), zeta):
            material = [float(x) for x in p[:K]]
            if abs(ball(z)).upper() > hard:
                return None
            old_key = (tuple(material), z)
            cache = self.old_cache.get(old_key)
            if cache is None:
                cache = {}
                self.old_cache[old_key] = cache
            y = {}
            for delta in dict.fromkeys((0., *action.deltas)):
                if delta not in cache:
                    H, B, R, *_ = self.model.kernels(float(delta))
                    T = acb_mat([[(1 / ball(material[i]) if i == j else 0) - ball(H[i, j])
                                  for j in range(K)] for i in range(K)])
                    try:
                        cache[delta] = ((row(np.ravel(R)) * T.solve(mat(B))) * (1 + ball(z))).entries()
                    except (ZeroDivisionError, ValueError):
                        return None
                    self.verified_solves += 1
                y[delta] = cache[delta]
            residual = [ball(obs) - f for obs, f in zip(self.case.old, y[0.])]
            score = sqnorm(residual) / eps**2 + sqnorm([ball(z)]) / gamma**2
            if not score < 1:
                return None
            old_scores.append(upper_float(score))
            calibrated.append(y)
        diff_square = sum((arb(float(E)) * sqnorm([x - y for x, y in zip(calibrated[0][d], calibrated[1][d])])
                           for d, E in zip(action.deltas, action.energies)), arb(0))
        if not diff_square < 4 * r**2:
            return None
        dist = (sum(((arb(float(p)) - arb(float(q)))**2 for p, q in zip(a[:K], b[:K])), arb(0)) / K).sqrt()
        common = [[complex_pair(.5 * np.sqrt(E) * (complex(float(x.real.mid()), float(x.imag.mid()))
                                                  + complex(float(y.real.mid()), float(y.imag.mid()))))
                   for x, y in zip(calibrated[0][d], calibrated[1][d])]
                  for d, E in zip(action.deltas, action.energies)]
        return dict(materials=[a[:K].tolist(), b[:K].tolist()], gains_zeta_ri=[complex_pair(z) for z in zeta],
                    rms_lower=max(0., lower_float(dist)), squared_old_contract_upper=old_scores,
                    new_mean_gap_upper=upper_float(diff_square.sqrt()), new_overlap_ratio_upper=upper_float(diff_square.sqrt() / (2*r)),
                    common_new_read_weighted_center_approx_ri=common,
                    strict_vie_factors=self.verified_solves-initial,
                    strict_vie_rhs=S * (self.verified_solves-initial))

    def action_search(self, points, action, maxiter):
        start, wall, calls = process_time(), perf_counter(), self.model.factorizations
        vstart = self.verified_solves
        radius = self.case.new_noise_radius(action)
        candidates = []
        paired = []
        invalid = []
        for a, b in combinations(points, 2):
            if self.inside(a, b, action, radius):
                candidates.append((rms(a[:6] - b[:6]), a, b))
                paired.append((a, b))
            else:
                invalid.append((rms(a[:6] - b[:6]), a, b))
        lines = []
        for x in points:
            T, _ = self.transform(x, action, radius)
            for j in (-1, -2):
                a, b = self.line(x, T[:, j], action, radius)
                d = rms(a[:6] - b[:6])
                if d > 0:
                    candidates.append((d, a, b))
                    lines.append((d, a, b))
        starts = sorted(lines, key=lambda x: -x[0])[:2]
        if paired:
            a, b = max(paired, key=lambda p: rms(p[0][:6] - p[1][:6]))
            starts.append((rms(a[:6] - b[:6]), a, b))
        if invalid:
            starts.append(max(invalid, key=lambda v: v[0]))
        details = []
        for _, a, b in starts:
            accepted, rec = self.optimize(a, b, action, radius, maxiter)
            candidates.extend(accepted)
            details.append(rec)
        result = None
        chosen = None
        verifications = 0
        for d, a, b in sorted(candidates, key=lambda p: -p[0]):
            verifications += 1
            result = self.verify(a, b, action, radius)
            if result is not None and result['rms_lower'] > 0:
                chosen = (a, b)
                break
            result = None
        return dict(action=action.name, band_offsets=action.deltas, energies=action.energies,
                    new_noise_radius=radius, witness=result, opt_runs=details, attempted_distinct_anchor_pairs=len(paired)+len(invalid),
                    initially_feasible_anchor_pairs=len(paired), float_line_candidates=len(lines),
                    strict_candidate_trials=verifications, strict_full_forward_factors=self.verified_solves-vstart,
                    strict_full_forward_rhs=len(self.case.old)*(self.verified_solves-vstart),
                    full_forward_with_derivative_calls=self.model.factorizations-calls,
                    cpu_s=process_time()-start, wall_s=perf_counter()-wall), chosen


def run(job):
    case, choices, maxiter = job
    start, wall = process_time(), perf_counter()
    search = PairSearch(case)
    points, scores = search.anchors(choices)
    work = dict(input_old_starts=len(choices['old_starts']), central_roots=len(choices['exact_central_roots']),
                eligible_anchors=len(points), attempted_anchor_scores=scores,
                full_forward_with_derivative_calls=search.model.factorizations,
                cpu_s=process_time()-start, measure_choices_precomputed_costs=choices['costs'])
    explored = {name: search.action_search(points, ACTIONS[name], maxiter) for name in ARMS}
    out = {name: result for name, (result, _) in explored.items()}
    candidates = [(name, result['witness']['rms_lower'], pair)
                  for name, (result, pair) in explored.items() if pair is not None]
    for name, record in out.items():
        start_cross, solves = process_time(), search.verified_solves
        baseline = record['witness']['rms_lower'] if record['witness'] is not None else 0.
        attempted = []
        for source, bound, pair in sorted(candidates, key=lambda x: -x[1]):
            if source == name or bound <= baseline:
                continue
            candidate = search.verify(*pair, ACTIONS[name], record['new_noise_radius'])
            attempted.append(dict(source=source, passed=candidate is not None))
            if candidate is not None:
                record['witness'] = candidate
                break
        record['cross_action_verification'] = dict(
            tried=attempted, witness_source=(attempted[-1]['source'] if attempted and attempted[-1]['passed'] else name),
            strict_full_forward_factors=search.verified_solves-solves,
            strict_full_forward_rhs=len(case.old)*(search.verified_solves-solves),
            cpu_s=process_time()-start_cross)
    return case.case_id, dict(cohort=case.cohort, anchor_search=work, actions=out,
                    full_forward_with_derivative_calls=search.model.factorizations,
                    strict_full_forward_factors=search.verified_solves,
                    strict_full_forward_rhs=len(case.old)*search.verified_solves,
                    cpu_s=process_time()-start, wall_s=perf_counter()-wall)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case', action='append')
    p.add_argument('--all', action='store_true')
    p.add_argument('--workers', type=int, default=3)
    p.add_argument('--maxiter', type=int, default=65)
    p.add_argument('--tag', default='pivot_diameter_witnesses')
    a = p.parse_args()
    choices_file = ROOT / 'results/measure_choices.json'
    choices = {x['case_id']: x for x in json.loads(choices_file.read_text())['cases']}
    ids = set(a.case or CASES)
    jobs = [(c, choices[c.case_id], a.maxiter) for c in load_cases() if a.all or c.case_id in ids]
    all_cases = {}
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for name, record in pool.map(run, jobs):
            all_cases[name] = record
            print(name, record['anchor_search']['eligible_anchors'], {k: None if v['witness'] is None else
                  round(v['witness']['rms_lower'], 8) for k, v in record['actions'].items()}, flush=True)
    report = dict(input_scope='Old observed rows, known geometry, noise radii, and old-only starts in measure_choices.json.',
                  certificate='Two materials and one explicit gain per material; strict old ellipsoid and optional hard gain disk; weighted new mean gap < 2r.',
                  maxiter_per_slsqp=a.maxiter, actions=ARMS, cases=all_cases)
    (ROOT / 'results' / f'{a.tag}.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')


if __name__ == '__main__':
    main()
