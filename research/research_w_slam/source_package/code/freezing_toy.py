#!/usr/bin/env python3
"""Algebra checks for W-SLAM freezing certificates; no Maxwell PDE is solved.

Reproduce with:
    python3 deliverables/code/freezing_toy.py

Output: deliverables/results/freezing_toy.json
Only NumPy and the Python standard library are required.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np


SEED = 20261025
TRIALS = 100


def relative_error(a, b):
    a, b = np.asarray(a), np.asarray(b)
    return float(np.linalg.norm(a - b) / max(1.0, np.linalg.norm(a), np.linalg.norm(b)))


def main():
    rng = np.random.default_rng(SEED)
    maxerr = dict(quadratic_gap=0.0, task_regret=0.0, task_covariance=0.0,
                  task_mean=0.0, dual_exact=0.0, resolvent_identity=0.0,
                  port_identity=0.0)
    intervals_pass = 0
    remainder_bounds_pass = 0
    port_bounds_pass = 0
    largest_condition = 0.0
    interval_example = None

    for trial in range(TRIALS):
        n_f, n_u, n_z = 4, 7, 3
        n = n_f + n_u
        raw = rng.normal(size=(n, n))
        h = raw.T @ raw + 0.3 * np.eye(n)
        g = rng.normal(size=n)
        c = rng.normal(size=(n_z, n))
        mroot = rng.normal(size=(n_z, n_z))
        m = mroot.T @ mroot + 0.1 * np.eye(n_z)
        # Symmetric square root avoids transpose conventions in task whitening.
        evals, evecs = np.linalg.eigh(m)
        mhalf = (evecs * np.sqrt(evals)) @ evecs.T
        hff, hfu = h[:n_f, :n_f], h[:n_f, n_f:]
        huf, huu = h[n_f:, :n_f], h[n_f:, n_f:]
        cf, cu = c[:, :n_f], c[:, n_f:]
        gf, gu = g[:n_f], g[n_f:]
        inv_huu_huf = np.linalg.solve(huu, huf)
        delta0_u = np.linalg.solve(huu, gu)
        s = hff - hfu @ inv_huu_huf
        r = gf - hfu @ delta0_u
        t = cf - cu @ inv_huu_huf
        delta = np.linalg.solve(h, g)
        delta0 = np.r_[np.zeros(n_f), delta0_u]
        q = lambda d: 0.5 * d @ h @ d - g @ d

        qgap_direct = q(delta0) - q(delta)
        qgap_schur = 0.5 * r @ np.linalg.solve(s, r)
        maxerr['quadratic_gap'] = max(maxerr['quadratic_gap'], relative_error(qgap_direct, qgap_schur))

        action_diff = c @ (delta - delta0)
        formula_diff = t @ np.linalg.solve(s, r)
        maxerr['task_mean'] = max(maxerr['task_mean'], relative_error(action_diff, formula_diff))
        regret_direct = action_diff @ m @ action_diff
        regret_formula = np.linalg.norm(mhalf @ formula_diff) ** 2
        maxerr['task_regret'] = max(maxerr['task_regret'], relative_error(regret_direct, regret_formula))

        cov_direct = c @ np.linalg.solve(h, c.T)
        cov_formula = cu @ np.linalg.solve(huu, cu.T) + t @ np.linalg.solve(s, t.T)
        maxerr['task_covariance'] = max(maxerr['task_covariance'], relative_error(cov_direct, cov_formula))

        d = mhalf @ t
        dual = np.linalg.solve(s, d.T)
        error_true = mhalf @ action_diff
        maxerr['dual_exact'] = max(maxerr['dual_exact'], relative_error(dual.T @ r, error_true))
        dual_approx = dual + 0.01 * rng.normal(size=dual.shape)
        xi = d.T - s @ dual_approx
        alpha = 0.999 * float(np.linalg.eigvalsh(s)[0])  # Certified here by exact small dense eigenvalue.
        beta = np.linalg.norm(xi, axis=0) * np.linalg.norm(r) / alpha
        radius = np.linalg.norm(beta)
        ehat = dual_approx.T @ r
        lower = max(0.0, np.linalg.norm(ehat) - radius) ** 2
        upper = (np.linalg.norm(ehat) + radius) ** 2
        valid = lower - 1e-10 <= regret_direct <= upper + 1e-10
        intervals_pass += int(valid)
        assert valid, (lower, regret_direct, upper)
        if interval_example is None:
            interval_example = dict(lower=float(lower), truth=float(regret_direct), upper=float(upper))
        largest_condition = max(largest_condition, float(np.linalg.cond(h)))

        # General finite-dimensional resolvent algebra, not a physical PDE model.
        n_field = 6
        matrix_t = rng.normal(size=(n_field, n_field))
        eta = rng.uniform(0.15, 0.7)
        matrix_t *= eta / np.linalg.norm(matrix_t, 2)
        e0 = rng.normal(size=n_field)
        obs = rng.normal(size=(3, n_field))
        e1 = np.linalg.solve(np.eye(n_field) + matrix_t, e0)
        remainder = obs @ (e1 - e0 + matrix_t @ e0)
        identity = obs @ matrix_t @ matrix_t @ np.linalg.solve(np.eye(n_field) + matrix_t, e0)
        maxerr['resolvent_identity'] = max(maxerr['resolvent_identity'], relative_error(remainder, identity))
        bound = np.linalg.norm(obs, 2) * eta**2 / (1 - eta) * np.linalg.norm(e0)
        valid = np.linalg.norm(remainder) <= bound + 1e-10
        remainder_bounds_pass += int(valid)
        assert valid

        # Four-port multiple scattering identity and uniform difference bound.
        p, receivers = 4, 3
        target = rng.normal(size=(p, p)) * 0.2
        ga = rng.normal(size=(p, p)) * 0.08
        gb = ga + rng.normal(size=(p, p)) * 0.01
        a1 = rng.normal(size=p)
        a2 = a1 + rng.normal(size=p) * 0.01
        rr1 = rng.normal(size=(receivers, p))
        rr2 = rr1 + rng.normal(size=(receivers, p)) * 0.01
        direct1 = rng.normal(size=receivers)
        direct2 = direct1 + rng.normal(size=receivers) * 0.01
        q1 = np.linalg.inv(np.eye(p) - ga @ target)
        q2 = np.linalg.inv(np.eye(p) - gb @ target)
        f1 = direct1 + rr1 @ target @ q1 @ a1
        f2 = direct2 + rr2 @ target @ q2 @ a2
        current1 = np.linalg.solve(np.eye(p) - target @ ga, target @ a1)
        maxerr['port_identity'] = max(maxerr['port_identity'], relative_error(f1, direct1 + rr1 @ current1))
        tau = np.linalg.norm(target, 2)
        kappa = max(np.linalg.norm(q1, 2), np.linalg.norm(q2, 2))
        amax = max(np.linalg.norm(a1), np.linalg.norm(a2))
        rmax = max(np.linalg.norm(rr1, 2), np.linalg.norm(rr2, 2))
        ed, er = np.linalg.norm(direct1 - direct2), np.linalg.norm(rr1 - rr2, 2)
        ea, eg = np.linalg.norm(a1 - a2), np.linalg.norm(ga - gb, 2)
        port_bound = ed + tau * kappa * (amax * er + rmax * ea) + rmax * amax * tau**2 * kappa**2 * eg
        valid = np.linalg.norm(f1 - f2) <= port_bound + 1e-10
        port_bounds_pass += int(valid)
        assert valid

    # A strictly improved objective can have worse task action regret.
    h = np.eye(2)
    mu = np.array([1.0, -1.0])
    g = mu.copy()
    c = np.array([[1.0, 1.0]])
    before, after = np.zeros(2), np.array([1.0, 0.0])
    q = lambda z: float(0.5 * z @ h @ z - g @ z)
    cancellation = dict(objective_before=q(before), objective_after=q(after),
                        task_regret_before=float(np.linalg.norm(c @ (before - mu))**2),
                        task_regret_after=float(np.linalg.norm(c @ (after - mu))**2))

    # Infinite-noise-average limit of two designs with an incorrectly frozen wall.
    wtrue, xtrue, bias = 0.8, 0.2, 0.1
    design = np.array([[1.0, 1.0], [2.0, 1.0]])
    data = design @ np.array([wtrue, xtrue])
    joint = np.linalg.solve(design, data)
    wfrozen = wtrue + bias
    xfrozen = float(np.mean(data - design[:, 0] * wfrozen))
    bias_case = dict(wall_truth=wtrue, target_truth=xtrue, wall_frozen=wfrozen,
                     target_frozen=xfrozen, target_bias=xfrozen-xtrue,
                     expected_target_bias=-1.5*bias, squared_bias=(xfrozen-xtrue)**2,
                     joint_solution=joint.tolist(), joint_design_rank=int(np.linalg.matrix_rank(design)))

    marginals_case = dict(posterior_a_support=[[0, 0, 0], [1, 1, 0]],
                          posterior_b_support=[[0, 1, 1], [1, 0, 1]],
                          within_posterior_weights=[0.5, 0.5],
                          voxel_marginals_both=[0.5, 0.5],
                          target_mean_a=0.0, target_mean_b=1.0,
                          covariance_w1_w2_a=0.25, covariance_w1_w2_b=-0.25,
                          support_columns=['wall_voxel_1', 'wall_voxel_2', 'target'])

    assert max(maxerr.values()) < 1e-10, maxerr
    results = dict(scope='Algebra verification only; no Maxwell PDE, inverse reconstruction, or speedup experiment.',
                   seed=SEED, trials=TRIALS, status='PASS', maximum_relative_errors=maxerr,
                   largest_joint_hessian_condition=largest_condition,
                   dual_intervals_valid=intervals_pass, dual_interval_example=interval_example,
                   resolvent_bounds_valid=remainder_bounds_pass,
                   four_port_bounds_valid=port_bounds_pass,
                   nonmonotone_task_counterexample=cancellation,
                   permanent_freezing_bias= bias_case,
                   equal_voxel_marginals_counterexample=marginals_case)
    out = Path(__file__).resolve().parents[1] / 'results' / 'freezing_toy.json'
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(results, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(dict(status=results['status'], trials=TRIALS,
                          maximum_relative_error=max(maxerr.values()), output=str(out)), indent=2))


if __name__ == '__main__':
    main()
