"""Distinguish scalar from shared constraints and exercise exact chart degeneracies."""
import json
from itertools import product
from time import process_time

import numpy as np
from flint import acb_mat, arb
from scipy.optimize import brentq
from scipy.stats import qmc

from data14 import Fit, FrequencyTM, ROOT, load_cases
from pivot_bounds import OldCase, pick_cols, sqnorm
from pivot_controls import Cofactors, profile_witness
from pivot_experiment import make_case


def diagnostics(public):
    case = make_case(public)
    cofactor = Cofactors(case)
    model = FrequencyTM(public.geometry)
    eps, gamma = public.additive_radius, public.gain_radius
    fit = Fit(model, [0.], [public.old], gain_calibration=((eps / gamma)**2, 1.))

    def score(m):
        return float(np.linalg.norm(fit.evaluate(m)[0] / eps)**2)

    lo, hi = public.box
    starts = [(case.upper + case.lower) / 2, *qmc.Sobol(case.k, scramble=True, seed=414).random_base2(2)
              * (hi - lo) + lo]
    center = min((fit.solve(start, 140)[0] for start in starts), key=score)
    result = dict(center_score=score(center), model_factorizations=model.factorizations, point=None)
    if score(center) >= 1:
        return result
    for target in (1.06, 1.15, 1.3):
        for j, sign in product(range(case.k), (-1, 1)):
            direction = sign * np.eye(case.k)[j]
            distance = min(1., center[j] - lo if sign < 0 else hi - center[j])
            if distance <= 0 or score(center + distance * direction) < target:
                continue
            t = brentq(lambda step: score(center + step * direction) - target, 0., distance)
            m = center + t * direction
            G, r = cofactor.evaluate([arb(float(x)) for x in m])
            scores = [sqnorm([a]) / sqnorm(g) for a, g in zip(r, G)]
            if not all(v < 1 for v in scores) or not cofactor.joint_bound(m, m):
                continue
            certified = None
            for radius in (1e-4, 1e-5, 1e-6, 1e-7, 1e-8):
                a, b = np.maximum(m - radius, lo), np.minimum(m + radius, hi)
                g_box, r_box = cofactor.evaluate([arb(float(x)).union(arb(float(y))) for x, y in zip(a, b)])
                ratios = [sqnorm([c]).upper() / sqnorm(g).lower() for c, g in zip(r_box, g_box)
                          if sqnorm(g).lower() > 0]
                if len(ratios) == len(r_box) and all(v < 1 for v in ratios) and cofactor.joint_bound(a, b):
                    certified = dict(region=[a.tolist(), b.tolist()], nominal_linf_radius=radius,
                                     scalar_score_upper=[float(v.upper()) for v in ratios])
                    break
            result.update(point=dict(material=m.tolist(), axis=j, sign=sign,
                                     scalar_scores=[float(v.upper()) for v in scores],
                                     direct_common_error_score=str(profile_witness(case, m)),
                                     shared_cofactor_exclusion=True, box=certified),
                          model_factorizations=model.factorizations)
            return result
    result['model_factorizations'] = model.factorizations
    return result


def degenerate_controls():
    k, s = 6, 3
    H = np.zeros((k, k))
    B = np.r_[np.eye(s), np.zeros((s, s))]
    R = np.array([1., 0., 2., 3., 4., 5.])
    host = np.zeros(s)
    m = np.ones(k)
    lower, upper = m - .5, m + .5
    dark = OldCase(H, R, B, host, R @ B, .001, .001, lower, upper)
    (a, b), status, count, q = dark.cover(lower, upper)
    assert status == 'dark_open' and a[1] == lower[1] and b[1] == upper[1]
    setup = dict(status=status, dark_axes=count, q=q, dark_axis=1,
                 retained_axis_bounds=[a[1], b[1]], v_axis=0., w_axis=0.)
    H[0, 3] = H[3, 0] = .25
    H[3, 3] = 1.
    D = R @ np.linalg.solve(np.eye(k) - H, B)
    singular = OldCase(H, R, B, host, D, .001, .001, lower, upper)
    chosen = singular.pivots
    singular.pivots = (3, 4, 5)
    singular.tail = (0, 1, 2)
    singular.NHI = pick_cols(singular.NH, singular.pivots)
    singular.HI = pick_cols(singular.H, singular.pivots)
    singular.NI = pick_cols(singular.N, singular.pivots)
    (a, b), status, _, _ = singular.cover(lower, upper)
    det_full = (singular.eye - singular.H).det()
    det_chart = (singular.NI - singular.NHI).det()
    assert status == 'center_singular' and np.array_equal(a, lower) and np.array_equal(b, upper)
    assert det_chart == 0 and not det_full.contains(0)
    return dict(exact_dark=setup, exact_bad_chart_regular_forward=dict(status=status,
                tested_pivots=singular.pivots, public_midpoint_choice=chosen, determinant_chart=str(det_chart),
                determinant_forward=str(det_full), whole_box_retained=True))


def main():
    tic = process_time()
    names = ('spectral-00', 'spectral-18', 'fresh6-00', 'fresh6-09', 'confirm6-00', 'confirm6-04')
    rows = {}
    for public in load_cases():
        if public.case_id in names:
            rows[public.case_id] = diagnostics(public)
            print(public.case_id, rows[public.case_id], flush=True)
    data = dict(input_scope='old observations and public model only', cases=rows,
                algebraic_controls=degenerate_controls(), cpu_s=process_time() - tic)
    (ROOT / 'results/pivot_diagnostics.json').write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')


if __name__ == '__main__':
    main()
