"""One-sided replay for the exact binary64 dense nominal observation contract.

Finite released coefficients are interpreted as exact dyadics. The argument
uses the inverse-defect test, with normwise forward-error allowances for
round-to-nearest binary64 products and inner products. It certifies the
finite released arrays. It does not enclose a continuous kernel or unmodeled
calibration. The large-grid original FFT residual check remains numerical.
"""
from __future__ import annotations

import argparse
import json

from shared import Model, OUT, elapsed, now, write_json

import numpy as np
from scipy.linalg import lu_factor, lu_solve
from threadpoolctl import threadpool_limits

EPS = np.finfo(float).eps
TINY_GUARD = 2.**-900


def up(value):
    value = float(value)
    if not np.isfinite(value) or value < 0:
        raise ArithmeticError("Nonfinite or negative quantity in an upper bound")
    return float(np.nextafter(value * (1 + 8*EPS) + TINY_GUARD, np.inf))


def gamma(operations):
    q = operations * EPS
    if q >= .01:
        raise ArithmeticError("Outside the rounding allowance contract")
    return up(q / (1 - q))


def norm(z):
    operations = 8 * z.size + 32
    return up(np.linalg.norm(z) * (1 + gamma(operations)))


def enclosed_witness(model, x, roi=None):
    t = now()
    if model.green.dense is None:
        raise ValueError("Use original-operator numerical replay for native FFT geometry")
    if np.any(~np.isfinite(x)) or np.any(x < 0) or np.any(x > 4.5):
        raise ValueError("Illegal contrast")
    if roi is not None and np.any(x[roi] != 0):
        raise ValueError("Zero-region candidate violates its exact box")
    active = np.flatnonzero(x)
    xo = x[active, None]
    # Inactive current rows are exactly zero, so only the active pivot is
    # needed for the uniqueness and residual test.
    g = model.green.dense[np.ix_(active, active)]
    b = model.B[active]
    s = model.S[:, active]
    y = model.Y
    n = len(active)
    eye = np.eye(n, dtype=complex)
    a = eye - xo * g
    f = xo * b
    fac = lu_factor(a, check_finite=False)
    inverse = lu_solve(fac, eye, check_finite=False)
    j = lu_solve(fac, f, check_finite=False)
    state_defect_observed = a @ j - f
    inverse_defect_observed = eye - inverse @ a
    data_defect_observed = s @ j - y
    xmax = float(np.max(xo))
    a_n, g_n, b_n = norm(a), norm(g), norm(b)
    inv_n, j_n, f_n = norm(inverse), norm(j), norm(f)
    eye_n, s_n, y_n = norm(eye), norm(s), norm(y)
    a_err = up(gamma(12) * (eye_n + xmax*g_n))
    f_err = up(gamma(8) * xmax*b_n)
    dot_err = gamma(16*n + 32)
    inv_defect = up(norm(inverse_defect_observed) + inv_n*a_err +
                    dot_err*(eye_n + inv_n*a_n))
    result = dict(contract="Exact released binary64 G,S,E0,Y and radius, nominal calibrated operator",
                  material_active=n, inverse_defect_upper=inv_defect)
    if inv_defect >= 1.:
        return dict(result, verified=False, reason="inverse_defect_not_below_one", **elapsed(t))
    inv_bound = up(inv_n / np.nextafter(1-inv_defect, 0.))
    state_upper = up(norm(state_defect_observed) + a_err*j_n + f_err +
                     dot_err*(a_n*j_n + f_n))
    data_upper = up(norm(data_defect_observed) + dot_err*(s_n*j_n+y_n) +
                    s_n*inv_bound*state_upper)
    return dict(result, verified=data_upper < model.radius,
                residual_float=float(np.linalg.norm(data_defect_observed)),
                residual_upper=data_upper, radius=model.radius,
                residual_upper_over_radius=up(data_upper/model.radius),
                exact_inverse_norm_upper=inv_bound, original_state_defect_upper=state_upper,
                arithmetic_dimensions=dict(active=n, observed=int(y.size)), **elapsed(t))


def case_replay(case, stationary=False):
    model = Model.read(case)
    p = OUT / case
    if stationary:
        p /= "stationary"
    with np.load(p / "baseline.npz", allow_pickle=False) as z:
        rois = z["rois"].reshape(-1, model.n)
    rec = {}
    with np.load(p / "b1_predictions.npz", allow_pickle=False) as z:
        report = json.loads((p / "b1.json").read_text())
        for k, roi in enumerate(rois):
            name = f"roi_{k+1:02d}"
            if report["witness_fits"][name]["residual_over_radius"] < 1:
                rec[name] = enclosed_witness(model, z[f"x_{name}"].ravel(), roi)
        choice = report["selected_roi"]
        if choice:
            rec["selected"] = enclosed_witness(model, z["chi"].ravel(),
                                                rois[int(choice.rsplit("_", 1)[1])-1])
    with np.load(p / "b1_controls.npz", allow_pickle=False) as z:
        x = z["physical_chi"].ravel()
        if model.audit(x)["residual_over_radius"] < 1:
            rec["unconstrained"] = enclosed_witness(model, x)
            rec["unconstrained"]["roi_means"] = [float(x[roi].mean()) for roi in rois]
    write_json(p / "one_sided_witness.json", rec)
    print(case, {k:(v["verified"],v.get("residual_upper_over_radius")) for k,v in rec.items()},
          flush=True)


if __name__ == "__main__":
    parse = argparse.ArgumentParser()
    parse.add_argument("cases", nargs="+")
    parse.add_argument("--stationary", action="store_true")
    args = parse.parse_args()
    with threadpool_limits(limits=1,user_api="blas"):
        for case in args.cases:
            case_replay(case, args.stationary)
