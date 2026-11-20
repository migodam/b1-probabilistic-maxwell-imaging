"""Independent full-wave checks for the material readout and transported receiver prediction."""

import json
from time import process_time

import numpy as np
from scipy.linalg import lu_solve
from threadpoolctl import threadpool_limits

from cohorts import load_cohort
from current_cycle import OUT, finite_line, norm, physical, voxel_ratio


def rel(lhs, rhs):
    return norm(lhs - rhs) / max(norm(lhs), norm(rhs), 1.0e-20)


def main():
    rng = np.random.default_rng(463)
    start = process_time()
    z = lambda shape: rng.normal(size=shape) + 1j * rng.normal(size=shape)
    j, e, d, gd = (z((32, 12)) for _ in range(4))
    prior = np.ones(len(j))
    bounds = (0.02, 20.0)
    e[0] = 0
    m = voxel_ratio(j, e, prior, bounds)
    a = np.sum(e.conj() * j, axis=1)
    power = np.sum(abs(e) ** 2, axis=1)
    rho = np.divide(a.real, power, out=prior.copy(), where=power > 0)
    minors = sum(
        abs(j[:, s] * e[:, t] - j[:, t] * e[:, s]) ** 2
        for s in range(12)
        for t in range(s + 1, 12)
    )
    left = power * np.sum(abs(j - m[:, None] * e) ** 2, axis=1)
    right = minors + a.imag**2 + power**2 * (rho - m) ** 2
    wedge = rel(left, right)
    j, e, d, gd = j[1:], e[1:], d[1:], gd[1:]
    u, v = z((48, 12)), z((48, 12))
    w = 0.32
    roots = finite_line(j, e, d, gd, u, v, w)

    def proxy(t):
        et = e + t * gd
        jt = j + t * d
        ee = np.sum(abs(et) ** 2, axis=1)
        jj = np.sum(abs(jt) ** 2, axis=1)
        inner = np.sum(et.conj() * jt, axis=1).real
        return norm(u + t * v) ** 2 + w * sum(
            (ee * jj - inner**2) / np.sum(abs(e) ** 2, axis=1)
        )

    h = 1.0e-5
    stationarity = max(
        abs(proxy(t + h) - proxy(t - h)) / (2 * h * proxy(t)) for t in roots
    )
    model, (case, *_) = load_cohort("r10_twod")
    G, S = model["G"], model["S"]
    m = 0.5 + rng.random(len(G))
    ja, ya, fa = physical(model, m)
    direction = rng.uniform(-1.0, 1.0, len(G))
    dj = z(ja.shape) * 0.1
    sensitivity = lu_solve(fa, S.T, trans=1, check_finite=False).T
    entries = []
    for strength in (0.2, 0.1, 0.05, 0.025):
        mp = m + strength * direction
        trial = ja + strength * dj
        ef = model["E0"] + G @ trial
        defect = trial - mp[:, None] * ef
        jp, yp, fp = physical(model, mp)
        transported = ya + sensitivity @ ((mp - m)[:, None] * ef)
        corrected = S @ trial - sensitivity @ defect
        remainder = -sensitivity @ (
            (mp - m)[:, None] * (G @ lu_solve(fp, defect, check_finite=False))
        )
        entries.append(
            dict(
                step=strength,
                prediction_err=norm(yp - transported),
                remainder_rel=rel(yp - transported, remainder),
                corrected_rel=rel(transported, corrected),
                fixed_residual_rel=rel(S @ trial - yp, S @ lu_solve(fp, defect)),
            )
        )
    print(
        "minors",
        wedge,
        "cubic FD",
        stationarity,
        "remainder",
        [x["remainder_rel"] for x in entries],
    )
    assert wedge < 1.0e-12 and m.size == 576 and stationarity < 1.0e-7
    assert max(x["remainder_rel"] for x in entries) < 1.0e-7
    assert max(x["fixed_residual_rel"] for x in entries) < 1.0e-10
    (OUT / "current_identities.json").write_text(
        json.dumps(
            dict(
                wedge_relative=wedge,
                cubic_roots=roots,
                cubic_relative_fd=stationarity,
                zero_field_keeps_previous=bool(
                    voxel_ratio(j[:1], e[:1] * 0, prior[:1], bounds)[0] == 1.0
                ),
                transport=entries,
                measured_cpu_seconds=process_time() - start,
            ),
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    with threadpool_limits(1):
        main()
