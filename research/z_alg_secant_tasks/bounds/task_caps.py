"""Exclude complete near-extremal task regions at a specified physical reading."""
import json
from time import perf_counter, process_time

from box_modulus import (ROOT, abs_upper, acb, arb, ball_array, ctx, dressed, np,
                         positive_resolvent, real_trial, realify, up)
from scipy.optimize import lsq_linear
from threadpoolctl import threadpool_limits


def sets(arm):
    base = [(0, 0), (1, 0), (0, 1)]
    s = 1 / arb(2).sqrt()
    return {"base": [(0, base, arb(1))],
            "angle": [(0, base + [(0, 2), (1, 2), (0, 3)], s)],
            "frequency": [(0, base, s), (1, base, (arb(6)/5).sqrt()*s)],
            "rich": [(0, [(a, b) for b in range(4) for a in range(6)], s)]}[arm]


def shift(dress, h):
    from box_modulus import acb_mat, eye, from_mat
    p, g, e = dress
    d = h[:, None]
    inv = from_mat(acb_mat((eye(len(h)) - d * g).tolist()).inv())
    x = inv @ (d * e)
    return p @ x, (p @ inv, g @ inv, e + g @ x)


def direction(j, diff, mask, radius):
    # Pick one separating linear functional, then rigorously evaluate it.
    mat, y = real_trial(realify(j)), real_trial(realify(diff))
    a = mat[:, ~mask] * real_trial(radius[~mask])
    fit = lsq_linear(a, -y, bounds=(-1, 1), method="bvls", tol=1e-12)
    z = y + a @ fit.x
    if not np.linalg.norm(z):
        return np.zeros(len(z))
    return z / np.linalg.norm(z)


def near_face(refs, mask, radii, sign, arm, eta, delta):
    spec = sets(arm)
    diff, j = [], []
    for f, channels, scale in spec:
        df, (p, g, e) = refs[f]
        for a, b in channels:
            diff.append(scale * df[a, b])
            j.append(scale * p[a] * e[:, b])
    diff, j = np.array(diff), np.array(j)
    n = len(radii)
    weight = [arb(float(x)) for x in direction(j, diff, mask, radii)]
    nc = len(diff)
    # Conjugated complex weights pair as Re sum weights[j] * data[j].
    w = [acb(weight[i], -weight[nc+i]) for i in range(nc)]
    wn = sum(x*x for x in weight).sqrt()
    if not wn > 0:
        return "no_direction", 0.
    qinv = int(mask.sum())
    rr = np.array([min(2 * r, qinv * eta) if m else r
                   for r, m in zip(radii, mask)], dtype=object)
    wd = sum(wi * v for wi, v in zip(w, diff)).real
    wj = np.array([(sum(wi * row[i] for wi, row in zip(w, j))).real
                   for i in range(n)], dtype=object)
    linear_out = sum(abs_upper(wj[~mask]) * rr[~mask])
    pos = [max(arb(0), (sign * x).upper()) for x in wj[mask]]
    linear_in = min(sum(v * r for v, r in zip(pos, rr[mask])), qinv * eta * max(pos))
    # Contract all port Hessians before taking any absolute value.
    Q = np.array([[arb(0) for _ in range(n)] for _ in range(n)], dtype=object)
    off, tail = 0, arb(0)
    for f, channels, scale in spec:
        _, (p, g, e) = refs[f]
        k = abs_upper(g) * rr[None, :]
        ip = positive_resolvent(k)
        if ip is None:
            return "no_majorant", 0.
        higher = (abs_upper(p) * rr[None, :]) @ k @ k @ ip @ abs_upper(e)
        for a, b in channels:
            c = scale * w[off]
            coeff = (c*p[a][:, None]) * g * e[:, b][None, :]
            Q += np.frompyfunc(lambda x: x.real, 1, 1)(coeff)
            tail += abs(c) * higher[a, b]
            off += 1
    quadratic_lo = arb(0)
    for i, ri in enumerate(rr):
        quadratic_lo += min(arb(0), Q[i, i].lower()) * ri**2
        for k in range(i):
            coeff = Q[i, k] + Q[k, i]
            lo = min(arb(0), coeff.lower()) if mask[i] and mask[k] else -abs(coeff)
            quadratic_lo += lo * ri * rr[k]
    lower = (wd - linear_in - linear_out + quadratic_lo - tail) / wn
    return ("excluded" if lower > delta else "unknown",
            max(0., float(np.nextafter(float(lower.lower()), -np.inf))))


if __name__ == "__main__":
    start, cpu = perf_counter(), process_time()
    data = json.loads((ROOT/"array/results/summary.json").read_text())
    raw = np.load(ROOT/"array/results/physical_arrays.npz")
    rows = []
    with threadpool_limits(limits=1):
        for case in data["cases"]:
            n = case["n"]
            center = np.array(case["material_center"])
            mask = raw[case["id"]+"__roi_mask"].astype(bool)
            nominal = [dressed(raw,n,center,f) for f in ("f0","f1")]
            for width in (.05,.25,.75):
                r = np.array([min(arb(width),arb(float(c))-arb(.5),arb(30)-arb(float(c)))
                              for c in center],dtype=object)
                width_task = 2*sum(r[mask])/int(mask.sum())
                for sign in (-1,1):
                    h = np.array([sign*radius if inside else arb(0)
                                  for radius,inside in zip(r,mask)],dtype=object)
                    refs = [shift(x,h) for x in nominal]
                    for arm in ("base","angle","frequency","rich"):
                        delta = arb(next(row["radius"] for row in data["pairs"]
                                        if row["case"]==case["id"] and row["arm"]==arm))
                        for portion in (0.,.05,.125,.25):
                            eta=arb(portion)*width_task
                            state, lb=near_face(refs,mask,r,sign,arm,eta,delta)
                            rows.append(dict(case=case["id"],width=width,sign=sign,arm=arm,
                                             prior_range=up(width_task),fraction=portion,
                                             removed_task_if_proved=float(eta),status=state,
                                             residual_lower=lb,delta=float(delta)))
    result={"scope":"archived discrete arrays, specified F(center), calibrated deterministic 10% budget",
            "precision_bits":ctx.prec,"wall_s":perf_counter()-start,"cpu_s":process_time()-cpu,"rows":rows}
    (ROOT/"bounds/caps.json").write_text(json.dumps(result,indent=2,allow_nan=False))
    print(len(rows),sum(z["status"]=="excluded" for z in rows),result["wall_s"])
