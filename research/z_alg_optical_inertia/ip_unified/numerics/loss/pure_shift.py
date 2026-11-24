"""Bound every scalar multiplier on the loss-only ellipsoid relaxation."""
import argparse
from collections import deque
import json

from flint import acb, acb_mat, arb, ctx
import numpy as np
from threadpoolctl import threadpool_limits

from finite_states import spectral_enclosure
from run import HERE, clock, elapsed, finite_geometry, write_json
from verify import (FiniteExact, ball, blocks, exact_number, eye, gram_interval,
                    hi, lo, norm, star)


def center_ball(center, error):
    radius = (error / 2).upper()
    return acb_mat(center.nrows(),center.ncols(),[
        value + acb(arb(0,radius),arb(0,radius)) for value in center.entries()
    ])


def certify_all_scalar(geo, exact, y, rho, cap, tangent, target_ratio="0.999"):
    start = clock()
    trial, error, _, setup = gram_interval(geo,exact,exact_number(cap),exact_number(tangent))
    r, c, center = blocks(geo,trial,arb(0))
    d = ball(y) - center_ball(center,error)
    if not spectral_enclosure(r)[0][0] > 0:
        raise ArithmeticError("Positive receiver Gram gap was not enclosed.")
    c_min = spectral_enclosure(c)[0][0].lower()
    if not c_min > 0:
        raise ArithmeticError("Positive source Gram gap was not enclosed.")
    # At tau >= ||D||_2^2/lambda_min(C), A(tau) <= 0.
    tau_max = hi(norm(d)**2/c_min)
    target = arb(rho)*arb(target_ratio)
    target2 = target**2
    evaluations = 0

    def trace_upper(a):
        nonlocal evaluations
        evaluations += 1
        coeff = r + eye(geo.p)*arb(a)
        solved = coeff.solve(d)
        a_mat = star(d)*solved-c
        return sum((max(arb(0),v.upper()) for v in spectral_enclosure(a_mat)[1]),arb(0)).upper()

    trace_at_0 = trace_upper(0.)
    a0 = max(np.finfo(float).tiny, lo(target2 / (2*trace_at_0)))
    ends = [a0]
    while ends[-1] < tau_max:
        ends.append(min(tau_max, 1.5*ends[-1]))
    queue = deque(zip(ends[:-1], ends[1:]))
    cache = {0.:trace_at_0}
    highest = (arb(a0)*trace_at_0).upper()
    rigorous_intervals = 1
    while queue:
        a,b = queue.popleft()
        if a not in cache:
            cache[a] = trace_upper(a)
        # A(tau) decreases in Loewner order. Therefore on [a,b],
        # tr[H_tau]_+ = tau tr[A(tau)]_+ <= b tr[A(a)]_+.
        bound = (arb(b)*cache[a]).upper()
        if bound < target2:
            highest = max(highest, bound)
            rigorous_intervals += 1
            continue
        if b/a < 1 + 2e-6:
            raise ArithmeticError(f"Scalar shift upper goal unresolved on [{a:.6e},{b:.6e}].")
        m = float(np.sqrt(a*b))
        queue.extend(((a,m),(m,b)))
    return dict(all_tau_distance_upper=hi(highest.sqrt()), tau_max_above_which_exact_zero=tau_max,
                alpha_setup=setup, PSD_cb_gap_lower=lo(c_min), intervals=rigorous_intervals,
                eigensolves=evaluations, **elapsed(start))


def main():
    data = json.loads((HERE/"finite_0.json").read_text())
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--weak-loss", action="store_true", help="Bound all common prices at t=0.05, u=0.75.")
    args = parser.parse_args()
    name = "loss_tm_01__fixed_t0.05__pair2pct" if args.weak_loss else "loss_tm_01__fixed_t0.2__pair2pct"
    row = next(o for o in data["observations"] if o["id"]==name)
    cap = .75 if args.weak_loss else 1.
    with np.load(HERE/"finite_0_inputs.npz",allow_pickle=False) as saved:
        y = saved[name]
    geo = finite_geometry(.6)
    bound = certify_all_scalar(geo,FiniteExact(geo,.6),y,row["rho"],cap,row["tangent"],
                               "0.81" if args.weak_loss else "0.999")
    filename = "finite_0_weak_loss_all_prices.json" if args.weak_loss else "finite_0_pure_shift_verified.json"
    write_json(HERE/filename,dict(id=name,tangent=row["tangent"],cap=cap,rho=row["rho"],**bound))
    print(bound,flush=True)


if __name__ == "__main__":
    ctx.prec, ctx.threads = 128,1
    with threadpool_limits(1):
        main()
