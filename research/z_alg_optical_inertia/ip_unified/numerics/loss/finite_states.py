"""Exact relative audit of selected finite-grid state and shift profiles."""
import json

from flint import acb, arb, ctx
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import minimize_scalar
from threadpoolctl import threadpool_limits

from run import HERE, adj, clock, elapsed, finite_geometry, tail, write_json
from verify import (FiniteExact, ball, blocks, exact_number, eye, gram_interval, hi, lo, norm,
                    num, star)


def spectral_enclosure(mat):
    spectrum, vectors = eigh(num(mat))
    vec, diag = ball(vectors), ball(np.diag(spectrum))
    unitary_error = norm(star(vec)*vec-eye(len(spectrum)))
    if not unitary_error < 1:
        raise ArithmeticError("Spectral trial columns were not verified independent.")
    # Let U be the polar unitary of vec. The perturbation from vec D vec*
    # to U D U* is bounded by ||D|| * ||vec*vec-I||.
    bound = (norm(mat - vec * diag * star(vec)) +
             max(abs(arb(float(v))) for v in spectrum) * unitary_error).upper()
    return [arb(float(v))-bound for v in spectrum], [arb(float(v))+bound for v in spectrum], hi(bound)


def verified_profile(geo, exact, y, cap, tangent, mode, profile):
    started = clock()
    beta = 1 / ((1+tangent*tangent)*cap)
    alpha = tangent * beta
    opt = arb(profile["weight"])
    mu = arb(profile["mu"])
    a_re = 1-opt*tangent if mode == "fixed_phase" else arb(1)
    delta = beta + opt*alpha if mode == "shift" else beta
    z = acb(a_re, opt) / 2
    yb = ball(y)
    k = star(exact.s) * exact.s / mu - a_re * (exact.g + star(exact.g)) / 2 + opt * exact.r + eye(geo.n)*delta
    l = star(exact.s) * yb / mu + z * exact.b
    k_lower,k_upper,k_error = spectral_enclosure(k)
    nu = sum(v < 0 for v in k_upper)
    if any(not (a > 0 or b < 0) for a,b in zip(k_lower,k_upper,strict=True)):
        return dict(status="inertia_unresolved", state_spectral_error=k_error, **elapsed(started))
    gap = min((-b if b < 0 else a).lower() for a,b in zip(k_lower,k_upper,strict=True))
    if nu >= geo.inputs:
        return dict(distance_lower=0.,distance_upper=0.,nu=nu,spectral_gap_lower=lo(gap), **elapsed(started))
    vals, vec = eigh(num(k))
    j_trial = ball(vec @ ((adj(vec) @ num(l))/vals[:,None]))
    kj = k*j_trial
    defect = norm(kj-l)
    gram = star(yb)*yb + mu*(star(j_trial)*kj - star(j_trial)*l - star(l)*j_trial)
    sterr = (mu*defect**2/gap).upper()
    low_est,high_est,h_error = spectral_enclosure(gram)
    lower = sum((max(arb(0),(v-sterr).lower()) for v in low_est[:geo.inputs-nu]),arb(0)).sqrt()
    upper = sum((max(arb(0),(v+sterr).upper()) for v in high_est[:geo.inputs-nu]),arb(0)).sqrt()
    return dict(distance_lower=max(0.,lo(lower)),distance_upper=hi(upper),nu=nu,
                spectral_gap_lower=lo(gap),state_defect_upper=hi(defect),
                state_spectral_error=k_error,source_spectral_error=h_error, **elapsed(started))


def pure_shift_score(geo, y, trial):
    started = clock()
    rs, cb, fc = map(num, blocks(geo, trial, arb(0)))
    vals, v = eigh(rs)
    z = adj(v) @ (y-fc)

    def distance(theta):
        sigma = np.exp(theta)
        return tail(sigma * (adj(z) @ (z/(vals[:,None] + sigma)) - cb))

    ts = np.linspace(-28, 4, 250)
    values = [distance(t) for t in ts]
    i = int(np.argmax(values))
    best = values[i], ts[i]
    if 0 < i < len(ts) - 1:
        opt = minimize_scalar(lambda v: -distance(v), bounds=(ts[i-1],ts[i+1]),method="bounded")
        if -opt.fun > best[0]:
            best = -float(opt.fun),float(opt.x)
    return dict(tau=float(np.exp(best[1])), distance=best[0], **elapsed(started))


def replay(study, width):
    start = clock()
    geo = finite_geometry(width)
    exact = FiniteExact(geo,width)
    res = json.loads((HERE/f"{study}.json").read_text())
    rows = []
    target = HERE/f"{study}_state_verified.json"
    with np.load(HERE/f"{study}_inputs.npz",allow_pickle=False) as data:
        for arm in res["receivers"]:
            cap,t = arm["cap"],arm["tangent"]
            ar_u,ar_t = exact_number(cap),exact_number(t)
            trial, _, _, times = gram_interval(geo,exact,ar_u,ar_t)
            for r in arm["rows"]:
                if "state_screen" not in r:
                    continue
                ident = r["id"]
                y, rho = data[ident],r["rho"]
                old = r["state_screen"]["best"]
                best = dict(id=ident,cap=cap,tangent=t,rho=rho,
                            pure_loss_shift_screen=pure_shift_score(geo,y,trial),
                            state={}, gram_setup=times)
                include = cap in (.75,1.,1.25,1.5) or max(old[m]["distance"] for m in ("old","shift","fixed_phase")) > rho
                if include:
                    for m in ("old","shift","fixed_phase"):
                        rec = verified_profile(geo,exact,y,ar_u,ar_t,m,old[m])
                        best["state"][m] = {**rec,"profile":old[m]}
                rows.append(best)
                if include:
                    print(study, cap,t,ident,[(m,round(v.get("distance_lower",0.)/rho,4),
                                round(v.get("distance_upper",0.)/rho,4)) for m,v in best["state"].items()],
                          round(best["pure_loss_shift_screen"]["distance"]/rho,4),flush=True)
            write_json(target,dict(study=study,rows=rows,**elapsed(start)))


if __name__ == "__main__":
    import sys
    ctx.prec,ctx.threads=128,1
    with threadpool_limits(1):
        for study in (sys.argv[1:] or ["finite_0","finite_1"]):
            replay(study,{"finite_0":.6,"finite_1":1.25}[study])
