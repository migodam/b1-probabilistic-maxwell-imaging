"""Fixed geometric state weights and full native versus old negative-span replay."""
import argparse
from pathlib import Path
import sys
from time import perf_counter, process_time

import numpy as np
from flint import ctx
from scipy.linalg import eigh, ldl
from threadpoolctl import threadpool_limits


HERE = Path(__file__).resolve().parent
BASE = HERE.parent / "native"
sys.path.insert(0, str(BASE))
from replay import (FourierEnvelope, Geometry, NativeState, analyze_row, configure, herm, load_json,
                    save_json, solve_stationarity, times)
from native import Weights
from verify_native import reusable_stationarity


MULTIPLIERS = (1.5, 2.5)
CAP = 1.0
MASK = tuple(i in (0,3) or j in (0,3) for i in range(4) for j in range(4))
FLOORS = (.03125, .0078125, .001953125)


def clock():
    return perf_counter(), process_time()


def profile(scale, optical, factor, rebase=1.):
    real_weights = tuple(scale * (factor if edge else 1.) for edge in MASK)
    return Weights(CAP,scale*rebase,real_weights,(optical,)*len(MASK))


def spectrum_report(a, radius, indexes):
    h = eigh(herm(a),eigvals_only=True)[::-1]
    rows = {}
    for nu in sorted(set(indexes)):
        retained = np.maximum(h[nu:],0.)
        rows[str(nu)] = float(np.sqrt(retained.sum())/radius)
    return dict(source_eigenvalues_desc=h.tolist(),distance_over_radius_by_removed_count=rows)


def row_spectra(native, observations, rhs, stationary, nu, v):
    start = clock()
    scale = native.w.scale
    matrix = herm(v.conj().T@native.numeric_action(v))
    ep, ev = eigh(matrix)
    projected_negative = int((ep<0).sum())
    q = ev.conj().T @ (v.conj().T @ rhs)
    records = []
    for meta, y in observations:
        coef = np.vstack((y,np.eye(y.shape[1])))
        l = rhs@coef
        jj = stationary@coef
        d = y.conj().T@y - scale*herm(l.conj().T@jj)
        small = q@coef
        p_schur = y.conj().T@y - scale*herm(small.conj().T@(small/ep[:,None]))
        radius = meta["padded_radius"]
        rows = {
            "id":meta["id"],"radius":radius,
            "full_native_H":spectrum_report(d,radius,[nu,projected_negative]),
            "old_13_Galerkin_H_diagnostic":spectrum_report(p_schur,radius,[projected_negative]),
        }
        records.append(rows)
    return records,dict(old_13_compressed_eigenvalues=ep.tolist(),old_13_compressed_negatives=projected_negative,
                        old_13_nearest_abs_eigenvalue=float(min(abs(ep))),**times(*start))


def numerical():
    start = clock()
    original, base_state, data, prior = configure("u100_best")
    g = Geometry(original.spec, original.sensor, original.incident,block_side=4)
    del original
    mu, optical = base_state.w.scale,base_state.w.optical[0]
    with np.load(BASE/"u100_best_conditioned.npz",allow_pickle=False) as saved:
        original_span = np.linalg.qr(saved["negative_vec"])[0]
    result = dict(geometry="4 x 4 tiles; outer 12 are edge, inner 4 are core; fixed grid coordinates",
                  contrast_cap=CAP, original_profile="scalar_grid_601",scale=mu,
                  constant_optical_weight=optical, cases=[])
    for multiplier in MULTIPLIERS:
        wall = clock()
        before = g.cost.copy()
        weights = profile(mu,optical,multiplier)
        engine = NativeState(g,weights)
        diag = engine.modes(21,equilibrate=True)
        print("native numerical inertia",multiplier,engine.nu,flush=True)
        n_negative = engine.nu
        rhs,stationary,solve_cost = solve_stationarity(engine)
        rows, projection = row_spectra(engine,data,rhs,stationary,n_negative,original_span)
        proof = engine.ensure_negative(n_negative,return_gap=True) if n_negative else None
        gd = dict(g.cost-before)
        out = dict(edge_multiplier=multiplier, weights=weights.as_dict(),floating=diag,
                   floating_negative_count=n_negative,
                   outward_negative_count_at_least=n_negative,
                   outward_negative_congruence_upper=proof[0] if proof else None,
                   outward_negative_gap_lower=proof[1] if proof else None,
                   stationary_solver=solve_cost, old_13_projection=projection,
                   rows=rows,operations=gd,**times(*wall))
        result["cases"].append(out)
        print("numerics",multiplier,n_negative,projection["old_13_compressed_negatives"],
              max(r["full_native_H"]["distance_over_radius_by_removed_count"][str(n_negative)] for r in rows),
              out["wall_seconds"],flush=True)
        save_json(HERE/"edge_core.json",{**result,**times(*start)})

def block_inertia(a):
    """Numerical proposal only; D includes scalar and 2-by-2 pivots."""
    _, piv, _ = ldl(herm(a), lower=True, hermitian=True, check_finite=False)
    diag = piv.diagonal().real
    pair = np.flatnonzero(np.diagonal(piv,offset=1) != 0)
    keep = np.ones(len(diag),dtype=bool)
    keep[pair] = False
    keep[pair+1] = False
    middle = .5*(diag[pair]+diag[pair+1])
    radius = np.hypot(.5*(diag[pair]-diag[pair+1]),abs(piv[pair,pair+1]))
    spectrum = np.r_[diag[keep],middle-radius,middle+radius]
    return int((spectrum < -1.e-9).sum())


def envelope_proposal(envelope, scale, optical, factor, n_negative, *, wide):
    start = clock()
    best = None
    trials, floor_failures = 0, 0
    norm_scales = (1.,) if wide else (1.,(1.+factor)/2.,factor**.5)
    cuts = (10,12) if wide else (4,6,8)
    choices = ((8.,.125),(5.,.25),(3.,.125),(2.,.125),(1.5,0.)) if wide else (
                (2.,.25),(4.,.25),(8.,.25),(16.,.25),(4.,.375),(8.,.375))
    floors = (.001953125,) if wide else FLOORS
    for norm_scale in norm_scales:
        w = profile(scale,optical,factor,norm_scale)
        for cut in cuts:
            for beta,shift in choices:
                for gap in floors:
                    trials += 1
                    value = envelope.numerical_matrix(w,cut,beta,shift,gap)
                    if value is None:
                        floor_failures += 1
                        continue
                    matrix, negative, modes, iso, lower_envelope = value
                    removed = block_inertia(matrix)
                    upper = max(0, negative-removed)
                    if upper < n_negative:
                        # Do not accept a proposed upper contradicting the
                        # outward negative-space count from the full matrix.
                        continue
                    score = upper,-(gap*w.scale)
                    if best is None or score < best[0]:
                        best = (score,w,(lower_envelope,beta,shift,gap,modes,iso,removed),
                                dict(rebase=norm_scale,cut=cut,beta=beta,shift=shift,gap=gap,
                                     shape=list(matrix.shape),retained_negative_blocks=negative,
                                     floating_inertia_upper=upper,compressed_negative_demonstrated=removed))
                    if upper == n_negative:
                        return best,dict(trials=trials,floor_failures=floor_failures,**times(*start))
    return best,dict(trials=trials,floor_failures=floor_failures,**times(*start))


def verify(factors, *, wide):
    total = clock()
    numeric = load_json(HERE/"edge_core.json")
    original, _, data, _ = configure("u100_best")
    g = Geometry(original.spec,original.sensor,original.incident,block_side=4)
    del original
    env = FourierEnvelope(g)
    scale, optical = numeric["scale"],numeric["constant_optical_weight"]
    path = HERE/("edge_core_fourier_wide.json" if wide else "edge_core_verified.json")
    scope = "Outward Fourier verification uses the same physical coefficients; rebase changes normalization only."
    out = load_json(path) if path.is_file() and factors else dict(cases=[],scope=scope)
    elapsed_before = {key:out.get(key,0.) for key in ("wall_seconds","cpu_seconds")}
    if factors:
        out["cases"] = [row for row in out["cases"] if row["edge_multiplier"] not in factors]
    for case in numeric["cases"]:
        edge = case["edge_multiplier"]
        if factors and edge not in factors:
            continue
        start = clock()
        before = g.cost.copy()
        count_lower = case["outward_negative_count_at_least"]
        best, searched = envelope_proposal(env,scale,optical,edge,count_lower,wide=wide)
        result = dict(edge_multiplier=edge,lower_count=count_lower,bounded_search=searched)
        if best is None:
            result["blocker"] = "No bounded-trial Fourier identity floor and usable compressed negative space."
        else:
            score,w,trial,candidate = best
            result["proposal"] = candidate
            env.select = lambda _w, _known: trial
            try:
                proof = env.verify(w,score[0])
                result["outward"] = proof
            except ArithmeticError as exc:
                result["blocker"] = str(exc)
            else:
                if score[0]==count_lower:
                    neg_gap = case["outward_negative_gap_lower"]
                    # The negative gap was certified at the original scale.
                    neg_gap = None if neg_gap is None else float(np.nextafter(
                        neg_gap*(scale/w.scale),-np.inf))
                    state = dict(exact_nu=count_lower,positive_gap_lower=proof["positive_gap_lower"],
                                 negative_gap_lower=neg_gap)
                    native = NativeState(g,w)
                    stat,gram,indices,stationary = reusable_stationarity(native,{32:g.incident})
                    result["stationary_interval"] = stationary
                    result["outward_rows"] = [
                        dict(id=metadata["id"],radius=metadata["padded_radius"],
                             **analyze_row(y,state,w.scale,stat,gram,indices[32],metadata["padded_radius"]))
                        for metadata,y in data
                    ]
                else:
                    result["blocker"] = "Fourier bound leaves an uncertain middle spectrum; no source Schur interval is asserted."
        result.update(operations=dict(g.cost-before),**times(*start))
        out["cases"].append(result)
        print("interval",edge,count_lower,result.get("outward",{}).get("upper_count"),
              result.get("blocker"), result["wall_seconds"],flush=True)
        cost = {key:value+elapsed_before[key] for key,value in times(*total).items()}
        save_json(path,{**out,**cost})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--interval", action="store_true")
    parser.add_argument("--factor", action="append", type=float, choices=MULTIPLIERS)
    parser.add_argument("--wide", action="store_true",help="Run the limited extended-cutoff Fourier proof.")
    args = parser.parse_args()
    ctx.prec,ctx.threads=96,1
    with threadpool_limits(1):
        if args.interval:
            verify(args.factor,wide=args.wide)
        else:
            numerical()
