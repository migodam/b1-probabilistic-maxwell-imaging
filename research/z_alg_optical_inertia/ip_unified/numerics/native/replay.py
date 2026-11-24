"""Full-grid frozen-array state, complete centers and compensation screens."""
import argparse
from pathlib import Path
import sys
from time import perf_counter, process_time

import numpy as np
from flint import acb_mat, arb, ctx
from scipy.linalg import eigh, qr
from scipy.sparse.linalg import LinearOperator, cg, gcrotmk, gmres
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
STUDY = HERE.parents[2]
ARRAY = STUDY / "physical_revision/array"
sys.path.insert(0, str(ARRAY))
from screen import inputs, load_json, save_json, times
from verify_native import (FourierEnvelope, Geometry, NativeState, analyze_row, ball,
                           digest_profile, frobenius_squared, gershgorin,
                           high, low, mid, outward_rhs, star)

CASES = {
    "u100_best": ("high_grid_p32", 1.0, "scalar_grid_601"),
    "u075": ("native_p32", 0.75, "stable_scalar_optical"),
    "u125": ("high_grid_p32", 1.25, "scalar_grid_566"),
    "u100_stable": ("native_p32", 1.0, "stable_scalar_optical"),
}
MATERIAL_FRACTIONS = np.linspace(0., 1., 121)


def herm(a):
    return (a+a.conj().T)*0.5


def configure(name):
    file, cap, profile = CASES[name]
    history = load_json(ARRAY / f"{file}.json")
    old = next(r for r in history["profiles"] if r["name"] == profile and r["cap"] == cap)
    spec, sensors, sources, observations = inputs()
    geo = Geometry(spec, sensors[32], sources[32])
    data = [(m, y) for m, y in observations[32][32] if m["noise_pct"] == 2]
    return geo, NativeState(geo, digest_profile(old["profile"], cap)), data, old


def solve_stationarity(native):
    g = native.g
    rhs = np.column_stack((g.sensor_star/native.w.scale, native.z[:, None]*g.incident))
    prop = np.empty_like(rhs)
    pre = LinearOperator((g.n,g.n), matvec=native.precondition, dtype=np.complex128)
    costs, recycled = [], []
    for j, b in enumerate(rhs.T):
        steps = []
        x, info = gcrotmk(native.operator, b, M=pre, CU=recycled, m=25, k=12,
                          maxiter=16, rtol=8e-12, atol=0., callback=lambda q: steps.append(1))
        if info:
            x, info = gmres(native.operator, b, x0=x, M=pre, restart=100, maxiter=12,
                            rtol=8e-12, atol=0.)
        error = np.linalg.norm(native.numeric_action(x)-b)
        if info or not np.isfinite(error):
            raise ArithmeticError(f"Stationary solve failed at column {j} with status {info}")
        prop[:, j] = x
        costs.append(dict(cycles=len(steps), status=int(info), abs_residual=float(error),
                          rel_residual=float(error/np.linalg.norm(b))))
    return rhs, prop, costs


def normalized_modes(native, vec):
    basis, _ = qr(vec, mode="economic", check_finite=False)
    vals, rot = eigh(herm(basis.conj().T@native.numeric_action(basis)),check_finite=False)
    if vals.max()>=0:
        raise ArithmeticError("Ritz span is not negative under Rayleigh rotation.")
    return vals,basis@rot


def frame_weight_upper(native, vec, eig):
    g = native.g
    norm2 = np.sum(abs(vec)**2, axis=1)
    symbol = native.numeric_symbol.real
    diagonal = 1 + symbol.mean() + np.sum(abs(g.sensor)**2, axis=0)/native.w.scale
    # Analytical upper uses K_ii + [(-K)_+]_ii over (1 - P_-[ii])**2.
    k_pos_diag = diagonal + (abs(vec)**2) @ (-eig)
    return k_pos_diag / (1-norm2)**2


def global_screen(native, n, gn, weight):
    a0 = herm(n.conj().T@(weight[:, None]*n))
    a1 = n.conj().T@(weight[:, None]*gn)
    a1 = a1+a1.conj().T
    a2 = herm(gn.conj().T@(weight[:, None]*gn))
    caps = native.w.cap*MATERIAL_FRACTIONS
    spectra = []
    for x in caps:
        spectra.append(np.linalg.eigvalsh(a0-x*a1+x*x*a2))
    spectra = np.array(spectra)
    # Choose the largest certified-dimension proposal, then the largest spectral slack.
    num = (spectra < 1).sum(axis=1)
    best = np.flatnonzero(num == num.max())
    i = int(best[np.argmin(spectra[best, max(0, num.max()-1)])])
    vals, sub = eigh(a0-caps[i]*a1+caps[i]**2*a2)
    return dict(x=float(caps[i]), contraction_dimension=int(num[i]),
                upper_compensated_dimension=int(n.shape[1]-num[i]),
                contraction_eigenvalues=vals.tolist()), caps, spectra, sub


def select_rows(n, gn, cap, weak):
    # Pivot rows in the weakest remaining coordinate directions and at both endpoints.
    stack = np.column_stack((n@weak, (n-cap*gn)@weak))
    p = qr(stack.T, pivoting=True, mode="r", check_finite=False)[1]
    return p[:4*n.shape[1]]


def projected_gram(native, v, pixels):
    g, u = native.g, native.w.cap
    eye = np.zeros((g.n, len(pixels)), dtype=complex)
    eye[pixels, np.arange(len(pixels))] = 1.
    endpoints = np.column_stack((eye, eye-u*g.numerical_green(eye, adjoint=True)))
    project = lambda t: t-v@(v.conj().T@t)
    # Add a unit block on the deflated span to leave a uniformly positive operator.
    def action(t):
        return project(native.numeric_action(project(t))) + v@(v.conj().T@t)
    def precondition(t):
        return project(native.precondition(project(t))) + v@(v.conj().T@t)
    p = LinearOperator((g.n,g.n), matvec=precondition, dtype=np.complex128)
    k = LinearOperator((g.n,g.n), matvec=action, dtype=np.complex128)
    trial, costs = np.empty_like(endpoints), []
    for j, r in enumerate(endpoints.T):
        b = project(r)
        steps = []
        value, status = cg(k, b, M=p, rtol=5e-11, atol=0., maxiter=240,
                           callback=lambda q: steps.append(1))
        trial[:, j] = project(value)
        error = np.linalg.norm(b-action(trial[:, j]))
        costs.append(dict(iterations=len(steps), status=int(status), abs_residual=float(error)))
    return endpoints, trial, herm(endpoints.conj().T@trial), costs


def outward_sources(native, state, data, response, saved_center, saved_field, saved_h):
    g = native.g
    rhs = outward_rhs(native, g.incident)
    x = ball(response)
    ax = native.exact_action(x)
    error = ax-rhs
    stat = star(x)*ax-star(x)*rhs-star(rhs)*x
    gram = star(error)*error
    e0 = g.exact_green(x)
    result = []
    for (meta, y), js, es, hs in zip(data, saved_center, saved_field, saved_h, strict=True):
        yb = ball(y)
        comb = acb_mat(2*g.s, g.s, yb.entries()+acb_mat(
            g.s,g.s,[1 if i == j else 0 for i in range(g.s) for j in range(g.s)]).entries())
        residual = error*comb
        gap = min(state["positive_gap_lower"], state["negative_gap_lower"])
        center_error = high(frobenius_squared(residual).sqrt()/arb(gap))
        center = x*comb
        e = g.incident_ball + e0*comb
        stat_trial = star(yb)*yb + arb(native.w.scale)*(star(comb)*stat*comb)
        row = analyze_row(y, state, native.w.scale, stat, gram,
                          list(range(2*g.s)), meta["padded_radius"])
        row.update(center_trial_defect_frobenius_upper=high(frobenius_squared(residual).sqrt()),
                   center_error_frobenius_upper=center_error,
                   source_form_saved_to_trial_F_error_upper=high(frobenius_squared(
                       ball(hs)-stat_trial).sqrt()),
                   center_saved_to_trial_F_error_upper=high(frobenius_squared(ball(js)-center).sqrt()),
                   total_center_field_saved_to_trial_F_error_upper=high(frobenius_squared(ball(es)-e).sqrt()))
        result.append(row)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=CASES)
    args = parser.parse_args()
    t, c = perf_counter(), process_time()
    geo, native, data, old = configure(args.case)
    path = HERE / args.case
    expected = old["state"]["exact_nu"]
    print(args.case, "state", flush=True)
    modes = native.modes(expected)
    if native.nu != expected or modes["numerical_space"] != "native K/mu":
        raise ArithmeticError("Ritz proposal does not match the inherited native negative count.")
    eig, vec = normalized_modes(native, native.eig_basis[:, :expected])
    native.eig_basis[:, :expected] = vec
    n = vec/np.sqrt(-eig)
    gn = geo.numerical_green(n)
    wupper = frame_weight_upper(native, vec, eig)
    global_probe, homogeneous_x, homogeneous_spectra, basis = global_screen(native, n, gn, wupper)
    print(args.case, "global probe", global_probe, flush=True)
    pixels = select_rows(n, gn, native.w.cap, basis[:, global_probe["contraction_dimension"]:])
    if len(pixels) == 0:
        pixels = select_rows(n, gn, native.w.cap, basis[:, -1:])
    print(args.case, "stationarity", flush=True)
    stage = perf_counter(), process_time()
    rhs, response, source_cost = solve_stationarity(native)
    charge = dict(stationarity=dict(**times(*stage), columns=len(source_cost), details=source_cost))
    y = np.stack([y for _, y in data])
    center = response[:, :32]@y+response[:, 32:]
    e0 = geo.incident + geo.numerical_green(response[:, :32])@y+geo.numerical_green(response[:, 32:])
    right = rhs[:, :32]@y+rhs[:, 32:]
    hs = np.stack([herm(yi.conj().T@yi-native.w.scale*ri.conj().T@ji)
                   for yi, ri, ji in zip(y, right, center, strict=True)])
    rows = []
    ordered = np.linalg.eigvalsh(hs)[:, ::-1]
    for (meta, _), vals in zip(data, ordered, strict=True):
        r = meta["padded_radius"]
        ds = np.array([np.sqrt(np.maximum(vals[expected-c:], 0).sum())
                       for c in range(expected+1)])
        feasible = np.flatnonzero(ds > r)
        rows.append(dict(**meta, numerical_H_eigenvalues=vals.tolist(), baseline_ratio=float(ds[0]/r),
                         optimistic_psd_charge_ratio_by_compensated_dimension=(ds/r).tolist(),
                         min_optimistic_dimension_for_exclusion=int(feasible[0]) if feasible.size else None))
    save_json(path.with_suffix(".json"), dict(case=args.case, profile=native.w.as_dict(), numerics=modes,
              rows=rows, global_probe=global_probe, stages=charge, total=times(t,c)))
    print(args.case, "endpoint Gram", flush=True)
    stage = perf_counter(), process_time()
    endpoint, trial, gamma, g_cost = projected_gram(native, vec, pixels)
    gram_val = np.linalg.eigvalsh(gamma)
    charge["projected_gram"] = dict(**times(*stage), columns=len(g_cost), details=g_cost,
                                   eigenvalues=gram_val.tolist())
    np.savez_compressed(path.with_suffix(".npz"), negative_vec=vec, negative_eig=eig, negative_gn=gn,
                        center=center, center_field=e0, hs=hs, source_response=response,
                        pixel_index=pixels, endpoint_trial=trial, endpoint_gram=gamma,
                        homogeneous_x=homogeneous_x, homogeneous_spectra=homogeneous_spectra,
                        homogeneous_subspace=basis, frame_weight_upper_numeric=wupper)
    print(args.case, "Arb state and source", flush=True)
    stage = perf_counter(), process_time()
    env = FourierEnvelope(geo)
    negative, negative_gap = native.ensure_negative(expected, return_gap=True)
    proof = env.verify(native.w, expected)
    state = dict(exact_nu=expected, negative_gap_lower=negative_gap,
                 negative_congruence_upper=negative, **proof)
    charge["outward_state"] = times(*stage)
    stage = perf_counter(), process_time()
    strict = outward_sources(native, state, data, response, center, e0, hs)
    charge["outward_source"] = times(*stage)
    for row, interval in zip(rows, strict, strict=True):
        row["outward_baseline"] = interval
    result = dict(case=args.case, profile=native.w.as_dict(), numerics=modes, verified_state=state,
                  rows=rows, global_probe=global_probe, stages=charge,
                  operations=dict(geo.cost), total=times(t,c))
    save_json(path.with_suffix(".json"), result)
    print(args.case, "done", result["total"], flush=True)


if __name__ == "__main__":
    ctx.prec, ctx.threads = 112, 2
    with threadpool_limits(1):
        main()
