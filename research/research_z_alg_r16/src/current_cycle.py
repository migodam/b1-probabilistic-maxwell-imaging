"""Alternate an explicit FFT current step and an unrestricted real material readout on MAX-STRUCT."""

import argparse
import csv
import json
from pathlib import Path
from time import perf_counter, process_time

import numpy as np
from scipy.linalg import cho_factor, cho_solve, lu_factor, lu_solve, qr, svd
from scipy.optimize import minimize
from threadpoolctl import threadpool_limits

from cohorts import historical_study, load_cohort
from state_rom import fft_dictionary, regularizer

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
BACKGROUND = 0.1
MODES = 48
STEPS = 12
POLISH = 24
DEFAULT_COHORTS = (
    "r10_twod",
    "r10_conditional4",
    "r11_independent48",
    "r11_crossgrid8",
    "r11_bothblocks16",
    "r12_confirm24",
)
ARMS = ("voxel_only", "lagged48", "shared48", "shared36_plus12")


def norm(x):
    return float(np.linalg.norm(x))


def physical(model, m):
    matrix = np.asfortranarray(-m[:, None] * model["G"])
    matrix.flat[:: len(m) + 1] += 1.0
    factor = lu_factor(matrix, overwrite_a=True, check_finite=False)
    current = lu_solve(factor, m[:, None] * model["E0"], check_finite=False)
    return current, model["S"] @ current, factor


def use_archived_objective(model):
    import sys

    sys.path.insert(0, str(historical_study(10).parents[1]))
    from research.research_z_alg_r4_diversity.src.physics import Regularizer

    model["original_penalty"] = Regularizer(model["shape"]).evaluate
    model["optimizer_scale"] = 1.0e6
    model["optimizer_options"] = dict(maxcor=10, ftol=1.0e-12, gtol=0.01, maxls=30)


def material_penalty(model, m):
    if "original_penalty" in model:
        return model["original_penalty"](m)[:2]
    return regularizer(m, model["shape"])


def objective(model, m, data, scale):
    current, y, fac = physical(model, m)
    residual = y - data
    r, _ = material_penalty(model, m)
    return 0.5 * norm(residual) ** 2 / scale + r, current, residual, fac


def true_gradient(model, m, current, residual, factor, scale):
    adj = lu_solve(factor, model["S"].conj().T @ residual, trans=2, check_finite=False)
    field = model["E0"] + model["G"] @ current
    return ((adj.conj() * field).sum(axis=1)).real / scale + material_penalty(model, m)[
        1
    ]


def voxel_ratio(current, field, previous, bounds):
    power = np.sum(abs(field) ** 2, axis=1)
    numerator = np.sum(field.conj() * current, axis=1).real
    best = np.divide(numerator, power, out=previous.copy(), where=power > 0)
    return np.clip(best, *bounds)


def field_projection(z, field, power):
    fitted = np.sum(field.conj() * z, axis=1).real / power
    return field * fitted[:, None]


def chosen_fourier_atoms(shape, frequencies):
    """User-supplied (ky,kx) DFT modes of the complex current."""
    frequencies = np.asarray(frequencies)
    if (
        frequencies.ndim != 2
        or frequencies.shape[1] != 2
        or not np.all(frequencies == np.round(frequencies))
    ):
        raise ValueError("give two integer indices per current mode: ky kx")
    idx = frequencies.astype(int) % np.asarray(shape)
    if len(np.unique(idx, axis=0)) != len(idx):
        raise ValueError("the selected list repeats a DFT mode")
    modes = np.zeros((*shape, len(idx)), complex)
    modes[idx[:, 0], idx[:, 1], np.arange(len(idx))] = 1.0
    return np.fft.ifft2(modes, axes=(0, 1), norm="ortho").reshape(-1, len(idx))


def prepare(model, current_modes=MODES, bright_ratio=1.0e-2, frequencies=None):
    start = perf_counter()
    S = model["S"]
    G = model["G"]
    ul, s, vh = svd(S, full_matrices=False, check_finite=False)
    # The bright current is anchored only in the numerically well measured modes.
    if not 0 < bright_ratio <= 1:
        raise ValueError("bright-ratio must be in (0,1]")
    bright = int(np.sum(s >= bright_ratio * s[0]))
    if frequencies is not None:
        current_modes = len(frequencies)
    if not 1 <= current_modes <= S.shape[1] - bright:
        raise ValueError("current-rank exceeds the available receiver-dark dimension")
    v = np.ascontiguousarray(vh[:bright].conj().T)
    seq = (
        fft_dictionary(model["shape"], current_modes)
        if frequencies is None
        else chosen_fourier_atoms(model["shape"], frequencies)
    )
    q, triangular = qr(
        seq - v @ (v.conj().T @ seq), mode="economic", check_finite=False
    )
    if min(abs(np.diag(triangular))) < 1.0e-12:
        raise ValueError("requested Fourier basis has dependent projected columns")
    return dict(
        left=ul[:, :bright],
        sigma=s[:bright],
        v=v,
        q=q,
        gq=G @ q,
        sq=S @ q,
        radiation_scale=float(s[0]),
        setup_s=perf_counter() - start,
    )


def real_dot(x, y):
    return float(np.vdot(x, y).real)


def finite_line(current, field, delta, gdelta, data_res, sdelta, weight):
    """Stationary points of the frozen-denominator shared-real consistency quartic."""

    def real_rows(x, y):
        return np.sum(x.conj() * y, axis=1).real

    def squared(v, d):
        return np.stack((real_rows(v, v), 2 * real_rows(v, d), real_rows(d, d)))

    jf = squared(current, delta)
    ef = squared(field, gdelta)
    inner = np.stack(
        (
            real_rows(field, current),
            real_rows(field, delta) + real_rows(gdelta, current),
            real_rows(gdelta, delta),
        )
    )
    den = np.maximum(ef[0], 1.0e-28)
    coeff = np.zeros(5)
    for i in range(3):
        for j in range(3):
            coeff[i + j] += weight * np.sum((jf[i] * ef[j] - inner[i] * inner[j]) / den)
    coeff[:3] += np.array(
        (
            real_dot(data_res, data_res),
            2 * real_dot(data_res, sdelta),
            real_dot(sdelta, sdelta),
        )
    )
    deriv = np.arange(1, 5) * coeff[1:]
    if not np.max(abs(deriv)):
        return []
    nz = np.flatnonzero(abs(deriv) > 1.0e-13 * np.max(abs(deriv)))
    roots = np.polynomial.polynomial.polyroots(deriv[: nz[-1] + 1])
    return sorted(
        set(
            float(t.real)
            for t in roots
            if abs(t.imag) < 1.0e-8 * (1 + abs(t.real)) and np.isfinite(t.real)
        )
    )


def gradient_material_direction(model, m, base, residual, factor, data):
    gradient = true_gradient(
        model, m, base, residual, factor, max(norm(data) ** 2, 1.0e-30)
    )
    lo, hi = model["bounds"]
    direction = -gradient
    direction[(m <= lo) & (direction < 0)] = 0
    direction[(m >= hi) & (direction > 0)] = 0
    return direction / max(norm(direction), 1.0e-30)


def tangent_columns(model, m, field, direction, factor):
    return lu_solve(factor, direction[:, None] * field, check_finite=False)


def compatible_blend(base, physical_field, trial, trial_field, proposed):
    """Fit the interpolating current to the constitutive equation at the proposed voxel material."""
    e0 = base - proposed[:, None] * physical_field
    change = (trial - proposed[:, None] * trial_field) - e0
    size = real_dot(change, change)
    if size == 0:
        return 0.0
    return float(np.clip(-real_dot(change, e0) / size, 0.0, 1.0))


def schur_current(
    model, base, m, data, bank, arm, noise_tolerance, factor, observe=None
):
    # Auxiliary current is a measured subspace innovation plus a freely optimized FFT update.
    r = data - model["S"] @ base
    obs = bank["left"].conj().T @ r
    # The noise-free equation and the noisy discrepancy tube use the same stable bright space.
    fraction = max(0.0, 1.0 - noise_tolerance / max(norm(obs), 1.0e-30))
    j0 = base + bank["v"] @ (fraction * obs / bank["sigma"][:, None])
    field = model["E0"] + model["G"] @ j0
    defect = j0 - m[:, None] * field
    data_res = model["S"] @ j0 - data
    power = np.maximum(np.sum(abs(field) ** 2, axis=1), 1.0e-28)
    coupled = arm.startswith("shared")
    # Squared signal-space scale explicitly accounts for the units of current and data.
    lam = bank["radiation_scale"] * 0.35
    weight = lam * lam

    def proj(z):
        return z - field_projection(z, field, power) if coupled else z

    tangent = arm in (
        "shared36_tangent_soft",
        "shared36_tangent_dark",
        "shared36_tangent_exactdark",
    )
    total = bank["q"].shape[1]
    extra_budget = min(data.shape[1], total - 1)
    k = total - extra_budget if arm == "shared36_plus12" or tangent else total
    q = bank["q"][:, :k]
    gq = bank["gq"][:, :k]
    sq = bank["sq"][:, :k]
    kept = 0
    if arm == "shared36_plus12" or tangent:
        if tangent:
            direction = gradient_material_direction(model, m, base, -r, factor, data)
            if arm == "shared36_tangent_exactdark":
                # Pull the receiver-bright constraint through the inverse constitutive derivative.
                forbidden = lu_solve(factor, bank["v"], trans=2, check_finite=False)
                forbidden, _ = qr(forbidden, mode="economic", check_finite=False)
                target = direction[:, None] * field
                g = lu_solve(
                    factor,
                    target - forbidden @ (forbidden.conj().T @ target),
                    check_finite=False,
                )
            else:
                g = tangent_columns(model, m, field, direction, factor)
        else:
            h = proj(defect)
            g = model["S"].conj().T @ data_res + weight * (
                h - model["G"].conj().T @ (m[:, None] * h)
            )
        missed = g - q @ (q.conj().T @ g)
        if arm != "shared36_tangent_soft":
            missed -= bank["v"] @ (bank["v"].conj().T @ g)
        extra, values, _ = svd(missed, full_matrices=False, check_finite=False)
        kept = int(np.sum(values > max(values[0] * 1.0e-9, 1.0e-15)))
        if kept and extra_budget:
            extra = extra[:, : min(extra_budget, kept)]
            extra -= q @ (q.conj().T @ extra)
            if arm != "shared36_tangent_soft":
                extra -= bank["v"] @ (bank["v"].conj().T @ extra)
            extra, _ = qr(extra, mode="economic", check_finite=False)
            q = np.column_stack((q, extra))
            gq = np.column_stack((gq, model["G"] @ extra))
            sq = np.column_stack((sq, model["S"] @ extra))
    a = q - m[:, None] * gq
    ah = np.ascontiguousarray(a.conj().T)
    reg = weight * 0.002
    h0 = sq.conj().T @ sq + weight * (ah @ a) + reg * np.eye(q.shape[1])
    fac = cho_factor(h0, check_finite=False)
    rhs = -(sq.conj().T @ data_res + weight * ah @ proj(defect))
    if coupled:
        # Real self-adjoint per-voxel projection couples ALL input illuminations.
        def hessian(x):
            return h0 @ x - weight * (ah @ field_projection(a @ x, field, power))

        x = np.zeros_like(rhs)
        res = rhs.copy()
        pre = cho_solve(fac, res, check_finite=False)
        direction = pre.copy()
        rho = real_dot(res, pre)
        rho0 = rho
        cg = 0
        for step in range(22):
            if rho <= 1.0e-12 * max(rho0, 1.0e-30):
                break
            hd = hessian(direction)
            den = real_dot(direction, hd)
            if den <= 0:
                break
            alpha = rho / den
            x += alpha * direction
            res -= alpha * hd
            cg += 1
            pre = cho_solve(fac, res, check_finite=False)
            new = real_dot(res, pre)
            direction = pre + (new / rho) * direction
            rho = new
    else:
        x = cho_solve(fac, rhs, check_finite=False)
        cg = 1
    delta = q @ x
    gd = gq @ x
    # Check the finite, shared-real-material readout on the actual quadratic/bilinear equations.
    # The same choice of line parameter is used in every cell and in every source.
    best = None
    line_candidates = [1.0, 0.5, 0.25, 0.0]
    cubic = arm in ("shared_cubic", "shared_transport", "shared_linear", "shared_blend")
    roots = finite_line(j0, field, delta, gd, data_res, sq @ x, weight) if cubic else []
    line_candidates.extend(roots)
    transport = arm in ("shared_transport", "shared_linear", "shared_blend")
    if transport:
        map_to_receiver = lu_solve(factor, model["S"].T, trans=1, check_finite=False).T
        anchor_res = model["S"] @ base - data
        physical_field = model["E0"] + model["G"] @ base
    for t in line_candidates:
        jp = j0 + t * delta
        ep = field + t * gd
        candidate = voxel_ratio(jp, ep, m, model["bounds"])
        con = jp - candidate[:, None] * ep
        if observe is not None:
            observe(t, candidate, jp, ep, con)
        merit = norm(data_res + t * sq @ x) ** 2 + weight * norm(con) ** 2
        if transport:
            mix = (
                compatible_blend(base, physical_field, jp, ep, candidate)
                if arm == "shared_blend"
                else float(arm == "shared_transport")
            )
            candidate_field = physical_field + mix * (ep - physical_field)
            r_pred = anchor_res + map_to_receiver @ (
                (candidate - m)[:, None] * candidate_field
            )
            scale = max(norm(data) ** 2, 1.0e-30)
            merit = (
                0.5 * norm(r_pred) ** 2 / scale + material_penalty(model, candidate)[0]
            )
        if best is None or merit < best[0]:
            best = (merit, candidate, t, con, jp)
    _, proposal, chosen, con, jp = best
    ratio = norm(field_projection(defect, field, power)) / max(norm(defect), 1.0e-30)
    return proposal, dict(
        cg=cg,
        fft_rank=int(q.shape[1]),
        added=min(extra_budget, kept),
        line_t=chosen,
        innovation_factor=float(fraction),
        aux_current_step=norm(jp - base) / max(norm(base), 1.0e-30),
        dark_current_step=chosen * norm(delta) / max(norm(base), 1.0e-30),
        voxel_step=norm(proposal - m) / np.sqrt(len(m)),
        explainable_squared=ratio**2,
        aux_data_rel=norm(model["S"] @ jp - data) / max(norm(data), 1.0e-20),
        defect_rel=norm(con) / max(norm(jp), 1.0e-20),
        stationary_roots=roots,
        tangent_rhs=(
            2 * data.shape[1]
            + (bank["v"].shape[1] if arm == "shared36_tangent_exactdark" else 0)
        )
        if tangent
        else 0,
        extra_adjoint_rhs=model["S"].shape[0] if transport else 0,
    )


def run_case(model, case, bank, arm, outer, polish, budget=0, initial=None):
    public_arm = arm
    aliases = {
        "lagged_current": "lagged48",
        "shared_current": "shared48",
        "tangent_soft": "shared36_tangent_soft",
        "tangent_dark": "shared36_tangent_exactdark",
    }
    arm = aliases.get(arm, arm)
    t0, c0 = perf_counter(), process_time()
    m = (
        np.full(len(model["G"]), BACKGROUND)
        if initial is None
        else np.array(initial, copy=True)
    )
    y = case["Y"]
    scale = norm(y) ** 2
    strong_noise = case["sigma"] * np.sqrt(len(bank["sigma"]) * y.shape[1])
    full_solves = 0
    gradient_calls = 0

    def call(mat):
        nonlocal full_solves
        full_solves += 1
        return objective(model, mat, y, scale)

    val, current, residual, factor = call(m)
    trace = []
    if arm != "voxel_only":
        for k in range(outer):
            p, info = schur_current(
                model, current, m, y, bank, arm, strong_noise, factor
            )
            initial = val
            accepted = 0.0
            for z in (1.0, 0.5, 0.25, 0.125):
                if budget and full_solves >= budget:
                    break
                trial = m + z * (p - m)
                vv, jj, rr, ff = call(trial)
                if vv < initial - 1.0e-12:
                    m = trial
                    val, current, residual, factor = vv, jj, rr, ff
                    accepted = z
                    break
            trace.append(
                dict(
                    cycle=k + 1,
                    mixing=accepted,
                    value_before=initial,
                    value_after=val,
                    full_rel=norm(residual) / np.sqrt(scale),
                    **info,
                )
            )
            if accepted == 0 or norm(residual) <= case["sigma"] * np.sqrt(y.size):
                break
    prefix = m.copy()
    pre_rel = norm(residual) / np.sqrt(scale)
    pre_eval = full_solves
    if polish > 0 and (budget == 0 or full_solves < budget):
        incumbent = [val, m.copy(), norm(residual) / np.sqrt(scale)]

        def fg(x):
            nonlocal gradient_calls
            if budget and full_solves >= budget:
                raise StopIteration
            gradient_calls += 1
            value, j, r, fac = call(x)
            if value < incumbent[0]:
                incumbent[:] = [value, x.copy(), norm(r) / np.sqrt(scale)]
            multiplier = model.get("optimizer_scale", 1.0)
            return multiplier * value, multiplier * true_gradient(
                model, x, j, r, fac, scale
            )

        try:
            minimize(
                fg,
                m,
                jac=True,
                method="L-BFGS-B",
                bounds=[model["bounds"]] * m.size,
                options=dict(
                    maxiter=polish,
                    **model.get(
                        "optimizer_options",
                        dict(maxcor=12, ftol=1.0e-13, gtol=1.0e-10, maxls=16),
                    ),
                ),
            )
        except StopIteration:
            pass
        val, m, final_rel = incumbent
    else:
        final_rel = pre_rel
    return [prefix, m], dict(
        arm=public_arm,
        selected_current_rank=bank["q"].shape[1],
        native_objective="original_penalty" in model,
        case=case["id"],
        full_rel=final_rel,
        initial_norm=norm(y),
        prefix_rel=pre_rel,
        prefix_full_solves=pre_eval,
        total_full_solves=full_solves,
        forward_rhs=full_solves * y.shape[1],
        gradient_rhs=gradient_calls * y.shape[1],
        tangent_rhs=sum(x["tangent_rhs"] for x in trace),
        transport_rhs=sum(x["extra_adjoint_rhs"] for x in trace),
        cycles=len(trace),
        accepted=sum(x["mixing"] > 0 for x in trace),
        trace=trace,
        cpu_s=process_time() - c0,
        wall_s=perf_counter() - t0,
    )


def main(args):
    rows = []
    images = []
    keys = []
    tables = []
    for name in args.cohorts:
        model, cases = load_cohort(name)
        if args.native_objective:
            use_archived_objective(model)
        frequencies = (
            None
            if args.current_modes_file is None
            else np.loadtxt(args.current_modes_file, ndmin=2)
        )
        params = prepare(model, args.current_rank, args.bright_ratio, frequencies)
        tables.append(
            dict(
                cohort=name,
                bright=len(params["sigma"]),
                current_rank=params["q"].shape[1],
                selected_modes=None
                if frequencies is None
                else frequencies.astype(int).tolist(),
                bright_ratio=args.bright_ratio,
                setup_s=params["setup_s"],
            )
        )
        for case in cases[: args.limit or None]:
            for arm in args.arms:
                out, row = run_case(
                    model, case, params, arm, args.cycles, args.polish, args.full_budget
                )
                keys.append([name, case["id"], arm])
                images.append(out)
                row["cohort"] = name
                rows.append(row)
        print(
            name,
            len(cases[: args.limit or None]),
            "last",
            rows[-1]["full_rel"],
            flush=True,
        )
    np.savez_compressed(
        OUT / ("current_cycle" + args.tag + ".npz"), key=keys, images=images
    )
    (OUT / ("current_cycle" + args.tag + ".json")).write_text(
        json.dumps(dict(shared_preparation=tables, rows=rows), indent=1) + "\n"
    )


def score(args):
    from cohorts import load_truth
    from score_materials import material_metrics

    rows = []
    with np.load(OUT / ("current_cycle" + args.tag + ".npz")) as archive:
        for (name, cid, arm), imgs in zip(archive["key"], archive["images"]):
            target = load_truth(str(name), str(cid))
            for stage, img in zip(("alternating", "after_voxel_polish"), imgs):
                rows.append(
                    dict(
                        cohort=name,
                        case=cid,
                        arm=arm,
                        stage=stage,
                        **material_metrics(img, target),
                    )
                )
    with (OUT / ("current_scores" + args.tag + ".csv")).open("w") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--cohorts", nargs="+", default=DEFAULT_COHORTS)
    p.add_argument(
        "--arms",
        nargs="+",
        default=ARMS,
        choices=(
            *ARMS,
            "shared_cubic",
            "shared_linear",
            "shared_transport",
            "shared_blend",
            "shared36_tangent_soft",
            "shared36_tangent_dark",
            "shared36_tangent_exactdark",
            "lagged_current",
            "shared_current",
            "tangent_soft",
            "tangent_dark",
        ),
    )
    p.add_argument("--cycles", type=int, default=STEPS)
    p.add_argument("--polish", type=int, default=POLISH)
    p.add_argument("--limit", type=int)
    p.add_argument("--tag", default="")
    p.add_argument("--native-objective", action="store_true")
    p.add_argument("--full-budget", type=int, default=0)
    p.add_argument(
        "--current-modes-file",
        type=Path,
        help="two integer Fourier indices ky kx per line; sets the current rank",
    )
    p.add_argument("--current-rank", type=int, default=MODES)
    p.add_argument("--bright-ratio", type=float, default=0.01)
    p.add_argument("--score", action="store_true")
    args = p.parse_args()
    with threadpool_limits(1):
        if args.score:
            score(args)
        else:
            main(args)
