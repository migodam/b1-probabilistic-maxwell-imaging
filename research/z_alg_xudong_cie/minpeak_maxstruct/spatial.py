#!/usr/bin/env python3
"""One-second numerical multiplier proposals and independent native certification."""

from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.linalg import eigh
from threadpoolctl import threadpool_limits

from minpeak import HERE, TRANSFER, clock, load, read, timing, write
from core.certificate import certify
from core.inertia import Prepared, adjoint, hermitian


def inherited(proofs, case, delta):
    passes = []
    complete = 0
    for row in proofs["receipts"].values():
        score = row.get("per_case", {}).get(case)
        if row["status"] != "geometry_verified" or score is None or "joint_interval" not in score:
            continue
        complete += 1
        if score["joint_interval"][0] > delta:
            passes.append(dict(profile=row["profile"], cap=row["cap"], mu=row["mu"],
                               lam=row["lambda_"], lower=score["joint_interval"][0],
                               nu=row["proof_nu"], fixed_R10=row["historical_r10_profile"]))
    if not passes:
        raise ArithmeticError("No inherited positive cap on this archived observation.")
    best = max(passes, key=lambda row: (row["cap"], row["lower"]))
    return dict(best=best, verified_profiles_reused=complete, strictly_passed=len(passes),
                frozen_R10=next((row for row in passes if row["fixed_R10"]), None))


def block_map(centers):
    # Every 3×3 group of native cells forms one block on a 24² regular square.
    axes = [np.unique(centers[:, i], return_inverse=True) for i in range(2)]
    if any(len(coords) != 24 for coords, index in axes):
        raise ValueError("Expected 24×24 regular saved cell centers.")
    return (axes[1][1] // 3) * 8 + (axes[0][1] // 3)


def evaluate(p, offsets, cap, w, lam):
    return p.evaluate_ldl(cap, w[offsets], lam[offsets], return_solution=True)


def summarize(row):
    return {key: value for key, value in row.items()
            if key not in ("solution_current", "linear_term")}


def key(row):
    lower = row.get("lower_bound_numeric") or 0.
    eig = row.get("H_eigenvalues_descending")
    signed = float(eig[row["nu"]]) if eig else -1.e30
    nu = row.get("nu")
    return lower, signed, -nu if nu is not None else -1.e30


def gradient(p, offsets, cap, row):
    X = row.get("solution_current")
    if X is None:
        return None
    L = row["linear_term"]
    H = hermitian(p.Y2 - adjoint(L) @ X)
    values, views = eigh(H, check_finite=False)
    values, views = values[::-1], views[:, ::-1]
    keep = np.flatnonzero((np.arange(p.s) >= row["nu"]) & (values > 0))
    if not len(keep):
        keep = [row["nu"]]  # Optimize its signed value until a positive mode survives.
    projected = X @ views[:, keep]
    field = (p.B + p.G @ X) @ views[:, keep]
    interaction = (projected.conj() * field).sum(axis=1)
    mu_gradient = np.bincount(offsets, (abs(projected)**2).sum(axis=1) -
                              cap * interaction.real, minlength=64)
    lam_gradient = np.bincount(offsets, interaction.imag, minlength=64)
    p.operations["source_weight_gradient_eigenproblems"] += 1
    p.operations["gradient_green_rhs"] += p.s
    return np.sign(mu_gradient), np.sign(lam_gradient)


def seeds(grid, index, cap, historical):
    available = []
    nearest = sorted(set(row["cap"] for row in grid),
                     key=lambda value: (abs(value - cap), -value))[:2]
    for row in grid:
        if row["cap"] not in nearest or not row["case_matrices"]:
            continue
        nu = row["nu"]
        if nu is None or nu >= 12:
            continue
        eig = np.asarray(row["case_matrices"][index]["eig"])
        numeric = float(np.sqrt(np.maximum(0., eig[nu:]).sum()))
        available.append((numeric, eig[nu], row["mu"], row["lambda_over_mu"]))
    take = [(historical["mu"], historical["lam"]/historical["mu"])]
    for _, _, mu, ratio in sorted(available, reverse=True):
        if (mu, ratio) not in take:
            take.append((mu, ratio))
        if len(take) >= 4:
            break
    return take


def explore(p, offsets, delta, plans, budget=1.):
    started = clock()
    start_time = perf_counter()
    deadline = start_time + budget
    records = []
    rows = []
    estimated_cost = .018

    def fits():
        return perf_counter() + estimated_cost * 1.35 < deadline

    def evaluate_costed(cap, w, l, stage):
        nonlocal estimated_cost
        before = perf_counter()
        result = evaluate(p, offsets, cap, w, l)
        measured = perf_counter() - before
        estimated_cost = max(.004, .7 * estimated_cost + .3 * measured)
        item = dict(cap=cap, stage=stage, sector=w.tolist(), optical=l.tolist(),
                    **summarize(result))
        rows.append(item)
        if key(result)[0] > delta * (1+1.e-7):
            records.append(item)
        return result

    for cap, proposals in plans:
        if not fits():
            break
        baseline = [(w, l, evaluate_costed(cap, w, l, "mixed_global_seed"))
                    for mu, ratio in proposals if fits()
                    for w, l in [(np.full(64, mu), np.full(64, mu*ratio))]]
        if not baseline:
            break
        w, l, current = max(baseline, key=lambda triple: key(triple[2]))
        if current.get("nu") is None or current["nu"] >= p.s:
            continue
        stage_started = perf_counter()
        stage_limit = min(deadline, stage_started + .20)
        nominal = float(np.mean(w))
        for repeat in range(4):
            if not fits() or perf_counter() + estimated_cost > stage_limit:
                break
            step = gradient(p, offsets, cap, current)
            if step is None:
                break
            dw, dl = step
            candidates = []
            for fraction in (.1, .35, .7, .025):
                if not fits() or perf_counter() + estimated_cost > stage_limit:
                    break
                w_next = np.maximum(.05 * nominal, w + fraction * nominal * dw)
                l_next = l + fraction * max(nominal, abs(l).max()) * dl
                trial = evaluate_costed(cap, w_next, l_next, "block_64_gradient")
                candidates.append((w_next, l_next, trial))
            best = max(candidates, key=lambda triple: key(triple[2]), default=(w, l, current))
            if key(best[2]) <= key(current):
                break
            w, l, current = best
    return records, rows, timing(started)


def main():
    batch = clock()
    entries, G, S, B, Y, delta, centers, _ = load(TRANSFER / "public_inputs.npz")
    archive_proofs = read(TRANSFER / "verified.json")
    archive_proposals = read(TRANSFER / "numerical_grid.json")
    grid = archive_proposals["values"]
    offsets = block_map(centers)
    if archive_proposals["cases"] != [entry["case"] for entry in entries]:
        raise ValueError("Historical ranking cannot be matched to the saved observations.")
    results = dict(scope="Matched archived 24² physical box [0, cap] and original CN99 radii",
                   online_wall_target=1., block_partition=[8, 8],
                   setup=timing(batch), cases={})
    for i, entry in enumerate(entries):
        case = entry["case"]
        t0 = clock()
        old = inherited(archive_proofs, case, delta[i])
        floor = old["best"]["cap"]
        opts = [floor + step for step in (.25, .5, .75, 1., 1.5, 2.)]
        plans = [(cap, seeds(grid, i, cap, old["best"])) for cap in opts]
        p = Prepared(G, S, B, Y[i])
        preparation = timing(t0)
        proposals, traces, cost = explore(p, offsets, float(delta[i]), plans)
        ranked = sorted(proposals, key=lambda row: (-row["cap"], -row["lower_bound_numeric"]))
        audited = []
        selected = None
        for item in ranked[:3]:
            w = np.asarray(item["sector"])[offsets]
            l = np.asarray(item["optical"])[offsets]
            verified = certify(p, item["cap"], w, l, radii={"original_CN99": delta[i]},
                               precision=128, threads=1)
            audited.append(dict(proposal=item, proof=verified))
            if verified.get("strictly_excludes", {}).get("original_CN99"):
                selected = audited[-1]
                break
        final_cap = selected["proposal"]["cap"] if selected else floor
        print(case, "published", floor, "final", final_cap, "elapsed", cost["wall_seconds"],
              "attempts", len(traces), flush=True)
        results["cases"][case] = dict(original_radius=float(delta[i]), historical=old,
                                       final_cap=final_cap, selected=selected, verified_proposals=audited,
                                       proposed_caps=opts, all_numeric_trials=traces,
                                       preparation=preparation, online_multiplier_search=cost,
                                       proposal_operations=dict(p.operations), **timing(t0))
        write(HERE / "lower_certificates.json", {**results, **timing(batch)})


if __name__ == "__main__":
    with threadpool_limits(1):
        main()
