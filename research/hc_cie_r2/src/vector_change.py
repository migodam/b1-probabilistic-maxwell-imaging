"""Paired two-voxel 3D monitoring using original SRAP dyadic Green actions."""
import argparse
import csv
import json
import os
from pathlib import Path
import time

import spectra as sp  # native linear algebra stays on one thread
import numpy as np
from scipy.linalg import svd
from scipy.optimize import least_squares

from spectra_inputs import _load_module


HERE = Path(__file__).resolve().parents[1]
CONFIG_ENV = "HC_CIE_VECTOR_ONLINE_CONFIG"
DESIGN = json.loads(Path(os.environ.get(CONFIG_ENV, HERE / "configs/vector_change.json")).read_text())
EPS = np.finfo(float).eps


def local_columns(support):
    return (3 * np.asarray(support)[:, None] + np.arange(3)[None, :]).ravel()


def prepare_change(parent, background, n):
    rng = np.random.default_rng(np.random.SeedSequence([DESIGN["generation_seed"], parent]))
    grid = np.indices((n,) * 3).reshape(3, -1).T
    central = np.logical_and(grid >= 2, grid < n - 2).all(axis=1)
    highest = np.flatnonzero((background >= np.quantile(background, 0.97)) & central)
    remove = int(rng.choice(highest))
    add = np.flatnonzero(central & (background <= 0.25) & (np.linalg.norm(grid - grid[remove], axis=1) >= 4))
    insert = int(rng.choice(add))
    material = background.copy()
    material[remove] *= 0.55
    material[insert] = min(4.0, material[insert] + 1.1)
    return material, tuple(sorted((remove, insert)))


def forward(physics, model, material):
    m = physics.t(material).real
    current, error = model.solve(m, m.repeat_interleave(3)[:, None] * model.E,
                                 tol=DESIGN["state_tolerance"])
    return physics.npy(model.S @ current), current, error


def get_cache(physics, model, background, cache, support):
    """Solve only selected unit-source RHS and reuse columns within a cache."""
    missing = [i for i in support if i not in cache]
    if missing:
        rows = local_columns(missing)
        rhs = np.zeros((model.M, len(rows)), dtype=complex)
        rhs[rows, np.arange(len(rows))] = 1
        state, error = model.solve(physics.t(background).real, physics.t(rhs), tol=DESIGN["state_tolerance"])
        propagated = physics.npy(model.G(state))
        states = physics.npy(state)
        for i, index in enumerate(missing):
            sel = slice(3 * i, 3 * (i + 1))
            cache[index] = (states[:, sel], propagated[:, sel], float(error))
    columns = np.concatenate([cache[i][0] for i in support], axis=1)
    block = np.concatenate([cache[i][1][local_columns(support)] for i in support], axis=1)
    return columns, block


def build_dictionary(physics, model, background, reference_current):
    receiver_adjoint, error = model.solve(physics.t(background).real, model.S.mH, adj=True,
                                         tol=DESIGN["state_tolerance"])
    observation_map = physics.npy(receiver_adjoint).conj().T
    background_field = physics.npy(model.E + model.G(reference_current))
    receiver_blocks = observation_map.T.reshape(model.N, 3, -1)
    interior_blocks = background_field.reshape(model.N, 3, -1)
    sensitivities = np.einsum("ncr,ncs->nrs", receiver_blocks, interior_blocks, optimize=True).reshape(model.N, -1)
    real = np.concatenate([sensitivities.real, sensitivities.imag], axis=1)
    return observation_map, background_field, np.ascontiguousarray(real), float(error)


def _real_vector(data):
    return np.r_[data.real.ravel(), data.imag.ravel()]


def bounded_pair_candidates(dictionary, observation, lower, upper):
    """Screen a whole-grid real Born dictionary with box-constrained pair fits."""
    y = _real_vector(observation)
    norm2 = np.einsum("ij,ij->i", dictionary, dictionary)
    dot = dictionary @ y
    single = np.clip(dot / norm2, lower, upper)
    score = 2 * dot * single - norm2 * single**2
    seeds = np.argsort(-score, kind="stable")[:DESIGN["beam_seeds"]]
    candidates = {}
    indices = np.arange(dictionary.shape[0])
    for i in seeds:
        h = dictionary @ dictionary[i]
        det = norm2[i] * norm2 - h**2
        safe = det > 32 * EPS * norm2[i] * norm2
        left = np.divide(norm2 * dot[i] - h * dot, det, where=safe, out=np.zeros_like(dot))
        right = np.divide(norm2[i] * dot - h * dot[i], det, where=safe, out=np.zeros_like(dot))
        interior = safe & (left >= lower[i]) & (left <= upper[i]) & (right >= lower) & (right <= upper)
        lx = [left]
        rx = [right]
        for first in (lower[i], upper[i]):
            lx.append(np.full_like(dot, first))
            rx.append(np.clip((dot - h * first) / norm2, lower, upper))
        for second in (lower, upper):
            lx.append(np.clip((dot[i] - h * second) / norm2[i], lower[i], upper[i]))
            rx.append(second)
        lx = np.asarray(lx)
        rx = np.asarray(rx)
        evaluations = 2 * (dot[i] * lx + dot * rx) - (norm2[i] * lx**2 + norm2 * rx**2 + 2 * h * lx * rx)
        evaluations[0, ~interior] = -np.inf
        evaluations[:, i] = -np.inf
        best = np.argmax(evaluations, axis=0)
        best_scores = evaluations[best, indices]
        for j in np.argsort(-best_scores, kind="stable")[:DESIGN["beam_extensions"]]:
            pair = (int(i), int(j))
            values = (float(lx[best[j], j]), float(rx[best[j], j]))
            if i > j:
                pair, values = pair[::-1], values[::-1]
            if pair not in candidates or candidates[pair][0] < best_scores[j]:
                candidates[pair] = (float(best_scores[j]), values)
    return sorted(candidates.items(), key=lambda line: (-line[1][0], line[0]))[:DESIGN["local_finalists"]]


def closed_readout(observation, o, u, k, lower, upper):
    left, singular, vh = svd(o, full_matrices=False, check_finite=False)
    keep = singular > DESIGN["relative_svd_floor"] * singular[0]
    v = (vh[keep].conj().T / singular[keep]) @ (left[:, keep].conj().T @ observation)
    total = (u + k @ v).reshape(len(lower), 3, -1)
    local = v.reshape(total.shape)
    field_energy = np.sum(abs(total)**2, axis=(1, 2))
    readout = np.sum((total.conj() * local).real, axis=(1, 2)) / field_energy
    clipped = np.clip(readout, lower, upper)
    return clipped, dict(singular_values=singular.tolist(), rank=int(np.count_nonzero(keep)),
                         raw_readout=readout.tolist(), clipped_readout=clipped.tolist(),
                         port_projection_relative=float(np.linalg.norm(o @ v - observation) / np.linalg.norm(observation)),
                         current_closure_relative=float(np.linalg.norm(local - clipped[:, None, None] * total) /
                                                        max(np.linalg.norm(local), 1e-300)),
                         field_energy=field_energy.tolist())


def local_prediction(o, u, k, d):
    diagonal = np.repeat(np.asarray(d), 3)[:, None]
    return o @ np.linalg.solve(np.eye(len(k)) - diagonal * k, diagonal * u)


def bounded_refit(o, u, k, y, born, closed, lower, upper):
    magnitude = np.linalg.norm(y)
    evaluations = 0

    def residual(d):
        nonlocal evaluations
        evaluations += 1
        return _real_vector(local_prediction(o, u, k, d) - y) / magnitude

    initial = min((born, closed), key=lambda d: np.linalg.norm(residual(d)))
    fit = least_squares(residual, initial, bounds=(lower, upper), max_nfev=DESIGN["bounded_local_maxfun"],
                        ftol=1e-10, xtol=1e-10, gtol=1e-10)
    return fit.x, float(np.linalg.norm(fit.fun) * magnitude), dict(local_evaluations=evaluations, stop=int(fit.status))


def independently_audit(physics, model, background, support, change, y0, observation):
    material = background.copy()
    material[list(support)] += change
    pred, _, state_error = forward(physics, model, material)
    differential = pred - y0
    return differential, float(np.linalg.norm(differential - observation)), state_error


def score_change(support, change, true_support, true_change, clean, reconstructed, observation, radius):
    recovered = dict(zip(support, change, strict=True))
    correct = dict(zip(true_support, true_change, strict=True))
    union = set((*recovered, *correct))
    error = np.linalg.norm([recovered.get(i, 0) - correct.get(i, 0) for i in union]) / np.linalg.norm(true_change)
    active = {i for i in recovered if abs(recovered[i]) > 1e-6}
    return dict(support_exact=active == set(true_support), exact_search_indices=set(support) == set(true_support),
                recovered_sites=len(active), support_overlap=len(active.intersection(true_support)),
                material_rel_l2=float(error),
                measured_residual_over_clean_delta=float(np.linalg.norm(reconstructed - observation) / np.linalg.norm(clean)),
                measured_residual_over_noise_radius=(float(np.linalg.norm(reconstructed - observation) / radius) if radius else None),
                clean_port_rel_l2=float(np.linalg.norm(reconstructed - clean) / np.linalg.norm(clean)))


def search_dictionary(physics, model, background, reference, measured, observation_map, u, dictionary):
    """Search unknown support using only the known background and observed differential."""
    lower, upper = -background, 4.0 - background
    shortlist = bounded_pair_candidates(dictionary, measured, lower, upper)
    if len(shortlist) < DESIGN["full_audit_finalists"]:
        raise ValueError("The whole-grid search produced too few candidate pairs")
    finalists = []
    search_cache = {}
    for support, (score, born) in shortlist:
        indices = list(support)
        local = local_columns(indices)
        _, k_pair = get_cache(physics, model, background, search_cache, indices)
        o_pair, u_pair = observation_map[:, local], u[local]
        closed, report = closed_readout(measured, o_pair, u_pair, k_pair, lower[indices], upper[indices])
        fitted, fitted_residual, small = bounded_refit(o_pair, u_pair, k_pair, measured, born, closed,
                                                        lower[indices], upper[indices])
        finalists.append(dict(support=list(support), change=fitted.tolist(), born=list(born), born_score=score,
                              closed=report, local_fitted_residual=fitted_residual, **small))
    audit_pool = sorted(finalists, key=lambda c: c["local_fitted_residual"])[:DESIGN["full_audit_finalists"]]
    inspected = []
    for item in audit_pool:
        differential, loss, state_error = independently_audit(physics, model, background, item["support"],
                                                               item["change"], reference, measured)
        item.update(independent_absolute_audit_loss=loss, independent_state_relative_error=state_error)
        inspected.append((loss, item, differential, state_error))
    loss, chosen, reconstructed, state_error = min(inspected, key=lambda q: q[0])
    return chosen, reconstructed, finalists, audit_pool, float(loss), state_error


def evaluate_case(physics, model, background, target, sites, parent, array):
    start, cpu = time.perf_counter(), time.process_time()
    counters = model.count.copy()
    reference, current, error_ref = forward(physics, model, background)
    perturbed, _, error_target = forward(physics, model, target)
    clean = perturbed - reference
    observation_map, u, dictionary, adjoint_error = build_dictionary(physics, model, background, current)
    receiver = physics.npy(model.S)
    true_change = target[list(sites)] - background[list(sites)]
    rng = np.random.default_rng(np.random.SeedSequence([DESIGN["noise_seed"], parent]))
    direction = rng.standard_normal(clean.shape) + 1j * rng.standard_normal(clean.shape)
    direction /= np.linalg.norm(direction)
    noise_radii = {"clean": 0.0, "diff10": 0.1 * np.linalg.norm(clean), "full10": 0.1 * np.linalg.norm(perturbed)}
    geometry_counts = dict(model.count - counters)
    setup = dict(wall_s=time.perf_counter() - start, cpu_s=time.process_time() - cpu,
                 counters=geometry_counts, reference_equation_residual=error_ref,
                 target_equation_residual=error_target, adjoint_equation_residual=adjoint_error,
                 all_voxel_count=model.N)
    rows, detailed = [], []
    saved = dict(background=background, target=target, reference_y=reference, target_y=perturbed,
                 true_sites=np.asarray(sites), true_change=true_change, clean_differential=clean,
                 clean_port_singular_values=svd(clean, compute_uv=False, check_finite=False), noise_direction=direction)
    oracle_cache = {}
    true_columns = local_columns(sites)
    lower, upper = -background, 4.0 - background
    for condition in DESIGN["conditions"]:
        clock, c0 = time.perf_counter(), time.process_time()
        offset = model.count.copy()
        radius = float(noise_radii[condition])
        measured = clean + radius * direction
        oracle_state, k = get_cache(physics, model, background, oracle_cache, sites)
        o, local_u = observation_map[:, true_columns], u[true_columns]
        equality_error = np.linalg.norm(receiver @ oracle_state - o) / np.linalg.norm(o)
        true_local_prediction = local_prediction(o, local_u, k, true_change)
        exact_difference_error = np.linalg.norm(true_local_prediction - clean) / np.linalg.norm(clean)
        oracle, readout = closed_readout(measured, o, local_u, k, lower[list(sites)], upper[list(sites)])
        local_oracle = local_prediction(o, local_u, k, oracle)
        audit, physical_loss, state_error = independently_audit(physics, model, background, sites, oracle, reference, measured)
        detail = dict(parent=parent, array=array, condition=condition,
                      differential_truth_to_full=float(np.linalg.norm(clean) / np.linalg.norm(perturbed)),
                      noise_over_true_differential=radius / float(np.linalg.norm(clean)),
                      truth_support=list(sites), truth_change=true_change.tolist())
        oracle_record = dict(**detail, method="oracle_support_closed", proposed_sites=list(sites),
                             proposed_changes=oracle.tolist(), shared_readout=readout,
                             direct_port_map_relative_error=float(equality_error),
                             true_two_voxel_resolvent_relative_error=float(exact_difference_error),
                             final_state_relative_error=state_error,
                             local_full_audit_relative_error=float(np.linalg.norm(local_oracle - audit) / np.linalg.norm(clean)),
                             wall_s=time.perf_counter() - clock, cpu_s=time.process_time() - c0,
                             counters=dict(model.count - offset),
                             **score_change(sites, oracle, sites, true_change, clean, audit, measured, radius))
        rows.append(oracle_record)
        saved.update({f"{condition}_observation": measured, f"{condition}_oracle_prediction": audit,
                      f"{condition}_oracle_change": oracle, "O_I": o, "U_I": local_u, "K_I": k})

        clock, c0 = time.perf_counter(), time.process_time()
        offset = model.count.copy()
        chosen, reconstructed, finalists, audit_pool, loss, state_error = search_dictionary(
            physics, model, background, reference, measured, observation_map, u, dictionary)
        rows.append(dict(**detail, method="unknown_support_dictionary", proposed_sites=chosen["support"],
                         proposed_changes=chosen["change"], final_state_relative_error=state_error,
                         audit_vs_exact_local_error_relative=float(abs(loss - chosen["local_fitted_residual"]) /
                                                                  np.linalg.norm(clean)),
                         true_pair_in_shortlist=any(set(c["support"]) == set(sites) for c in finalists),
                         true_pair_in_full_audit=any(set(c["support"]) == set(sites) for c in audit_pool),
                         shortlisted_pairs=len(finalists), full_audits=len(audit_pool),
                         wall_s=time.perf_counter() - clock, cpu_s=time.process_time() - c0,
                         counters=dict(model.count - offset),
                         **score_change(chosen["support"], chosen["change"], sites, true_change,
                                        clean, reconstructed, measured, radius)))
        detailed.append(dict(parent=parent, array=array, condition=condition, candidates=finalists))
        saved.update({f"{condition}_unknown_prediction": reconstructed,
                      f"{condition}_unknown_sites": np.asarray(chosen["support"]),
                      f"{condition}_unknown_change": np.asarray(chosen["change"])})
    return rows, saved, detailed, setup


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--workspace", type=Path, default=Path.cwd())
    cli.add_argument("--output", type=Path, default=HERE / "results/vector_change")
    args = cli.parse_args()
    root = args.workspace.resolve()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    srap = root / "research/srap_cie_array_study_20261121/vector"
    os.environ.setdefault("SRAP_DEVICE", "cpu")
    physics = _load_module("hc_cie_vector_change_physics", srap / "physics.py")
    physics.CD = physics.torch.complex128
    all_rows, all_details, all_setups = [], [], []
    start, cpu = time.perf_counter(), time.process_time()
    models = {array: physics.Model(12, array) for array in DESIGN["arrays"]}
    startup = time.perf_counter() - start
    for parent in DESIGN["parents"]:
        history = srap / f"predictions/hemi_{parent:02d}_Dual/materials.npz"
        with np.load(history, allow_pickle=False) as payload:
            background = payload["light_180"].astype(float)
        target, sites = prepare_change(parent, background, 12)
        for array, model in models.items():
            rows, saved, detailed, setup = evaluate_case(physics, model, background, target, sites, parent, array)
            all_rows.extend(rows)
            all_details.extend(detailed)
            all_setups.append(dict(parent=parent, array=array, known_background_source=str(history.relative_to(root)), **setup))
            np.savez_compressed(output / f"{array}_{parent:02d}.npz", **saved)
            print(array, parent, [(r["condition"], r["support_overlap"], round(r["material_rel_l2"], 4))
                                  for r in rows if r["method"] == "unknown_support_dictionary"], flush=True)
    receipt = dict(rows=all_rows, setups=all_setups, candidate_audits=all_details, config=DESIGN,
                   original_operator=str((srap / "physics.py").relative_to(root)),
                   wall_s=time.perf_counter() - start, cpu_s=time.process_time() - cpu,
                   model_startup_wall_s=startup, native_precision="CPU complex128 single-thread")
    (output / "report.json").write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n")
    fields = list(dict.fromkeys(k for r in all_rows for k in r if not isinstance(r[k], dict)))
    with (output / "scores.csv").open("w", newline="", encoding="utf8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(all_rows)


if __name__ == "__main__":
    main()
