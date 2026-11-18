"""Finite weak-material alternatives and source selection on archived GN–NN R4 TEST observations."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from time import perf_counter

import numpy as np
from scipy.linalg import cho_factor, cho_solve, eigh, svd
from scipy.ndimage import gaussian_filter
from scipy.special import hankel1
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parents[1]
DEFAULT_ARCHIVE = HERE.parents[3] / "z_alg_r2_publication" / "repository"
PARENT_IDS = (
    "test_p0000", "test_p0034", "test_p0058",
    "test_p0064", "test_p0066", "test_p0072",
    "test_p0087", "test_p0097", "test_p0107",
)
THRESHOLD = .0011
LOWER, UPPER = .02, 20.
R4 = Path("research/research_gn_nn_preservation_r4")


def read_npz(path):
    with np.load(path, allow_pickle=False) as z:
        return {key: z[key] for key in z.files}


def rms(x):
    return float(np.sqrt(np.mean(x**2)))


def ratio(x, scale):
    return float(np.linalg.norm(x) / scale)


def new_sources(ops, radius, offset):
    angles = (np.arange(12) + offset) * (np.pi / 6)
    xy = radius * np.stack((np.cos(angles), np.sin(angles)), axis=1)
    distances = np.linalg.norm(ops["centers"][:, None] - xy[None], axis=2)
    E0 = .25j * hankel1(0, float(ops["k"]) * distances)
    return xy, E0


def weak_candidates(Model, D, ops, Y, initial, amplitude):
    start = perf_counter()
    model = Model(ops)
    n = initial.size
    scale = np.linalg.norm(Y)
    base = model.forward(initial, jac=True)
    J = base["J"] / scale
    _, s, vh = svd(J, full_matrices=False, check_finite=False)
    weak = vh[s <= 2e-4 * s[0]].T
    # Ask for the smoothest component among Jacobian-weak material directions.
    # Break the common rotational near-degeneracy with a fixed spatial envelope,
    # using only the existing geometry and no new-source information.
    centers = ops["centers"]
    envelope = 1 + .25 * centers[:, 0] + .15 * centers[:, 1]
    DW = D @ weak
    roughness = DW.T @ DW
    roughness += 2e-5 * (weak.T @ (envelope[:, None] * weak))
    _, smooth_mix = eigh(roughness, subset_by_index=(0, 0), check_finite=False)
    direction = (weak @ smooth_mix[:, 0]).reshape(-1)
    direction *= np.sign(direction[np.argmax(np.abs(direction))])
    direction /= rms(direction)
    bright = vh[:np.count_nonzero(s > 1e-3 * s[0])].T

    def candidate(sign, magnitude):
        m = np.clip(initial + sign * magnitude * direction, LOWER, UPPER)
        state = model.forward(m, jac=True)
        res = ratio(state["y"] - Y, scale)
        # Correct only in the old bright directions so a small residual
        # cannot be obtained by stepping back along the tested weak direction.
        for _ in range(3):
            if res <= .96 * THRESHOLD:
                break
            Jb = state["J"] @ bright / scale
            residual = Model.pack(state["y"] - Y) / scale
            H = Jb.T @ Jb
            step = bright @ cho_solve(cho_factor(H + 1e-11 * np.eye(H.shape[0])),
                                      -Jb.T @ residual)
            improved = False
            for fraction in (1., .5, .25):
                trial = model.forward(np.clip(state["m"] + fraction * step, LOWER, UPPER))
                new_res = ratio(trial["y"] - Y, scale)
                if new_res < res:
                    res, state, improved = new_res, trial, True
                    break
            if not improved:
                break
            model.jacobian(state)
        return state, res

    trials = []
    for mag in (amplitude, .75 * amplitude, .5 * amplitude, .25 * amplitude, .125 * amplitude):
        plus, r_plus = candidate(1, mag)
        minus, r_minus = candidate(-1, mag)
        trials.append([mag, r_plus, r_minus])
        if max(r_plus, r_minus) <= THRESHOLD:
            break
    return dict(
        plus=plus, minus=minus, base=base, singular_values=s,
        weak_basis=weak, direction=direction, trials=trials,
        residuals=[ratio(base["y"]-Y, scale), r_plus, r_minus],
        separated=bool(max(r_plus, r_minus) <= THRESHOLD),
        distance_rmse=rms(plus["m"] - minus["m"]),
        distance_linf=float(np.max(np.abs(plus["m"] - minus["m"]))),
        weak_direction_relative_gain=ratio(J @ direction, np.linalg.norm(direction) * s[0]),
        weak_dimension=weak.shape[1], bright_dimension=bright.shape[1],
        counts=dict(model.counts), seconds=perf_counter()-start,
    )


def source_analysis(Model, ops, cand, lib, noise_per_source, budget, rng):
    start = perf_counter()
    E0 = np.column_stack([row[2] for row in lib])
    model = Model(dict(ops, E0=E0))
    up = model.forward(cand["plus"]["m"])
    down = model.forward(cand["minus"]["m"])
    stage_cost = {"finite_predictions": dict(counts=dict(model.counts), seconds=perf_counter()-start)}
    delta = up["y"] - down["y"]
    old_delta = cand["plus"]["y"] - cand["minus"]["y"]
    weights = np.linalg.lstsq(ops["E0"], E0, rcond=None)[0]
    field_perp = E0 - ops["E0"] @ weights
    combo = old_delta @ weights
    # Validate the remainder by solving the VIE for the actual projected field
    # on both finite materials. For fixed m the VIE is linear in E0.
    tick = perf_counter()
    perp_model = Model(dict(ops, E0=field_perp))
    out = perp_model.forward(cand["plus"]["m"])["y"] - perp_model.forward(cand["minus"]["m"])["y"]
    closure = ratio(delta - combo - out, np.linalg.norm(delta))
    stage_cost["projected_field_diagnostic"] = dict(counts=dict(perp_model.counts), seconds=perf_counter()-tick)
    tick, old_counts = perf_counter(), model.counts.copy()
    probe = model.forward(cand["base"]["m"], jac=True)
    weak = cand["weak_basis"]
    material_kernel = probe["J"].reshape(len(lib), 2 * model.nr, model.n)
    sensitivity = np.linalg.norm(material_kernel @ weak, axis=(1, 2)) / noise_per_source
    stage_cost["weak_local_control"] = dict(counts=dict(model.counts-old_counts), seconds=perf_counter()-tick)
    gap = np.linalg.norm(delta, axis=0) / noise_per_source

    def top(scores, ids):
        return [int(i) for i in sorted(ids, key=lambda i: (-scores[i], i))[:budget]]

    selections = {}
    for name, lo in (("same_radius", 0), ("near_field", 12)):
        ids = np.arange(lo, lo+12)
        clearances = [min(np.linalg.norm(ops["centers"] - lib[i][1], axis=1)) for i in ids]
        closest = np.argsort(np.round(np.asarray(clearances), 12), kind="stable")[:budget]
        selections[name] = dict(
            candidate_gap=top(gap, ids),
            weak_local=top(sensitivity, ids),
            geometric=[lo + k for k in (0, 3, 6, 9)][:budget],
            geometric_closest_box=[lo+int(k) for k in closest],
            random=sorted(rng.choice(ids, size=budget, replace=False).tolist()),
        )
    selections["repeat"] = {"equally_spaced": [24+k for k in (0, 3, 6, 9)][:budget]}
    decomposition = []
    incident_norms = np.linalg.norm(E0, axis=0)
    field_match_scale = np.ones(len(lib))
    field_match_scale[12:24] = incident_norms[:12] / incident_norms[12:24]
    for i, (name, xy, _) in enumerate(lib):
        scale = field_match_scale[i]
        decomposition.append(dict(
            i=i, library=name, x=float(xy[0]), y=float(xy[1]),
            closest_disk_clearance=float(np.linalg.norm(ops["centers"]-xy, axis=1).min()
                                         - ops["radius"]),
            incident_field_norm=incident_norms[i],
            field_perp_fraction=ratio(field_perp[:, i], incident_norms[i]),
            old_combo_over_noise=ratio(combo[:, i], noise_per_source),
            perpendicular_over_noise=ratio(out[:, i], noise_per_source),
            total_over_noise=gap[i], local_sensitivity=sensitivity[i],
            field_match_scale=scale, field_match_emitted_energy=scale**2,
            field_match_old_combo_over_noise=scale*ratio(combo[:, i], noise_per_source),
            field_match_perpendicular_over_noise=scale*ratio(out[:, i], noise_per_source),
            field_match_total_over_noise=scale*gap[i],
        ))
    counts = model.counts + perp_model.counts
    return dict(E0=E0, selections=selections, decomposition=decomposition,
                closure=closure, counts=dict(counts), seconds=perf_counter()-start, stage_cost=stage_cost,
                delta=delta, combo=combo, perpendicular=out, field_match_scale=field_match_scale)


def optimize(Model, reg, ops, Y, noise_norm, initial, stages, smooth, label):
    start = perf_counter()
    model = Model(ops)
    scale = np.linalg.norm(Y)
    target = 1.10 * noise_norm / scale
    history = []
    m = initial.copy()
    D = reg.D
    for B, iterations in stages:
        x = m.copy() if B is None else np.linalg.lstsq(B, m, rcond=None)[0]
        x = np.clip(x, LOWER, UPPER)
        expand = (lambda u: u) if B is None else (lambda u: B @ u)
        DB = D if B is None else D @ B
        state = model.forward(expand(x), label=label)
        def objective(z):
            data = ratio(z["y"] - Y, scale)
            r = reg.evaluate(z["m"])[0] if smooth else 0.
            return .5*data**2 + r, data
        value, data = objective(state)
        damping = 3e-9 if smooth else 1e-8
        for _ in range(iterations):
            history.append(data)
            if data <= target:
                break
            J = model.jacobian(state) / scale
            if B is not None:
                J = J @ B
            r = Model.pack(state["y"] - Y) / scale
            g = J.T @ r
            H = J.T @ J
            if smooth:
                _, rg, rw = reg.evaluate(state["m"])
                g += rg if B is None else B.T @ rg
                priorH = (D.T @ D.multiply(rw[:, None])).toarray() if B is None else DB.T @ (rw[:, None] * DB)
                H += priorH
            accepted = False
            for _ in range(4):
                step = cho_solve(cho_factor(H + damping*np.eye(x.size), check_finite=False),
                                 -g, check_finite=False)
                for fraction in (1., .5, .25, .125):
                    proposal = np.clip(x + fraction * step, LOWER, UPPER)
                    trial = model.forward(expand(proposal), label=label)
                    trial_value, trial_data = objective(trial)
                    if trial_value < value:
                        x, state, value, data = proposal, trial, trial_value, trial_data
                        accepted = True
                        damping = max(1e-12, damping * .35)
                        break
                if accepted:
                    break
                damping *= 10
            if not accepted:
                break
        m = state["m"]
    return dict(material=m, relative_residual=data, discrepancy=target,
                noise_multiple=data * scale / noise_norm, history=history,
                discrepancy_reached=bool(data <= target),
                seconds=perf_counter()-start, counts=dict(model.counts))


def archived_material_score(m, truth):
    from research.research_gn_nn_completion_r3.src.score import material_metrics

    result = material_metrics(m, truth["m"], truth["amplitude"], truth["alpha"])
    return dict(rmse=result["RMSE"], nrmse=result["NRMSE"], iou=result["IoU_Otsu"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, default=DEFAULT_ARCHIVE)
    parser.add_argument("--cases", nargs="+", default=PARENT_IDS)
    parser.add_argument("--budget", type=int, default=4, choices=range(1, 5))
    parser.add_argument("--output", type=Path, default=HERE/"results"/"weak_r4")
    parser.add_argument("--candidates-only", action="store_true")
    args = parser.parse_args()
    sys.path.insert(0, str(args.archive.resolve()))
    from research.research_z_alg_r4_diversity.src.physics import Model, Regularizer, material_basis

    r4 = args.archive / R4
    ops = read_npz(r4/"data/operators.npz")
    records = json.loads((r4/"data/manifests/evaluation.json").read_text())["records"]
    records = {row["parent_id"]: row for row in records if row["split"] == "test" and row["noise"] == .001}
    reg = Regularizer((24, 24), smooth=1e-5, tv=1e-5, epsilon=.01)
    assert args.cases and all(i in records for i in args.cases)
    near_r = .72
    assert near_r > np.linalg.norm(ops["centers"],axis=1).max() + float(ops["radius"])
    lib = []
    old_radius = np.linalg.norm(ops["sources"][0])
    for name, radius in (("same_radius", old_radius), ("near_field", near_r)):
        xy, E0 = new_sources(ops, radius, .5)
        lib.extend((name, xy[i], E0[:, i]) for i in range(E0.shape[1]))
    lib.extend(("repeat", ops["sources"][i], ops["E0"][:, i]) for i in range(ops["E0"].shape[1]))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result = {"config": dict(parent_ids=args.cases, supplemental_sources=args.budget, near_radius=near_r,
                             same_radius=old_radius, old_residual_threshold=THRESHOLD, source_angle_offset_deg=15,
                             old_noise_relative=.001, source_energy="unit scalar monopole coefficient squared per exposure; same k and pulse duration",
                             field_match_ablation="scale each near source by the incident-field norm ratio between corresponding far and near angles; equal grid-mean incident intensity",
                             field_match_noise="same absolute receiver variance; unequal emitted energy is an ablation only",
                             initial_target_rms=.05, old_jacobian_weak_cutoff=.0002,
                             old_case_source=str(R4 / "data"), noise_seeds=[20261119, 702],
                             new_noise_calibration="equal absolute receiver variance estimated from old data norm",
                             material_scores="Archived GN–NN convention: RMSE and amplitude-normalized RMSE, 256-bin Otsu-IoU against alpha >= .5",
                             source_design_input="public noisy Y, public GN material, known operators"),
              "cases": []}
    material_arrays = {}
    with threadpool_limits(1):
        for case_id in args.cases:
            row = records[case_id]
            seed = int(case_id.split("p")[-1])
            obs = read_npz(r4/"data"/row["public_path"])
            Y = obs["Y"]
            anchor = read_npz(r4/"anchors"/f'{row["obs_id"]}_discrepancy_allscales.npz')["m0"]
            cand = weak_candidates(Model, reg.D, ops, Y, anchor, .05)
            noise_norm = float(obs["noise"] * np.linalg.norm(Y))
            noise_per_source = noise_norm/np.sqrt(Y.shape[1])
            sources = source_analysis(Model, ops, cand, lib, noise_per_source, args.budget,
                                      np.random.default_rng([20261119, seed]))
            save = {field: cand[field] for field in ("trials", "residuals", "separated", "distance_rmse", "distance_linf",
                                                      "weak_direction_relative_gain", "weak_dimension",
                                                      "bright_dimension", "counts", "seconds")}
            save.update(parent_id=case_id, obs_id=row["obs_id"], group=row["group"],
                        family=row["family"], digit=row.get("digit"), amplitude=row["amplitude"],
                        selection=sources["selections"], source_decomposition=sources["decomposition"],
                        decomposition_closure=sources["closure"],
                        distance_rmse_over_archived_amplitude=cand["distance_rmse"]/row["amplitude"],
                        old_pair_gap_over_noise=ratio(cand["plus"]["y"]-cand["minus"]["y"], noise_norm),
                        old_pair_relative_data_gap=ratio(cand["plus"]["y"]-cand["minus"]["y"], np.linalg.norm(Y)),
                        source_cost=dict(counts=sources["counts"], seconds=sources["seconds"],
                                         stages=sources["stage_cost"]))
            print(case_id, "old raw residuals", save["residuals"], "separation", save["distance_rmse"],
                  "max far/near/repeat gap over noise", [max(d["total_over_noise"] for d in sources["decomposition"]
                                                            if d["library"]==x) for x in ("same_radius","near_field","repeat")],
                  flush=True)
            material_arrays.update({case_id+"_anchor": anchor, case_id+"_plus": cand["plus"]["m"],
                                    case_id+"_minus": cand["minus"]["m"], case_id+"_s": cand["singular_values"],
                                    case_id+"_weak_direction": cand["direction"]})
            result["cases"].append(save)
            if not args.candidates_only:
                finish_case(Model, ops, reg, material_basis, row, r4, Y, anchor, cand, sources, noise_norm,
                            noise_per_source, seed, save, material_arrays)
            args.output.with_suffix(".json").write_text(json.dumps(result, indent=2, default=float)+"\n")
            np.savez_compressed(args.output.with_suffix(".npz"), **material_arrays)


def finish_case(Model, ops, reg, material_basis, row, r4, Y, anchor, cand, sources, noise_norm,
                noise_per_source, seed, save, arrays):
    case_id = row["parent_id"]
    setup_start = perf_counter()
    truth_record = read_npz(r4/"data"/row["truth_path"])
    gt = truth_record["m"]
    gain_model = Model(dict(ops, E0=sources["E0"]))
    clean_all = gain_model.forward(gt)["y"]
    noise_rng = np.random.default_rng([702, seed])
    noise = (noise_rng.normal(size=clean_all.shape) + 1j * noise_rng.normal(size=clean_all.shape))
    noise *= noise_per_source/np.linalg.norm(noise, axis=0)
    noisy_all = clean_all + noise
    arrays[case_id+"_truth"] = gt
    arrays[case_id+"_old_Y"] = Y
    arrays[case_id+"_new_Y_all"] = noisy_all
    save["truth_simulation_cost"] = dict(counts=dict(gain_model.counts), seconds=perf_counter()-setup_start)
    save["legacy"] = []
    for name, m in (
        ("allscales_discrepancy", anchor),
        ("fulltrf", read_npz(r4/"anchors"/f'{row["obs_id"]}.npz')["m0"]),
        ("voxel_discrepancy", read_npz(r4/"anchors"/f'{row["obs_id"]}_discrepancy.npz')["m0"]),
        ("weakreg1e-5", read_npz(r4/"anchors"/f'{row["obs_id"]}_weakreg1e-5.npz')["m0"]),
        ("finite_plus", cand["plus"]["m"]), ("finite_minus", cand["minus"]["m"]),
    ):
        model = Model(ops)
        state = model.forward(m)
        save["legacy"].append(dict(name=name, **archived_material_score(m, truth_record),
                                   relative_residual=ratio(state["y"]-Y, np.linalg.norm(Y)),
                                   noise_multiple=ratio(state["y"]-Y, noise_norm)))
    budgets = dict(warm=((None, 14),),
                   cold=((material_basis((24, 24), 8), 8), (material_basis((24, 24), 12), 8), (None, 14)))
    old_fit = optimize(Model, reg, ops, Y, noise_norm, np.full(576, .1), budgets["cold"],
                       True, "old_cold_reg")
    old_reg_m = old_fit.pop("material")
    save["old_cold_reg"] = dict(**old_fit, **archived_material_score(old_reg_m, truth_record))
    arrays[case_id+"_old_cold_reg"] = old_reg_m
    save["arms"] = []
    cached = {}
    for library, controls in sources["selections"].items():
        for selector, ids in controls.items():
            if tuple(ids) not in cached:
                cached[tuple(ids)] = reconstruct_arm(Model, reg, ops, Y, anchor, cand, sources, ids, noisy_all,
                                                     noise_norm, noise_per_source, budgets, truth_record)
            result = cached[tuple(ids)]
            save["arms"].append(dict(library=library, selector=selector, **{k: v for k, v in result.items()
                                                                           if k != "materials"}))
            for k, m in result["materials"].items():
                arrays[f"{case_id}_{library}_{selector}_{k}"] = m
    print(case_id, "allscales old RMSE", save["legacy"][0]["rmse"],
          [(arm["library"], arm["selector"], round(arm["gap_over_noise"], 3),
            round(arm["solutions"]["warm_reg"]["rmse"], 4))
           for arm in save["arms"]], flush=True)


def reconstruct_arm(Model, reg, ops, Y, anchor, cand, sources, ids, noisy_all, noise_norm, noise_per_source,
                    budgets, truth_record):
    begin = perf_counter()
    new_Y = noisy_all[:, ids]
    joint_E0 = np.column_stack((ops["E0"], sources["E0"][:, ids]))
    joint_Y = np.column_stack((Y, new_Y))
    joint_noise = np.hypot(noise_norm, noise_per_source*np.sqrt(len(ids)))
    arm_ops = dict(ops, E0=joint_E0)
    model = Model(arm_ops)
    candidate_ms = [cand["plus"]["m"], cand["minus"]["m"]]
    predictions = [model.forward(m)["y"] for m in candidate_ms]
    new_noise = noise_per_source * np.sqrt(len(ids))
    selection_err = [ratio(p[:, Y.shape[1]:] - new_Y, new_noise) for p in predictions]
    joint_errors = [ratio(p - joint_Y, joint_noise) for p in predictions]
    chosen = int(np.argmin(selection_err))
    initializations = (
        ("winner_reg", candidate_ms[chosen], budgets["warm"], True),
        ("warm_reg", anchor, budgets["warm"], True),
        ("cold_reg", np.full(576, .1), budgets["cold"], True),
        ("voxel_fullJ", anchor, budgets["warm"], False),
    )
    selected_m = candidate_ms[chosen]
    solutions = dict(winner_raw=dict(**archived_material_score(selected_m, truth_record),
                                     seconds=perf_counter()-begin, counts=dict(model.counts),
                                     relative_residual=ratio(predictions[chosen]-joint_Y,np.linalg.norm(joint_Y)),
                                     noise_multiple=joint_errors[chosen],
                                     discrepancy_reached=bool(joint_errors[chosen] <= 1.10)))
    materials = dict(winner_raw=selected_m)
    for key, m, stages, smooth in initializations:
        fit = optimize(Model, reg, arm_ops, joint_Y, joint_noise, m, stages, smooth, key)
        materials[key] = fit.pop("material")
        solutions[key] = dict(**fit, **archived_material_score(materials[key], truth_record))
    return dict(source_ids=ids,
                gap_over_noise=ratio(sources["delta"][:, ids], new_noise),
                old_combo_over_joint_new_noise=ratio(sources["combo"][:, ids], new_noise),
                perpendicular_over_joint_new_noise=ratio(sources["perpendicular"][:, ids], new_noise),
                field_match_gap_over_noise=ratio(sources["delta"][:, ids]*sources["field_match_scale"][ids], new_noise),
                field_match_perpendicular_over_joint_new_noise=ratio(
                    sources["perpendicular"][:, ids]*sources["field_match_scale"][ids], new_noise),
                selected_branch=["plus","minus"][chosen], selection_new_noise_multiples=selection_err,
                candidate_joint_noise_multiples=joint_errors,
                joint_noise_norm=float(joint_noise), solutions=solutions, materials=materials)


if __name__ == "__main__":
    main()
