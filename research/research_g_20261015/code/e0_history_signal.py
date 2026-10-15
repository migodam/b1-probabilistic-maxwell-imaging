"""Zero-solver, parent-grouped predictive audit of frozen Research A histories.

This program imports only the Python standard library and NumPy. Labels are
read AFTER feature construction, and no truth or offline metrics are loaded.
The old radial structural predicate is an exposed diagnostic, not G's metric.
"""
import argparse
import csv
import json
import os
from pathlib import Path
import time

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parents[1]
OLD = ROOT / "research/research_a_20261015/results"
ARMS = ("finite_static", "multistart", "population", "mdn", "adaptive")
SEED = 261015410


def dump(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False) + "\n")


def structural_difference(a, b):
    """Literal frozen Research A predicate; no import of a physics module."""
    return abs(a[0] - b[0]) >= .05 or abs(a[1] - b[1]) >= .025


def material_distance(a, b):
    """Frozen radial volume-weighted physical-field distance, using physical p."""
    edges = {0., 1.}
    for p in (a, b):
        if p[0] > 0:
            edges.add(float(p[0]))
        if p[1] > 0:
            edges.update([.65-float(p[1]), .65+float(p[1])])
    edges = sorted(edges)
    def eps(p, r):
        return 1.+0j if r < p[0] or (p[1] > 0 and .65-p[1] <= r < .65+p[1]) else complex(p[2], p[3])
    v = sum((hi**3-lo**3)*abs(eps(a, (hi+lo)/2)-eps(b, (hi+lo)/2))**2
            for lo, hi in zip(edges[:-1], edges[1:]))
    return float(np.sqrt(v)/abs(1.5+.1j))


def csv_write(path, records):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def base_features(parent, method, attempt, step):
    h = attempt["history"]
    y = np.asarray(parent["observation"], float)
    sigma = np.asarray(parent["sigma"], float)
    cur = np.asarray(h["parent_p"], float)
    prop = np.asarray(h["proposal_p"], float)
    box = h["corrector"].get("protected_geometry")
    # The BOX is a declared action constraint available before correction;
    # nfev, convergence and corrected residual are deliberately excluded.
    names = [f"observed_y_{i:02d}" for i in range(len(y))]
    x = list(y)
    names += ["noise_mean", "noise_min", "noise_max", "observation_l2"]
    x += [float(sigma.mean()), float(sigma.min()), float(sigma.max()), float(np.linalg.norm(y))]
    for group, values in (("current_material", cur), ("raw_proposal", prop), ("raw_delta", prop-cur)):
        names += [f"{group}_{j}" for j in range(4)]
        x += values.tolist()
    names += [f"policy_{a}" for a in ARMS]
    x += [float(method["arm"] == a) for a in ARMS]
    names += ["step_fraction", "remaining_endpoint_fraction", "protected_action", "edit_id_normalized", "raw_proxy_residual"]
    x += [step / len(method["attempts"]), (len(method["attempts"])-step)/len(method["attempts"]), float(box is not None), (h["edit_id"] if h["edit_id"] is not None else -1)/8., h["raw_proxy_residual"]]
    names += [f"protected_bound_{j}" for j in range(4)]
    x += list(box) if box is not None else [0.] * 4
    return names, x


def history_features(previous, edit_id, structural_count):
    """All arguments contain outcomes strictly before the current proposal."""
    n = len(previous)
    same = [a for a in previous if a["history"]["edit_id"] == edit_id]
    n_same = len(same)
    def summary(items):
        if not items:
            return [0., 0., 0., 0., 0.]
        return [float(np.mean([a["valid"] for a in items])),
                float(np.mean([a["discovery_status"] == "NEW_MATERIAL" for a in items])),
                float(np.mean([a["discovery_status"] == "DUPLICATE_MATERIAL" for a in items])),
                float(np.mean([a["residual"] for a in items])),
                float(np.min([a["residual"] for a in items]))]
    names = ["prefix_attempts", "prefix_same_action_attempts", "prefix_structural_representatives"]
    values = [float(n), float(n_same), float(structural_count)]
    for label, items in (("all", previous), ("recent2", previous[-2:]), ("same_action", same)):
        names += [f"prefix_{label}_{s}" for s in ("valid_rate", "discovery_rate", "duplicate_rate", "mean_residual", "best_residual")]
        values += summary(items)
    return names, values


def build_rows(parents):
    rows, names_a, names_h = [], None, None
    for parent in parents:
        for method in parent["methods"]:
            previous = []
            reps = [parent["base"]["p"]] if parent["base"]["valid"] else []
            full_valid_ledger = [parent["base"]["p"]] if parent["base"]["valid"] else []
            for step, attempt in enumerate(method["attempts"]):
                names_a, xa = base_features(parent, method, attempt, step)
                names_h, xh = history_features(previous, attempt["history"]["edit_id"], len(reps))
                # First access to CURRENT validation/endpoint comes here.
                archive_structural_novel = bool(attempt["valid"] and attempt["discovery_status"] == "NEW_MATERIAL"
                                        and all(structural_difference(attempt["p"], p) for p in reps))
                structural_novel = bool(attempt["valid"] and all(structural_difference(attempt["p"], p) for p in full_valid_ledger))
                material_novel = bool(attempt["valid"] and all(material_distance(attempt["p"], p) >= .03 for p in full_valid_ledger))
                rows.append(dict(parent=parent["id"], arm=method["arm"], step=step,
                                 features_a=xa, features_h=xh, valid=int(attempt["valid"]),
                                 corrected_structural_novelty=int(structural_novel),
                                 corrected_material_novelty=int(material_novel),
                                 legacy_archive_structural_discovery=int(archive_structural_novel),
                                 discovery_status=attempt["discovery_status"]))
                if archive_structural_novel:
                    reps.append(attempt["p"])
                if attempt["valid"]:
                    full_valid_ledger.append(attempt["p"])
                previous.append(attempt)
    return rows, names_a, names_h


def fit_predict(x, y, train, test):
    # Fixed ridge strength, no validation-selected hyperparameters.
    mean, scale = x[train].mean(axis=0), x[train].std(axis=0)
    scale[scale < 1e-10] = 1.
    z = np.column_stack([np.ones(len(x)), (x-mean)/scale])
    xt, yt = z[train], y[train]
    rate = (float(yt.sum()) + .5) / (len(yt)+1.)
    beta = np.zeros(z.shape[1]); beta[0] = np.log(rate/(1-rate))
    penalty = np.eye(z.shape[1]); penalty[0, 0] = 1e-6
    objective = lambda b: float(np.logaddexp(0., xt@b).sum()-yt@(xt@b)+.5*b@penalty@b)
    for _ in range(35):
        p = 1./(1.+np.exp(-np.clip(xt@beta, -35, 35)))
        w = np.maximum(p*(1-p), 1e-8)
        gradient = xt.T@(p-yt) + penalty@beta
        hessian = xt.T@(w[:, None]*xt) + penalty
        delta = np.linalg.solve(hessian, gradient)
        old_obj = objective(beta)
        alpha = 1.
        while objective(beta-alpha*delta) > old_obj-1e-4*alpha*float(gradient@delta):
            alpha *= .5
            if alpha < 2**-24:
                break
        beta -= alpha*delta
        assert objective(beta) <= old_obj+1e-6
        if np.max(np.abs(alpha*delta)) < 1e-8:
            break
    return np.clip(1./(1.+np.exp(-np.clip(z[test]@beta, -35, 35))), 1e-8, 1-1e-8)


def auc(y, p):
    pos, neg = p[y == 1], p[y == 0]
    if not len(pos) or not len(neg):
        return None
    return float(np.mean((pos[:, None] > neg[None, :]) + .5*(pos[:, None] == neg[None, :])))


def metrics(y, p):
    ece = 0.
    bins = []
    for i in range(10):
        keep = (p >= i/10) & (p < (i+1)/10 if i < 9 else p <= 1.)
        if keep.any():
            c, r = float(p[keep].mean()), float(y[keep].mean())
            ece += float(keep.mean())*abs(c-r)
            bins.append(dict(bin=i, count=int(keep.sum()), predicted=c, observed=r))
    return dict(log_loss=float(np.mean(-y*np.log(p)-(1-y)*np.log(1-p))),
                brier=float(np.mean((p-y)**2)), roc_auc=auc(y, p), ece_10bins=ece,
                positive_rate=float(y.mean()), calibration_bins=bins)


def matched_derangement(ids, observed, training_ids):
    """Nearest-observation pairs/triples; a permutation, within split only."""
    obs = np.stack([observed[p] for p in ids])
    train = np.stack([observed[p] for p in training_ids])
    scale = train.std(axis=0); scale[scale < 1e-10] = 1.
    dist = np.linalg.norm((obs[:, None, :]-obs[None, :, :])/scale, axis=2)
    remaining = list(range(len(ids))); donor = {}
    while len(remaining) > 3:
        a, b = min(((a, b) for i, a in enumerate(remaining) for b in remaining[i+1:]), key=lambda ab: dist[ab])
        donor[ids[a]], donor[ids[b]] = ids[b], ids[a]
        remaining.remove(a); remaining.remove(b)
    if len(remaining) in (2, 3):
        for a, b in zip(remaining, remaining[1:]+remaining[:1]):
            donor[ids[a]] = ids[b]
    if len(remaining) == 1:
        raise ValueError("A singleton test fold cannot support matched-history permutation")
    assert sorted(donor) == sorted(donor.values())
    assert all(k != v for k, v in donor.items())
    return donor


def confidence_delta(rows, labels, pred, comparator, seed=SEED+2):
    groups = sorted(set(r["parent"] for r in rows))
    loss = lambda p: -labels*np.log(p)-(1-labels)*np.log(1-p)
    delta = loss(pred)-loss(comparator)
    parent_delta = np.array([delta[np.array([r["parent"] == g for r in rows])].mean() for g in groups])
    rng = np.random.default_rng(seed)
    bootstrap = parent_delta[rng.integers(len(groups), size=(2000, len(groups)))].mean(axis=1)
    return dict(mean_log_loss_delta=float(parent_delta.mean()), parent_bootstrap95=np.quantile(bootstrap, [.025, .975]).tolist(),
                negative_is_better=True, parents=len(groups), resamples=2000,
                scope="descriptive cross-fitted predictions; models not refit in bootstrap; not confirmatory policy benefit")


def main():
    parser = argparse.ArgumentParser(); parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    if not args.execute:
        raise SystemExit("Use --execute for the frozen zero-physics analysis")
    os.nice(10)
    start, cpu = time.perf_counter(), time.process_time()
    parents = json.loads((OLD/"CONFIRM_CAMPAIGN.json").read_text())
    rows, na, nh = build_rows(parents)
    assert len(parents) == 32 and len(rows) == 1920
    x_a = np.asarray([r["features_a"] for r in rows]); x_h = np.asarray([r["features_h"] for r in rows])
    assert np.isfinite(x_a).all() and np.isfinite(x_h).all()
    # Direct feature-leakage regression: mutating current/future outcomes must
    # leave features of the current decision unchanged.
    import copy
    altered = copy.deepcopy(parents[:1])
    altered[0]["methods"][0]["attempts"][3]["valid"] = not altered[0]["methods"][0]["attempts"][3]["valid"]
    altered[0]["methods"][0]["attempts"][3]["residual"] += 999.
    altered_rows, _, _ = build_rows(altered)
    assert rows[3]["features_a"] == altered_rows[3]["features_a"]
    assert rows[3]["features_h"] == altered_rows[3]["features_h"]
    groups = sorted(p["id"] for p in parents)
    rng = np.random.default_rng(SEED)
    order = rng.permutation(groups)
    folds = {str(g): int(i % 5) for i, g in enumerate(order)}
    observed = {p["id"]: np.asarray(p["observation"], float) for p in parents}
    row_lookup = {(r["parent"], r["arm"], r["step"]): i for i, r in enumerate(rows)}
    donor_rows = np.empty(len(rows), int); permutations = []
    split = np.array([folds[r["parent"]] for r in rows])
    predictions = {label: {model: np.empty(len(rows)) for model in ("intercept", "A", "B", "C", "D")} for label in ("valid", "corrected_structural_novelty", "corrected_material_novelty", "legacy_archive_structural_discovery")}
    for fold in range(5):
        train, test = np.flatnonzero(split != fold), np.flatnonzero(split == fold)
        train_ids = [g for g in groups if folds[g] != fold]
        test_ids = [g for g in groups if folds[g] == fold]
        assert not set(train_ids).intersection(test_ids)
        perm = {}
        for ids, split_name in ((train_ids, "train"), (test_ids, "test")):
            mapping = matched_derangement(ids, observed, train_ids)
            perm.update(mapping)
            permutations += [dict(fold=fold, split=split_name, recipient=k, history_donor=v) for k, v in mapping.items()]
        donor = np.array([row_lookup[(perm[r["parent"]], r["arm"], r["step"])] for r in rows])
        donor_rows[test] = donor[test]
        matrices = dict(A=x_a, B=np.column_stack([x_a, x_h]), C=np.column_stack([x_a, x_h[donor]]), D=x_h)
        for label, models in predictions.items():
            y = np.asarray([r[label] for r in rows], float)
            models["intercept"][test] = (y[train].sum()+.5)/(len(train)+1.)
            for model, matrix in matrices.items():
                models[model][test] = fit_predict(matrix, y, train, test)
        if time.perf_counter()-start > 100:
            raise RuntimeError("Analysis pre-emption before the 120-second task wall bound")
    result = dict(status="EXECUTED_EXPOSED_GROUPED_PREDICTIVE_DIAGNOSTIC", seed=SEED,
                  source_revision="f1b674eb4893e39126701ed48f3d3dcef3098c3c",
                  source="research/research_a_20261015/results/CONFIRM_CAMPAIGN.json", parents=32,
                  rows=len(rows), steps_per_policy=12, folds=folds,
                  feature_count=dict(A=len(na), history=len(nh), B=len(na)+len(nh), C=len(na)+len(nh), D=len(nh)),
                  fit="logistic Newton/IRLS with monotone backtracking; L2=1 (intercept near-unpenalized), max35 iterations; fixed without target-based selection",
                  C_definition="B with prefix history permuted between near-observation parents; same arm/step; train and held-out folds permuted separately; no self donor",
                  labels={}, comparisons={}, truth_read=False, solver_instantiated=False,
                  logged_propensities="NOT_AVAILABLE", causal_policy_gain="NOT_IDENTIFIED",
                  labels_definition=dict(corrected_structural_novelty="valid corrected physical p, separated from EVERY prefix full-valid endpoint including duplicates by delta-r>=.05R OR delta-h>=.025R; base included when valid",
                                         corrected_material_novelty="valid corrected physical p, distance>=.03 from EVERY prefix full-valid endpoint under exact shell-volume material L2 normalized by abs(1.5+.1j)",
                                         legacy_archive_structural_discovery="secondary: original NEW_MATERIAL plus historical greedy geometry representative criterion"),
                  topology_scope="Old radial geometry predicate delta-r>=.05 or delta-halfwidth>=.025; no claim of new 2D G topology discovery")
    records = []
    for label, models in predictions.items():
        y = np.asarray([r[label] for r in rows], float)
        result["labels"][label] = {m: metrics(y, p) for m, p in models.items()}
        result["comparisons"][label] = {"B_minus_A": confidence_delta(rows, y, models["B"], models["A"]),
                                             "B_minus_C": confidence_delta(rows, y, models["B"], models["C"])}
        for i, r in enumerate(rows):
            record = {k: r[k] for k in ("parent", "arm", "step", "discovery_status")}
            record.update(fold=folds[r["parent"]], label=label, outcome=int(y[i]), history_donor=rows[donor_rows[i]]["parent"])
            record.update({f"prediction_{m}": float(p[i]) for m, p in models.items()})
            records.append(record)
    # Mechanistic descriptive counts: do not read offline MAE/truth arrays.
    policy = {}
    for arm in ARMS:
        subset = [r for r in rows if r["arm"] == arm]
        methods = [next(m for m in p["methods"] if m["arm"] == arm) for p in parents]
        policy[arm] = dict(valid=sum(r["valid"] for r in subset), new_material=sum(r["discovery_status"] == "NEW_MATERIAL" for r in subset),
                           duplicate=sum(r["discovery_status"] == "DUPLICATE_MATERIAL" for r in subset),
                           invalid=sum(not r["valid"] for r in subset), new_structural=sum(r["corrected_structural_novelty"] for r in subset),
                           corrected_material_novelty=sum(r["corrected_material_novelty"] for r in subset),
                           legacy_archive_structural_discovery=sum(r["legacy_archive_structural_discovery"] for r in subset),
                           proxy_evaluations=sum(m["proxy_evaluations"] for m in methods),
                           summed_search_wall_seconds=sum(m["wall_seconds"] for m in methods))
    result["descriptive_policy_outcomes"] = policy
    e2 = json.loads((OLD/"E2/RESULT_SUMMARY.json").read_text())
    result["A30_E2_history_limit"] = dict(parents=4, conditions=8, new_structural_family_beyond_multistart=e2["conclusions"]["new_structural_family_beyond_multistart"],
                                           validity_cv="NOT_FIT_ON_NONMATCHING_TARGET", novelty_cv="UNIDENTIFIABLE_NO_POSITIVE_NEW_STRUCTURE", reason="Do not pool scalar-TM trajectories with radial Mie; only four parents and zero new structure labels")
    provenance = []
    for name in na:
        origin = "observation/sigma" if name.startswith(("observ", "noise")) else "history.parent_p/proposal_p/source/edit_id plus declared method horizon"
        if name.startswith("protected_"):
            origin = "history.corrector.protected_geometry: selected action box, available before correction; optimizer results excluded"
        if name == "raw_proxy_residual":
            origin = "history.raw_proxy_residual, computed by original search before correct()"
        provenance.append(dict(feature=name, group="A", source=origin, availability="PRE_CORRECTION", truth_used=False))
    for name in nh:
        provenance.append(dict(feature=name, group="history", source="same-arm steps strictly less than t; source outcomes only after their past validation", availability="PREFIX_ONLY", truth_used=False))
    csv_write(OUT/"results/E0_HISTORY_GROUPED_PREDICTIONS.csv", records)
    csv_write(OUT/"results/E0_HISTORY_FEATURE_PROVENANCE.csv", provenance)
    csv_write(OUT/"results/E0_HISTORY_SHUFFLE_MAPPING.csv", permutations)
    result["resources"] = dict(controlled_wall_seconds=time.perf_counter()-start, cpu_seconds=time.process_time()-cpu,
                                 new_full_state_attempts=0, complex_full_solver_RHS=0, GPU_seconds=0,
                                 priority="nice+10", blas_threads=os.environ.get("OPENBLAS_NUM_THREADS"), hard_wall_limit=120)
    dump(OUT/"results/E0_HISTORY_SIGNAL.json", result)
    print(json.dumps(dict(resources=result["resources"], comparisons=result["comparisons"], policies=policy), indent=2))


if __name__ == "__main__":
    main()
