"""Acquire a fourth physical read for each naturally generated six-cell case."""

import json

from branch_feedback import cnormal, design, library, power_basis, unique_roots
from branch_discovery import Experiment, HI, LO, ROOT, rms
from homotopy_control import run as find_central_roots
from noise_sets import EPS, GAIN, noisy_vector
import numpy as np

SOURCE_SEED = 892019
NOISE_SEED = 481119
ARMS = ("active", "random", "fisher")


def select_feeds(model, gram, old_count, central, seed):
    if len(central) >= 2:
        feeds, _ = design(model, gram, old_count, central, np.random.default_rng(seed),
                          list(range(len(central))))
        return dict(active=feeds["group maximin"], random=feeds["fixed random new"],
                    fisher=feeds["local Fisher"]), "central roots maximin"

    basis = power_basis(gram)
    draw = cnormal(np.random.default_rng(seed), basis.shape[1])
    random = basis @ (draw / np.linalg.norm(draw))
    if not len(central):
        return dict.fromkeys(ARMS, random), "random fallback, zero central roots"

    _, j = model.response(central[0])
    weakest = np.linalg.svd(j[:2*old_count], full_matrices=True)[2][-1]
    sensitivity = ((j[::2] + 1j*j[1::2]) @ weakest) @ basis
    feed = basis @ (sensitivity.conj() / np.linalg.norm(sensitivity))
    return dict(active=feed, random=random, fisher=feed), "Fisher fallback, one central root"


def as_complex_pair(value):
    return [float(value.real), float(value.imag)]


def report_feed(feed, gram, predictions):
    differences = [abs(a-b) for i,a in enumerate(predictions) for b in predictions[i+1:]]
    return {
        "radiated_power": float(np.vdot(feed, gram @ feed).real),
        "feed_l2": float(np.linalg.norm(feed)),
        "central_pair_min_gap": float(min(differences)) if differences else None,
    }


def acquire_old(geometry, truth, noise, gain_error):
    model = Experiment(geometry).model
    clean = model.forward(truth, independent=True)["y"][0]
    return (1 + gain_error) * clean + noise[:len(clean)]


def acquire_extra(model, truth, feeds, noise, gain_error):
    response = model.model.forward(truth, independent=True)["y"][0]
    return {name: (1 + gain_error) * (response @ feed) + noise[-1] for name,feed in feeds.items()}


def main():
    path = ROOT / "results" / "branch_holdout"
    metadata = json.loads(path.with_suffix(".json").read_text())["cases"]
    names = sorted(name for name in metadata if name.startswith("fresh6-"))
    data = dict(case_ids=names, arms=ARMS, eps=EPS, gamma=GAIN, box=[LO, HI])
    settings = dict(additive_radius=float(EPS), gain_radius=float(GAIN),
                    additive_radius_fraction=.55, gain_radius_fraction=.8,
                    joint_radius_used=float(np.hypot(.55,.8)),
                    source_seed_base=SOURCE_SEED, noise_seed_base=NOISE_SEED,
                    pivot_index_origin=0, truth_use="independent full VIE acquisitions and post-design scoring")
    report = {"settings": settings, "cases": {}}
    with np.load(path.with_suffix(".npz")) as cases:
        for row_id, name in enumerate(names):
            geometry = metadata[name]["geometry"]
            old = Experiment(geometry)
            rng = np.random.default_rng(NOISE_SEED + row_id)
            epsilon = noisy_vector(rng, 4) * EPS * .55
            gain_error = GAIN * .8 * np.exp(2j*np.pi*rng.uniform())

            y_old = acquire_old(geometry, cases[f"{name}/truth"], epsilon, gain_error)
            hits = find_central_roots(geometry, y_old)
            candidates = np.asarray(hits["validated_roots"]).reshape(-1, old.n)
            if len(candidates):
                midpoint = (LO+HI)/2
                candidates = candidates[np.argsort(np.linalg.norm(candidates-midpoint, axis=1), kind="stable")]
            model, gram, old_count, source_positions = library(geometry)
            feeds, rule = select_feeds(model, gram, old_count, candidates, SOURCE_SEED+row_id)
            pivot = np.argsort(-np.ptp(candidates, axis=0), kind="stable")[:2] if len(candidates)>=2 else np.array([0,1])
            # Source choices and pivots above use only the noisy old read and known geometry.
            truth = cases[f"{name}/truth"]
            new = acquire_extra(model, truth, feeds, epsilon, gain_error)
            clean_roots, labels = unique_roots(name, metadata[name], cases)
            clean_roots = clean_roots.reshape(-1,old.n)
            record = {
                "rule": rule,
                "old_shared_gain_error": as_complex_pair(gain_error),
                "additive_noise": [as_complex_pair(z) for z in epsilon],
                "pivots": pivot.tolist(),
                "central_roots_count": len(candidates),
                "central_roots": candidates.tolist(),
                "central_truth_RMS": [float(rms(c-truth)) for c in candidates],
                "clean_witness_count": len(clean_roots),
                "clean_witness_methods": labels,
                "central_path_classifications": hits["path_classifications"],
                "central_paths_attempted": hits["paths_attempted"],
                "central_solve_status": hits["solve_status"],
                "central_total_cpu_s": hits["total_cpu_seconds"],
                "old_read": [as_complex_pair(z) for z in y_old],
                "arms": {},
            }
            payload = dict(B0=old.b, H=old.h, R=old.r, host=np.zeros(old.n), y_old=y_old,
                           central_roots=candidates, clean_witnesses=clean_roots, truth=truth,
                           pivots=pivot, source_positions=source_positions, power_gram=gram,
                           additive_noise=epsilon, gain_error=gain_error)
            data.update({f"{name}/{k}":v for k,v in payload.items()})
            predictions = [model.response(c)[0] for c in candidates]
            for arm, feed in feeds.items():
                b = model.b @ feed
                B = np.column_stack((old.b,b))
                y = np.r_[y_old,new[arm]]
                responses = np.asarray([v @ feed for v in predictions])
                record["arms"][arm] = dict(report_feed(feed, gram, responses),
                                            reading=as_complex_pair(new[arm]),
                                            role="random fallback" if not len(candidates) else arm,
                                            rank=int(np.linalg.matrix_rank(B)))
                for k,v in dict(b=b,B=B,y=y,y_host=np.zeros(4,dtype=complex),feed=feed,pivots=pivot).items():
                    data[f"{name}/{arm}/{k}"] = v
            report["cases"][name] = record
            print(name, len(candidates), rule, pivot.tolist(), flush=True)
    np.savez_compressed(ROOT/"results"/"natural_inputs.npz", **data)
    (ROOT/"results"/"natural_acquisition.json").write_text(json.dumps(report,indent=2,allow_nan=False))


if __name__ == "__main__":
    main()
