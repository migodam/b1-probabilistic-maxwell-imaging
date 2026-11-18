"""Closed-loop discrimination of independently discovered branches of the scalar TM VIE."""

from itertools import combinations
from pathlib import Path
from time import process_time
import json

from branch_discovery import Experiment, LO, HI, SPAN, MATERIAL_TOL, pack, rms
import numpy as np
from scipy.linalg import lu_factor, lu_solve
from scipy.optimize import least_squares, minimize
from scipy.special import j0
from scipy.stats import chi2, qmc

ROOT = Path(__file__).resolve().parents[1]
COHORTS = ("branch_holdout", "branch", "branch_archive")
SIGMA = 2e-4             # RMS amplitude of each complex receiver noise sample
EXACT_SCALE = 2e-8       # Numerical scale for fitting exact old values
EXACT_LIMIT = 1e-7       # Absolute old-block L2 tolerance
POWER_FLOOR = 1 / 16     # Unit-power source feeds have coefficient norm <= 4
START_LIMIT = 100
NOISE_REPEATS = 3
MENU_SIZE = 64


def power_basis(gram):
    s, u = np.linalg.eigh(gram)
    use = s >= POWER_FLOOR
    return u[:, use] / np.sqrt(s[use])


def cnormal(rng, shape):
    return (rng.standard_normal(shape) + 1j * rng.standard_normal(shape)) / np.sqrt(2)


def library(geometry):
    old = np.asarray(geometry["sources"])
    source_radius = float(np.median(np.linalg.norm(old, axis=1)))
    angles = 2 * np.pi * (np.arange(20) + .371) / 20
    extra = source_radius * np.column_stack((np.cos(angles), np.sin(angles)))
    sources = np.vstack((old, extra))
    model = Experiment(dict(geometry, sources=sources.tolist()))
    distances = np.linalg.norm(sources[:, None] - sources[None, :], axis=-1)
    gram = j0(model.model.k * distances)
    return model, gram, len(old), sources


def unique_roots(name, records, arrays):
    rows, labels = [], []
    ordered = [("anchor", "anchor"), *((m, m) for m in records["methods"])]
    for method, key in ordered:
        item = f"{name}/{key}"
        if item not in arrays:
            continue
        sample = np.atleast_2d(arrays[item])
        for candidate in sample:
            match = next((i for i, row in enumerate(rows) if rms(row-candidate) <= MATERIAL_TOL), None)
            if match is None:
                rows.append(candidate)
                labels.append([method])
            elif method not in labels[match]:
                labels[match].append(method)
    return np.asarray(rows), labels


def maximin_direction(difference, seeds):
    """Optimize the worst pair gap in the row span of the trained pair differences."""
    _, s, vh = np.linalg.svd(difference, full_matrices=False)
    basis = vh[s > s[0] * 1e-11].conj().T
    reduced = difference @ basis / s[0]
    directions = seeds @ basis.conj()
    directions /= np.maximum(np.linalg.norm(directions, axis=1, keepdims=True), 1e-30)
    start_gap = np.min(abs(directions @ reduced.T), axis=1)
    order = np.argsort(start_gap)[::-1][:min(12, len(seeds))]
    best = directions[order[0]]
    best_gap = start_gap[order[0]]
    rank = basis.shape[1]
    if rank == 1:
        return basis[:, 0]

    def complex_vector(x):
        return x[:rank] + 1j * x[rank:2*rank]

    def constraints(x):
        a = complex_vector(x)
        return np.r_[1. - np.vdot(a, a).real, abs(reduced @ a)**2 - x[-1]]

    def jacobian(x):
        a = complex_vector(x)
        product = (reduced @ a).conj()[:, None] * reduced
        return np.vstack((np.r_[-2*x[:-1], 0.],
                          np.column_stack((2*product.real, -2*product.imag, -np.ones(len(reduced))))))

    objective_jacobian = np.r_[np.zeros(2*rank), -1.]
    for i in order:
        a = directions[i]
        x0 = np.r_[a.real, a.imag, start_gap[i]**2]
        result = minimize(lambda x: -x[-1], x0, jac=lambda x: objective_jacobian,
                          constraints=[{"type": "ineq", "fun": constraints, "jac": jacobian}],
                          method="SLSQP", options={"maxiter": 100, "ftol": 1e-12})
        a = complex_vector(result.x)
        a /= max(np.linalg.norm(a), 1e-30)
        gap = np.min(abs(reduced @ a))
        if gap > best_gap:
            best, best_gap = a, gap
    return basis @ best


def design(model, gram, old_count, roots, rng, training):
    setup_start = process_time()
    basis = power_basis(gram)
    state = [model.response(row) for row in roots[training]]
    y = np.array([row[0] for row in state])
    pairs = list(combinations(range(len(training)), 2))
    difference = np.array([(y[i]-y[j]) @ basis for i, j in pairs])
    norms = np.linalg.norm(difference, axis=1)
    matched = difference.conj() / np.maximum(norms[:, None], 1e-30)
    menu = cnormal(rng, (MENU_SIZE, basis.shape[1]))
    menu /= np.linalg.norm(menu, axis=1, keepdims=True)
    menu_gaps = np.min(abs(menu @ difference.T), axis=1)
    pair_gaps = np.min(abs(matched @ difference.T), axis=1)
    best_pair = int(np.argmax(pair_gaps))
    # Anchor is always included. Fisher sensitivity is local to this fitted root.
    anchor_jac = state[0][1]
    local = np.linalg.svd(anchor_jac[:2*old_count], full_matrices=True)[2][-1]
    complex_jac = anchor_jac[::2] + 1j * anchor_jac[1::2]
    fisher_direction = (complex_jac @ local) @ basis
    fisher = fisher_direction.conj() / max(np.linalg.norm(fisher_direction), 1e-30)
    old_basis = power_basis(gram[:old_count, :old_count])
    old_draw = cnormal(rng, old_basis.shape[1])
    old_feed = old_basis @ (old_draw / np.linalg.norm(old_draw))
    feeds = {
        "group maximin": basis @ maximin_direction(difference, np.vstack((matched, menu))),
        "fixed random new": basis @ menu[0],
        "local Fisher": basis @ fisher,
        "best pair match": basis @ matched[best_pair],
        "best random menu": basis @ menu[np.argmax(menu_gaps)],
        "old recombination": np.r_[old_feed, np.zeros(len(gram)-old_count)],
    }
    records = {
        "source_modes": basis.shape[1],
        "pairs": [[training[i], training[j]] for i, j in pairs],
        "training_roots": training,
        "pair_match_training_min_gap": pair_gaps.tolist(),
        "pair_match_feed_pair_gap_upper_bound": float(norms.min()),
        "random_menu_training_min_gap": menu_gaps.tolist(),
        "feeds": {},
        "select_cpu_s": process_time()-setup_start,
    }
    for rule, feed in feeds.items():
        pred = y @ feed
        records["feeds"][rule] = {
            "training_min_gap": float(min(abs(pred[i]-pred[j]) for i,j in pairs)),
            "radiated_power": float(np.vdot(feed, gram @ feed).real),
            "feed_l2": float(np.linalg.norm(feed)),
            "coefficient_real": feed.real.tolist(),
            "coefficient_imag": feed.imag.tolist(),
        }
    return feeds, records


class JointData:
    def __init__(self, model, old_count, feed, observed, noisy_old):
        self.b = np.column_stack((model.b[:, :old_count], model.b @ feed))
        self.h, self.r, self.eye = model.h, model.r, model.eye
        self.observed = observed
        self.scales = np.r_[np.full(old_count, SIGMA if noisy_old else EXACT_SCALE), SIGMA]
        self.evals = 0

    def response(self, material):
        lu = lu_factor(self.eye - material[:, None] * self.h, check_finite=False)
        current = lu_solve(lu, material[:, None] * self.b, check_finite=False)
        field = self.b + self.h @ current
        transfer = lu_solve(lu, self.r, trans=1, check_finite=False)
        self.evals += 1
        return self.r @ current, (transfer[:, None] * field).T

    def fit(self, initial):
        t0 = process_time()
        calls_before = self.evals
        previous, res, derivative = None, None, None

        def evaluate(x):
            nonlocal previous, res, derivative
            if previous is None or not np.array_equal(x, previous):
                previous = x.copy()
                y, j = self.response(LO + SPAN*x)
                res = pack((y-self.observed)/self.scales)
                j = SPAN * j / self.scales[:, None]
                derivative = np.stack((j.real, j.imag), axis=1).reshape(-1, j.shape[1])

        def fun(x):
            evaluate(x)
            return res

        def jac(x):
            evaluate(x)
            return derivative

        res = least_squares(fun, np.clip((initial-LO)/SPAN, 0, 1), jac=jac, bounds=(0,1),
                            max_nfev=START_LIMIT, ftol=1e-12, xtol=1e-12, gtol=1e-10)
        material = LO + SPAN*res.x
        y, _ = self.response(material)
        return material, y-self.observed, self.evals-calls_before, process_time()-t0


def acceptance(error, noisy_old):
    old_l2, new_abs = np.linalg.norm(error[:-1]), abs(error[-1])
    # Two real components with variance SIGMA**2/2 per complex reading.
    old_chi2, new_chi2 = 2*(old_l2/SIGMA)**2, 2*(new_abs/SIGMA)**2
    joint_chi2 = old_chi2 + new_chi2
    old_pass = old_chi2 <= chi2.ppf(.95, 2*(len(error)-1)) if noisy_old else old_l2 <= EXACT_LIMIT
    joint_pass = joint_chi2 <= chi2.ppf(.95, 2*len(error)) if noisy_old else old_pass and new_chi2 <= chi2.ppf(.95, 2)
    return dict(old_l2=float(old_l2), new_abs=float(new_abs), joint_accept=bool(joint_pass),
                old_accept=bool(old_pass), new_accept=bool(new_chi2 <= chi2.ppf(.95, 2)))


def record_fit(joint, initial, truth, roots, noisy_old, label):
    material, err, evaluations, seconds = joint.fit(initial)
    distance = [float(rms(material-root)) for root in roots]
    return {
        "init": label, "contrast": material.tolist(), "truth_RMS": float(rms(material-truth)),
        "distance_to_old_roots_RMS": distance, "nearest_old_root": int(np.argmin(distance)),
        "weighted_sum_sq": float(np.linalg.norm(err/joint.scales)**2),
        "forward_and_J_calls": evaluations, "cpu_s": seconds, **acceptance(err, noisy_old),
    }


def report_fit(fit):
    picked = min(fit, key=lambda r: r["weighted_sum_sq"])
    accepted = [r for r in fit if r["joint_accept"]]
    distinct = []
    for candidate in sorted(accepted, key=lambda r: r["weighted_sum_sq"]):
        m = np.asarray(candidate["contrast"])
        if all(rms(m-np.asarray(row["contrast"])) > MATERIAL_TOL for row in distinct):
            distinct.append(candidate)
    diameter = max((rms(np.asarray(a["contrast"])-b["contrast"]) for a,b in combinations(distinct,2)), default=0.)
    return dict(selected=picked, all_endpoints=fit, feasible_distinct_endpoints=distinct,
                feasible_diameter_RMS=float(diameter), calls=sum(r["forward_and_J_calls"] for r in fit),
                cpu_s=sum(r["cpu_s"] for r in fit))


def acquire(model, roots, truth, old_count, feeds, training, rng, repeat_count):
    reference = model.model.forward(truth, independent=True)["y"][0]
    exact = np.array([model.model.forward(m, independent=True)["y"][0] for m in roots])
    pairs = list(combinations(range(len(roots)), 2))
    held = next((i for i in range(len(roots)) if i not in training), None)
    n_start = max(4, len(training))
    ordinary = LO + SPAN*qmc.Sobol(len(truth), scramble=True, seed=rng.integers(2**31).item()).random_base2(3)
    ordinary[0] = (LO+HI)/2
    ordinary = ordinary[:n_start]
    warm = np.vstack((roots[training], ordinary[:n_start-len(training)]))
    epsilon = SIGMA * cnormal(rng, (repeat_count, old_count+1))
    output = {}
    for rule, feed in feeds.items():
        branch_prediction = exact @ feed
        truth_prediction = reference @ feed
        stats = {
            "all_pair_min_gap": float(min(abs(branch_prediction[i]-branch_prediction[j]) for i,j in pairs)),
            "held_pair_min_gap": None if held is None else float(min(abs(branch_prediction[held]-branch_prediction[i]) for i in training)),
            "independent_reading_real_imag": [float(truth_prediction.real), float(truth_prediction.imag)],
            "scores": {},
        }
        predictions = np.column_stack((exact[:, :old_count], branch_prediction))
        for noisy_old in (False, True):
            trials = []
            for draw in epsilon:
                old = reference[:old_count] + noisy_old*draw[:-1]
                observed = np.r_[old, truth_prediction+draw[-1]]
                joint = JointData(model, old_count, feed, observed, noisy_old)
                base_trf = [record_fit(joint, m, truth, roots, noisy_old, f"restart {i}")
                            for i,m in enumerate(ordinary)]
                branch_trf = [record_fit(joint, m, truth, roots, noisy_old,
                                         f"old root {training[i]}" if i<len(training) else f"restart {i-len(training)}")
                              for i,m in enumerate(warm)]
                residuals = [acceptance(pred-observed, noisy_old) for pred in predictions]
                trials.append(dict(
                    observation_real=observed.real.tolist(), observation_imag=observed.imag.tolist(),
                    uncorrected_old_roots=residuals,
                    uncorrected_joint_accepted=[i for i,r in enumerate(residuals) if r["joint_accept"]],
                    uncorrected_new_excluded=[i for i,r in enumerate(residuals) if r["old_accept"] and not r["new_accept"]],
                    restart=report_fit(base_trf), warm=report_fit(branch_trf),
                ))
            stats["scores"]["old noisy" if noisy_old else "old exact"] = trials
        output[rule] = stats
    return output


def plot_feedback(result):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    cases = [r for r in result["cases"].values() if r["cohort"] == "branch_holdout"]
    scores = [r["variants"]["all roots"] for r in cases]
    rules = ["group maximin", "best pair match", "local Fisher", "best random menu", "fixed random new"]
    labels = ["Group max–min", "Best pair match", "Local Fisher", "Best of 64 random", "Fixed random"]
    colors = ["#153d57", "#447487", "#7a4c81", "#9a572f", "#78937c"]
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.5), layout="constrained")
    for i, (rule, color) in enumerate(zip(rules, colors)):
        gaps = [row["natural_truth"][rule]["all_pair_min_gap"] / SIGMA for row in scores]
        jitter = np.linspace(-.14, .14, len(gaps))
        axes[0, 0].scatter(gaps, i+jitter, s=19, color=color, alpha=.60)
        axes[0, 0].scatter(np.median(gaps), i, s=43, marker="D", color=color, edgecolor="white", zorder=5)
    axes[0, 0].axvline(1, color=".52", linestyle="--", linewidth=1)
    axes[0, 0].set(xscale="log", yticks=range(len(labels)), yticklabels=labels,
                   xlabel="Minimum pair separation / complex noise RMS", title="A   Candidate separation at equal power")
    axes[0, 0].invert_yaxis()

    arms = ["restart", "warm"]
    arm_colors = ["#a56436", "#217476"]
    for j, (style, color) in enumerate(zip(arms, arm_colors)):
        passed = []
        for rule in [*rules, "old recombination"]:
            trials = [t for row in scores for t in row["natural_truth"][rule]["scores"]["old exact"]]
            passed.append(sum(t[style]["selected"]["joint_accept"] for t in trials))
        axes[0, 1].barh(np.arange(len(passed))+(.18 if j else -.18), passed, height=.30, color=color)
    axes[0, 1].set(xlim=(0,39), xlabel="Jointly accepted selections (39 trials)",
                   yticks=range(len(labels)+1), yticklabels=[*labels, "Old recombination"],
                   title="B   Exact old reads, noisy new read")
    axes[0, 1].set_ylim(6.0, -.5)
    axes[0, 1].legend(["TRF restart", "TRF warm start"], loc="lower right", frameon=False, fontsize=8)

    for i, (rule, color) in enumerate(zip(rules, colors)):
        trials = [t for row in scores for t in row["natural_truth"][rule]["scores"]["old noisy"]]
        for j, style in enumerate(arms):
            m = [t[style]["selected"]["truth_RMS"] for t in trials]
            dot = -.13 if j == 0 else .13
            axes[1, 0].scatter(m, i+dot+np.linspace(-.045, .045, len(m)), s=10,
                               color=arm_colors[j], alpha=.31)
            axes[1, 0].scatter(np.median(m), i+dot, s=33, marker="D", color=arm_colors[j],
                               edgecolor="white", zorder=5)
    axes[1, 0].set(yticks=range(len(labels)), yticklabels=labels, xlabel="Selected material RMS error",
                   title="C   Both old and new reads noisy")
    axes[1, 0].invert_yaxis()

    held = [(r["case"], r["variants"]["hold last discovery"]["excluded_root_as_target"][rules[0]]) for r in cases
            if "hold last discovery" in r["variants"]]
    for i, regime in enumerate(("old exact", "old noisy")):
        for _, row in held:
            q = row["scores"][regime][0]["warm"]["selected"]
            axes[1, 1].scatter(row["held_pair_min_gap"]/SIGMA, q["truth_RMS"], s=32,
                               facecolors=arm_colors[i] if q["joint_accept"] else "none",
                               edgecolors=arm_colors[i], marker=("o", "s")[i])
    axes[1, 1].axvline(1, color=".52", linestyle="--", linewidth=1)
    axes[1, 1].set(xscale="log", yscale="log", xlabel="Held-out pair separation / complex noise RMS",
                   ylabel="Warm-start material RMS error", title="D   Withheld branch as the target")
    keys = [Line2D([],[],marker=marker,color=arm_colors[i],linestyle="None",label=label)
            for i,(marker,label) in enumerate((("o","Exact old"),("s","Noisy old")))]
    keys.append(Line2D([],[],marker="o",color=".45",linestyle="None",markerfacecolor="none",label="Fails data test"))
    axes[1, 1].legend(handles=keys, loc="lower left", frameon=False, fontsize=8)
    for ax in axes.flat:
        ax.spines[["right", "top"]].set_visible(False)
        ax.grid(True, axis="x", alpha=.14)
        ax.set_axisbelow(True)
        ax.title.set_fontweight("bold")
    for ext in ("pdf", "png"):
        fig.savefig(ROOT / "figures" / f"feedback.{ext}", dpi=220)
    plt.close(fig)


def main():
    plate = {}
    experiments = []
    rng = np.random.default_rng(22110619)
    for cohort in COHORTS:
        path = ROOT / "results" / cohort
        records = json.loads(path.with_suffix(".json").read_text())["cases"]
        with np.load(path.with_suffix(".npz")) as arrays:
            for name, record in records.items():
                if name == "R5-T0" or (cohort == "branch" and name.startswith("R5")):
                    continue
                roots, labels = unique_roots(name, record, arrays)
                if len(roots) < 2:
                    continue
                model, gram, old_count, sources = library(record["geometry"])
                key = f"{cohort}/{name}"
                training = list(range(len(roots)))
                protocols = {"all roots": training}
                if len(roots) >= 3:
                    protocols["hold last discovery"] = training[:-1]
                y_old = arrays[f"{name}/y"]
                # A frozen design file is written before any ground truth arrays are accessed.
                info = dict(case=name, cohort=cohort, roots=roots.tolist(), method_labels=labels,
                            old_data_relative_errors=[float(np.linalg.norm(model.response(m)[0][:old_count]-y_old)/np.linalg.norm(y_old)) for m in roots],
                            old_channels=old_count, source_positions=sources.tolist(),
                            source_power_eigenvalues=np.linalg.eigvalsh(gram).tolist(), variants={})
                for variant, t in protocols.items():
                    # Compare full and held-out design with the exact same random menu.
                    feed_rng_seed = rng.integers(2**31).item() if variant=="all roots" else feed_rng_seed
                    feeds, row = design(model, gram, old_count, roots, np.random.default_rng(feed_rng_seed), t)
                    info["variants"][variant] = row
                    experiments.append((key, variant, path, name, model, roots, old_count, feeds, t))
                plate[key] = info
    (ROOT / "results" / "feedback_design.json").write_text(json.dumps(plate, indent=2, allow_nan=False))
    for key, variant, path, name, model, roots, old_count, feeds, training in experiments:
        with np.load(path.with_suffix(".npz")) as arrays:
            truth = arrays[f"{name}/truth"]
        row = plate[key]["variants"][variant]
        row["truth_distance_to_old_roots_RMS"] = [float(rms(m-truth)) for m in roots]
        row["natural_truth"] = acquire(model, roots, truth, old_count, feeds, training, rng, NOISE_REPEATS)
        if variant != "all roots":
            row["excluded_root_as_target"] = acquire(model, roots, roots[-1], old_count, feeds, training, rng, 1)
        print(key, variant, "roots", len(roots),
              "maximin min gap", f"{row['feeds']['group maximin']['training_min_gap']:.4g}", flush=True)
    result = dict(settings=dict(rng_seed=22110619, complex_noise_RMS=SIGMA, repeats=NOISE_REPEATS,
                                exact_old_fit_scale=EXACT_SCALE, exact_old_L2_limit=EXACT_LIMIT,
                                power_floor=POWER_FLOOR, menu_count=MENU_SIZE, max_TRF_nfev=START_LIMIT,
                                endpoint_dedup_RMS=MATERIAL_TOL), cases=plate)
    (ROOT / "results" / "feedback.json").write_text(json.dumps(result, separators=(",",":"), allow_nan=False))
    plot_feedback(result)


if __name__ == "__main__":
    main()
