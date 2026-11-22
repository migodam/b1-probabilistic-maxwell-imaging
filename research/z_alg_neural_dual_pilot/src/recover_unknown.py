"""Bounded second Fourier search with a gap proposed from native geometry."""

import argparse
from itertools import product
import json
from pathlib import Path
import sys
from collections import Counter
from time import sleep

sys.dont_write_bytecode = True

import numpy as np
from flint import acb_mat, arb
from scipy.linalg import eigh

from native import BoundedEnvelope, Native, PAPER, ball, gershgorin, hermitian, high, star, _paper
from score_batch import (append_receipt, atomic_json, clock, elapsed, find_open, inputs,
                         load_receipts, profile_charge, read_new_lines, shard_pair)


PILOT = Path(__file__).resolve().parents[1]
MAX_EIGENSOLVES = 14
MAX_CUTOFF = 13
BETA = (.5, .7, .85, 1., 1.15, 1.3, 1.45, 1.6, 1.85, 2.2, 3., 4., 6., 8., 12., 24., 48.)
SHIFTS = (0., .1875, .375, .5625)


class BorderedProofFinished(ArithmeticError):
    pass


def strict_bordered_upper(action, y, nu, radius, proposal):
    """Upper screen by exact congruence on [M, X; 0, I], before any K-gap proof."""
    started = clock()
    trial, _, source_eigenvalues, numerical = proposal
    if numerical["distance"] >= radius:
        return dict(status="INCONCLUSIVE", upper=None, reason="numeric_distance_at_or_above_radius",
                    numerical_distance=numerical["distance"], **elapsed(started))
    g, w = action.g, action.w
    scale = arb(w.scale)
    yb = ball(y)
    x = ball(trial)
    ax = action.exact_action(x)
    _, z = action.exact_scaled_coefficients()
    rhs = star(g.sensor_ball) * yb / scale + _paper.mult_rows(g.incident_ball, z)
    h0 = star(yb) * yb + scale * (star(x) * ax - star(rhs) * x - star(x) * rhs)
    s = y.shape[1]
    matrix = acb_mat(nu + s, nu + s)
    for j in range(s):
        for k in range(s):
            matrix[nu + j, nu + k] = h0[j, k]
    if nu:
        state = ball(action.eig_basis[:, :nu] / np.sqrt(w.scale))
        negative = scale * (star(state) * action.exact_action(state))
        cross = scale * (star(state) * (ax - rhs))
        for j in range(nu):
            for k in range(nu):
                matrix[j, k] = negative[j, k]
            for k in range(s):
                matrix[j, nu + k] = cross[j, k]
                matrix[nu + k, j] = cross[j, k].conjugate()
    mid = hermitian(_paper.mid(matrix))
    checks = []

    def negative_at(a):
        center = mid.copy()
        center[nu:, nu:] -= a * np.eye(s)
        values, vectors = eigh(center, subset_by_index=(0, s - 1), check_finite=False)
        count = int(np.count_nonzero(values < -1.e-13))
        if not count:
            return 0
        v = ball(vectors[:, :count])
        sources = ball(vectors[nu:, :count])
        vv = star(v) * (matrix * v) - arb(a) * (star(sources) * sources)
        while count:
            folded = acb_mat(count, count, [vv[j, k] for j in range(count) for k in range(count)])
            edge = gershgorin(folded, upper=True)
            if edge < 0:
                checks.append(dict(a=a, negative_at_least=count, upper=high(edge)))
                return count
            count -= 1
        return 0

    q = s - negative_at(0.)
    upper = None
    if q == 0:
        upper = 0.
    else:
        a = float(np.nextafter(radius**2 / q * (1 - 1.e-11), 0.))
        if a > 0. and negative_at(a) == s:
            candidate = high((arb(q) * arb(a)).sqrt())
            if candidate <= radius:
                upper = candidate
        if upper is None:
            deficit = radius**2 - numerical["distance"]**2
            padding = deficit / (2 * q)
            bins = [float(max(0., v) + padding) for v in source_eigenvalues[nu:nu + q]]
            if (len(bins) == q and all(a > 0. and negative_at(a) >= s - i
                                      for i, a in enumerate(bins))):
                candidate = high(sum((arb(a) for a in bins), arb(0)).sqrt())
                if candidate <= radius:
                    upper = candidate
    status = "PROVED_ZERO" if upper == 0. else "PROVED_BELOW_RADIUS" if upper is not None else "INCONCLUSIVE"
    return dict(status=status, upper=upper, numerical_distance=numerical["distance"],
                trials=checks, state_columns=nu, native_stationary_columns=s, **elapsed(started))


class RecoveryEnvelope(BoundedEnvelope):
    """Trial gap and candidate priority use only the state operator."""

    def select(self, weights, nu):
        self.trial_eigensolves = 0
        edge = float(self.action.eigs[nu])
        self.proposal_gap, self.proposed_edge = None, edge
        if not edge > 0.:
            raise ArithmeticError("No positive native numerical edge to propose a gap.")
        # This number chooses a proposal only. The inherited Arb check supplies
        # spectral evidence. In particular the numerical edge is never a bound.
        gap = float(2. ** np.floor(np.log2(min(1 / 1024, edge / 4))))
        if gap <= 2. ** -25:
            raise ArithmeticError("Numerical gap is below the fixed recovery proposal floor.")
        sector = np.asarray(weights.real) / weights.scale
        optical = np.asarray(weights.optical) / weights.scale
        variation = np.hypot(weights.cap * (sector - 1), optical - optical.mean()).max() / 2
        if weights.is_global:
            schedule = ((8, 1., 0.), (10, 1., 0.), (12, 1., 0.), (13, 1., 0.))
        else:
            cuts = (8, 10, 12, 13) if variation < .075 else (12, 13, 10)
            candidates = {}
            for cut in cuts:
                profiles = []
                for beta, shift in product(BETA, SHIFTS):
                    proposal = self._proposal(weights, cut, beta, shift, gap)
                    if proposal is not None:
                        profiles.append((float(np.min(proposal[0])), beta, shift))
                # Include the strongest diagonal-floor proposal by shift so one
                # successful shift family cannot consume the entire allowance.
                winners = [max((row for row in profiles if row[2] == shift),
                               default=None) for shift in SHIFTS]
                candidates[cut] = sorted((row for row in winners if row is not None), reverse=True)
            schedule = [(cut, beta, shift) for cut in cuts
                        for _, beta, shift in candidates[cut]]
            # Try the physically unshifted multiplier cross term near its useful
            # Young-inequality scales before revisiting narrow original settings.
            schedule += [(13, beta, 0.) for beta in (1.15, 1.3, 1.45, .85)]
        self.proposal_gap = gap
        self.proposed_edge = edge
        for cut, beta, shift in dict.fromkeys(schedule):
            proposal = self._proposal(weights, cut, beta, shift, gap)
            if proposal is None:
                continue
            need = max(0, int(np.count_nonzero(proposal[1] < 0)) - nu)
            if need:
                if self.trial_eigensolves >= self.max_eigensolves:
                    break
                compact = self._matrix(cut, proposal, weights, beta)
                self.trial_eigensolves += 1
                edge_value = eigh(hermitian(compact), eigvals_only=True, check_finite=False,
                                  subset_by_index=(need - 1, need - 1))[0]
                if edge_value >= -1.e-8:
                    continue
            iy, ix, *_ = self._mode_data(cut)
            return proposal[-1], beta, shift, gap, tuple(zip(iy.tolist(), ix.tolist())), weights.is_global, need
        raise ArithmeticError("Bounded secondary Fourier search did not prove a positive spacing.")

    def verify(self, weights, nu):
        self.trial_eigensolves = 0
        self.proposal_gap, self.proposed_edge = None, None
        if self.runner.bordered and self.runner.radius is not None:
            started = clock()
            try:
                proposal = self.action.numerical_distance(self.runner.y, nu)
                self.action.numerical_distance = lambda y, nu: proposal
                detail = strict_bordered_upper(self.action, self.runner.y, nu, self.runner.radius, proposal)
            except (ArithmeticError, ValueError) as error:
                detail = dict(status="INCONCLUSIVE", upper=None, explanation=str(error))
            self.runner.last_screen = detail | {"total": elapsed(started)}
            if detail["upper"] is not None:
                raise BorderedProofFinished("The exact bordered upper already resolves the threshold.")
        return super().verify(weights, nu)


class RecoveryNative(Native):
    def __init__(self, *, bordered=True):
        super().__init__(max_cutoff=MAX_CUTOFF, max_eigensolves=MAX_EIGENSOLVES)
        self.bordered = bordered

    def _profile(self, w, geometry, action, stages, nu_hint):
        side = geometry.block_side
        if side not in self.envelopes:
            self.envelopes[side] = RecoveryEnvelope(geometry, self.max_cutoff, self.max_eigensolves)
        envelope = self.envelopes[side]
        envelope.action, envelope.runner = action, self
        proof = super()._profile(w, geometry, action, stages, nu_hint)
        if proof.get("envelope_proposal_eigensolves") is not None and not proof.get("cache_hit"):
            proof["secondary_trial_gap"] = getattr(envelope, "proposal_gap", None)
            proof["secondary_numeric_edge"] = getattr(envelope, "proposed_edge", None)
        return proof

    def score(self, y, u, mu, lam, *, radius=None, nu_hint=None):
        self.y, self.radius, self.last_screen = y, radius, None
        result = super().score(y, u, mu, lam, radius=radius, nu_hint=nu_hint)
        if self.last_screen is not None:
            result["bordered_upper"] = self.last_screen
            if self.last_screen["upper"] is not None:
                # The state-only profile cache cannot contain a Y-dependent upper.
                self.profiles.pop(self._weights(u, mu, lam), None)
                result.update(status=self.last_screen["status"], upper=self.last_screen["upper"])
                result.pop("explanation", None)
                result["proof"]["method"] = "bordered_source_congruence"
                result["proof"].pop("explanation", None)
        return result


def training_check(bordered=True):
    p = PILOT / "data/synthetic_observations.npz"
    with np.load(p, allow_pickle=False) as f:
        ix = f["object_id"].tolist().index("dualpilot_train_000")
        iy = f["noise_fraction"].tolist().index(.02)
        assert f["split"][ix] == "train"
        y, radius = f["Y"][ix, iy], float(f["radius"][ix, iy])
    profiles = json.loads((PAPER / "submission/single/coherent_general_profiles.json").read_text())
    base, recovery = Native(), RecoveryNative(bordered=bordered)
    records = []
    for ix in (0, 2):
        p = profiles[ix]
        mu = p["mu"] * np.asarray(p["sector"])
        lam = p["mu"] * np.asarray(p["optical"])
        primary = base.score(y, .75, mu, lam, radius=radius)
        print("training control", ix, primary["status"], primary["timing"]["total"]["wall_seconds"], flush=True)
        supplement = recovery.score(y, .75, mu, lam, radius=radius)
        print("training recovery", ix, supplement["status"], supplement["lower"], supplement["upper"],
              supplement.get("explanation"), supplement["timing"]["total"]["wall_seconds"], flush=True)
        records.append(dict(profile=ix, primary=primary, supplement=supplement))
        result = dict(scope="Only training observation dualpilot_train_000, 2% paired noise", rows=records,
                      first_setup=base.cache_summary(), supplement_setup=recovery.cache_summary())
        name = "native_recovery_bordered_training_check.json" if bordered else "native_recovery_training_check.json"
        atomic_json(PILOT / "results" / name, result)


def eligible(row, *, fixed=False):
    if not fixed and row["status"] != "UNKNOWN":
        return False
    radius = row["padded_radius" if fixed else "radius"]
    return row["lower"] <= 0. and (row["upper"] is None or row["upper"] > radius)


def examine(proposal, primary, object_name, observations, engine, cache, profile_costs):
    time0 = clock()
    try:
        y, u, mu, lam, radius, key = inputs(proposal, object_name, observations)
        previous = cache.get(key)
        if previous is not None:
            result = {
                label: previous[label]
                for label in ("lower", "upper", "nu", "nu_lower", "status", "passed", "native")
            }
            result["reused_recovery_id"] = previous["id"]
            charge = {unit: previous[f"charged_proof_{unit}"] for unit in ("wall", "cpu")}
        else:
            native = engine.score(y, u, mu, lam, radius=radius)
            result = {label: native[label] for label in ("lower", "upper", "nu", "nu_lower", "status")}
            result.update(passed=native["strictly_excludes"], native=native, reused_recovery_id=None)
            extra = profile_charge(result, profile_costs)
            spent = elapsed(time0)
            charge = {unit: spent[unit] + extra[unit] for unit in ("wall", "cpu")}
    except (KeyError, ValueError, TypeError, ArithmeticError, RuntimeError) as error:
        key = None
        result = dict(status="UNKNOWN", lower=0., upper=None, nu=None, nu_lower=0, passed=False,
                      native=None, reused_recovery_id=None, explanation=str(error))
        charge = elapsed(time0)
    spent = elapsed(time0)
    result.update(id=primary["id"], oid=object_name, arm=proposal["arm"], u=proposal["u"],
                  noise=proposal["noise"], side=proposal["side"], radius=proposal["radius"],
                  initial_status=primary["status"])
    for unit in ("wall", "cpu"):
        initial = primary[f"proof_{unit}"]
        initial_charge = primary.get(f"charged_proof_{unit}", initial)
        result.update({
            f"proof_{unit}": spent[unit],
            f"charged_proof_{unit}": charge[unit],
            f"initial_proof_{unit}": initial,
            f"total_proof_{unit}": initial + spent[unit],
            f"charged_total_proof_{unit}": initial_charge + charge[unit],
        })
    if key is not None:
        cache.setdefault(key, result)
    return result


def run(args):
    find_open(args.watch, 2.)
    from holdout import observations

    started = clock()
    results = PILOT / "results"
    if args.supplement:
        results = results / "supplement"
    target = results / "recovery"
    target.mkdir(parents=True, exist_ok=True)
    cases = {}
    for oid, noise, radius, y, caps in observations(supplement=args.supplement):
        cases.setdefault(oid, {})[noise] = (radius, np.ascontiguousarray(y, dtype=np.complex128), caps)
    rank, workers = args.shard
    names = sorted(cases)[rank::workers]
    if args.object:
        if not set(args.object) <= set(names):
            raise ValueError("Requested recovery object falls outside the shard.")
        names = [name for name in names if name in set(args.object)]
    scopes = names[0] if args.object and len(args.object) == 1 else f"{rank}-of-{workers}"
    progress_path = target / f"worker-{scopes}.json"
    receipts = {name: load_receipts(target / f"{name}.jsonl") for name in names}
    engine = RecoveryNative()
    profile_costs, reuse = {}, {}
    complete = set()
    while True:
        updated = 0
        counts = Counter()
        for name in names:
            if name in complete:
                continue
            proposals, _ = read_new_lines(results / "proposals" / f"{name}.jsonl", 0)
            primary, _ = read_new_lines(results / "native" / f"{name}.jsonl", 0)
            proposal_records = {row["id"]: row for row in proposals if "id" in row}
            for row in primary:
                if not eligible(row):
                    continue
                if row["id"] in receipts[name]:
                    continue
                if row["id"] not in proposal_records:
                    raise ValueError(f"Missing sealed weights for {row['id']}.")
                proposal = proposal_records[row["id"]]
                receipt = examine(proposal, row, name, cases[name], engine, reuse, profile_costs)
                append_receipt(target / f"{name}.jsonl", receipt)
                receipts[name][receipt["id"]] = receipt
                updated += 1
                print(receipt["id"], receipt["status"], f"{receipt['proof_wall']:.2f}s", flush=True)
            if (any(row.get("_complete") is True for row in proposals)
                    and {row["id"] for row in primary} == set(proposal_records)):
                complete.add(name)
        fixed_file = results / "fixed_baselines.json"
        fixed = json.loads(fixed_file.read_text()) if fixed_file.is_file() else {"status": "pending", "rows": []}
        fixed_rows = [row for row in fixed["rows"] if row["object_id"] in names]
        fixed_eligible = [row for row in fixed_rows if eligible(row, fixed=True)]
        # The exact normalized archived control coefficients are kept with their
        # original scorer. In this run every such row has a radius bracket or a
        # positive lower bound. Any future ambiguous row is explicitly visible.
        all_receipts = [row for name in names for row in receipts[name].values()]
        costs = {f"{prefix}proof_{unit}": sum(row[f"{prefix}proof_{unit}"] for row in all_receipts)
                 for prefix in ("", "charged_", "initial_", "total_", "charged_total_")
                 for unit in ("wall", "cpu")}
        counts.update(row["status"] for row in all_receipts)
        finished = len(complete) == len(names) and fixed["status"] == "complete" and not fixed_eligible
        atomic_json(progress_path, dict(supplement=args.supplement, shard=[rank, workers], objects=names,
                                        finished=finished, complete_objects=sorted(complete),
                                        attempted=len(all_receipts), status_counts=dict(counts),
                                        fixed_rows_under_same_rule=len(fixed_rows), fixed_remaining=fixed_eligible,
                                        timing=costs, setup=engine.cache_summary(), elapsed=elapsed(started)))
        if finished or not args.watch:
            return 0
        if len(complete) == len(names) and fixed["status"] == "complete":
            return 2
        if not updated:
            sleep(2.)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    stage = parser.add_mutually_exclusive_group(required=True)
    stage.add_argument("--training", action="store_true")
    stage.add_argument("--training-fourier-only", action="store_true")
    stage.add_argument("--run", action="store_true")
    parser.add_argument("--shard", type=shard_pair, default=(0, 1))
    parser.add_argument("--object", action="append")
    parser.add_argument("--supplement", action="store_true")
    parser.add_argument("--watch", action="store_true")
    args = parser.parse_args()
    if args.run:
        raise SystemExit(run(args))
    training_check(bordered=args.training)
