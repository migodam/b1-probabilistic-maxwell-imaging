"""Verify high-contrast pixel counts on two saved scalar-TM observations."""

import argparse
import json
from pathlib import Path
import sys

import numpy as np
from flint import arb, ctx
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parents[1]
OPTICAL = RESEARCH / "z_alg_optical_inertia"
ARCHIVE = RESEARCH / "z_alg_secant_tasks/highdim"
for path in (ARCHIVE, OPTICAL / "core", OPTICAL / "fine"):
    sys.path.insert(0, str(path))

from certificate import ball, certify, frobenius_squared, safe_float
from inertia import Prepared, clock, timing
from replay import OBSERVATION_RADII, OPERATORS, SOURCE, load_operators
from verify_fourier import ExactGeometry, grid, one_y, profile_proof

RING_ID = "r10_twod_0001"
NATIVE_ID = "optfresh_01__pw360_10x20__n02"
RING_THRESHOLD = 4.0
NATIVE_KEY = ("PW360_10x20", 12, 0.5, 0.0, 0.001, 0.78125)
RING_PROFILES = (
    {"id": "archived_peak_scalar", "mu": 0.01, "lambda": -0.01},
    {"id": "half_percent_scalar", "mu": 5e-6, "lambda": -2e-5},
    {"id": "cardinality1_scalar", "mu": 4.571871609296927e-5, "lambda": -6.845644029318437e-5},
)
PRECISION, ARB_THREADS = 120, 2


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def save_json(path, data):
    temporary = path.with_suffix(".partial")
    temporary.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def as_path(path):
    return str(path.relative_to(RESEARCH))


def lower(x):
    return max(0.0, safe_float(x, upper=False))


def upper(x):
    return 0.0 if x == 0 else safe_float(x, upper=True)


def summarize_tails(sources, profiles, budgets):
    """Envelope over verified proposals. Each sum and comparison is outward rounded."""
    rows = []
    for q in range(sources + 1):
        bounds = {
            profile["id"]: sum((arb(v["lower"]) for v in profile["eigenvalue_lowers"]
                                if v["index_descending_one_based"] > profile["nu"] + q), arb(0))
            for profile in profiles
            if profile["status"] == "Arb_verified"
        }
        if not bounds:
            raise ArithmeticError("No rigorously verified profile was available.")
        identity, best_square = max(bounds.items(), key=lambda item: item[1].lower())
        decisions = {}
        for label, radius in budgets.items():
            if not (radius >= 0 and radius.is_finite()):
                raise ValueError("Every noise radius must be finite and nonnegative.")
            decisions[label] = bool(best_square > radius**2)
        rows.append(dict(
            maximum_high_pixels_hypothesized=q, squared_residual_lower=lower(best_square),
            residual_lower=lower(best_square.sqrt()), winning_profile=identity,
            strictly_excluded=decisions,
        ))
    return rows


def result_metrics(n, s, profiles, budgets):
    rows = summarize_tails(s, profiles, budgets)
    area_target_count = (n + 199) // 200
    results = {}
    for name, radius in budgets.items():
        count = max((q["maximum_high_pixels_hypothesized"] + 1 for q in rows
                     if q["strictly_excluded"][name]), default=0)
        results[name] = dict(
            radius_lower=lower(radius), radius_upper=upper(radius),
            proved_high_pixels=count, proved_pixel_fraction=count / n,
            strongest_direct_rank_ceiling=s, host_grid_pixels=n,
            half_percent_area_target_pixels=area_target_count,
            proves_half_percent_area=bool(count >= area_target_count),
            direct_rank_method_can_reach_half_percent_area=bool(s >= area_target_count),
        )
    return dict(budgets=results, tails=rows)


def ring_data():
    reader = load_operators(RING_ID)
    path = SOURCE / f"{RING_ID}.npz"
    with np.load(path, allow_pickle=False) as stored:
        return Prepared(reader["G"], reader["S"], reader["E0"], stored["Y_archive"]), stored["positions"]


def solve_ring(extra_profiles):
    started = clock()
    model, positions = ring_data()
    observed_norm = frobenius_squared(ball(model.Y)).sqrt()
    radius = read_json(OBSERVATION_RADII)[RING_ID]["frozen_observation_ball"]
    budgets = {
        "archived_cn99": arb(radius),
        "hypothetical_0.5pct_of_saved_observation_norm": observed_norm / 200,
        "hypothetical_noiseless_at_same_center": arb(0),
    }

    profiles = []
    for proposal in (*RING_PROFILES, *extra_profiles):
        mu, lam = proposal["mu"], proposal["lambda"]
        if "partition_side" in proposal:
            side = proposal["partition_side"]
            frac = (positions - positions.min(axis=0)) / np.ptp(positions, axis=0)
            block = np.minimum((frac * side).astype(int), side - 1)
            order = side * block[:, 1] + block[:, 0]
            mu, lam = np.asarray(mu)[order], np.asarray(lam)[order]
        proof = certify(model, RING_THRESHOLD, mu, lam, precision=PRECISION, threads=ARB_THREADS)
        profiles.append(dict(
            id=proposal["id"], proposed_multipliers=proposal,
            nu=proof.get("verified_negative_eigenvalues"),
            eigenvalue_lowers=proof.get("retained_H_eigenvalue_lower_bounds", []),
            **proof,
        ))
    metrics = result_metrics(model.n, model.s, profiles, budgets)
    return dict(
        model="original dense 24×24 matched scalar TM finite arrays",
        provenance="MAX-STRUCT R10 archive, exposed task case",
        receiver_count=model.r, illumination_count=model.s, threshold=RING_THRESHOLD,
        geometry_path=as_path(OPERATORS), observation_path=as_path(SOURCE / f"{RING_ID}.npz"),
        radius_path=as_path(OBSERVATION_RADII),
        observed_residual_Frobenius_normalization_upper=upper(observed_norm),
        profiles=profiles, **metrics, **timing(started),
    )


def solve_native():
    started = clock()
    geometry, cut, threshold, optical, mu, floor = NATIVE_KEY
    public = OPTICAL / "synthetic/public"
    geo = public / "fine_geometry.npz"
    obs = public / f"{NATIVE_ID}.npz"
    feasibility = OPTICAL / "results/theory/fine_feasibility.json"
    with np.load(geo, allow_pickle=False) as data:
        spec, sensor, incident = (data[key] for key in ("spectrum", f"{geometry}__S", f"{geometry}__B"))
    with np.load(obs, allow_pickle=False) as data:
        y = data["Y"]
    radii = next(row for row in read_json(feasibility)["observations"] if row["id"] == NATIVE_ID)
    budgets = {
        "archived_2pct_with_generation_padding": arb(radii["inflated_radius"]),
        "archived_2pct_raw": arb(radii["sample_radius_binary64"]),
        "hypothetical_noiseless_at_same_center": arb(0),
    }

    geom_started = clock()
    q = grid(cut, len(spec))
    exact = ExactGeometry(sensor, incident, q, len(spec), cut)
    setup = timing(geom_started)
    u, nu, proof = profile_proof(exact, spec, q, len(sensor), incident.shape[1],
                                 threshold, optical, mu, floor)
    lower_started = clock()
    residual_lower, spectrum = one_y(u, nu, y, threshold, optical, mu)
    profiles = [dict(
        id="archived_native_proof", status="Arb_verified", nu=nu,
        profile_parameters=dict(geometry=geometry, cutoff=cut, threshold=threshold,
                                optical_over_mu=optical, mu=mu, symbol_floor=floor),
        eigenvalue_lowers=[dict(index_descending_one_based=i, lower=v)
                           for i, v in enumerate(spectrum, start=nu + 1)],
        peak_residual_lower=residual_lower, proof=proof, source_proof_cost=timing(lower_started),
    )]
    metrics = result_metrics(incident.shape[0], incident.shape[1], profiles, budgets)
    return dict(
        model="original native 96×96 scalar TM, exact saved Fourier multiplier",
        provenance="Codex/Sol archived optfresh forked-block shape; exposed follow-up",
        receiver_count=sensor.shape[0], illumination_count=incident.shape[1], threshold=threshold,
        geometry_path=as_path(geo), observation_path=as_path(obs), radius_path=as_path(feasibility),
        setup=setup, profiles=profiles, **metrics, **timing(started),
    )


def run(output, extra):
    started = clock()
    ctx.prec, ctx.threads = PRECISION, ARB_THREADS
    ring = solve_ring(extra)
    print(RING_ID, {k: v["proved_high_pixels"] for k, v in ring["budgets"].items()}, flush=True)
    native = solve_native()
    print(NATIVE_ID, {k: v["proved_high_pixels"] for k, v in native["budgets"].items()}, flush=True)
    receipt = dict(
        scope="strict lower counts of cells above the threshold; public observed arrays and physical metadata",
        precision_bits=PRECISION, interval_threads=ARB_THREADS, numeric_threads=1,
        material_forward_solves=0,
        cases={RING_ID: ring, NATIVE_ID: native}, **timing(started),
    )
    save_json(output, receipt)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE / "area_receipt.json")
    parser.add_argument("--extra-ring-profiles", type=Path)
    args = parser.parse_args()
    with threadpool_limits(1):
        run(args.output, read_json(args.extra_ring_profiles)["profiles"] if args.extra_ring_profiles else ())
