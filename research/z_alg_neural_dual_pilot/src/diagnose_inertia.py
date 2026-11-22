"""Certify negative state directions for the sealed optfresh multiplier profiles."""

import json
from collections import Counter
from itertools import product
from pathlib import Path
from statistics import fmean
from time import perf_counter, process_time

import numpy as np
from flint import acb_mat, arb
from scipy.linalg import eigh
from threadpoolctl import threadpool_info

from native import Native, RegionalWitness, Weights, ball, gershgorin, hermitian, high, low
from native import _paper, star, timed
from score_batch import atomic_json, read_new_lines


PILOT = Path(__file__).resolve().parents[1]
RESULTS = PILOT / "results"
REFERENCES = "../z_alg_optical_inertia/submission/single/coherent_general_profiles.json"
OBJECTS = tuple(f"optfresh_{number:02}" for number in range(1, 7))
ARMS = ("neural", "no_feature", "hybrid", "gradient1", "gradient5")
CAPS = (1.0, 0.75)
NOISES = (1, 2, 10)
SIDES = (4, 8)


def profile_weights(cap, mu, lam):
    mu = np.asarray(mu, dtype=float).ravel()
    lam = np.asarray(lam, dtype=float).ravel()
    if mu.shape != lam.shape or not np.isfinite(mu).all() or not np.isfinite(lam).all():
        raise ValueError("Expected finite sector and optical weights with the same shape.")
    if np.min(mu) <= 0.0:
        raise ValueError("Sector weights must be positive.")
    return Weights(cap, float(np.mean(mu)), tuple(map(float, mu)), tuple(map(float, lam)))


def fixed_basis(geometry):
    """Build the same 24 directions before reading any of the sealed proposals."""
    archive = json.loads((PILOT / REFERENCES).read_text())
    references = [
        ("uniform_coarse_seed", 14, profile_weights(1.0, [0.001 / 9] * 16, [0.0005 / 9] * 16)),
    ]
    for index in (0, 2):
        old = archive[index]
        scale = old["mu"]
        weights = profile_weights(
            1.0, scale * np.asarray(old["sector"]), scale * np.asarray(old["optical"]),
        )
        references.append((f"archived_fixed_16_index_{index}", 5, weights))

    directions, receipts = [], []
    for name, count, weights in references:
        action = _paper.Native(geometry, weights)
        search = action.modes(count)
        if not np.all(action.eigs[:count] < 0):
            raise ArithmeticError(f"Reference {name} did not suggest {count} negative directions.")
        directions.append(action.eig_basis[:, :count])
        receipts.append(dict(name=name, candidate_columns=count, profile=weights.as_dict(), **search))

    # QR normalizes the joint space. Fix each column's phase via its largest entry.
    basis, triangular = np.linalg.qr(np.column_stack(directions), mode="reduced")
    if min(abs(triangular.diagonal())) < 1.0e-10:
        raise ArithmeticError("Fixed reference vectors are numerically dependent.")
    pivots = basis[np.argmax(abs(basis), axis=0), np.arange(basis.shape[1])]
    basis *= pivots.conj() / abs(pivots)
    return basis, receipts


class GeometryWitness(RegionalWitness):
    """Exact fine-grid compression, with 4×4 cells assembled from the 8×8 cells."""

    def __init__(self, geometry, basis):
        super().__init__(geometry, basis, "fixed_geometry_24")
        self.regions = {8: self.regional}
        self.regions[4] = []
        for row, column in np.ndindex(4, 4):
            start = 16 * row + 2 * column
            a, b, c, d = (self.regional[start + offset] for offset in (0, 1, 8, 9))
            self.regions[4].append(tuple(x + y + z + w for x, y, z, w in zip(a, b, c, d)))

    def certify(self, weights, side):
        cap, scale = arb(weights.cap), arb(weights.scale)
        projected = self.receiver
        for mu, lam, (gram, sector, optical) in zip(
            weights.real, weights.optical, self.regions[side], strict=True,
        ):
            projected = projected + arb(mu) * (gram - cap * sector) + arb(lam) * optical
        projected = projected / scale

        values, rotation = eigh(hermitian(_paper.mid(projected)), check_finite=False)
        q = ball(rotation)
        rotated = star(q) * (projected * q)
        for count in range(len(values), 0, -1):
            block = acb_mat(count, count, [rotated[i, j] for i in range(count) for j in range(count)])
            edge = gershgorin(block, upper=True)
            if edge < 0:
                selected = ball(rotation[:, :count])
                gram = star(selected) * (self.gram * selected)
                gram_upper = gershgorin(gram, upper=True)
                return dict(
                    nu_lower=count, negative_congruence_upper_scaled=high(edge),
                    negative_gap_lower_scaled=low(-edge / gram_upper),
                    projected_midpoint_negative_count=int(np.count_nonzero(values < 0)),
                    tenth_projected_midpoint_eigenvalue_scaled=float(values[9]),
                )
        return dict(
            nu_lower=0, negative_congruence_upper_scaled=None, negative_gap_lower_scaled=None,
            projected_midpoint_negative_count=int(np.count_nonzero(values < 0)),
            tenth_projected_midpoint_eigenvalue_scaled=float(values[9]),
        )


def proposal_rows():
    for oid in OBJECTS:
        path = RESULTS / "proposals" / f"{oid}.jsonl"
        rows, offset = read_new_lines(path, 0)
        if offset != path.stat().st_size or not rows or not rows[-1].get("_complete"):
            raise ValueError(f"Proposal file {path.name} is not complete.")
        indexed = {(row["u"], row["noise"], row["side"], row["arm"]): row for row in rows[:-1]}
        if len(indexed) != len(rows) - 1:
            raise ValueError(f"Proposal file {path.name} has duplicate profiles.")
        for cap, noise, side, arm in product(CAPS, NOISES, SIDES, ARMS):
            row = indexed[cap, noise, side, arm]
            row_id = f"{oid}__n{noise}__u{cap:g}__b{side}__{arm}"
            if row["oid"] != oid or row["id"] != row_id:
                raise ValueError(f"Proposal ID is inconsistent with {path.name}.")
            yield row


def attach_prior_counts(receipts):
    prior = {}
    for oid in OBJECTS:
        rows, _ = read_new_lines(RESULTS / "native" / f"{oid}.jsonl", 0)
        for row in rows:
            if "id" in row:
                prior[row["id"]] = row
    for row in receipts:
        old = prior.get(row["id"], {})
        tiled = (old.get("upper_check") or {}).get("nu_lower")
        native = (old.get("native") or {}).get("proof") or {}
        row.update(
            first_pass_available=bool(old), first_pass_tile_state_nu_lower=tiled,
            first_pass_full_state_nu_lower=native.get("nu_lower"),
            first_pass_full_state_nu_exact=native.get("nu"),
        )


def counts_summary(rows, source_count):
    lower = [row["nu_lower"] for row in rows]
    coarse = [row["coarse_nu"] for row in rows if row["coarse_nu"] is not None]
    tiled = [row["first_pass_tile_state_nu_lower"] for row in rows
             if row["first_pass_tile_state_nu_lower"] is not None]
    old = [max(row["first_pass_tile_state_nu_lower"] or 0, row["first_pass_full_state_nu_lower"] or 0)
           for row in rows if row["first_pass_available"]]
    return dict(
        profiles=len(rows), native_nu_lower_min=min(lower), native_nu_lower_mean=fmean(lower),
        native_nu_lower_max=max(lower), native_nu_lower_histogram=dict(sorted(Counter(lower).items())),
        native_proven_exhausted=sum(n >= source_count for n in lower),
        native_exhaustion_unresolved=sum(n < source_count for n in lower),
        coarse_count=len(coarse), coarse_nu_mean=fmean(coarse) if coarse else None,
        coarse_nu_histogram=dict(sorted(Counter(coarse).items())),
        coarse_count_at_least_sources=sum(n >= source_count for n in coarse),
        first_pass_receipts=len(old), first_pass_tile_count=len(tiled),
        first_pass_tile_nu_lower_mean=fmean(tiled) if tiled else None,
        first_pass_tile_nu_lower_histogram=dict(sorted(Counter(tiled).items())),
        first_pass_proven_exhausted=sum(n >= source_count for n in old),
    )


def record_profile(row, witness, source_count):
    side = row["side"]
    mu, lam = [np.asarray(row[name], dtype=np.float64) / 9.0 for name in ("mu", "lam")]
    if mu.shape != (side * side,) or lam.shape != mu.shape:
        raise ValueError(f"Malformed profile {row['id']}.")
    weight = profile_weights(row["u"], mu, lam)
    stages = {}
    with timed(stages, "native_projected_interval_proof"):
        proof = witness.certify(weight, side)
    return {
        **{name: row[name] for name in ("id", "oid", "noise", "u", "side", "arm")},
        "observation_id": f"{row['oid']}__pw360_10x20__dual{row['noise']}pct",
        "coarse_nu": row["nu"],
        "coarse_proposal_wall_seconds": row["proposal_wall"],
        "coarse_proposal_cpu_seconds": row["proposal_cpu"],
        "coarse_diagnostic_wall_seconds": row.get("diagnostic_wall"),
        "coarse_diagnostic_cpu_seconds": row.get("diagnostic_cpu"),
        "native_weights": weight.as_dict(),
        **proof,
        "proven_source_exhaustion": proof["nu_lower"] >= source_count,
        "timing": stages,
    }


def timing_summary(rows):
    timing = {}
    for unit in ("wall_seconds", "cpu_seconds"):
        native = [r["timing"]["native_projected_interval_proof"][unit] for r in rows]
        proposal = [r[f"coarse_proposal_{unit}"] for r in rows]
        diagnostic = [r[f"coarse_diagnostic_{unit}"] for r in rows
                      if r[f"coarse_diagnostic_{unit}"] is not None]
        timing[unit] = dict(
            native_projected_total=sum(native), native_projected_mean=fmean(native),
            archived_coarse_proposal_total=sum(proposal), archived_coarse_proposal_mean=fmean(proposal),
            archived_coarse_diagnostic_total=sum(diagnostic),
            archived_coarse_diagnostic_mean=fmean(diagnostic) if diagnostic else None,
            archived_coarse_diagnostic_count=len(diagnostic),
        )
    return timing


def main():
    started = perf_counter(), process_time()
    precision, threads = 112, 1
    verifier = Native(precision=precision, threads=threads)
    setup = verifier.setup
    with verifier._precision():
        with timed(setup, "native_partition_geometry"):
            geometry4 = verifier._geometry(4, {})
            geometry8 = verifier._geometry(8, {})
        with timed(setup, "fixed_geometry_basis"):
            basis, references = fixed_basis(geometry4)
        with timed(setup, "exact_fourier_and_regional_congruences"):
            witness = GeometryWitness(geometry8, basis)
        with timed(setup, "read_completed_proposal_files"):
            proposals = list(proposal_rows())
        receipts = [record_profile(row, witness, verifier.sources) for row in proposals]
        with timed(setup, "read_existing_state_counts"):
            attach_prior_counts(receipts)
        blas = [{"internal_api": lib["internal_api"], "num_threads": lib["num_threads"]}
                for lib in threadpool_info() if lib["user_api"] == "blas"]
        source_count = verifier.sources
        summaries = {}
        for cap in CAPS:
            group = [row for row in receipts if row["u"] == cap]
            summaries[str(cap)] = dict(
                overall=counts_summary(group, source_count),
                by_partition={str(side**2): counts_summary(
                    [row for row in group if row["side"] == side], source_count) for side in SIDES},
                by_arm={arm: counts_summary(
                    [row for row in group if row["arm"] == arm], source_count) for arm in ARMS},
                timing=timing_summary(group),
            )
        payload = dict(
            test="fixed-reference geometry-only native state inertia",
            operator="archived 96² PW360_10x20 binary64 sensor and complete 192² Green symbol",
            scope="optfresh 01–06; paired 1%, 2%, 10%; 16 and 64 blocks; five sealed dynamic arms",
            precision_bits=precision, arb_threads=threads, blas=blas, native_dim=geometry8.n,
            basis_dim=basis.shape[1], source_count=source_count,
            basis_selection="before reading proposals; no Y, target labels, or prior outcomes",
            archived_reference_source=REFERENCES, reference_candidates=references,
            reference_convention="uniform coarse seed divided by nine; archived 16-block profiles already native",
            assessed_weight_convention="both stored coarse weights divided by nine in binary64",
            certification="outward ACB regional compression; strict upper Gershgorin on rotated prefix",
            exhaustion_rule="at least ten strict native K negative directions imply distance zero",
            lower_bound_semantics="a count below ten leaves source exhaustion unresolved",
            historical_coarse_semantics="recorded 32² numerical inertia, not an interval proof",
            historical_native_semantics="first-pass state-only fields read after independent proof",
            optional_cap_semantics="0.75 reuses the 1.0 basis and only diagnoses possible state subspace deficiency",
            setup=setup, basis_search_operations=dict(geometry4.cost),
            exact_witness_operations=dict(geometry8.cost), summaries=summaries, rows=receipts,
            run_wall_seconds=perf_counter() - started[0], run_cpu_seconds=process_time() - started[1],
        )
    atomic_json(RESULTS / "geometry_inertia.json", payload)
    for cap, summary in summaries.items():
        print(cap, json.dumps(summary["overall"], allow_nan=False), flush=True)


if __name__ == "__main__":
    main()
