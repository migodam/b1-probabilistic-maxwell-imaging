"""Replay three frozen Paper 1 profiles against the 98 opened pilot conditions.

Run as a separate Python process with PYTHONDONTWRITEBYTECODE=1. This path
imports the original Paper 1 native verifier and uses the archived coefficient
convention, including exact products of the saved normalized multipliers.
"""

import argparse
import csv
from dataclasses import dataclass
from importlib.util import module_from_spec, spec_from_file_location
import json
from pathlib import Path
import platform
import sys
from time import perf_counter, process_time

import numpy as np
from flint import arb, ctx
from scipy.linalg import LinAlgError
from threadpoolctl import threadpool_limits

from holdout import observations


PILOT = Path(__file__).resolve().parents[1]
PAPER = PILOT.parent / "z_alg_optical_inertia"
SUBMISSION = PAPER / "submission"
CURVE = SUBMISSION / "curve"
SINGLE = SUBMISSION / "single"
FINE = PAPER / "synthetic" / "public" / "fine_geometry.npz"
ORIGINAL_ONLINE = PAPER / "validation" / "srap" / "online.json"
RESULTS = PILOT / "results"
CAPS = (0.5, 0.75, 1.0)

# Generic local module names in Paper 1 are intentionally isolated in this
# standalone process, so they do not resolve to the new pilot's native.py.
sys.dont_write_bytecode = True
sys.path.insert(0, str(CURVE))
from native import Geometry, Native, Weights
from envelope import FourierEnvelope
from fixed_blocks import SavedProfile


def _load_reference():
    spec = spec_from_file_location("_paper1_fixed_source", ORIGINAL_ONLINE.with_suffix(".py"))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


_reference = _load_reference()


def read_json(path):
    return json.loads(path.read_text())


def started():
    return perf_counter(), process_time()


def timing(start):
    return dict(zip(("wall_seconds", "cpu_seconds"), (x - y for x, y in zip(started(), start))))


class FrozenEnvelope(FourierEnvelope):
    """Use rounded weights to propose and exact normalized products to verify."""

    def select(self, weights, nu):
        trial = Weights(
            weights.cap, weights.scale, tuple(map(float, weights.real)), tuple(map(float, weights.optical))
        )
        return super().select(trial, nu)


@dataclass
class Profile:
    name: str
    cap: float
    specification: dict
    archived_inertia: dict
    index: int | None = None

    @property
    def block_side(self):
        return 4 if self.index is not None else 2

    def action(self, geometry):
        if self.index is not None:
            return SavedProfile(geometry, self.specification, self.cap)
        return Native(geometry, Weights.global_(self.cap, self.specification["scale"], 0.0))

    def verification_weights(self):
        if self.index is None:
            return Weights.global_(self.cap, self.specification["scale"], 0.0)
        scale = self.specification["mu"]
        real = tuple(arb(scale) * arb(v) for row in self.specification["sector"] for v in row)
        optical = tuple(arb(scale) * arb(v) for row in self.specification["optical"] for v in row)
        return Weights(self.cap, scale, real, optical)


def fixed_profiles():
    profiles = read_json(SINGLE / "coherent_general_profiles.json")
    old = {p["profile_index"]: p for p in read_json(SINGLE / "coherent_general_verified.json")["cases"].values()}
    highcap_rows = {
        p["profile_index"]: p
        for p in read_json(CURVE / "fixed_blocks_higher_caps.json")["cases"] if p["cap"] == 1.0
    }
    archived_global = {}
    for cap, name, scale in (
        (0.5, "global.5.json", 0.0025),
        (0.75, "global.75.json", 0.0004366418693307708),
    ):
        q = next(
            p for p in read_json(CURVE / name)["profiles"]
            if p["scale"] == scale and set(p["weights"]) == {scale} and set(p["optical"]) == {0.0}
        )
        if q["status"] != "native_gap_and_inertia_verified" or q["cap"] != cap:
            raise ValueError("Archived global profile does not match the specified fixed control.")
        archived_global[cap] = dict(
            method="inherited_global_Fourier_proof", source=str((CURVE / name).relative_to(PAPER)),
            nu=q["exact_nu"], nu_lower=q["exact_nu"], positive_gap=q["positive_gap_lower"],
            historical_wall_seconds=q["profile_wall_seconds"],
            historical_cpu_seconds=q["profile_cpu_seconds"],
        )

    global_negative = next(
        p for p in read_json(ORIGINAL_ONLINE)["profiles"] if p["label"] == "global_real_075"
    )
    key = global_negative["parameters"]
    scale = 0.0004366418693307708
    if key != Weights.global_(0.75, scale, 0.0).as_dict() or global_negative["nu"] != 4:
        raise ValueError("Archived negative gap cannot be transferred to the requested profile.")
    archived_global[0.75].update(
        negative_gap=global_negative["negative_proof"]["negative_gap_lower"],
        negative_gap_source=str(ORIGINAL_ONLINE.relative_to(PAPER)),
        additional_negative_historical_timing=global_negative["negative_proof"],
    )

    boxes = read_json(SUBMISSION / "highcaps" / "certificates.json")
    box = next(b for b in boxes["boxes"] if b["box"] == [0.5, 0.06, 0.0, 0.0])
    if boxes["cap_range"][0] != 1.0 or boxes["source_count"] != 10 or not (
        box["verified"] and box["max_projected_eigenvalue_upper"] < 0
    ):
        raise ValueError("The archived exact receiver-kernel box must include the unit global weights.")
    archived_global[1.0] = dict(
        method="inherited_global_member_of_receiver_kernel_box", nu=None, nu_lower=10,
        source="submission/highcaps/certificates.json", box=box["box"],
        negative_congruence_upper=box["max_projected_eigenvalue_upper"],
        historical_wall_seconds=boxes["wall_seconds"], historical_cpu_seconds=boxes["cpu_seconds"],
    )

    for cap in CAPS:
        scale = 0.0025 if cap == 0.5 else 0.0004366418693307708
        yield Profile(
            "global_real", cap, Weights.global_(cap, scale, 0.0).as_dict(), archived_global[cap],
        )
        for index in (0, 2):
            profile = profiles[index]
            proof = {}
            if cap == 0.75:
                q = old[index]
                if q["cap"] != cap or q["proved_K_scaled_gap"] <= 0:
                    raise ValueError("The frozen profile must carry the matching positive gap proof.")
                proof = dict(
                    method="inherited_fixed_16_Fourier_proof", nu=0, nu_lower=0,
                    positive_gap=q["proved_K_scaled_gap"], source="submission/single/coherent_general_verified.json",
                    historical_wall_seconds=q["wall_seconds"], historical_cpu_seconds=q["cpu_seconds"],
                )
            elif cap == 1.0:
                q = highcap_rows[index]
                if q["native_negative_lower_count"] < 10 or q["negative_congruence_upper"] >= 0:
                    raise ValueError("The archived fixed-profile high-cap proof is unresolved.")
                proof = dict(
                    method="inherited_fixed_16_negative_subspace", nu=None, nu_lower=10,
                    negative_congruence_upper=q["negative_congruence_upper"],
                    source="submission/curve/fixed_blocks_higher_caps.json",
                    historical_wall_seconds=q["wall_seconds"], historical_cpu_seconds=q["cpu_seconds"],
                )
            yield Profile(f"fixed_16_{index}", cap, profile, proof, index)


def native_geometry():
    with np.load(FINE, allow_pickle=False) as data:
        return data["spectrum"], data["PW360_10x20__S"], data["PW360_10x20__B"]


def condition_meta(object_id, noise, radius, cap):
    return dict(
        object_id=object_id, noise_pct=noise, padded_radius=radius, cap=cap,
        condition=f"{object_id}__n{noise:02d}__u{cap:.2f}",
    )


def apply_one(profile, proof, Y, meta, action=None, stationary=None):
    time0 = started()
    if proof.get("nu_lower", 0) >= 10:
        lower, upper, status, detail = 0.0, 0.0, "PROVED_ZERO", {}
    elif stationary is not None:
        m, defect = stationary
        detail = _reference.bound_source(
            Y, action, proof["nu"], proof["positive_gap"], proof.get("negative_gap"), m, defect, meta,
        )
        lower, upper = detail["joint_distance_lower"], detail["joint_distance_upper"]
        if not (0 <= lower <= upper and np.isfinite(upper)):
            raise ArithmeticError("The archived reusable proof produced inconsistent endpoints.")
        status = "CERTIFIED_POSITIVE" if lower > 0 else ("PROVED_ZERO" if upper == 0 else "UNKNOWN")
    else:
        lower, upper, status, detail = 0.0, None, "UNKNOWN", {}
    radius = meta["padded_radius"]
    return dict(
        **meta, profile_id=profile.name, status=status, lower=lower, upper=upper,
        strictly_excludes=lower > radius, verified_at_or_below_radius=upper is not None and upper <= radius,
        nu=proof.get("nu"), nu_lower=proof.get("nu_lower"), stationarity_source=detail, **timing(time0),
    )


def replay(profiles, geometry_arrays, data, progress):
    geometries, envelopes, shared, rows = {}, {}, [], []
    for p in profiles:
        time0 = started()
        proof = p.archived_inertia.copy()
        receipt = dict(profile_id=p.name, cap=p.cap, coefficient_convention=(
            "archived normalized binary64 products" if p.index is not None else "constant native binary64"
        ), specification=p.specification)
        action, stationary = None, None
        if proof.get("nu_lower", 0) < 10:
            side = p.block_side
            if side not in geometries:
                t = started()
                geometries[side] = Geometry(*geometry_arrays, block_side=side)
                receipt["partition_geometry"] = timing(t)
            geometry = geometries[side]
            action = p.action(geometry)
            before = geometry.cost.copy()
            try:
                # Timing a fresh frozen global proof gives a cold-start control.
                # The original proof remains a valid fallback if the new run is inconclusive.
                if not proof or (p.name == "global_real" and p.cap == 0.5):
                    t = started()
                    try:
                        if side not in envelopes:
                            envelopes[side] = FrozenEnvelope(geometry)
                        outward = envelopes[side].verify(p.verification_weights(), 0)
                        if proof:
                            receipt["inherited_inertia"] = proof
                        proof = dict(
                            method="same_machine_Fourier_proof", nu=0, nu_lower=0,
                            positive_gap=outward["positive_gap_lower"], full_proof=outward,
                        )
                    except (ArithmeticError, LinAlgError, ValueError) as error:
                        receipt["fresh_gap_failure"] = str(error)
                    finally:
                        receipt["fresh_gap_attempt"] = timing(t)
                if proof.get("positive_gap") is not None:
                    m, defect, work = _reference.exact_common_response(action, proof["nu"])
                    stationary = (m, defect)
                    receipt["shared_stationary_work"] = work
            except (ArithmeticError, LinAlgError, OverflowError, ValueError) as error:
                receipt["stationarity_failure"] = str(error)
            receipt["operations"] = dict(geometry.cost - before)
        receipt.update(inertia_proof=proof, shared_timing=timing(time0))
        local = []
        for object_id, noise, radius, Y, caps in data:
            if p.cap not in caps:
                continue
            meta = condition_meta(object_id, noise, radius, p.cap)
            try:
                local.append(apply_one(p, proof, Y, meta, action, stationary))
            except (ArithmeticError, LinAlgError, OverflowError, ValueError) as error:
                row = apply_one(p, proof, Y, meta)
                row["source_failure"] = str(error)
                local.append(row)
        rows.extend(local)
        shared.append(receipt)
        progress(shared, rows)
        print(
            p.name, p.cap, proof.get("method", "UNKNOWN"),
            sum(r["strictly_excludes"] for r in local), "/", len(local),
            "shared wall", round(receipt["shared_timing"]["wall_seconds"], 4),
            "online median", round(float(np.median([r["wall_seconds"] for r in local])), 6),
            flush=True,
        )
    return shared, rows


def portfolio(rows):
    groups = {}
    for row in rows:
        groups.setdefault(row["condition"], []).append(row)
    united = []
    for group in groups.values():
        exemplar = group[0]
        best = max(group, key=lambda r: r["lower"])
        upper = max(r["upper"] for r in group) if all(r["upper"] is not None for r in group) else None
        united.append(dict(
            **{k: exemplar[k] for k in ("condition", "object_id", "noise_pct", "cap", "padded_radius")},
            lower=best["lower"], upper=upper, best_profile=best["profile_id"],
            strictly_excludes=best["strictly_excludes"],
            verified_at_or_below_radius=upper is not None and upper <= exemplar["padded_radius"],
            online_wall_seconds=sum(r["wall_seconds"] for r in group),
            online_cpu_seconds=sum(r["cpu_seconds"] for r in group),
        ))
    return united


def crosscheck_history(rows):
    keyed = {(r["profile_id"], r["object_id"], r["noise_pct"], r["cap"]): r for r in rows}
    archived_sets = []
    for source in read_json(ORIGINAL_ONLINE)["profiles"]:
        labels = {
            "global_real_050": "global_real", "global_real_075": "global_real",
            "fixed_4x4_0_075": "fixed_16_0", "fixed_4x4_2_075": "fixed_16_2",
        }
        if source["label"] in labels:
            arch = [(r["object_id"], 10, source["cap"], r["padded_radius"], r["joint_distance_lower"],
                     r["joint_distance_upper"]) for r in source["rows"]]
            archived_sets.append((labels[source["label"]], arch))
    for source in read_json(CURVE / "fixed_blocks.75.json")["profiles"]:
        arch = [(r["object_id"], r["noise_pct"], r["cap"], r["radius"], r["distance_lower"],
                 r["distance_upper"]) for r in source["rows"] if r["noise_pct"] in (1, 2, 10)]
        archived_sets.append((f"fixed_16_{source['profile_index']}", arch))

    checks = []
    for label, archived in archived_sets:
        for object_id, noise, cap, radius, lo, hi in archived:
            new = keyed[label, object_id, noise, cap]
            checks.append(dict(
                profile_id=label, condition=new["condition"], same_radius=new["padded_radius"] == radius,
                intervals_overlap=new["upper"] is not None and new["lower"] <= hi and lo <= new["upper"],
                strict_outcome_agrees=new["strictly_excludes"] == (lo > radius),
            ))
    return dict(
        compared=len(checks), radii_agree=sum(r["same_radius"] for r in checks),
        outward_intervals_overlap=sum(r["intervals_overlap"] for r in checks),
        strict_outcomes_agree=sum(r["strict_outcome_agrees"] for r in checks),
        differences=[r for r in checks if not all(
            r[k] for k in ("same_radius", "intervals_overlap", "strict_outcome_agrees")
        )],
    )


def save(out, doc):
    tmp = out.with_suffix(".tmp")
    tmp.write_text(json.dumps(doc, indent=2, allow_nan=False) + "\n")
    tmp.replace(out)


def supplement_observations():
    bundle = PILOT / "data" / "srap_supplement.npz"
    if not bundle.with_suffix(".json").is_file():
        raise ValueError("The supplemental native-radius replay receipt is incomplete.")
    with np.load(bundle, allow_pickle=False) as data:
        if set(data.files) != {"object_id", "noise_pct", "Y", "radius"}:
            raise ValueError("Supplement may contain only ids, noisy fields and radii.")
        ids, noise, y, radii = (data[k] for k in ("object_id", "noise_pct", "Y", "radius"))
    expected = [f"s2d_{i:03d}" for i in _reference.SCOPE]
    if ids.tolist() != expected or noise.tolist() != [1, 2] or y.shape != (22, 2, 20, 10) or (
        radii.shape != (22, 2) or not np.isfinite(radii).all() or np.any(radii <= 0)
        or not np.isfinite(y).all()
    ):
        raise ValueError("Supplement does not match the SRAP22 paired noise contract.")
    return [
        (str(oid), int(pct), float(radii[i, j]), y[i, j], (0.75,))
        for i, oid in enumerate(ids) for j, pct in enumerate(noise)
    ]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--supplement", action="store_true", help="Replay only the supplemental SRAP22 bands at 0.75.")
    args = parser.parse_args()
    marker = RESULTS / "holdout_open.json"
    if not marker.is_file():
        raise SystemExit("Holdout marker is absent. No observations were opened.")
    t = started()
    ctx.prec, ctx.threads = 112, 2
    profiles = list(fixed_profiles())
    arrays = native_geometry()
    if args.supplement:
        profiles = [p for p in profiles if p.cap == 0.75]
        data = supplement_observations()
        expected = (22, 44, 44)
        destination = RESULTS / "supplement" / "fixed_baselines.json"
    else:
        data = list(observations())
        expected = (28, 40, 98)
        destination = RESULTS / "fixed_baselines.json"
    out_path = args.output or destination
    out_path.parent.mkdir(parents=True, exist_ok=True)
    case_count = sum(len(caps) for _, _, _, _, caps in data)
    object_count = len({item[0] for item in data})
    if (object_count, len(data), case_count) != expected:
        raise ValueError("The opened data do not match the specified fixed baseline conditions.")

    out = dict(
        contract="Three frozen archived profiles, no proposal selection from observation values",
        shared_ARb_precision_bits=ctx.prec, interval_threads=ctx.threads, BLAS_threads_requested=1,
        machine=platform.platform(), observations=len(data), condition_count=case_count,
        object_count=object_count, supplemental_noise_sensitivity=args.supplement,
        setup=timing(t), status="in_progress",
    )

    def update(shared, rows):
        out.update(profiles=shared, rows=rows, **timing(t))
        save(out_path, out)

    with threadpool_limits(limits=1):
        shared, rows = replay(profiles, arrays, data, update)
    aggregate = portfolio(rows)
    if len(aggregate) != case_count or len(rows) != 3 * case_count:
        raise ValueError("Three scores per pre-specified condition are required.")
    out.update(status="complete", portfolio=aggregate)
    if not args.supplement:
        out["archived_crosscheck"] = crosscheck_history(rows)
    update(shared, rows)
    with out_path.with_suffix(".csv").open("w", newline="") as file:
        fields = [
            "condition", "object_id", "noise_pct", "cap", "padded_radius", "profile_id",
            "nu", "nu_lower", "status", "lower", "upper", "strictly_excludes",
            "verified_at_or_below_radius", "wall_seconds", "cpu_seconds",
        ]
        writer = csv.DictWriter(file, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
