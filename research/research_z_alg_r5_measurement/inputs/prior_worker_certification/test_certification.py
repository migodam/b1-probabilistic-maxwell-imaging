#!/usr/bin/env python3
"""Physical 4/6-cell regression tests for the local interval certifier.

Run with:
    Gaussian/.venv_nn/bin/python research/delegated/z_alg_r5_certification/test_certification.py

The fixtures are read-only historical R2 files.  The test writes only
``tests.json`` beside this script.
"""
from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path
import sys
import time
from typing import Any

import mpmath as mp
import numpy as np


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[2]
R2 = (
    PROJECT
    / "research/delegated/z_alg_r2_publication/repository/research/research_z_alg_r2"
)
sys.path.insert(0, str(R2 / "src"))

from physics import Model  # noqa: E402
from certification import certify  # noqa: E402


def _json_counts(value: Any) -> dict[str, int]:
    return {str(key): int(item) for key, item in value.items()}


def _json_timings(value: Any) -> dict[str, float]:
    return {str(key): float(item) for key, item in value.items()}


def _run_case(
    *,
    name: str,
    contract_path: Path,
    observation_path: Path,
    result_path: Path,
    root_index: int,
) -> dict[str, Any]:
    wall_start, cpu_start = time.perf_counter(), time.process_time()
    contract = json.loads(contract_path.read_text())
    archived = json.loads(result_path.read_text())
    target = np.load(observation_path)["y"]
    root = np.asarray(archived["independently_validated"][root_index]["m"], dtype=np.float64)

    setup_start, setup_cpu_start = time.perf_counter(), time.process_time()
    model = Model(None, **contract["geometry"])
    setup_wall = time.perf_counter() - setup_start
    setup_cpu = time.process_time() - setup_cpu_start
    counts_before = dict(model.counts)
    setup_counts = _json_counts(counts_before)
    setup_timings = _json_timings(model.timings)

    forward_start, forward_cpu_start = time.perf_counter(), time.process_time()
    state = model.forward(root, jac=True, label=f"r5-cert-test-{name}")
    forward_wall = time.perf_counter() - forward_start
    forward_cpu = time.process_time() - forward_cpu_start
    forward_counts = {
        key: int(model.counts[key] - counts_before.get(key, 0))
        for key in model.counts
    }
    response_relative = float(
        np.linalg.norm(state["y"] - target)
        / max(np.linalg.norm(target), np.finfo(float).tiny)
    )

    certificate = certify(
        model.G,
        model.E0,
        model.S,
        target,
        root,
        box=tuple(contract["material_box"]),
        timeout=30.0,
    )
    json.dumps(certificate, allow_nan=False)
    success = certificate["status"] == "CERTIFIED_LOCAL_REAL_ROOT"
    strict = bool(
        success
        and certificate["attempts"][-1]["strict_inclusion"]
        and certificate["attempts"][-1]["contraction_strict"]
    )
    jacobian_contained = False
    if certificate["attempts"] and "J_interval" in certificate["attempts"][0]:
        J_enclosure = certificate["attempts"][0]["J_interval"]
        jacobian_contained = all(
            J_enclosure[row][col]["lo"]
            <= float(state["J"][row, col])
            <= J_enclosure[row][col]["hi"]
            for row in range(state["J"].shape[0])
            for col in range(state["J"].shape[1])
        )
    return {
        "name": name,
        "n": int(model.n),
        "s": int(model.ns),
        "contract": str(contract_path.relative_to(PROJECT)),
        "observation": str(observation_path.relative_to(PROJECT)),
        "archived_results": str(result_path.relative_to(PROJECT)),
        "archived_root_index": root_index,
        "archived_m": root.tolist(),
        "fixture_forward": {
            "response_relative": response_relative,
            "state_relative": float(state["state_relative_residual"]),
            "resolvent_sigma_min": float(state["resolvent_sigma_min"]),
            "counts": forward_counts,
            "cpu_seconds": forward_cpu,
            "wall_seconds": forward_wall,
        },
        "physical_setup": {
            "counts": setup_counts,
            "timings": setup_timings,
            "cpu_seconds": setup_cpu,
            "wall_seconds": setup_wall,
        },
        "model_total_counts_after_fixture_forward": _json_counts(model.counts),
        "assertions": {
            "fresh_forward_matches_archived_target": response_relative <= 1.0e-10,
            "certified": success,
            "strict_inclusion_and_contraction": strict,
            "R2_real_jacobian_contained_by_interval_J": jacobian_contained,
            "json_serializable": True,
        },
        "certificate": certificate,
        "cpu_seconds": time.process_time() - cpu_start,
        "wall_seconds": time.perf_counter() - wall_start,
    }


def _negative_control(
    *,
    reference_case: dict[str, Any],
    contract_path: Path,
    observation_path: Path,
    result_path: Path,
) -> dict[str, Any]:
    """Perturb the target so the archived material must not certify."""
    wall_start, cpu_start = time.perf_counter(), time.process_time()
    contract = json.loads(contract_path.read_text())
    archived = json.loads(result_path.read_text())
    target = np.load(observation_path)["y"].copy()
    root = np.asarray(archived["independently_validated"][0]["m"], dtype=np.float64)
    target[0, 0] += 1.0e-3 + 2.0e-3j
    model = Model(None, **contract["geometry"])
    certificate = certify(
        model.G,
        model.E0,
        model.S,
        target,
        root,
        box=tuple(contract["material_box"]),
        timeout=30.0,
    )
    json.dumps(certificate, allow_nan=False)
    success = certificate["status"] == "UNKNOWN"
    return {
        "name": "negative_shifted_target_6",
        "reference_case": reference_case["name"],
        "target_shift": {"real": 1.0e-3, "imag": 2.0e-3, "source_index": 1},
        "assertions": {
            "unexpected_root_rejected": success,
            "failed_attempts_preserved": len(certificate["attempts"]) > 0,
            "json_serializable": True,
        },
        "certificate": certificate,
        "cpu_seconds": time.process_time() - cpu_start,
        "wall_seconds": time.perf_counter() - wall_start,
    }


def main() -> int:
    wall_start, cpu_start = time.perf_counter(), time.process_time()
    mp_precision_before = mp.iv.prec
    output: dict[str, Any] = {
        "schema": "Z-ALG-R5-LOCAL-KRAWCZYK-REGRESSION-v1",
        "tested_utc": datetime.now(timezone.utc).isoformat(),
        "environment": {
            "python": sys.version,
            "mpmath": mp.__version__,
            "mpmath_iv_available": hasattr(mp, "iv"),
            "numpy": np.__version__,
            "precision_schedule_bits": [64, 128, 256],
            "radius_multipliers": [1.0e-7, 1.0e-6, 1.0e-5],
        },
        "inputs_read_only": True,
        "outputs_modified": [str((HERE / "tests.json").relative_to(PROJECT))],
        "cases": [],
        "negative_controls": [],
    }
    failures: list[str] = []
    try:
        specifications = [
            (
                "high4_root0",
                R2 / "experiments/HIGH4_CONTRACT.json",
                R2 / "experiments/high4_public/E2_HIGH4_00.npz",
                R2 / "results/high4_E2_HIGH4_00_FullJ8.json",
                0,
            ),
            (
                "high4_root1",
                R2 / "experiments/HIGH4_CONTRACT.json",
                R2 / "experiments/high4_public/E2_HIGH4_00.npz",
                R2 / "results/high4_E2_HIGH4_00_FullJ8.json",
                1,
            ),
            (
                "general6_root0",
                R2 / "experiments/GENERAL6_CONTRACT.json",
                R2 / "experiments/general6_public/GENERAL6_00.npz",
                R2 / "results/general6_GENERAL6_00_general_quotient20.json",
                0,
            ),
            (
                "general6_root1",
                R2 / "experiments/GENERAL6_CONTRACT.json",
                R2 / "experiments/general6_public/GENERAL6_00.npz",
                R2 / "results/general6_GENERAL6_00_general_quotient20.json",
                1,
            ),
        ]
        for specification in specifications:
            case = _run_case(
                name=specification[0],
                contract_path=specification[1],
                observation_path=specification[2],
                result_path=specification[3],
                root_index=specification[4],
            )
            output["cases"].append(case)
            for assertion, passed in case["assertions"].items():
                if not passed:
                    failures.append(f"{case['name']}:{assertion}")

        negative = _negative_control(
            reference_case=output["cases"][2],
            contract_path=specifications[2][1],
            observation_path=specifications[2][2],
            result_path=specifications[2][3],
        )
        output["negative_controls"].append(negative)
        for assertion, passed in negative["assertions"].items():
            if not passed:
                failures.append(f"{negative['name']}:{assertion}")

        if mp.iv.prec != mp_precision_before:
            failures.append("mpmath_interval_precision_not_restored")
    except Exception as exc:
        failures.append(f"EXCEPTION:{type(exc).__name__}:{exc}")
    finally:
        output["status"] = "PASS" if not failures else "FAIL"
        output["failures"] = failures
        output["mpmath_interval_precision_before"] = mp_precision_before
        output["mpmath_interval_precision_after"] = mp.iv.prec
        output["cpu_seconds_total"] = time.process_time() - cpu_start
        output["wall_seconds_total"] = time.perf_counter() - wall_start
        output["cost_note"] = (
            "CPU/wall includes historical fixture construction, one fresh "
            "forward per archived root, all certificate attempts, reports, and "
            "tests.json serialization."
        )
        (HERE / "tests.json").write_text(
            json.dumps(output, indent=2, allow_nan=False) + "\n"
        )

    print(
        json.dumps(
            {
                "status": output["status"],
                "cases": [
                    {
                        "name": case["name"],
                        "status": case["certificate"]["status"],
                        "precision_bits": case["certificate"].get("precision_bits"),
                        "radius": case["certificate"].get("radius"),
                        "attempts": len(case["certificate"]["attempts"]),
                    }
                    for case in output["cases"]
                ],
                "negative_control": output["negative_controls"][0]["certificate"][
                    "status"
                ]
                if output["negative_controls"]
                else None,
                "cpu_seconds_total": output["cpu_seconds_total"],
                "wall_seconds_total": output["wall_seconds_total"],
                "failures": failures,
            },
            indent=2,
        )
    )
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
