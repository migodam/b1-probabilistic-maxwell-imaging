"""Independently verify the frozen 96² materials against their public noise balls."""

import argparse
from pathlib import Path
import sys

import numpy as np
from flint import arb, ctx
from threadpoolctl import threadpool_limits

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from research.z_alg_xudong_cie.holdout.data import (
    DEFAULT_ROOT, NATIVE_SIDE, OBJECT_IDS, OPTICAL, elapsed, planned_cases, public_directory,
    public_geometry, require_freeze, timers, write_json,
)

sys.path.insert(0, str(OPTICAL / "synthetic"))
from flint import acb_mat
from nonvacuity import ExactFourierGeometry, StableActiveModel, ball, hi, lo, norm


def verified_row(name, measured, selection, corrected, correction):
    rows, cols = selection
    prediction = acb_mat(len(rows), len(cols),
                         [corrected[int(i), int(j)] for i in rows for j in cols])
    residual = norm(prediction - ball(measured["Y"]))
    lower_bound = max(0., lo(residual.lower() - correction))
    upper_bound = hi(residual.upper() + correction)
    radius = float(measured["delta"])
    raw = float(measured["raw_radius"])
    y_nrm = norm(ball(measured["Y"]))
    expected_padding = arb(raw) + (arb(1) / 10**6) * y_nrm
    # The stored finite radius is used in both generation and scoring. Report
    # the Arb discrepancy of its 1 ppm normalization separately from fit.
    return dict(observation_id=name, native_gt_residual_lower=lower_bound,
                native_gt_residual_upper=upper_bound,
                raw_radius=raw, padded_radius=radius,
                difference_from_nominal_padding_upper=hi(abs(arb(radius)-expected_padding)),
                true_material_in_padded_ball=bool(arb(upper_bound) <= arb(radius)),
                verified_margin_lower=lo(arb(radius)-arb(upper_bound)))


def audit(root, freeze_record, output, objects=None):
    freeze = require_freeze(freeze_record)
    started = timers()
    settings = (128, 1)
    ctx.prec, ctx.threads = settings
    ids = tuple(objects) if objects else OBJECT_IDS
    setup_start = timers()
    ref = ExactFourierGeometry(public_geometry(root, NATIVE_SIDE))
    setup = elapsed(setup_start)
    cases = list(planned_cases())
    report = dict(scope="offline truth feasibility, independent interval forward",
                  freeze_record=str(freeze), precision_bits=settings[0], Arb_threads=settings[1],
                  geometry_setup_cost=setup, objects={}, failed={}, observations=[])
    for oid in ids:
        object_start = timers()
        with np.load(Path(root) / "offline" / "objects" / f"{oid}.npz", allow_pickle=False) as private:
            material = private["m_native"]
        try:
            model = StableActiveModel(ref, material)
            inv, inverse = model.inverse_bound()
            corrected, correction, solve = model.enclose_field(inv, 128)
            results = []
            for row in (r for r in cases if r["object_id"] == oid):
                name = row["case"]
                path = public_directory(root) / "observations" / f"{name}.npz"
                with np.load(path, allow_pickle=False) as measured:
                    results.append(verified_row(name, measured, ref.inds[row["array"]],
                                                corrected, correction))
            report["objects"][oid] = dict(active_cells=model.n, physical_LU=1, inverse_rhs=model.n,
                                           inverse=inverse, native_primal_and_adjoint=solve,
                                           **elapsed(object_start))
            report["observations"].extend(results)
            print(oid, "verified", sum(r["true_material_in_padded_ball"] for r in results),
                  "of", len(results), flush=True)
        except (ArithmeticError, ValueError) as exc:
            report["failed"][oid] = dict(reason=str(exc), **elapsed(object_start))
            print(oid, "unresolved", str(exc), flush=True)
        report.update(elapsed(started))
        write_json(output, report)
    report["all_requested_fitting_balls_verified"] = (
        len(report["observations"]) == 4 * len(ids) and
        all(r["true_material_in_padded_ball"] for r in report["observations"])
    )
    write_json(output, report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--freeze-record", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--object", choices=OBJECT_IDS, action="append")
    args = parser.parse_args()
    with threadpool_limits(1):
        audit(args.root, args.freeze_record, args.output or args.root / "offline/truth_audit.json", args.object)
