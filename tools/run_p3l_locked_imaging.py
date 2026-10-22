"""Future serial locked online image panel. Authoring does not authorize a run.

Requires a parent-created METHODS_FROZEN_FOR_LOCKED.json with locked_imaging
fields documented in validate_freeze(). No truth reads or offline scoring.
Cheap16, independent native96, and transformed64 costs are kept separately.
"""
from pathlib import Path
import gc
import json
import sys
import time
import shutil

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from p3l.core import OUT, Book, Target, BudgetStop, write, save, plain
from p3l.inversion import reconstruct
from p3l.learned_images import (validation_selection, load_selected_density, context_features,
                                make_proposal, learned_image_actions)
from run_p3l_joint_gaussian import optimize as transformed_optimize


def validate_freeze(selection):
    path = OUT / "METHODS_FROZEN_FOR_LOCKED.json"
    if not path.exists():
        raise RuntimeError("locked methods freeze absent; no paid execution authorized")
    frozen = json.loads(path.read_text())
    spec = frozen.get("locked_imaging", {})
    expected = [f"SBI_locked_{i:04d}" for i in range(128)]
    requirements = dict(parents=expected, context_fresh_state_cap=16,
                        independent_native_fresh_state_cap=96,
                        transformed_fresh_state_cap=64, learned_draws=2048,
                        heldout_risk_requirement_registered=True,
                        checkpoint_selection="VALIDATION_NLL_ALL15_AND_SIMPLE_ONLY")
    for name, value in requirements.items():
        if spec.get(name) != value:
            raise RuntimeError(f"locked_imaging freeze missing or incompatible: {name}")
    if spec.get("selected_model") != selection["selected"]:
        raise RuntimeError("frozen selected model differs from all15+simple validation NLL selection")
    if not isinstance(spec.get("sampling_seed"), int) or spec["sampling_seed"] < 0:
        raise RuntimeError("parent must predeclare nonnegative learned sampling_seed")
    if not isinstance(spec.get("max_wall_seconds"), (int, float)) or not 0 < spec["max_wall_seconds"] <= 7200:
        raise RuntimeError("parent must predeclare positive panel wall cap ≤7200 seconds")
    return spec


def cost(book, before, started):
    after = book._latest()["stages"]["D"]
    return dict(RHS=int(after["RHS"]-before["RHS"]),
                fresh_forward_states=int(after["states"]-before["states"]),
                LU=int(after["LU"]-before["LU"]), wall_seconds=time.perf_counter()-started,
                resource_source="authoritative D-stage ledger delta inside serial paid worker")


def main():
    # Existence check precedes model report loading; it never writes a freeze.
    if not (OUT / "METHODS_FROZEN_FOR_LOCKED.json").exists():
        raise RuntimeError("METHODS_FROZEN_FOR_LOCKED.json required")
    selection = validation_selection(OUT)
    spec = validate_freeze(selection)
    manifest = json.loads((OUT / "configs/PARENTS.json").read_text())
    parents = [r["parent"] for r in manifest if r["group"] == "locked"]
    if parents != spec["parents"]:
        raise RuntimeError("original frozen 128 locked identities/order differ")
    generator = OUT / "AUDIT/locked_DATA_GENERATOR_EXIT_REQUIRED.json"
    if not generator.exists() or json.loads(generator.read_text()).get("generated") != 128:
        raise RuntimeError("separate locked generator must finish before online imaging")
    loaded = load_selected_density(selection["selected"], out=OUT)
    if loaded["best_validation_raw_target_nll"] != selection["selected"]["best_validation_raw_target_nll"]:
        raise RuntimeError("loaded checkpoint validation NLL differs from frozen selection")
    if (loaded["model_type"] != selection["selected"]["model_type"] or
            loaded["metadata"].get("seed") != selection["selected"]["seed"]):
        raise RuntimeError("loaded checkpoint family/seed differs from frozen selection")
    with np.load(OUT / "D/FIXED_FRAME.npz", allow_pickle=False) as archive:
        frame = {name: archive[name].copy() for name in ("Vw", "Vs")}
    frame_metadata = json.loads((OUT / "D/FIXED_FRAME.json").read_text())
    if frame_metadata.get("training_only") is not True or frame_metadata.get("GT_used") is not False:
        raise RuntimeError("training-only weak frame provenance absent")
    book = Book("D", "locked_imaging_frozen_panel", max_wall=spec["max_wall_seconds"])
    status, error = "COMPLETE", None
    rows = [dict(parent=p, status="NOT_RUN", phases={}) for p in parents]
    summary = OUT / "D/LOCKED_ONLINE_IMAGING.json"
    try:
        write(book.jobdir/'PROTOCOL.json',dict(method_freeze='METHODS_FROZEN_FOR_LOCKED.json',locked_imaging=spec,GT_read=False))
        for relative in ['tools/run_p3l_locked_imaging.py','tools/run_p3l_joint_gaussian.py','src/p3l/learned_images.py','src/p3l/learned_proposal.py','src/p3l/inversion.py','src/p3l/simple_density.py','src/p3l/core.py','src/p1_material/prior.py']:
            destination=book.jobdir/'SOURCE'/relative;destination.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(ROOT/relative,destination)
        target = Target(book, keep_fields=True)
        for index, row in enumerate(rows):
            parent = row["parent"]
            folder = OUT / f"D/locked_imaging/{parent}"
            cheap = strong = mode = proposal = None
            active_phase = None
            try:
                y = target.observation(parent)
                for phase, cap, polish in (("cheap16", 16, False), ("independent_native96", 96, True)):
                    active_phase = phase
                    before = dict(book._latest()["stages"]["D"]); started = time.perf_counter()
                    # Both calls independently start from zero z and midpoint
                    # properties, matching training contexts and fair baseline.
                    try:
                        point, report = reconstruct(target, y, parent+"_"+phase,
                            method="coarse_TRF", max_evaluations=cap, polish=polish,
                            initial_point=None)
                        target.persist_cache(point, folder / f"{phase}_CACHE.npz")
                        J = target.jacobian(point)
                        save(folder / f"{phase}_ESTIMATE.npz", q=point.q, chi=point.state.chi,
                             prediction=point.state.data, J=J)
                        row["phases"][phase] = dict(status="SAVED", source_endpoint=point.key,
                            source_cache=str((folder/f"{phase}_CACHE.npz").relative_to(OUT)), report=plain(report))
                        if phase == "cheap16":
                            cheap = point
                            scalar, image = context_features(target.prior, point.q, point.state.data, J, y, target.sigma)
                            save(folder / "CHEAP16_CONTEXT.npz", scalar=scalar, image=image,
                                 q=point.q, prediction=point.state.data, J=J,
                                 residual=target.white(point.state.data-y))
                            try:
                                proposal = make_proposal(loaded["model"], scalar, image, frame,
                                                         point.q, J, target.white(point.state.data-y))
                                arrays, diagnostics = learned_image_actions(proposal, target.prior,
                                    seed=spec["sampling_seed"]+index, count=2048, batch_size=512)
                                save(folder / "LEARNED_REFERENCE_IMAGES.npz", **arrays)
                                write(folder / "LEARNED_REFERENCE_IMAGES.json", diagnostics)
                                row["learned_images"] = dict(status="SAVED_UNQUALIFIED_MODEL_REFERENCE",
                                    selected_model=selection["selected"], diagnostics=diagnostics,
                                    context_source="CHEAP16_CONTEXT.npz; exactly training feature API",
                                    full_target_correction="NOT_RUN_BY_THIS_SCRIPT")
                            except BudgetStop:
                                raise
                            except Exception as exc:
                                row["learned_images"] = dict(status="FAILED_RETAINED_POINT_BASELINES_UNCHANGED", error=repr(exc))
                            arrays = proposal = cheap = None
                        else:
                            strong = point
                    except BudgetStop:
                        raise
                    except Exception as exc:
                        row["phases"][phase] = dict(status="FAILED_RETAINED", error=repr(exc))
                    finally:
                        row["phases"].setdefault(phase, {})["cost"] = cost(book, before, started)
                        point = None
                        if target.model is not None:
                            target.model._transfers.clear()
                        gc.collect()
                active_phase = "probit_refinement64"
                if strong is None:
                    row["phases"][active_phase] = dict(status="NOT_RUN_NATIVE_BASELINE_FAILED")
                else:
                    before = dict(book._latest()["stages"]["D"]); started = time.perf_counter()
                    try:
                        mode, x, H, stationarity = transformed_optimize(target, y, parent,
                                                                      folder/"probit64", strong)
                        save(folder / "PROBIT64_ESTIMATE.npz", q=mode.q, x=x,
                             chi=mode.state.chi, prediction=mode.state.data, H_GN=H)
                        row["phases"][active_phase] = dict(status="SAVED_PROBIT_COORDINATE_MODE_NOT_NATIVE_MAP",
                            source_endpoint=mode.key, stationarity=plain(stationarity),
                            source_cache=str((folder/"probit64/RAW_MODE_CACHE.npz").relative_to(OUT)))
                    except BudgetStop:
                        raise
                    except Exception as exc:
                        row["phases"][active_phase] = dict(status="FAILED_RETAINED", error=repr(exc))
                    finally:
                        row["phases"].setdefault(active_phase, {})["cost"] = cost(book, before, started)
                row["status"] = "ONLINE_ESTIMATES_RECORDED_FAILURES_RETAINED"
            except BudgetStop:
                row.update(status="BUDGET_STOP_PARTIAL_ARTIFACTS_RETAINED", interrupted_phase=active_phase)
                raise
            except Exception as exc:
                row.update(status="PARENT_FAILURE_RETAINED", error=repr(exc))
            finally:
                cheap = strong = mode = proposal = None
                if target.model is not None:
                    target.model._transfers.clear()
                gc.collect()
                write(folder / "ONLINE_RESULT.json", row)
                write(summary, dict(rows=rows, all128_retained=True, GT_read=False,
                    freeze="METHODS_FROZEN_FOR_LOCKED.json", model_selection=selection,
                    full_target_correction="SEPARATE_RUNNER_NOT_RUN_HERE",
                    interpretation="learned×Gaussian-reference MC actions, plug-in control and point estimators; no posterior calibration or teacher claim",
                    prospective_memory_policy="reconstruct: four-entry Point cache plus pinned best; small residual cache; evicted revisits paid afresh inside unchanged state cap; serial transfer cache limit4; transformed optimizer current+best",
                    nominal_max_estimation_RHS=901120,
                    nominal_bound_note="128*(16+96+64)*40; actual source-polish/receiver/cache costs retained separately in ledger; truth1024 and future IS separate"))
    except BaseException as exc:
        status, error = "ERROR", repr(exc)
        raise
    finally:
        book.close(status, error)


if __name__ == "__main__":
    main()
