"""Fresh scalar-TM blind-case generation, write-once prediction seal, and scoring.

Only `generate` imports the physical forward model or opens private material.
`seal` reads public case IDs and prediction JSON. `score` requires that seal.
Keep the private directory local and out of every solver upload.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
import secrets
import time
import uuid
import sys

import numpy as np

HERE = Path(__file__).resolve().parents[1]
REPO = HERE.parents[1]
R4 = REPO / "research/research_z_alg_r4_diversity"
sys.path.insert(0, str(HERE))
from src.illumination_design import incident_field, source_positions

NOISY_SHAPE_SLOTS = (0, 3, 6, 9)
SOURCE_LAYOUT_SEED = 20261031


def dump_new(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x") as f:
        json.dump(value, f, indent=2, allow_nan=False)
        f.write("\n")


def _shape(index: int, rng: np.random.Generator) -> np.ndarray:
    """Nine deformed new-font digits, three shapes, two contrast controls."""
    if 0 <= index <= 8:
        from PIL import Image, ImageDraw, ImageFilter, ImageFont
        from matplotlib import font_manager
        n = 384
        font_family = "DejaVu Sans" if index % 2 else "DejaVu Serif"
        font_path = font_manager.findfont(font_manager.FontProperties(
            family=font_family, weight="bold"), fallback_to_default=False)
        font = ImageFont.truetype(font_path, 340)
        digit = str(index + 1)
        canvas = Image.new("L", (n, n), 0)
        box = ImageDraw.Draw(canvas).textbbox((0, 0), digit, font=font)
        glyph = Image.new("L", (box[2]-box[0], box[3]-box[1]), 0)
        ImageDraw.Draw(glyph).text((-box[0], -box[1]), digit, font=font, fill=255)
        target_h = rng.uniform(.63, .72)*n
        scale = min(target_h/glyph.height, .62*n/glyph.width)
        glyph = glyph.resize((max(1, round(glyph.width*scale)),
                              max(1, round(glyph.height*scale))), Image.Resampling.LANCZOS)
        if rng.random() < .5:
            glyph = glyph.filter(ImageFilter.MaxFilter(int(rng.choice([3, 5]))))
        else:
            glyph = glyph.filter(ImageFilter.MinFilter(3))
        dx, dy = rng.integers(-17, 18, size=2)
        canvas.paste(glyph, ((n-glyph.width)//2 + int(dx),
                             (n-glyph.height)//2 + int(dy)))
        canvas = canvas.rotate(float(rng.uniform(-13, 13)),
                               resample=Image.Resampling.BICUBIC)
        alpha = np.asarray(canvas.resize((24, 24), Image.Resampling.BOX), float)/255.
        if not (.025 < np.mean(alpha >= .5) < .50):
            raise RuntimeError("Generated digit has implausible support")
        return alpha
    n = 192
    y, x = np.mgrid[:n, :n]
    x = (x + .5)/n*2 - 1
    y = (y + .5)/n*2 - 1
    angle = rng.uniform(-.25, .25)
    c, s = np.cos(angle), np.sin(angle)
    u, v = c*x-s*y, s*x+c*y
    u -= rng.uniform(-.07, .07)
    v -= rng.uniform(-.07, .07)
    family = {9: 0, 10: 1, 11: 2, 12: 0, 13: 2}.get(index)
    if family is None:
        raise ValueError("Only 14 frozen shape slots are supported")
    if family == 0:  # asymmetric elliptical annulus
        a, b = rng.uniform(.55, .71), rng.uniform(.45, .62)
        outside = (u/a)**2+(v/b)**2 <= 1
        inside = ((u-rng.uniform(-.12, .12))/(a*rng.uniform(.42, .58)))**2 + ((v-rng.uniform(-.12, .12))/(b*rng.uniform(.42, .58)))**2 < 1
        mask = outside & ~inside
    elif family == 1:  # unequal separated inclusions
        l = ((u+.38)/rng.uniform(.20, .29))**2+((v+.12)/rng.uniform(.27, .38))**2 <= 1
        r = ((u-.34)/rng.uniform(.23, .34))**2+((v-.21)/rng.uniform(.17, .29))**2 <= 1
        mask = l | r
    elif family == 2:  # interrupted curved bar
        radius = np.sqrt((u-.03)**2+(v+.01)**2)
        phi = np.mod(np.arctan2(v+.01, u-.03), 2*np.pi)
        mask = (np.abs(radius-rng.uniform(.49, .58)) < rng.uniform(.10, .15)) & (phi > .25) & (phi < 5.25)
    alpha = mask.reshape(24, 8, 24, 8).mean(axis=(1, 3))
    if not (.04 < np.mean(alpha >= .5) < .55):
        raise RuntimeError("Generated shape has implausible support")
    return alpha


def _model():
    spec = importlib.util.spec_from_file_location("r4_physics_blind_generation", R4 / "src/physics.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.Model


def generate(public_out: Path, private_out: Path, operators: Path,
             method_starts: dict[str, int]) -> None:
    n_shapes = 14
    if not method_starts or not all(isinstance(k, str) and k and isinstance(v, int) and v > 0
                                   for k, v in method_starts.items()):
        raise ValueError("Freeze a nonempty method -> planned start count mapping")
    public_out, private_out = public_out.resolve(), private_out.resolve()
    if public_out == private_out or public_out in private_out.parents or private_out in public_out.parents:
        raise ValueError("Private truth and public data must be separate directories")
    if (public_out.exists() and any(public_out.iterdir())) or (private_out.exists() and any(private_out.iterdir())):
        raise FileExistsError("Preserve existing realization; both output directories must be empty")
    public_out.mkdir(parents=True, exist_ok=True)
    private_out.mkdir(parents=True, exist_ok=True)
    with np.load(operators, allow_pickle=False) as f:
        ops = {k: f[k] for k in f.files}
    if (ops["G"].shape, ops["E0"].shape, ops["S"].shape) != ((576, 576), (576, 12), (48, 576)):
        raise ValueError("Expected original 576/12/48 scalar-TM operators")
    # Both physical 12-source acquisitions are fixed before any object truth.
    uniform_ops = dict(ops)
    random_sources = source_positions(12, "random", ops["sources"], SOURCE_LAYOUT_SEED)
    random_ops = dict(ops)
    random_ops["sources"] = random_sources
    random_ops["E0"] = incident_field(ops["centers"], random_sources, float(ops["k"]))
    random_ops["U0"] = random_ops["E0"]
    public_ops_dir = public_out / "operators"
    public_ops_dir.mkdir()
    np.savez_compressed(public_ops_dir / "uniform12.npz", **uniform_ops)
    np.savez_compressed(public_ops_dir / "random12.npz", **random_ops)
    rng_seed = secrets.randbits(128)
    rng = np.random.default_rng(rng_seed)
    Model = _model()
    models = {"uniform12": Model(uniform_ops), "random12": Model(random_ops)}
    public_records, private_records = [], []
    started = time.perf_counter()
    # The order and IDs are random. All shapes are selected before any inference.
    for i in range(n_shapes):
        alpha = _shape(i, rng)
        amplitude = .4 if i == 12 else 2.5 if i == 13 else 12.
        background = float(rng.uniform(.08, .15))
        material = background + amplitude*alpha.ravel()
        layout = "uniform12" if i % 2 == 0 else "random12"
        forward = models[layout].forward(material, label="fresh_blind_generation")
        clean = forward["y"]
        group_id = "case_" + uuid.uuid4().hex[:16]
        variants = [("clean", 0)]
        if i in NOISY_SHAPE_SLOTS:
            variants.extend((("noisy", 1), ("noisy", 2)))
        for noise, realization in variants:
            obs_id = "obs_" + uuid.uuid4().hex[:16]
            noise_seed = None
            if noise == "clean":
                Y = clean
            else:
                noise_seed = secrets.randbits(128)
                noise_rng = np.random.default_rng(noise_seed)
                direction = (noise_rng.normal(size=clean.shape)
                             + 1j*noise_rng.normal(size=clean.shape))
                direction /= np.linalg.norm(direction)
                Y = clean + 1e-3*np.linalg.norm(clean)*direction
            np.savez_compressed(public_out / f"{obs_id}.npz", Y=Y,
                                observation_id=np.array(obs_id))
            private_records.append({"observation_id": obs_id, "group_id": group_id,
                                    "shape_slot": i, "shape_family":
                                    f"digit{i+1}" if i < 9 else
                                    ("ring" if i in (9, 12) else "two_islands" if i == 10 else "curved_bar"),
                                    "amplitude": amplitude,
                                    "background": background, "noise": noise,
                                    "noise_realization": realization,
                                    "noise_seed_PRIVATE": None if noise_seed is None else str(noise_seed),
                                    "source_layout": layout,
                                    "truth_file": f"{obs_id}.npz"})
            np.savez_compressed(private_out / f"{obs_id}.npz", m=material,
                                alpha=alpha, mask=alpha >= .5, Y_clean=clean)
            public_records.append({"observation_id": obs_id,
                                   "public_path": f"{obs_id}.npz", "noise": noise,
                                   "operator_path": f"operators/{layout}.npz"})
    rng.shuffle(public_records)
    dump_new(public_out / "manifest.json", {
        "schema": "fbg-r2-fresh-blind-v1", "status": "FRESH_UNSCORED",
        "records": public_records,
        "shape_count": n_shapes, "high_contrast_shapes": 12,
        "control_shapes": 2, "observation_count": len(public_records),
        "noise_plan": "14 clean; four anonymous objects have two independently drawn 1e-3 relative noise realizations each",
        "source_layouts": ["operators/uniform12.npz", "operators/random12.npz"],
        "source_layout_seed_public": SOURCE_LAYOUT_SEED,
        "layout_assignment": "alternating by hidden shape slot; paired method comparisons within each observation; layout and shape effects confounded",
        "planned_methods_and_starts": method_starts,
        "public_payload_keys": ["Y", "observation_id"],
        "solver_truth_access": "NONE", "planned_prediction_format": "methods -> starts -> result.m",
    })
    dump_new(private_out / "mapping_PRIVATE.json", {
        "schema": "fbg-r2-private-v1", "entropy": str(rng_seed),
        "noisy_shape_slots_PRIVATE": list(NOISY_SHAPE_SLOTS),
        "records": private_records,
        "truth_policy": "LOCAL ONLY; score after every planned prediction disposition is sealed",
        "forward_solves": n_shapes,
        "generation_costcounts": {k: sum(model.counts.get(k, 0) for model in models.values())
                                  for k in set().union(*(m.counts for m in models.values()))},
        "generation_wall_seconds": time.perf_counter()-started,
    })
    print(json.dumps({"status": "GENERATED_PRIVATE_TRUTH_HELD",
                      "public_manifest": str(public_out / "manifest.json"),
                      "observations": len(public_records), "forward_solves": n_shapes}))


def _finite_scalar(value):
    try:
        x = float(value)
        return x if np.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def _start_id(value, fallback: int) -> str:
    if value is None:
        return f"s{fallback}"
    if isinstance(value, int) or (isinstance(value, str) and value.isdigit()):
        return f"s{value}"
    return str(value)


def _normalize_start(source: dict | None, start_id: str) -> dict:
    if source is None:
        return {"start_id": start_id, "disposition": "NOT_RUN", "prediction": None}
    result = source.get("result") or source
    material = result.get("m")
    try:
        material_array = np.asarray(material, float).ravel() if material is not None else None
        finite = material_array is not None and material_array.size == 576 and np.isfinite(material_array).all()
    except (TypeError, ValueError):
        material_array, finite = None, False
    top = result.get("top1") or {}
    ranking = _finite_scalar(result.get("ranking_objective",
                             result.get("ranking_score", top.get("ranking_score"))))
    residual = _finite_scalar(result.get("relative_residual", top.get("relative_residual")))
    disposition = "COMPLETE_FINITE" if finite and ranking is not None else "FAILED_OR_NONFINITE"
    return {"start_id": start_id, "disposition": disposition,
            "prediction": {"m": material_array.tolist() if finite else None,
                           "ranking_objective": ranking, "relative_residual": residual,
                           "counts": result.get("counts", result.get("costcounts")),
                           "wall_seconds": _finite_scalar(source.get("wall_s", result.get("task_wall_s"))),
                           "status": source.get("status", result.get("status")),
                           "stop_reason": result.get("wrapper_stop_reason", result.get("stop_reason")),
                           "task_id": (source.get("task") or {}).get("task_id", source.get("task_id"))}}


def seal(manifest_path: Path, predictions_dir: Path, output: Path) -> None:
    manifest = json.loads(manifest_path.read_text())
    records = manifest.get("records") or []
    if (len(records) != 22 or len({r["observation_id"] for r in records}) != 22
            or sum(r["noise"] == "clean" for r in records) != 14
            or sum(r["noise"] == "noisy" for r in records) != 8
            or {r.get("operator_path") for r in records} !=
            {"operators/uniform12.npz", "operators/random12.npz"}):
        raise ValueError("Stage 4 manifest must contain exactly 14 clean + 8 noisy observations and two layouts")
    planned = manifest.get("planned_methods_and_starts")
    if not isinstance(planned, dict) or not planned:
        raise ValueError("Manifest must freeze planned methods and starts before generation")
    rows = []
    for record in records:
        ident = record["observation_id"]
        path = predictions_dir / f"{ident}.json"
        aggregate = json.loads(path.read_text()) if path.exists() else None
        if aggregate is not None and aggregate.get("observation_id", ident) != ident:
            raise ValueError(f"Prediction ID mismatch: {path}")
        payloads = (aggregate or {}).get("methods") or {}
        if not isinstance(payloads, dict):
            raise ValueError(f"methods must be a mapping: {path}")
        method_rows = []
        for method, K in planned.items():
            payload = payloads.get(method)
            if payload is None:
                raw_starts = []
            elif isinstance(payload, dict) and "starts" in payload:
                value = payload["starts"]
                raw_starts = ([{"start_id": k, **v} for k, v in value.items()]
                              if isinstance(value, dict) else value)
            else:
                raw_starts = [payload]  # direct K=1 payload also accepted
            if not isinstance(raw_starts, list):
                raise ValueError(f"starts must be a list or mapping: {path} / {method}")
            by_id = {}
            for i, item in enumerate(raw_starts):
                if not isinstance(item, dict):
                    raise ValueError(f"start item must be a mapping: {path} / {method}")
                task = item.get("task") or {}
                sid = _start_id(item.get("start_id", item.get("start_index", task.get("start_index"))), i)
                if sid in by_id:
                    raise ValueError(f"Duplicate start ID {sid}: {path} / {method}")
                by_id[sid] = item
            starts = [_normalize_start(by_id.get(f"s{i}"), f"s{i}") for i in range(K)]
            unplanned = sorted(set(by_id) - {f"s{i}" for i in range(K)})
            eligible = [s for s in starts if s["disposition"] == "COMPLETE_FINITE"]
            selected = min(eligible, key=lambda s: (s["prediction"]["ranking_objective"], s["start_id"])) if eligible else None
            method_rows.append({"method": method, "planned_starts": K, "starts": starts,
                                "online_selected_start_id": selected["start_id"] if selected else None,
                                "online_selection": "minimum finite public physical ranking_objective; tie by start_id",
                                "unplanned_start_ids": unplanned,
                                "unplanned_starts_preserved": [by_id[k] for k in unplanned]})
        rows.append({"observation_id": ident, "noise": record["noise"],
                     "operator_path": record["operator_path"],
                     "prediction_file_present": aggregate is not None,
                     "methods": method_rows,
                     "unplanned_methods_preserved": {k: v for k, v in payloads.items() if k not in planned}})
    counts = {"COMPLETE_FINITE": 0, "FAILED_OR_NONFINITE": 0, "NOT_RUN": 0}
    for row in rows:
        for method in row["methods"]:
            for start in method["starts"]:
                counts[start["disposition"]] += 1
    dump_new(output, {"schema": "fbg-r2-prediction-seal-v2",
                      "status": "ALL_PLANNED_PREDICTIONS_CLOSED",
                      "manifest": str(manifest_path), "sealed_unix": time.time(),
                      "planned_observations": len(rows),
                      "planned_methods_and_starts": planned,
                      "planned_predictions": sum(planned.values())*len(rows),
                      "disposition_counts": counts, "rows": rows,
                      "truth_access": "NONE"})
    print(json.dumps({"sealed_observations": len(rows), "predictions": counts,
                      "output": str(output)}))


def otsu(x: np.ndarray) -> float:
    x = np.asarray(x, float).ravel()
    low, high = float(x.min()), float(x.max())
    if high <= low:
        return high
    hist, edges = np.histogram(x, bins=256, range=(low, high))
    centers = (edges[:-1]+edges[1:])/2
    weight = np.cumsum(hist)
    moment = np.cumsum(hist*centers)
    n, total = len(x), moment[-1]
    denominator = weight*(n-weight)
    scores = np.divide((total*weight-n*moment)**2, denominator,
                       out=np.zeros(256), where=denominator > 0)
    return float(centers[int(np.argmax(scores))])


def iou(a: np.ndarray, b: np.ndarray) -> float:
    union = np.count_nonzero(a | b)
    return float(np.count_nonzero(a & b)/union) if union else 1.


def _scored_start(sealed_start: dict, truth: np.ndarray, mask: np.ndarray,
                  amplitude: float, noise: str) -> dict:
    disposition = sealed_start["disposition"]
    row = {"start_id": sealed_start["start_id"], "disposition": disposition,
           "RMSE": None, "NRMSE": None, "MAE": None, "IoU_Otsu": None,
           "relative_data_residual": None,
           "data_gate": "NOT_RUN" if disposition == "NOT_RUN" else "FAIL",
           "material_gate": "NOT_RUN" if disposition == "NOT_RUN" else "FAIL",
           "joint_gate": "NOT_RUN" if disposition == "NOT_RUN" else "FAIL"}
    if disposition != "COMPLETE_FINITE":
        return row
    pred = sealed_start["prediction"]
    m = np.asarray(pred["m"], float).ravel()
    rmse = float(np.sqrt(np.mean((m-truth)**2)))
    threshold = otsu(m)
    jaccard = iou(m > threshold, mask)
    residual = pred.get("relative_residual")
    limit = .001 if noise == "clean" else .003
    data_ok = residual is not None and residual <= limit
    material_ok = rmse/amplitude <= .15 and jaccard >= .7
    row.update({"RMSE": rmse, "NRMSE": rmse/amplitude,
                "MAE": float(np.mean(np.abs(m-truth))),
                "IoU_Otsu": jaccard,
                "otsu_threshold_reconstruction_only": threshold,
                "relative_data_residual": residual,
                "data_gate": "PASS" if data_ok else "FAIL",
                "material_gate": "PASS" if material_ok else "FAIL",
                "joint_gate": "PASS" if data_ok and material_ok else "FAIL"})
    return row


def score(seal_path: Path, private_dir: Path, output: Path) -> None:
    closed = json.loads(seal_path.read_text())
    if closed.get("status") != "ALL_PLANNED_PREDICTIONS_CLOSED" or closed.get("schema") != "fbg-r2-prediction-seal-v2":
        raise ValueError("Complete v2 prediction closure required before truth access")
    mapping = json.loads((private_dir / "mapping_PRIVATE.json").read_text())
    lookup = {r["observation_id"]: r for r in mapping["records"]}
    sealed_ids = {r["observation_id"] for r in closed["rows"]}
    if sealed_ids != set(lookup) or len(closed["rows"]) != 22:
        raise ValueError("Sealed observation set must match all 22 private observations")
    groups = {}
    for item in mapping["records"]:
        groups.setdefault(item["group_id"], []).append(item)
    if (len(groups) != 14 or
            sorted(len(group) for group in groups.values()) != [1]*10 + [3]*4 or
            any(sum(item["noise"] == "clean" for item in group) != 1 for group in groups.values())):
        raise ValueError("Private pairing must be 14 objects with one clean each and two noisy on four")
    rows = []
    for sealed in closed["rows"]:
        ident = sealed["observation_id"]
        meta = lookup[ident]
        if sealed["operator_path"] != f"operators/{meta['source_layout']}.npz":
            raise ValueError("Sealed public source layout disagrees with private generation receipt")
        with np.load(private_dir / meta["truth_file"], allow_pickle=False) as f:
            truth, mask = f["m"].ravel(), f["mask"].ravel().astype(bool)
        for method in sealed["methods"]:
            starts = [_scored_start(s, truth, mask, meta["amplitude"], meta["noise"])
                      for s in method["starts"]]
            predictions = [s["prediction"] for s in method["starts"] if s["prediction"] is not None]
            known_lu = [p["counts"].get("LU") for p in predictions
                        if isinstance(p.get("counts"), dict) and p["counts"].get("LU") is not None]
            known_wall = [p["wall_seconds"] for p in predictions if p.get("wall_seconds") is not None]
            known_rhs = []
            for pred in predictions:
                counts = pred.get("counts")
                if isinstance(counts, dict):
                    if counts.get("actual_complex_RHS") is not None:
                        known_rhs.append(counts["actual_complex_RHS"])
                    elif counts.get("source_RHS") is not None:
                        known_rhs.append(sum(counts.get(k, 0) for k in
                                             ("source_RHS", "adjoint_RHS", "tangent_RHS", "receiver_RHS")))
            selected_id = method["online_selected_start_id"]
            online = next((s for s in starts if s["start_id"] == selected_id), None)
            finite = [s for s in starts if s["disposition"] == "COMPLETE_FINITE"]
            oracle = min(finite, key=lambda s: (s["RMSE"], s["start_id"])) if finite else None
            rows.append({"observation_id": ident, "group_id": meta["group_id"],
                         "shape_slot": meta["shape_slot"], "shape_family_OFFLINE": meta["shape_family"],
                         "amplitude": meta["amplitude"], "noise": meta["noise"],
                         "noise_realization_OFFLINE": meta["noise_realization"],
                         "source_layout": meta["source_layout"],
                         "method": method["method"], "planned_starts": method["planned_starts"],
                         "finite_starts": len(finite),
                         "complete_start_set": all(s["disposition"] != "NOT_RUN" for s in starts),
                         "measured_costs_across_starts": {
                             "LU_sum_known": sum(known_lu) if known_lu else None,
                             "LU_known_starts": len(known_lu),
                             "actual_complex_RHS_sum_known": sum(known_rhs) if known_rhs else None,
                             "RHS_known_starts": len(known_rhs),
                             "wall_seconds_sum_known": sum(known_wall) if known_wall else None,
                             "wall_known_starts": len(known_wall),
                             "wall_scope": "sum of independent per-start task walls; not parallel batch elapsed wall"},
                         "start_scores_OFFLINE": starts,
                         "online_selected_start_id": selected_id,
                         "online": online,
                         "offline_best_at_K_ORACLE": None if oracle is None else {
                             "start_id": oracle["start_id"], "RMSE": oracle["RMSE"],
                             "NRMSE": oracle["NRMSE"], "IoU_Otsu": oracle["IoU_Otsu"],
                             "joint_gate": oracle["joint_gate"],
                             "selection_uses_GT": True},
                         "offline_oracle_scope": "finite planned starts only; not independent roots or a deployable selector"})
    by_method = []
    for method in closed["planned_methods_and_starts"]:
        group = [r for r in rows if r["method"] == method]
        by_method.append({"method": method, "planned_observations": len(group),
                          "independent_objects": len({r["group_id"] for r in group}),
                          "online_joint_pass": sum(r["online"] is not None and r["online"]["joint_gate"] == "PASS" for r in group),
                          "online_joint_fail": sum(r["online"] is not None and r["online"]["joint_gate"] == "FAIL" for r in group),
                          "online_not_run_or_no_finite": sum(r["online"] is None for r in group),
                          "oracle_best_at_K_joint_pass_OFFLINE": sum(r["offline_best_at_K_ORACLE"] is not None and r["offline_best_at_K_ORACLE"]["joint_gate"] == "PASS" for r in group)})
    dump_new(output, {"schema": "fbg-r2-blind-scores-v2",
                      "seal": str(seal_path), "truth_disclosed_only_after_seal": True,
                      "new_physics_LU": 0, "rows": rows, "by_method": by_method,
                      "planned_observations": 22, "independent_objects": 14,
                      "scope": "12 high A=12 (digits 1-9 plus three shapes), 2 contrast controls; 14 clean plus two noisy variants on each of four objects",
                      "layout_scope": "uniform/random 12-source layouts alternate by hidden shape slot; layout and shape effects are confounded; method comparisons remain paired within observation",
                      "oracle_warning": "best@K is offline GT selection among finite starts; online top-1 was frozen by public physical score"})
    print(json.dumps({"method_rows": len(rows), "methods": len(by_method),
                      "output": str(output)}))


def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="command", required=True)
    g = sub.add_parser("generate")
    g.add_argument("--public-out", type=Path, required=True)
    g.add_argument("--private-out", type=Path, required=True)
    g.add_argument("--operators", type=Path, default=R4 / "public_data/operators.npz")
    g.add_argument("--method-starts", type=Path, required=True,
                   help="Frozen JSON mapping method name to planned number of starts")
    g.add_argument("--go", action="store_true")
    s = sub.add_parser("seal")
    s.add_argument("--manifest", type=Path, required=True)
    s.add_argument("--predictions", type=Path, required=True)
    s.add_argument("--output", type=Path, required=True)
    s.add_argument("--close", action="store_true")
    c = sub.add_parser("score")
    c.add_argument("--seal", type=Path, required=True)
    c.add_argument("--private-dir", type=Path, required=True)
    c.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    if a.command == "generate":
        if not a.go:
            ap.error("Fresh generation requires --go after design selection")
        method_starts = json.loads(a.method_starts.read_text())
        generate(a.public_out, a.private_out, a.operators, method_starts)
    elif a.command == "seal":
        if not a.close:
            ap.error("Parent must close all planned prediction dispositions with --close")
        seal(a.manifest, a.predictions, a.output)
    else:
        score(a.seal, a.private_dir, a.output)


if __name__ == "__main__":
    main()
