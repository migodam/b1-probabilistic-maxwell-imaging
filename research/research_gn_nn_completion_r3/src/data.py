"""Fresh grouped material generation; only this module opens generation truth.

FAEMIS procedural.py inspired the primitive/composite/island geometry choices.
Its scattering model, labels and reconstructions are deliberately not imported.
All observations and truth currents come from the retained R4 full nonlinear VIE.
"""
from __future__ import annotations

import json
import time
from collections import Counter
from pathlib import Path
import uuid
import numpy as np
from PIL import Image, ImageDraw, ImageFont
from scipy.ndimage import gaussian_filter
from threadpoolctl import threadpool_limits
from research.research_z_alg_r4_diversity.src.physics import Model

SHAPE = (24, 24)
FAMILIES = ("digits", "primitive_composites", "islands", "smooth")


def dump_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def read_npz(path):
    with np.load(path, allow_pickle=False) as source:
        return {key: source[key] for key in source.files}


def font_path(name):
    candidates = [Path("/System/Library/Fonts/Supplemental") / name,
                  Path("C:/Windows/Fonts") / name,
                  Path("/usr/share/fonts/truetype/msttcorefonts") / name]
    aliases = {"Arial Bold.ttf": "arialbd.ttf", "Times New Roman Bold.ttf": "timesbd.ttf",
               "Courier New Bold.ttf": "courbd.ttf", "Georgia Bold.ttf": "georgiab.ttf"}
    if name in aliases:
        candidates.append(Path("C:/Windows/Fonts") / aliases[name])
    for path in candidates:
        if path.exists():
            return path
    raise FileNotFoundError("Required font missing: " + name)


def material_alpha(rng, family, ood=False, index=0):
    """Return a fresh 24² alpha and evaluator-only geometry provenance."""
    size = 192
    image = Image.new("L", (size, size), 0)
    draw = ImageDraw.Draw(image)
    meta = {"family": family, "distribution": "OOD" if ood else "ID"}
    if family == "digits":
        digit = 1 + index % 9
        names = ("Courier New Bold.ttf", "Times New Roman Bold.ttf") if ood else ("Arial Bold.ttf", "Georgia Bold.ttf")
        name = names[(index // 9) % len(names)]
        font = ImageFont.truetype(str(font_path(name)), 170)
        box = draw.textbbox((0, 0), str(digit), font=font)
        glyph = Image.new("L", (box[2] - box[0], box[3] - box[1]))
        ImageDraw.Draw(glyph).text((-box[0], -box[1]), str(digit), font=font, fill=255)
        scale = min(rng.uniform(.57, .78)*size/glyph.height, .70*size/glyph.width)
        glyph = glyph.resize((max(1, round(glyph.width*scale*rng.uniform(.85, 1.1))), max(1, round(glyph.height*scale))), Image.Resampling.LANCZOS)
        glyph = glyph.rotate(float(rng.uniform(-24, 24)), resample=Image.Resampling.BICUBIC, expand=True)
        x = int((size-glyph.width)/2 + rng.uniform(-.07, .07)*size)
        y = int((size-glyph.height)/2 + rng.uniform(-.07, .07)*size)
        image.paste(glyph, (x, y))
        meta.update(digit=digit, font=name, generator="fresh_transformed_font")
    elif family in ("primitive_composites", "islands"):
        count = int(rng.integers(2, 6)) if family == "primitive_composites" else int(rng.integers(3, 9))
        for _ in range(count):
            cx, cy = rng.uniform(.18, .82, 2)*size
            rx, ry = rng.uniform(.045, .20 if family == "primitive_composites" else .11, 2)*size
            box = (cx-rx, cy-ry, cx+rx, cy+ry)
            fill = int(rng.uniform(.65, 1)*255)
            if ood:
                if family == "primitive_composites":
                    vertices = int(rng.integers(3, 7))
                    theta = np.linspace(0, 2*np.pi, vertices, endpoint=False)+rng.uniform(0, 2*np.pi)
                    draw.polygon([(cx+rx*np.cos(t), cy+ry*np.sin(t)) for t in theta], fill=fill)
                else:
                    draw.arc(box, start=int(rng.integers(0, 90)), end=int(rng.integers(230, 350)), fill=fill, width=max(2, int(min(rx, ry)*.65)))
            elif family == "islands" or rng.random() < .65:
                draw.ellipse(box, fill=fill)
            else:
                draw.rectangle(box, fill=fill)
        meta.update(generator=("polygons" if family == "primitive_composites" else "arc_islands") if ood else "ellipse_rectangle_union", object_count=count)
    elif family == "smooth":
        if not ood:
            field = gaussian_filter(rng.normal(size=SHAPE), sigma=rng.uniform(1.7, 3.5))
            meta["generator"] = "gaussian_random_field"
        else:
            yy, xx = np.mgrid[-1:1:complex(24), -1:1:complex(24)]
            field = np.zeros(SHAPE)
            for _ in range(7):
                a, b = rng.integers(1, 6, 2)
                field += rng.normal()*np.sin(a*xx+b*yy+rng.uniform(-np.pi, np.pi))/(a+b)
            meta["generator"] = "random_fourier_field"
        field = np.maximum(field-np.quantile(field, .40), 0)
        alpha = field/max(float(field.max()), 1e-12)
        return alpha, meta
    else:
        raise ValueError(family)
    alpha = np.asarray(image.resize(SHAPE[::-1], Image.Resampling.BOX), dtype=float)/255
    return alpha, meta


def prepare_dataset(output, operators_path, seed=20261107, train=256, val=64, test=64, dev=8):
    """Write evaluator and online manifests with disjoint parent geometries."""
    output = Path(output)
    if (output/"manifests/evaluation.json").exists():
        raise FileExistsError("Dataset already exists; preserve its grouped split")
    started = time.perf_counter()
    ops = read_npz(operators_path)
    ops = {key: ops[key] for key in ("G", "E0", "S", "grid_shape", "k", "centers", "sources", "receivers", "radius") if key in ops}
    if "E0" not in ops or "S" not in ops:
        original = read_npz(operators_path)
        ops.update(E0=original.get("E0", original.get("U0")), S=original.get("S", original.get("H")))
    ops["grid_shape"] = np.array(SHAPE)
    if (ops["G"].shape, ops["E0"].shape, ops["S"].shape) != ((576, 576), (576, 12), (48, 576)):
        raise ValueError("Require unchanged R4 24x24 / 12-source / 48-receiver operators")
    # One shared receiver right basis; no per-observation SVD.
    basis_started = time.perf_counter()
    cache = Path(__file__).resolve().parents[2]/"research_gn_green_error_20261106/results/fixed_spectra.npz"
    basis_source = "one_economy_S_SVD"
    if cache.exists():
        cached = read_npz(cache)
        vh = cached.get("ob_Vh")
        if vh is not None and vh.shape == (48, 576):
            ops["V48"] = vh.conj().T
            basis_source = str(cache)
    if "V48" not in ops:
        with threadpool_limits(limits=1):
            _, singular, vh = np.linalg.svd(ops["S"], full_matrices=False)
        ops["V48"] = vh.conj().T
        ops["S_singular_values"] = singular
    basis_seconds = time.perf_counter()-basis_started
    output.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output/"operators.npz", **ops)
    online, evaluation, parents, receipts = [], [], [], []
    model = Model(ops)
    with threadpool_limits(limits=1):
        for split, count in (("train", train), ("val", val), ("test", test), ("dev", dev)):
            split_code = {"train": 1, "val": 2, "test": 3, "dev": 4}[split]
            for index in range(count):
                # Each parent owns an independent stream; variants never change split.
                rng = np.random.default_rng(np.random.SeedSequence([seed, split_code, index]))
                ood = split == "test" and index >= test//2
                family = FAMILIES[index % len(FAMILIES)]
                alpha, meta = material_alpha(rng, family, ood, index//len(FAMILIES))
                background = float(rng.uniform(.08, .15))
                amplitude = (.1, .2, .4)[index % 3]
                parent_id = f"{split}_p{index:04d}"
                # Four consecutive families per eight parents: half of every
                # family and each distribution, fixed before any reconstruction.
                stress = split == "test" and index % 8 < 4
                parents.append(dict(parent_id=parent_id, split=split, seed_sequence=[seed, split_code, index], background=background, amplitude=amplitude, stress_designated=stress, **meta))
                direction = rng.normal(size=(48, 12))+1j*rng.normal(size=(48, 12))
                direction /= np.linalg.norm(direction)
                levels = [("low", amplitude)] + ([("stress", 2.5)] if stress else [])
                for level, amp in levels:
                    m = background+amp*alpha.reshape(-1)
                    tick = time.perf_counter()
                    state = model.forward(m, label="full_vie_generation")
                    truth_name = parent_id if level == "low" else parent_id+"_stress"
                    truth_path = f"truth/{split}/{truth_name}.npz"
                    (output/truth_path).parent.mkdir(parents=True, exist_ok=True)
                    np.savez_compressed(output/truth_path, m=m, j=state["j"], alpha=alpha, amplitude=np.array(amp))
                    receipts.append(dict(parent_id=parent_id, level=level, seconds=time.perf_counter()-tick, state_relative_residual=state["state_relative_residual"]))
                    for noise in (0., 1e-3):
                        obs_id = "obs_"+uuid.uuid5(uuid.NAMESPACE_URL, f"gn-nn-r3/{seed}/{parent_id}/{level}/{noise}").hex[:20]
                        public_path = f"public/{obs_id}.npz"
                        (output/public_path).parent.mkdir(parents=True, exist_ok=True)
                        Y = state["y"]+noise*np.linalg.norm(state["y"])*direction
                        np.savez_compressed(output/public_path, Y=Y, noise=np.array(noise))
                        public = dict(obs_id=obs_id, parent_id=parent_id, split=split, noise=noise, public_path=public_path)
                        online.append(public)
                        evaluation.append(dict(**public, truth_path=truth_path, amplitude=amp, level=level, background=background, stress_designated=stress, **meta))
    header = dict(schema="gn-nn-r3-grouped-data-v1", operators_path="operators.npz", root=".", parent_grouped_splits=True)
    dump_json(output/"manifests/online.json", dict(**header, records=online))
    dump_json(output/"manifests/evaluation.json", dict(**header, records=evaluation, parents=parents, seed=seed, geometry_source="FAEMIS inspired geometry only; no FAEMIS scattering"))
    receipt = dict(status="GENERATED_FULL_VIE", parents=len(parents), observations=len(online), split_parents={"train": train, "val": val, "test": test, "dev": dev}, stress_parents=sum(p["stress_designated"] for p in parents), wall_seconds=time.perf_counter()-started, counts=dict(model.counts), forward_receipts=receipts, operators_source=str(operators_path), public_keys=["Y", "noise"], truth_keys=["m", "j", "alpha", "amplitude"], gt_exposed_online=False, receiver_basis_source=basis_source, receiver_basis_seconds=basis_seconds)
    dump_json(output/"GENERATION_RECEIPT.json", receipt)
    return receipt


def expand_training_dataset(output, target_train):
    """Append fresh train parents only; existing files/records are never rewritten.

    This entry point is implemented for a later parent decision. It does not
    choose a sample using truth, validation scores, or the current reconstruction.
    """
    if target_train not in (1024, 2048):
        raise ValueError("Predesignated training parent targets are 1024 and 2048")
    output = Path(output)
    online_path, evaluator_path = output/"manifests/online.json", output/"manifests/evaluation.json"
    online = json.loads(online_path.read_text(encoding="utf-8"))
    evaluator = json.loads(evaluator_path.read_text(encoding="utf-8"))
    existing = [p for p in evaluator["parents"] if p["split"] == "train"]
    current = len(existing)
    if current >= target_train:
        raise ValueError(f"Training already has {current} parents; expansion must increase it")
    if {p["parent_id"] for p in existing} != {f"train_p{i:04d}" for i in range(current)}:
        raise ValueError("Training parent indices are not a contiguous fixed prefix")
    seed = int(evaluator["seed"])
    if any(p["seed_sequence"] != [seed, 1, int(p["parent_id"].split("p")[-1])] for p in existing):
        raise ValueError("Existing training seed streams differ from the fixed generator contract")
    # Preflight only checks file identity/collision, not material content.
    for index in range(current, target_train):
        parent_id = f"train_p{index:04d}"
        paths = [output/f"truth/train/{parent_id}.npz"]
        for noise in (0., 1e-3):
            obs_id = "obs_"+uuid.uuid5(uuid.NAMESPACE_URL, f"gn-nn-r3/{seed}/{parent_id}/low/{noise}").hex[:20]
            paths.append(output/f"public/{obs_id}.npz")
        if any(path.exists() for path in paths):
            raise FileExistsError(f"Expansion output already exists for {parent_id}; preserve partial artifacts")
    started = time.perf_counter()
    model = Model(read_npz(output/"operators.npz"))
    new_online, new_evaluation, new_parents, forward_receipts = [], [], [], []
    status, error = "GENERATING", None
    receipt_path = output/f"GENERATION_EXPAND_TRAIN_{current:04d}_TO_{target_train:04d}_RECEIPT.json"
    if receipt_path.exists():
        raise FileExistsError("Preserve the existing expansion receipt")
    try:
        with threadpool_limits(limits=1):
            for index in range(current, target_train):
                rng = np.random.default_rng(np.random.SeedSequence([seed, 1, index]))
                family = FAMILIES[index % len(FAMILIES)]
                alpha, metadata = material_alpha(rng, family, False, index//len(FAMILIES))
                background = float(rng.uniform(.08, .15))
                amplitude = (.1, .2, .4)[index % 3]
                parent_id = f"train_p{index:04d}"
                new_parents.append(dict(parent_id=parent_id, split="train", seed_sequence=[seed, 1, index],
                                        background=background, amplitude=amplitude, stress_designated=False, **metadata))
                direction = rng.normal(size=(48, 12))+1j*rng.normal(size=(48, 12))
                direction /= np.linalg.norm(direction)
                m = background+amplitude*alpha.reshape(-1)
                tick = time.perf_counter()
                state = model.forward(m, label="full_vie_generation_train_expansion")
                truth_path = f"truth/train/{parent_id}.npz"
                np.savez_compressed(output/truth_path, m=m, j=state["j"], alpha=alpha, amplitude=np.array(amplitude))
                forward_receipts.append(dict(parent_id=parent_id, level="low", seconds=time.perf_counter()-tick,
                                             state_relative_residual=state["state_relative_residual"]))
                for noise in (0., 1e-3):
                    obs_id = "obs_"+uuid.uuid5(uuid.NAMESPACE_URL, f"gn-nn-r3/{seed}/{parent_id}/low/{noise}").hex[:20]
                    public_path = f"public/{obs_id}.npz"
                    Y = state["y"]+noise*np.linalg.norm(state["y"])*direction
                    np.savez_compressed(output/public_path, Y=Y, noise=np.array(noise))
                    public = dict(obs_id=obs_id, parent_id=parent_id, split="train", noise=noise, public_path=public_path)
                    new_online.append(public)
                    new_evaluation.append(dict(**public, truth_path=truth_path, amplitude=amplitude, level="low",
                                               background=background, stress_designated=False, **metadata))
        # Old manifest entries are retained verbatim and in their original order.
        # Append at the end rather than regrouping existing val/test/DEV entries.
        online["records"] = online["records"]+new_online
        evaluator["records"] = evaluator["records"]+new_evaluation
        evaluator["parents"] = evaluator["parents"]+new_parents
        for path, value in ((online_path, online), (evaluator_path, evaluator)):
            temporary = path.with_suffix(".expansion.tmp")
            dump_json(temporary, value)
            temporary.replace(path)
        status = "APPENDED_TRAIN_FULL_VIE"
    except Exception as exc:
        status, error = "FAILED_PARTIAL_ARTIFACTS_PRESERVED", f"{type(exc).__name__}: {exc}"
        raise
    finally:
        receipt = dict(schema="gn-nn-r3-training-expansion-v1", status=status, error=error,
                       train_parents_before=current, train_parents_target=target_train,
                       new_parents_complete=len(forward_receipts), new_observations_complete=len(new_online),
                       seed=seed, new_seed_indices=[current, target_train-1], wall_seconds=time.perf_counter()-started,
                       counts=dict(model.counts), forward_receipts=forward_receipts,
                       original_generation_receipt_preserved=True, existing_files_rewritten=0,
                       old_manifest_entries_preserved=True, selection="fixed parent seed streams only; no GT or reconstruction selection",
                       manifest_commit_completed=status == "APPENDED_TRAIN_FULL_VIE")
        dump_json(receipt_path, receipt)
    return receipt


def summarize_anchor_receipts(data_root, anchor_root):
    """Read-only aggregation; never opens truth, NPZ data, or runs physics."""
    data_root, anchor_root = Path(data_root), Path(anchor_root)
    online = json.loads((data_root/"manifests/online.json").read_text(encoding="utf-8"))
    by_obs = {r["obs_id"]: r for r in online["records"]}
    groups, unknown, counted_files, base_ids, batches = {}, [], [], set(), []
    for path in sorted(anchor_root.rglob("*.json")):
        try:
            receipt = json.loads(path.read_text(encoding="utf-8"))
        except (ValueError, OSError) as exc:
            unknown.append(dict(path=str(path), reason=str(exc)))
            continue
        if path.name.startswith("BATCH_"):
            batches.append(dict(path=str(path), status=receipt.get("status"),
                                wall_seconds=receipt.get("wall_seconds"), launched=receipt.get("launched")))
            continue
        if receipt.get("schema") != "gn-nn-r3-anchor-v1":
            continue
        obs_id = path.stem.split("_weakreg")[0]
        record = by_obs.get(obs_id)
        if record is None:
            unknown.append(dict(path=str(path), reason="No public manifest record for anchor"))
            continue
        kind = "weakreg1e-5" if receipt.get("weakreg", 0) else "base"
        if kind == "base": base_ids.add(obs_id)
        backend = receipt.get("backend", "cpu")  # This anchor schema uses the CPU R4 Model.
        status = receipt.get("stop_reason") or receipt.get("status", "UNKNOWN")
        key = (record["split"], status, backend, kind)
        group = groups.setdefault(key, dict(split=key[0], status=key[1], backend=key[2], anchor_kind=key[3],
                                           anchors=0, counts=Counter(), feature_counts=Counter(), solve_wall_seconds_sum=0.,
                                           feature_wall_seconds_sum=0., input_io_seconds_sum=0.,
                                           process_task_wall_seconds_sum=0., all_full_24_stage=True))
        group["anchors"] += 1
        group["counts"].update(receipt.get("counts", {}))
        group["feature_counts"].update(receipt.get("feature_counts", {}))
        for target, source in (("solve_wall_seconds_sum", "solve_wall_seconds"),
                               ("feature_wall_seconds_sum", "feature_wall_seconds"),
                               ("input_io_seconds_sum", "input_io_seconds"),
                               ("process_task_wall_seconds_sum", "process_task_wall_seconds")):
            group[target] += float(receipt.get(source, 0.))
        group["all_full_24_stage"] &= bool(receipt.get("all_unknowns_released", False))
        counted_files.append(str(path))
    rows = []
    for key in sorted(groups):
        row = groups[key]
        row["counts"] = dict(row["counts"])
        row["feature_counts"] = dict(row["feature_counts"])
        rows.append(row)
    missing = Counter(r["split"] for obs_id, r in by_obs.items() if obs_id not in base_ids)
    reference_path = anchor_root.parent/"cache/reference_jacobian.json"
    reference = json.loads(reference_path.read_text(encoding="utf-8")) if reference_path.exists() else None
    return dict(schema="gn-nn-r3-anchor-cost-summary-v1", groups=rows,
                pending_base_observations_by_split=dict(missing), counted_receipts=len(counted_files),
                receipt_files=counted_files, unclassified_or_unreadable=unknown,
                shared_reference_setup=reference, batch_wall_records=batches,
                wall_accounting="Per-task sums are work totals; overlapping batch wall records must not be summed into elapsed time",
                truth_access="NONE_JSON_PUBLIC_MANIFEST_AND_SCALAR_RECEIPTS_ONLY", physics_operations=0)
