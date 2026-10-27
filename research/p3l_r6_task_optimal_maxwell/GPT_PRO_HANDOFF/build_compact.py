"""Offline, lossless main-D handoff. No model, decoder, physics or scoring imports.

Run with ../../.venv_nn/bin/python from the R6 repository root. Existing outputs
are never overwritten; --verify-only repeats all source/compact array checks.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import time
import numpy as np

HERE = Path(__file__).resolve().parent
P = HERE.parent
GROUPS = {"DLOW": 0.5, "DHIGH": 2.0}


def read(path):
    return json.loads(Path(path).read_text())


def identity(path):
    path = Path(path)
    digest = hashlib.blake2b()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return {"path": str(path.resolve()), "bytes": path.stat().st_size,
            "blake2b": digest.hexdigest()}


def checked(record):
    actual = identity(record["path"])
    if actual != record:
        raise ValueError("Sealed source identity mismatch: " + record["path"])
    return Path(record["path"])


def exact(left, right, label):
    if (left.dtype != right.dtype or left.shape != right.shape or
            left.tobytes(order="C") != right.tobytes(order="C")):
        raise ValueError("Array not byte-identical: " + label)


def load(path):
    with np.load(path, allow_pickle=False) as source:
        return {key: source[key].copy() for key in source.files}


def main(verify_only=False):
    wall = time.perf_counter()
    cpu = time.process_time()
    unlock_path = P / "track_d/GT_UNLOCK.json"
    unlock = read(unlock_path)
    if (unlock["status"] != "ROOT_UNLOCKED_AFTER_SEAL" or
            not unlock["full512_and_supplement128_predictions_sealed_before_truth"]):
        raise RuntimeError("GT access requires the existing root unlock.")
    global_path = checked(unlock["global_seal_identity"])
    global_seal = read(global_path)
    if global_seal["status"] != "ALL_512_PREDICTIONS_SEALED_BEFORE_D_GT":
        raise ValueError("Wrong main-D seal")
    freeze_path = checked(global_seal["freeze_identity"])
    freeze = read(freeze_path)
    seeds = freeze["seeds"]
    actions = list(freeze["estimator_arms"])
    if freeze["groups"] != list(GROUPS) or len(seeds) != 3 or len(actions) != 6:
        raise ValueError("Expected two separate groups, six actions and three seeds")
    dest = HERE / "compact"
    dest.mkdir(exist_ok=True)
    manifest_path = HERE / "COMPACT_NUMERICAL_MANIFEST.json"
    if manifest_path.exists() and not verify_only:
        raise FileExistsError("Manifest already exists; use --verify-only")
    results = []
    all_parents = set()
    for group, scale in GROUPS.items():
        path = dest / (group + "_COMPACT_NUMERICAL.npz")
        if path.exists() and not verify_only:
            raise FileExistsError(path)
        if verify_only and not path.exists():
            raise FileNotFoundError(path)
        group_identity = next(record for record in global_seal["group_seals"]
                              if Path(record["path"]).parent.name == group)
        seal_path = checked(group_identity)
        seal = read(seal_path)
        checked(seal["freeze_identity"])
        parent_manifest_path = checked(freeze["manifest_identities"][group])
        parent_manifest = read(parent_manifest_path)
        parents = [row["parent"] for row in parent_manifest["parents"]]
        if parents != seal["parents"] or len(parents) != 256 or len(set(parents)) != 256:
            raise ValueError("Parent manifest/order mismatch")
        if all_parents.intersection(parents):
            raise ValueError("Cohorts overlap")
        all_parents.update(parents)
        phase_path = checked(seal["phases_identity"])
        phase = load(phase_path)
        exact(phase["parent"], np.asarray(parents), "phase parent order")
        arrays = {"parent": phase["parent"], "seeds": np.asarray(seeds, dtype=np.int64),
                  "action_names": np.asarray(actions),
                  "sigma_scale": np.asarray(scale, dtype=np.float64),
                  "full_anchor_phase_code": phase["phase_code"],
                  "full_anchor_phase_fractions": phase["fractions"]}
        sources = [{"role": "group_seal", **group_identity},
                   {"role": "parent_manifest", **freeze["manifest_identities"][group]},
                   {"role": "full_anchor_phases", **seal["phases_identity"]}]
        # Each mapping says exactly which compact slice is compared with which
        # original array; it also makes the full provenance machine-readable.
        mappings = [dict(compact="parent", source=str(phase_path), array="parent"),
                    dict(compact="full_anchor_phase_code", source=str(phase_path), array="phase_code"),
                    dict(compact="full_anchor_phase_fractions", source=str(phase_path), array="fractions")]
        records = {(r["method"], r["seed"]): r for r in seal["records"]}
        if len(records) != 18:
            raise ValueError("Expected every frozen action/seed exactly once")
        for action in actions:
            predictions = []
            for seed_index, seed in enumerate(seeds):
                record = records[action, seed]
                source_path = checked(record["identity"])
                prediction = load(source_path)
                exact(prediction["parent"], arrays["parent"], "prediction parent order")
                predictions.append(prediction["prediction"])
                sources.append(dict(role="frozen_action_prediction", action=action,
                                    seed=seed, **record["identity"]))
                mappings.append(dict(compact="prediction_" + action,
                                     seed_index=seed_index, source=str(source_path), array="prediction"))
            if len({a.dtype.str for a in predictions}) != 1:
                raise ValueError("Stack would alter source dtype")
            arrays["prediction_" + action] = np.stack(predictions, axis=1)
        context_records = {(r["parent"], r["context"]): r for r in seal["context_artifacts"]}
        if len(context_records) != 512:
            raise ValueError("Expected both physical contexts for every parent")
        collected = {}
        def append(key, value, source, original_key, index):
            collected.setdefault(key, []).append(value)
            mappings.append(dict(compact=key, parent_index=index,
                                 source=str(source), array=original_key))
        for index, parent in enumerate(parents):
            truth_path = P / "raw/offline_truth" / (parent + ".npz")
            online_path = P / "raw/online" / (parent + ".npz")
            truth, online = load(truth_path), load(online_path)
            if str(truth["parent"]) != parent or str(online["parent"]) != parent:
                raise ValueError("GT/observation parent mismatch")
            sources.extend([dict(role="unlocked_offline_GT", **identity(truth_path)),
                            dict(role="original_online_observation", **identity(online_path))])
            for key, original in (("GT_material", "chi"), ("GT_native_q", "q"),
                                  ("GT_gaussianized_x", "x"), ("GT_noiseless_response", "mean"),
                                  ("noise_innovation", "noise_innovation"),
                                  ("prior_seed", "prior_seed"), ("noise_seed", "noise_seed")):
                append(key, truth[original], truth_path, original, index)
            for key in ("observation", "sigma"):
                append(key, online[key], online_path, key, index)
            # This is saved-data arithmetic only, not a new forward call.
            exact(truth["mean"] + online["sigma"] * truth["noise_innovation"],
                  online["observation"], "saved noiseless+noise=observation")
            for context in ("full", "low"):
                record = context_records[parent, context]
                source_path = checked(record["context_identity"])
                checked(record["report_identity"])
                value = load(source_path)
                exact(value["observation"], online["observation"], "context observation")
                if str(value["parent"]) != parent:
                    raise ValueError("Context parent mismatch")
                sources.extend([dict(role=context + "_context", **record["context_identity"]),
                                dict(role=context + "_context_report", **record["report_identity"])])
                for label, original in (("physics_native_q", "q0"), ("physics_material", "chi0"),
                                        ("native_stage_q", "native_q"), ("native_stage_material", "native_chi"),
                                        ("physics_response", "prediction"), ("physics_gaussianized_x", "x")):
                    append(context + "_" + label, value[original], source_path, original, index)
        for key, values in collected.items():
            if len({a.dtype.str for a in values}) != 1:
                raise ValueError("Stack would alter source dtype: " + key)
            arrays[key] = np.stack(values)
        if not np.all(arrays["sigma"] == arrays["sigma"][0]):
            raise ValueError("Unexpected varying sigma within fixed condition")
        schema = {key: {"shape": list(a.shape), "dtype": a.dtype.str}
                  for key, a in arrays.items()}
        if not verify_only:
            temporary = path.with_suffix(".tmp.npz")
            np.savez_compressed(temporary, **arrays)
            os.replace(temporary, path)
        # Verify every output including metadata, then independently reopen each
        # original source and compare every included slice. No tolerance/casting.
        compact = load(path)
        if compact.keys() != arrays.keys():
            raise ValueError("Compact schema mismatch")
        for key, value in arrays.items():
            exact(compact[key], value, "roundtrip " + key)
        source_cache = {}
        for mapping in mappings:
            source = mapping["source"]
            if source not in source_cache:
                source_cache[source] = load(source)
            expected = source_cache[source][mapping["array"]]
            actual = compact[mapping["compact"]]
            if "parent_index" in mapping:
                actual = actual[mapping["parent_index"]]
            if "seed_index" in mapping:
                actual = actual[:, mapping["seed_index"]]
            exact(actual, expected, str(mapping))
            # Release J/H immediately; do not retain all contexts in memory.
            del source_cache[source]
        results.append(dict(group=group, parents=256, sigma=float(arrays["sigma"][0]),
                            sigma_scale=scale, actions=actions, seeds=seeds,
                            identity=identity(path), schema=schema, sources=sources,
                            source_array_mappings=mappings,
                            validation=dict(status="PASS_BYTE_IDENTICAL", included_source_arrays=len(mappings),
                                            compact_arrays=len(arrays), dtype_shape_and_bytes=True,
                                            saved_noiseless_plus_original_noise_exact=True,
                                            full_low_observations_identical=True)))
        print(json.dumps(dict(group=group, bytes=path.stat().st_size,
                              verified_source_arrays=len(mappings))), flush=True)
    if verify_only:
        old = read(manifest_path)
        for actual, previous in zip(results, old["groups"], strict=True):
            if actual != previous:
                raise ValueError("Source/schema/compact differs from original compact manifest")
        print("VERIFY_ONLY_PASS; no output overwritten")
        return
    manifest = dict(status="OFFLINE_LOSSLESS_HANDOFF_VERIFIED", format_version=1,
                    GT_exposed=True, exposure_scope="All 512 main-D ground truths unlocked by root after prediction seals",
                    not_posterior_samples=True, conditions_statistically_pooled=False,
                    GT_best_seed_or_parent_selection=False, parent_selection="Frozen manifest order; all 256 per condition",
                    semantics={
                        "prediction_axis": "Each prediction_ACTION has [parent, seed, material_cell]; six actions are frozen decision outputs, not posterior samples",
                        "material": "Complex contrast chi, flattened 24x24 grid; preserves source precision",
                        "native_q": "66 native prior coordinates: Gaussian field64, Uniform property2; optimizer endpoints and truth, not posterior draws",
                        "gaussianized_x": "66 Gaussianized prior coordinates; retained separately from native q",
                        "full_physics": "Native96 then Probit64 final q0/chi0 physical endpoint",
                        "low_physics": "Native24 then Probit64 final q0/chi0 physical endpoint",
                        "native_stage": "Endpoint before the subsequent Probit stage, explicitly distinct from final physics endpoint",
                        "phase": "Full final-anchor decoder phase, not GT; code -1 outside ROI,0 host,1 transition,2 pure; fractions columns host/transition/pure",
                        "observation": "Original complex128 32x8 noisy observation; exact source sigma, saved mean and innovation retained",
                        "scope": "Two fixed noise cohorts only; no new layout; main512 distinct from exposed CDEV and supplement128"},
                    unlock_identity=identity(unlock_path), global_seal_identity=identity(global_path),
                    freeze_identity=identity(freeze_path), builder_identity=identity(__file__), groups=results,
                    cost=dict(new_physics_RHS=0, new_training=0, new_inference=0,
                              wall_seconds=time.perf_counter()-wall, CPU_seconds=time.process_time()-cpu))
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n")
    print(str(manifest_path))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()
    main(args.verify_only)
