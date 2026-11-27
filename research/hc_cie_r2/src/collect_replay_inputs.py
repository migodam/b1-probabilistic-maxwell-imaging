#!/usr/bin/env python3
"""Collect the frozen inputs required by HC-CIE-R2 using a file and array-key allowlist."""
import argparse
from io import BytesIO
import json
from pathlib import Path
from time import perf_counter
import zipfile

import numpy as np


HERE = Path(__file__).resolve().parents[1]
OLD = Path("research/delegated/max_struct_r1")
SRAP = Path("research/srap_cie_array_study_20261121")
R13 = Path("research/research_max_struct_r13")
R6 = Path("research/research_max_struct_r6")
R13_FIELDS = ("G", "E0", "S", "Y", "grid_shape", "delta", "sigma", "m_init",
              "lower", "upper", "forward_side", "inverse_side", "sources", "receivers")
MODEL_FIELDS = ("G", "S", "E0", "coords", "sources", "receivers", "grid_shape", "k", "radius")
COHORT_FIELDS = ("Y", "sigma", "material_bounds", "grid_shape", "noise_radius")
ENDPOINT_FIELDS = ("m", "Y_prediction")
BASE_COST_FIELDS = ("whole_process_CPU_seconds", "cpu_seconds", "wall_seconds")
TOY_NUMERICAL_FIELDS = {
    "tm_disjoint_elimination_counterexample.json": (
        "k", "area", "positions", "receiver_angles", "incident_angles", "B", "chi", "J_bad", "weight"),
    "tm_disjoint_exact_integral_counterexample_verified.json": (
        "k", "area", "positions", "receiver_angles", "incident_angles", "chi",
        "G", "B", "S", "Y", "J_bad", "J_true", "weight"),
    "tm_hardfiber_counterexample.json": (
        "k", "area", "positions", "receiver_angles", "incident_angles",
        "G", "B", "S", "Y", "N", "J0", "Z_bad", "chi"),
    "variables_exact_convex_certificate.json": ("strong_convexity_modulus_lower_exact",),
    "variables_exact_nonconvex_witness.json": ("second_directional_derivative_exact",),
    "exact_two_material_ambiguity.json": ("G", "B", "S"),
}


def load_json(path):
    return json.loads(path.read_text())


def select(record, keys):
    return {key: record[key] for key in keys}


def historical(workspace, worktree, nominal):
    direct = workspace/nominal
    return direct if direct.is_dir() else workspace/OLD/worktree/nominal


class Archive:
    def __init__(self, zipped, workspace):
        self.zipped = zipped
        self.workspace = workspace
        self.entries = {}

    def describe(self, source):
        try:
            return str(source.relative_to(self.workspace))
        except ValueError:
            return source.name

    def put(self, path, payload, role, source, **details):
        path = Path(path)
        if path.is_absolute() or ".." in path.parts:
            raise ValueError(path)
        name = path.as_posix()
        if name in self.entries:
            raise ValueError(f"Duplicate archive path {name}")
        entry = zipfile.ZipInfo(name, date_time=(2026, 11, 28, 0, 0, 0))
        compression = zipfile.ZIP_STORED if path.suffix == ".npz" else zipfile.ZIP_DEFLATED
        self.zipped.writestr(entry, payload, compress_type=compression, compresslevel=5)
        self.entries[name] = dict(role=role, original=self.describe(source), bytes=len(payload), **details)

    def npz(self, path, source, keys, role, partial=False):
        with np.load(source, allow_pickle=False) as original:
            available = tuple(key for key in keys if key in original)
            if not partial and set(available) != set(keys):
                raise ValueError(f"Missing input keys in {source}: {set(keys)-set(available)}")
            arrays = {key: original[key] for key in available}
        if not arrays or any(value.dtype.hasobject for value in arrays.values()):
            raise ValueError(f"Unexpected empty or object-valued subset {source}")
        staged = BytesIO()
        np.savez_compressed(staged, **arrays)
        self.put(path, staged.getvalue(), role, source, arrays={
            key: dict(dtype=str(arr.dtype), shape=list(arr.shape)) for key, arr in arrays.items()})

    def json(self, path, source, data, role):
        self.put(path, (json.dumps(data, ensure_ascii=False, indent=2)+"\n").encode(),
                 role, source)


def r13_assets(archive):
    source = historical(archive.workspace, "r13_worktree", R13)
    arms = ("cie", "dual_cie", "a_plus_arc_r7", "g19", "dr")
    for i in range(48):
        case = f"r13_fresh_{i:04d}"
        path = Path("data_public/fresh48")/(case+".npz")
        archive.npz(R13/path, source/path, R13_FIELDS, "online_observation")
        path = Path("sealed_truth/fresh48")/(case+".npz")
        archive.npz(R13/path, source/path, ("material_coarse",), "offline_label")
        for arm in arms:
            family = (Path("supplement_aplus48/stack_run/RESULTS/aplus_stack112")
                      if arm == "a_plus_arc_r7" else Path("RESULTS/confirmation_core48"))
            folder = family/f"{case}__{arm}"
            path = folder/"primary.npz"
            archive.npz(R13/path, source/path, ENDPOINT_FIELDS, "online_historical_endpoint")
            receipt_path = folder/"RECEIPT.json"
            receipt_source = source/receipt_path
            if receipt_source.is_file():
                costs = load_json(receipt_source)
                archive.json(R13/receipt_path, receipt_source,
                             {k:costs[k] for k in BASE_COST_FIELDS if k in costs},
                             "historical_cost")
        if i in (0, 1):
            path = Path("RESULTS/confirmation_core48")/f"{case}__cie"/"corrected.npz"
            archive.npz(R13/path, source/path, ("m",), "online_historical_endpoint")


def r6_assets(archive):
    source = historical(archive.workspace, "r6_worktree", R6)
    source_root = source.parents[1]
    cases = {"obs_b20d5a55677d4eca", "obs_7b8d5df517f1497d"}
    arms = ("cie_single", "cie_dual", "a_material", "global_dual")
    manifest = source/"configs/primary_manifest.json"
    keys = ("case_id", "arm", "public_files", "output_dir")
    tasks = [select(t, keys) for t in load_json(manifest)["tasks"]
             if t["case_id"] in cases and t["arm"] in arms]
    assert len(tasks) == 8
    archive.json(R6/"configs/primary_manifest.json", manifest, {"tasks":tasks}, "online_index")
    for path in sorted({Path(p) for t in tasks for p in t["public_files"]}):
        fields = (*R13_FIELDS, "E", "noise_radius")
        archive.npz(path, source_root/path, fields, "online_observation", partial=True)
    for t in tasks:
        stem = Path(t["output_dir"])
        path = stem/"prediction.npz"
        archive.npz(path, source_root/path, ENDPOINT_FIELDS, "online_historical_endpoint")
        path = stem/"RECEIPT.json"
        original = load_json(source_root/path)
        archive.json(path, source_root/path, {k:original[k] for k in BASE_COST_FIELDS if k in original},
                     "historical_cost")
    case_source = source/"data_public/cases.json"
    public_cases = [select(t, ("case_id", "observation_file", "amplitude"))
                    for t in load_json(case_source)["records"] if t["case_id"] in cases]
    archive.json(R6/"data_public/cases.json", case_source, {"records":public_cases}, "online_index")
    mapping_source = source/"sealed_truth/scoring/historical_mapping_PRIVATE.json"
    mapping = []
    for original in load_json(mapping_source):
        if original["case_id"] not in cases:
            continue
        path = Path(original["truth_path"])
        archive.npz(path, source_root/path, ("m", "mask"), "offline_label")
        mapping.append({k:original[k] for k in ("case_id", "truth_path", "object_id") if k in original})
    archive.json(R6/"sealed_truth/scoring/offline_replay_mapping.json", mapping_source, mapping,
                 "offline_index")


def scalar_legacy_assets(archive):
    r15 = Path("research/research_max_struct_r15")
    source = historical(archive.workspace, "r14_worktree", r15)
    for case in ("r14_new_0000", "r14_new_0002", "r14_new_0003"):
        path = Path("data_public/independent24")/(case+".npz")
        archive.npz(r15/path, source/path, R13_FIELDS, "online_observation")

    # R16 load_cohort validates every member before HC-CIE selects the submitted cases.
    for round_id, views, group, prefix, count in (
        ("r10", ("",), "twod", "r10_twod", 24),
        ("r11", ("cross48to24","matched24"), "crossgrid8_paired", "r11_crossgrid8", 8)):
        study = Path("research")/f"research_max_struct_{round_id}"
        source = historical(archive.workspace, f"{round_id}_worktree", study)
        directory = Path("data_public")/group
        path = directory/"operators.npz"
        archive.npz(study/path, source/path, MODEL_FIELDS, "online_operator")
        for view in views:
            for i in range(count):
                path = directory/view/f"{prefix}_{i:04d}.npz"
                archive.npz(study/path, source/path, COHORT_FIELDS, "online_observation")


def srap_assets(archive):
    study = archive.workspace/SRAP
    spectral = load_json(HERE/"configs/spectra.json")
    obs = {f"{obj}__{layout}_10x20": ("srap" if (obj,layout)==("s2d_010","ls_arc") else "dual")
           for obj in spectral["srap_objects"] for layout in spectral["srap_layouts"]}
    # The first object of each frozen family is regenerated separately by model_gap.
    gap = {f"s2d_{i:03d}__{layout}_10x20" for i in (1, 2, 3, 4, 7, 15)
           for layout in ("pw360","ls_line_same")}
    obs.update({case:"dual" for case in sorted(gap)})
    for case, arm in sorted(obs.items()):
        path = Path("data/public")/(case+".npz")
        archive.npz(SRAP/path, study/path, ("case_json","E","S","Y","centers","k","radius","delta"),
                    "online_observation")
        path = Path("results/predictions")/f"{case}__{arm}.npz"
        archive.npz(SRAP/path, study/path, ("m_30_light",), "online_historical_endpoint")
        if case in gap:
            path = Path("data/truth")/(case+".npz")
            archive.npz(SRAP/path, study/path, ("m","family"), "offline_label")

    for layout in ("contiguous", "interleaved"):
        case = f"fresnel_twin_8GHz__{layout}"
        path = Path("data/public")/(case+".npz")
        archive.npz(SRAP/path, study/path, ("case_json","E","S","Y","centers","k","radius",
                    "rx_ids","valid_mask","delta"), "online_public_processed_Fresnel")
        path = Path("results/lazy_fresnel")/f"{case}__dual.npz"
        archive.npz(SRAP/path, study/path, ("m_30_light",), "online_historical_endpoint")

    for array in ("hemi", "plane_grid"):
        for i in range(12):
            case = f"{array}_{i:02d}"
            path = Path("vector/data")/(case+".npz")
            archive.npz(SRAP/path, study/path, ("Y","delta"), "online_observation")
            path = Path("vector/predictions")/(case+"_Dual")/"materials.npz"
            archive.npz(SRAP/path, study/path, ("light_180",), "online_historical_endpoint")


def classical_assets(archive):
    top = Path("research/top_venue_killtest_20261120")
    runtime = top/".runs"
    config = load_json(HERE/"configs/classical_v1.json")
    observed = ("E", "S", "Y", "centers", "k", "radius", "delta", "case_json")
    labels = {method+config["parent_suffix"] for method in config["parent_methods"]}
    cost_keys = ("key", "case", "suite", "method", "budget", "view", "wall_s", "cpu_s",
                 "status", "iterations", "proposal_status", "proposal_iterations", "LU", "FFT",
                 "complex_RHS", "phi", "relative_data_residual", "candidate_wall",
                 "candidate_step", "candidate_stage")
    for case in config["observation_ids"]:
        if case.startswith("austria"):
            path = runtime/"public"/(case+".npz")
            archive.npz(path, archive.workspace/path, observed, "online_observation")
        # The separately licensed original Foam experimental measurements are supplied by the reader.
        path = runtime/"pred"/(case+".npz")
        source = archive.workspace/path
        with np.load(source, allow_pickle=False) as saved:
            subset = {key: saved[key] for label in sorted(labels) for key in (label,label+"_Y")}
            historical = json.loads(str(saved["rows_json"]))
        rows = [{key:row[key] for key in cost_keys if key in row}
                for row in historical if row["key"] in labels]
        if len(rows) != len(labels):
            raise ValueError(f"Missing historical receipts in {source}")
        subset["rows_json"] = np.asarray(json.dumps(rows, ensure_ascii=False))
        staged = BytesIO()
        np.savez_compressed(staged, **subset)
        archive.put(path, staged.getvalue(), "online_historical_endpoint", source, arrays={
            key:dict(dtype=str(arr.dtype),shape=list(arr.shape)) for key,arr in subset.items()})


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--workspace",type=Path,required=True)
    p.add_argument("--toy-source",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args = p.parse_args()
    workspace = args.workspace.resolve()
    dest = args.out.resolve()
    dest.parent.mkdir(parents=True,exist_ok=True)
    stage = dest.with_suffix(".zip.part")
    t0 = perf_counter()
    with zipfile.ZipFile(stage,"w",allowZip64=True) as zipped:
        archive = Archive(zipped,workspace)
        r13_assets(archive)
        r6_assets(archive)
        scalar_legacy_assets(archive)
        srap_assets(archive)
        classical_assets(archive)
        for reference in ("r13_fresh_0000", "r13_fresh_0001"):
            source = HERE/"results/real_weights"/f"r13_{reference}.npz"
            member = Path("research/hc_cie_r2/replay_inputs/real_frontend")/(reference+".npz")
            archive.npz(member, source, ("minimum_real_weight",), "online_frozen_weight")
        for name,keys in TOY_NUMERICAL_FIELDS.items():
            source = args.toy_source.resolve()/"research_work"/name
            archive.json(Path("research/hc_cie_r2/replay_inputs/toy/research_work")/name,
                         source,select(load_json(source),keys),"toy_numeric_fixture")
        manifest = dict(instructions="research/hc_cie_r2/INPUTS.md",
                        extraction="repository root; use unzip -n to preserve files already tracked",
                        contents=archive.entries)
        zipped.writestr("research/hc_cie_r2/replay_inputs/manifest.json",
                        json.dumps(manifest,ensure_ascii=False,indent=2)+"\n",
                        compress_type=zipfile.ZIP_DEFLATED)
    stage.replace(dest)
    print(f"{len(archive.entries)} inputs, {dest.stat().st_size:,} bytes, {perf_counter()-t0:.2f} s")


if __name__ == "__main__":
    main()
