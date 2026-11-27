"""Load only the public inputs named by the frozen spectral design."""
from dataclasses import dataclass
import importlib.util
import json
import os
from pathlib import Path
import sys

import numpy as np


@dataclass
class Observation:
    key: str
    cohort: str
    physical_id: str
    g: np.ndarray
    s: np.ndarray
    b: np.ndarray
    y: np.ndarray
    delta: float
    sigma: float | None
    lower: float
    upper: float
    base: np.ndarray
    provenance: dict


def _relative(path, root):
    return str(path.relative_to(root))


def _load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _study(root, fallback, name):
    direct = root / "research" / name
    return direct if direct.is_dir() else root / fallback / "research" / name


def load_observations(root, config):
    historical = Path("research/delegated/max_struct_r1")
    for tag, round_dir, directory in (
        ("r13", "r13_worktree", "fresh48"),
        ("r15", "r14_worktree", "independent24"),
    ):
        study = _study(root, historical / round_dir, f"research_max_struct_{tag}")
        for case in config[f"{tag}_cases"]:
            path = study / "data_public" / directory / f"{case}.npz"
            with np.load(path, allow_pickle=False) as archive:
                yield Observation(
                    f"{tag}_{case}", tag, case, archive["G"], archive["S"], archive["E0"],
                    archive["Y"], float(archive["delta"]), float(archive["sigma"]),
                    float(archive["lower"]), float(archive["upper"]), archive["m_init"],
                    dict(observation=_relative(path, root), operator=_relative(path, root),
                         base="public m_init", forward_side=int(archive["forward_side"]),
                         inverse_side=int(archive["inverse_side"]),
                         source_receiver_min_distance=float(np.min(np.linalg.norm(
                             archive["sources"][:, None] - archive["receivers"][None, :], axis=2)))),
                )

    r16 = _study(root, Path("research/delegated/z_alg_r5_delivery/r16_publication"),
                 "research_z_alg_r16")
    loader_path = r16 / "src" / "cohorts.py"
    loader = _load_module("hc_cie_r16_cohorts", loader_path)
    r16_input = root if (root / "research/research_max_struct_r10").is_dir() else root / historical
    os.environ.setdefault("R16_MAX_STRUCT_ROOT", str(r16_input))
    for name, requested in config["r16_cases"].items():
        model, cases = loader.load_cohort(name)
        study, subgroup, _ = loader._cohort_paths(name)
        operator_path = study / "data_public" / subgroup / "operators.npz"
        view = {"r11_crossgrid8": "cross48to24", "r11_crossgrid8_matched": "matched24"}.get(name, "")
        for case in cases:
            if case["id"] not in requested:
                continue
            path = operator_path.parent / view / f"{case['id']}.npz"
            with np.load(path, allow_pickle=False) as archive:
                radius = float(archive["noise_radius"])
            base = 0.1 + 3.0 * (np.linalg.norm(model["coords"], axis=1) <= 0.20)
            yield Observation(
                f"r16_{name}_{case['id']}", "r16", case["id"], model["G"], model["S"],
                model["E0"], case["Y"], radius * case["sigma"], case["sigma"],
                *model["bounds"], base,
                dict(observation=_relative(path, root), operator=_relative(operator_path, root),
                     loader=_relative(loader_path, root), base="0.1 plus centered disk r=0.20 peak=3.1",
                     forward_side=48 if name == "r11_crossgrid8" else 24, inverse_side=24,
                     source_receiver_min_distance=float(np.min(np.linalg.norm(
                         model["sources"][:, None] - model["receivers"][None, :], axis=2)))),
            )

    srap = root / "research/srap_cie_array_study_20261121"
    os.environ.setdefault("LA_WORKSPACE", str(root))
    legacy = _load_module("hc_cie_legacy_la", srap / "legacy_la.py")
    dense_green = None
    for obj in config["srap_objects"]:
        for layout in config["srap_layouts"]:
            case = f"{obj}__{layout}_10x20"
            path = srap / "data/public" / f"{case}.npz"
            arm = "srap" if (obj, layout) == ("s2d_010", "ls_arc") else "dual"
            prediction = srap / "results/predictions" / f"{case}__{arm}.npz"
            with np.load(prediction, allow_pickle=False) as saved:
                base = saved["m_30_light"]
            with np.load(path, allow_pickle=False) as archive:
                meta = json.loads(str(archive["case_json"]))
                if dense_green is None:
                    public = {k: archive[k] for k in ("centers", "k", "radius", "E", "S")}
                    model = legacy.PaperModel(public, meta, legacy.Clock(3600.0))
                    dense_green = model.G
                    probe = archive["E"][:, :2]
                    if np.linalg.norm(dense_green @ probe - model.green.apply(probe)) > 1e-12 * np.linalg.norm(probe):
                        raise ValueError("Dense SRAP Green and the historical FFT action disagree")
                tx = meta.get("tx_line_source_m")
                receiver_distance = None if tx is None else float(np.min(np.linalg.norm(
                    np.asarray(tx)[:, None] - np.asarray(meta["rx_m"])[None, :], axis=2)))
                yield Observation(
                    f"srap_{case}", "srap", obj, dense_green, archive["S"], archive["E"],
                    archive["Y"], float(archive["delta"]), None, 0.0, 20.0, base,
                    dict(observation=_relative(path, root), operator="legacy_la.PaperModel.G, same 32² Green across layouts",
                         base=f"historical {arm} 30 s light", base_path=_relative(prediction, root),
                         forward_side=meta["fine"], inverse_side=meta["side"], layout=meta["array"],
                         tx_kind=meta["tx_kind"], source_receiver_min_distance=receiver_distance),
                )
