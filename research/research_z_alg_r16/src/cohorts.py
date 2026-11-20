"""Load original scalar-TM observations. Material is available only via load_truth."""

from os import environ
from pathlib import Path

import numpy as np


# (historical round, subgroup directory, case prefix, number of observations)
COHORTS = {
    "r10_twod": (10, "twod", "r10_twod", 24),
    "r10_conditional4": (10, "conditional_confirmation", "r10_conditional", 4),
    "r11_independent48": (11, "independent48", "r11_independent", 48),
    "r11_bothblocks16": (11, "bothblocks16", "r11_bothblocks16", 16),
    "r11_crossgrid8": (11, "crossgrid8_paired", "r11_crossgrid8", 8),
    "r11_crossgrid8_matched": (11, "crossgrid8_paired", "r11_crossgrid8", 8),
    "r12_dev6": (12, "dev6", "r12_dev6", 6),
    "r12_confirm24": (12, "confirm24", "r12_confirm24", 24),
}
R3_CURVED_BAR_ID = "obs_b20d5a55677d4eca"


def historical_study(round_number):
    """Resolve a checkout or the directory containing the original worktrees."""
    source = environ.get("R16_MAX_STRUCT_ROOT")
    roots = (
        [Path(source).expanduser()]
        if source
        else [
            parent / "max_struct_r1"
            for parent in Path(__file__).resolve().parents
            if (parent / "max_struct_r1").is_dir()
        ]
    )
    folder = f"research_max_struct_r{round_number}"
    if not source:
        local = next(
            (
                p
                for p in Path(__file__).resolve().parents
                if (p / "research" / folder).is_dir()
            ),
            None,
        )
        if local is not None:
            roots.append(local)
    for root in roots:
        for path in (
            root / f"r{round_number}_worktree" / "research" / folder,
            root / "research" / folder,
            root / folder,
        ):
            if path.is_dir():
                return path
    raise FileNotFoundError(
        f"Cannot find {folder}. Set R16_MAX_STRUCT_ROOT to the historical checkout "
        "or the directory containing the MAX-STRUCT worktrees."
    )


def _cohort_paths(name):
    round_number, subgroup, prefix, count = COHORTS[name]
    study = historical_study(round_number)
    if round_number == 12:
        study = study / "delegated" / "confirmation"
    return study, subgroup, tuple(f"{prefix}_{i:04d}" for i in range(count))


def _read_model(archive):
    model = {
        key: archive[key] for key in ("G", "S", "E0", "coords", "sources", "receivers")
    }
    model["shape"] = tuple(int(i) for i in archive["grid_shape"])
    model["k"] = float(archive["k"])
    model["radius"] = float(archive["radius"])
    return model


def load_cohort(name):
    """Return (shared model, cases); each case is only an ID, raw noisy Y, and sigma."""
    study, subgroup, case_ids = _cohort_paths(name)
    operators = study / "data_public" / subgroup
    observations = operators
    if name == "r11_crossgrid8":
        observations = operators / "cross48to24"
    elif name == "r11_crossgrid8_matched":
        observations = operators / "matched24"

    with np.load(operators / "operators.npz", allow_pickle=False) as archive:
        model = _read_model(archive)

    cases = []
    bounds = None
    expected_shape = (model["S"].shape[0], model["E0"].shape[1])
    for case_id in case_ids:
        path = observations / f"{case_id}.npz"
        with np.load(path, allow_pickle=False) as archive:
            y = archive["Y"]
            sigma = float(archive["sigma"])
            case_bounds = tuple(float(x) for x in archive["material_bounds"])
            if bounds is None:
                bounds = case_bounds
            if (
                y.shape != expected_shape
                or not np.isfinite(y).all()
                or not np.isfinite(sigma)
                or sigma <= 0
                or case_bounds != bounds
                or tuple(archive["grid_shape"]) != model["shape"]
            ):
                raise ValueError(f"Inconsistent observation schema in {path}")
        cases.append({"id": case_id, "Y": y, "sigma": sigma})

    model["bounds"] = bounds
    return model, cases


def load_truth(name, case_id):
    """Offline scoring only. Returns flat real contrast; ignores other archive fields."""
    study, subgroup, case_ids = _cohort_paths(name)
    if case_id not in case_ids:
        raise KeyError(f"{case_id!r} is not a member of {name!r}")
    with np.load(
        study / "sealed_truth" / subgroup / f"{case_id}.npz", allow_pickle=False
    ) as archive:
        return archive["m"]


def load_r3_curved_bar():
    """Return the exposed MAX-STRUCT R3 stress case with its distinct R2 geometry."""
    path = historical_study(2) / "data_public" / f"{R3_CURVED_BAR_ID}.npz"
    with np.load(path, allow_pickle=False) as archive:
        model = _read_model(archive)
        model["bounds"] = tuple(float(x) for x in archive["material_bounds"])
        case = {
            "id": R3_CURVED_BAR_ID,
            "Y": archive["Y"],
            "sigma": float(archive["sigma"]),
        }
    return model, [case]


def load_r3_curved_bar_truth():
    """Offline only. Return material and the saved area-fraction support mask."""
    path = historical_study(2) / "sealed_truth" / f"{R3_CURVED_BAR_ID}.npz"
    with np.load(path, allow_pickle=False) as archive:
        return archive["m"], archive["mask"].ravel()
