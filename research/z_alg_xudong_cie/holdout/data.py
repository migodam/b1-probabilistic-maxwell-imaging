"""Public-only loader and blind identifiers for the B3 paired acquisition."""

import argparse
import json
from pathlib import Path
from time import perf_counter, process_time

import numpy as np

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parents[1]
WORKSPACE = RESEARCH.parent
OPTICAL = RESEARCH / "z_alg_optical_inertia"
ARCHIVE = OPTICAL / "synthetic/public"
DEFAULT_ROOT = WORKSPACE / ".artifacts/z_alg_xudong_cie/b3_holdout"
FREQUENCY = 400e6
LIGHT_SPEED = 3e8
WIDTH = 2.0
INVERSE_SIDE, NATIVE_SIDE = 32, 96
K = 2 * np.pi * FREQUENCY / LIGHT_SPEED
OBJECT_IDS = tuple(f"b3m_{i:04d}" for i in range(1, 13))
APERTURES = ("PW360_10x20", "PW180_10x20")
NOISE_TAGS = {0.02: "02", 0.10: "10"}


def observation_id(oid, aperture, noise):
    if oid not in OBJECT_IDS or aperture not in APERTURES:
        raise ValueError("Unknown planned material ID or aperture.")
    try:
        tag = NOISE_TAGS[float(noise)]
    except KeyError as exc:
        raise ValueError("Use noise fractions 0.02 or 0.10.") from exc
    return f"{oid}__{aperture.lower()}__n{tag}"


def planned_cases():
    for oid in OBJECT_IDS:
        for aperture in APERTURES:
            for noise in NOISE_TAGS:
                name = observation_id(oid, aperture, noise)
                yield dict(case=name, object_id=oid, suite="LA_SYNTH", cohort="confirmation",
                           array=aperture, side=INVERSE_SIDE, width=WIDTH, fine=NATIVE_SIDE,
                           frequency=FREQUENCY, noise_relative=noise,
                           tx_count=10, rx_count=20, train_channels=200,
                           nominal_aperture_deg=int(aperture[2:5]))


def public_directory(root):
    return Path(root) / "public"


def manifest_path(root):
    return public_directory(root) / "manifest.json"


def public_geometry(root, side):
    if side not in (INVERSE_SIDE, NATIVE_SIDE):
        raise ValueError("This acquisition has 32² inverse and 96² native grids.")
    return public_directory(root) / f"geometry{side}.npz"


def get_model(aperture, root=DEFAULT_ROOT, *, side=INVERSE_SIDE):
    """Load the public operator once per aperture/grid; no measurements are read."""
    if aperture not in APERTURES:
        raise ValueError("Unknown frozen fitting aperture.")
    geometry = public_geometry(root, side)
    with np.load(geometry, allow_pickle=False) as src:
        out = {key: src[key] for key in ("centers", "radius", "k", "spectrum")}
        out["E"], out["S"] = (src[f"{aperture}__{key}"] for key in ("B", "S"))
    out.update(B=out["E"], E0=out["E"], grid_shape=(side, side),
               fft_geometry=dict(coords=out["centers"], k=float(out["k"]), radius=float(out["radius"])))
    return out


def get_native_model(aperture, root=DEFAULT_ROOT):
    return get_model(aperture, root=root, side=NATIVE_SIDE)


def get_case(name, root=DEFAULT_ROOT, *, side=INVERSE_SIDE):
    """Return engine-compatible metadata and fitting operators, with no private reads."""
    row = next((x for x in planned_cases() if x["case"] == name), None)
    if row is None:
        raise KeyError(f"Unknown planned B3 observation {name}.")
    case = dict(row, side=side)
    geometry = public_geometry(root, side)
    manifest = manifest_path(root)
    if not manifest.is_file():
        raise FileNotFoundError(f"No sealed public manifest exists at {manifest}.")
    aperture = row["array"]
    observation = public_directory(root) / "observations" / f"{name}.npz"
    out = get_model(aperture, root, side=side)
    with np.load(observation, allow_pickle=False) as src:
        if set(src.files) != {"Y", "delta", "raw_radius"}:
            raise ValueError("A fitting observation must contain only the public measurement and radii.")
        out["Y"], out["delta"], out["raw_radius"] = (src[key] for key in ("Y", "delta", "raw_radius"))
    if (out["Y"].shape != (20, 10) or out["E"].shape != (side**2, 10) or
            out["S"].shape != (20, side**2) or not np.isfinite(out["Y"]).all() or
            not 0 < out["raw_radius"] <= out["delta"]):
        raise ValueError("Public B3 arrays do not match the frozen measurement contract.")
    case["measurement_files"] = [str(p.relative_to(Path(root))) for p in (observation, geometry, manifest)]
    return case, out


def get_native_case(name, root=DEFAULT_ROOT):
    """Share the public fitting measurement with the exact archived 96² geometry."""
    return get_case(name, root=root, side=NATIVE_SIDE)


def load_public(oid, aperture, noise, root=DEFAULT_ROOT):
    return get_case(observation_id(oid, aperture, noise), root=root)


def require_freeze(path):
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"A prior algorithm freeze record is required: {path}")
    return path.resolve()


def timers():
    return perf_counter(), process_time()


def elapsed(start):
    return dict(zip(("wall_seconds", "cpu_seconds"),
                    (x - y for x, y in zip(timers(), start, strict=True))))


def write_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    stage = path.with_suffix(".tmp")
    stage.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    stage.replace(path)


def write_npz(path, **arrays):
    path.parent.mkdir(parents=True, exist_ok=True)
    stage = path.with_suffix(".tmp")
    with stage.open("wb") as out:
        np.savez_compressed(out, **arrays)
    stage.replace(path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--list", action="store_true", required=True)
    parser.parse_args()
    print(json.dumps(list(planned_cases()), indent=2))
