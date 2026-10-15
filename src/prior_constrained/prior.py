"""Original A25 TRAIN-distribution materials; no physics dependency.

P0_EXACT denotes a replayable generator endpoint, never just a phase test.
The original deterministic rasterizer is loaded as its AST function, without
importing the dataset module (and therefore without any Maxwell entry point).
"""
from __future__ import annotations
import ast
import json
from pathlib import Path
from functools import lru_cache
import numpy as np
from scipy import ndimage

FAMILIES = ("shell", "asymmetric_piecewise", "smooth", "mixed_scale")
V = complex(1, .2) / abs(complex(1, .2))
AMPLITUDES = (.3, 1.2)
THRESHOLDS = (.1, .25, .5, .75, .9)
SCALES = (0., 1., 2.)
ROOT = Path(__file__).resolve().parents[2]
A26 = ROOT.parents[1] / "A26/a26-om-feedback-prior-completion"
A30 = ROOT.parents[1] / "A30/a30-maxwell-material-manifold"
INPUTS = A30 / "inputs/a29_release/inputs/a25"
SNAPSHOTS = ROOT / "research/prior_constrained_20261015/source_snapshots"


def config():
    path = A26 / "configs/a25/campaign.json"
    return json.loads((path if path.exists() else SNAPSHOTS / "campaign.json").read_text())


def geometry_points():
    """Public fixed grid; no material or locked geometry is read."""
    path = INPUTS / "geometry.json"
    return np.asarray(json.loads((path if path.exists() else SNAPSHOTS / "geometry.json").read_text())["points"], float)


def background(name, points=None):
    p = geometry_points() if points is None else np.asarray(points)
    c = config()
    z = np.full(len(p), complex(*c["background_weak"]), complex)
    if name == "ring":
        r = np.linalg.norm(p, axis=1)
        z[(r >= c["ring_radius_inner"]) & (r <= c["ring_radius_outer"])] = complex(*c["background_ring"])
    elif name != "weak":
        raise ValueError("unknown frozen background")
    return z


def draw_train_row(master_seed, family_index, local_index, *, parent_id=None):
    """Exact TRAIN branch of generate_manifest; excludes every locked branch.

    Original RNG draw order and all parameter supports are unchanged. A new
    master seed creates fresh materials with the same TRAIN geometry law.
    """
    fi, local = int(family_index), int(local_index)
    if not (0 <= fi < 4) or local < 0:
        raise ValueError("TRAIN family/local required")
    seed = [int(master_seed), 250, 0, fi, local]
    rng = np.random.default_rng(np.random.SeedSequence(seed))
    angle = float(rng.uniform(-np.pi, np.pi))
    center = rng.uniform(-.12, .12, 2).tolist()
    scale = float(rng.uniform(.78, 1.04))
    ratio = float(rng.uniform(.82, 1.18))
    d = {"center": center, "rotation": angle, "scale": scale, "axis_ratio": ratio,
         "shell_radius": float(rng.uniform(.14, .19)), "shell_width": float(rng.uniform(.035, .055)),
         "disk_radius": float(rng.uniform(.065, .10)), "disk_weight": float(rng.uniform(.45, .8)),
         "second_shell_weight": float(rng.uniform(.45, .75)), "rectangle_width": float(rng.uniform(.17, .25)),
         "rectangle_height": float(rng.uniform(.19, .29)), "subtype": int(rng.integers(0, 2)),
         "gaussian_centers": rng.uniform(-.13, .13, (3, 2)).tolist(),
         "gaussian_widths": rng.uniform(.05, .095, 3).tolist(),
         "gaussian_weights": rng.uniform(.4, .9, 3).tolist(), "segment_offset": float(rng.uniform(-.06, .06))}
    return {"parent_id": parent_id, "split": "train_distribution_material_only", "family": FAMILIES[fi],
            "parent_seed": seed, "geometry": d, "ood_type": "ID", "amplitudes": list(AMPLITUDES)}


@lru_cache(maxsize=1)
def _original_rasterizer():
    path = A26 / "src/a25/dataset.py"
    if not path.exists(): path = SNAPSHOTS / "A25_dataset.py"
    tree = ast.parse(path.read_text())
    fn = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "parent_material")
    module = ast.Module(body=[fn], type_ignores=[])
    ns = {"np": np}
    exec(compile(ast.fix_missing_locations(module), str(path), "exec"), ns)
    return ns["parent_material"]


def decode(row, *, background_id="weak", amplitude=1.2, points=None):
    if row.get("ood_type") != "ID":
        raise ValueError("only TRAIN-distribution ID geometry is authorized")
    points = geometry_points() if points is None else np.asarray(points)
    pattern = _original_rasterizer()(points, row, config())
    return background(background_id, points) + float(amplitude) * pattern


def rho_field(chi, chi_b):
    return np.real(np.conj(V) * (np.asarray(chi) - np.asarray(chi_b)))


def phase_audit(chi, chi_b, amplitude=1.2, tolerance=1e-10):
    """Necessary PHASE conditions, explicitly insufficient for P0_EXACT."""
    x = np.asarray(chi); b = np.asarray(chi_b)
    delta = x - b
    aligned = np.conj(V) * delta
    rho, orth = aligned.real, aligned.imag
    peak = float(np.max(rho))
    scale = max(float(np.linalg.norm(delta)), np.finfo(float).tiny)
    box = bool(np.min(x.real) >= -.5-tolerance and np.max(x.real) <= 4+tolerance and
               np.min(x.imag) >= -tolerance and np.max(x.imag) <= .5+tolerance)
    phase = bool(np.max(abs(orth)) <= tolerance and np.min(rho) >= -tolerance and
                 peak <= float(amplitude)+tolerance and abs(peak-float(amplitude)) <= tolerance)
    return {"box_valid": box, "phase_valid": phase,
            "prior_status": "PHASE" if phase else "PRIOR_UNKNOWN",
            "phase_orth_rms": float(np.sqrt(np.mean(orth**2))),
            "phase_orth_relative": float(np.linalg.norm(orth)/scale),
            "phase_orth_max": float(np.max(abs(orth))),
            "rho_min": float(np.min(rho)), "rho_max": peak,
            "negative_rho_cells": int(np.sum(rho < -tolerance)),
            "negative_rho_l2": float(np.linalg.norm(np.minimum(rho, 0))),
            "over_amplitude_cells": int(np.sum(rho > float(amplitude)+tolerance)),
            "peak_normalization_defect": abs(peak-float(amplitude)),
            "real_background_under_cells": int(np.sum(x.real < b.real-tolerance)),
            "imag_background_under_cells": int(np.sum(x.imag < b.imag-tolerance)),
            "p0_membership": "NOT_ESTABLISHED_BY_PHASE_TEST"}


def parameter_support_audit(row):
    """Support of the original continuous TRAIN parameter law (not a density)."""
    d = row.get("geometry", {})
    scalar = {"rotation": (-np.pi,np.pi), "scale": (.78,1.04), "axis_ratio": (.82,1.18),
              "shell_radius": (.14,.19), "shell_width": (.035,.055), "disk_radius": (.065,.10),
              "disk_weight": (.45,.8), "second_shell_weight": (.45,.75),
              "rectangle_width": (.17,.25), "rectangle_height": (.19,.29),
              "segment_offset": (-.06,.06)}
    array = {"center": ((2,),-.12,.12), "gaussian_centers": ((3,2),-.13,.13),
             "gaussian_widths": ((3,),.05,.095), "gaussian_weights": ((3,),.4,.9)}
    failures=[]
    for name,(lo,hi) in scalar.items():
        v=np.asarray(d.get(name,np.nan))
        if v.shape != () or not np.isfinite(v) or not lo <= float(v) <= hi:
            failures.append(name)
    for name,(shape,lo,hi) in array.items():
        v=np.asarray(d.get(name,[]),float)
        if v.shape != shape or not np.all(np.isfinite(v)) or np.any(v<lo) or np.any(v>hi):
            failures.append(name)
    if d.get("subtype") not in (0,1): failures.append("subtype")
    if row.get("family") not in FAMILIES: failures.append("family")
    if row.get("ood_type") != "ID": failures.append("ood_type")
    return {"parameter_support_valid":not failures,"failed_parameters":failures,
            "measure":"original TRAIN equal four-family mixture and declared independent continuous parameter draws"}


def verify_exact_endpoint(chi, row, *, background_id="weak", amplitude=1.2, tolerance=1e-12):
    support = parameter_support_audit(row)
    if not support["parameter_support_valid"] or float(amplitude) not in AMPLITUDES:
        return {"prior_status":"PRIOR_UNKNOWN","generator_endpoint_max_error":None,**support}
    replay = decode(row, background_id=background_id, amplitude=amplitude)
    defect = float(np.max(abs(np.asarray(chi)-replay)))
    return {"prior_status": "P0_EXACT" if defect <= tolerance else "RETRACTION_NOT_FOUND",
            "generator_endpoint_max_error": defect, "generator_source": "src/a25/dataset.py:parent_material",
            "parent_seed": row["parent_seed"], "parameters": row["geometry"],
            "family": row["family"], "amplitude": float(amplitude), "background_id": background_id, **support}


def structure_descriptor(rho, *, amplitude=1.2):
    """Increment structures, never absolute-air labels; primary grid is 24x24."""
    x = np.asarray(rho, float).reshape(24, 24) / float(amplitude)
    smooth = np.stack([ndimage.gaussian_filter(x, s, mode="nearest") if s else x for s in SCALES])
    masks = smooth[:, None] >= np.asarray(THRESHOLDS)[None, :, None, None]
    boundaries = []
    components, holes = [], []
    fg = ndimage.generate_binary_structure(2, 1)
    bg = ndimage.generate_binary_structure(2, 2)
    for mask in masks.reshape(-1, 24, 24):
        boundaries.append(mask ^ ndimage.binary_erosion(mask, structure=fg))
        labels, n = ndimage.label(mask, structure=fg)
        components.append(sum(np.sum(labels == k) >= 4 for k in range(1, n+1)))
        inv, m = ndimage.label(~mask, structure=bg)
        edge = set(np.r_[inv[0], inv[-1], inv[:, 0], inv[:, -1]])
        holes.append(sum(k not in edge and np.sum(inv == k) >= 4 for k in range(1, m+1)))
    return {"normalized_field": x, "smoothed": smooth, "masks": masks,
            "boundaries": np.asarray(boundaries).reshape(len(SCALES), len(THRESHOLDS),24,24),
            "components": np.asarray(components).reshape(len(SCALES),len(THRESHOLDS)),
            "holes": np.asarray(holes).reshape(len(SCALES),len(THRESHOLDS)),
            "dc": float(np.mean(x))}


def structural_distance(a, b, *, amplitude=1.2):
    da, db = structure_descriptor(a, amplitude=amplitude), structure_descriptor(b, amplitude=amplitude)
    disagreement = np.mean(da["masks"] != db["masks"], axis=(-2,-1))
    return {"continuous_multiscale_rms": float(np.sqrt(np.mean((da["smoothed"]-db["smoothed"])**2))),
            "shape_multiscale_rms": float(np.sqrt(np.mean(((da["smoothed"]-np.mean(da["smoothed"],axis=(-2,-1),keepdims=True))-(db["smoothed"]-np.mean(db["smoothed"],axis=(-2,-1),keepdims=True)))**2))),
            "boundary_multithreshold_mean": float(np.mean(da["boundaries"] != db["boundaries"])),
            "occupancy_multithreshold_mean": float(np.mean(disagreement)),
            "threshold_disagreement": disagreement.tolist(),
            "persistent_thresholds_ge_005": int(np.sum(np.mean(disagreement,axis=0)>=.05)),
            "topology_differing_levels": int(np.sum((da["components"] != db["components"]) | (da["holes"] != db["holes"]))),
            "dc_difference": abs(da["dc"]-db["dc"]),
            "historical_normalized_l2": float(np.sqrt(np.mean((np.asarray(a)-np.asarray(b))**2)) / float(amplitude))}
