"""Changed-cell rank necessary conditions on stored 3D and Fresnel data."""
import argparse
import csv
import itertools
import json
import os
from pathlib import Path
import time

import spectra as sp  # pins native BLAS threads before loading NumPy
import numpy as np
from scipy.linalg import svdvals

from spectra_inputs import _load_module


def lower_bound(y, prediction, radius, components, padding=0.0):
    residual = np.asarray(y, dtype=complex) - np.asarray(prediction, dtype=complex)
    singular_values = svdvals(residual, check_finite=False)
    tail = np.r_[np.sqrt(np.cumsum(singular_values[::-1]**2))[::-1], 0.0]
    rounded_radius = radius + padding + 32 * np.finfo(float).eps * np.linalg.norm(residual)
    minimum_rank = int(np.count_nonzero(tail > rounded_radius))
    minimum_cells = (minimum_rank + components - 1) // components
    return dict(minimum_changed_cells=minimum_cells, minimum_response_rank=minimum_rank,
                residual_over_radius=float(tail[0] / radius), radius=radius,
                numerical_padding=padding, residual_singular_values=singular_values.tolist(),
                best_rank_residual_norms=tail.tolist(), receive_channels=y.shape[0], sources=y.shape[1])


def vector_audits(root, srap, rows, arrays):
    os.environ.setdefault("SRAP_DEVICE", "cpu")
    model_path = srap / "vector/physics.py"
    physics = _load_module("hc_cie_srap_vector_physics", model_path)
    physics.CD = physics.torch.complex128
    for array in ("hemi", "plane_grid"):
        model_started = time.perf_counter()
        model = physics.Model(12, array)
        model_setup = time.perf_counter() - model_started
        for index in range(12):
            case = f"{array}_{index:02d}"
            source = srap / f"vector/data/{case}.npz"
            endpoint = srap / f"vector/predictions/{case}_Dual/materials.npz"
            with np.load(source, allow_pickle=False) as saved:
                y, radius = saved["Y"], float(saved["delta"])
            with np.load(endpoint, allow_pickle=False) as saved:
                material = saved["light_180"]
            start, cpu = time.perf_counter(), time.process_time()
            d = physics.t(material).real.repeat_interleave(3)[:, None]
            counters = model.count.copy()
            j_early, coarse_error = model.solve(physics.t(material).real, d * model.E, tol=2e-5)
            j, error = model.solve(physics.t(material).real, d * model.E, tol=1e-10)
            predict = physics.npy(model.S @ j)
            epsilon = 2 * np.linalg.norm(predict - physics.npy(model.S @ j_early))
            rows.append(dict(case=case, array=array, anchor="air", dimension=3, model_size="12^3",
                             observation=str(source.relative_to(root)), endpoint=None,
                             **lower_bound(y, np.zeros_like(y), radius, components=3)))
            rows.append(dict(case=case, array=array, anchor="Dual light_180", dimension=3, model_size="12^3",
                             observation=str(source.relative_to(root)), endpoint=str(endpoint.relative_to(root)),
                             operator=str(model_path.relative_to(root)),
                             model_setup_wall_s=model_setup if index == 0 else 0.0,
                             wall_s=time.perf_counter() - start, cpu_s=time.process_time() - cpu,
                             loose_original_equation_relative=float(coarse_error),
                             tight_original_equation_relative=float(error),
                             original_counters=dict(model.count - counters),
                             **lower_bound(y, predict, radius, components=3, padding=float(epsilon))))
            arrays[f"{case}_Y"] = y
            arrays[f"{case}_Dual_recomputed"] = predict
            print(case, rows[-2]["minimum_changed_cells"], rows[-1]["minimum_changed_cells"], flush=True)


def common_rectangle(receive_sets):
    """Most measured entries with a common receiver set, using only the mask."""
    best = ()
    best_channels = ()
    for count in range(2, len(receive_sets) + 1):
        for sources in itertools.combinations(range(len(receive_sets)), count):
            receivers = tuple(sorted(set.intersection(*(receive_sets[i] for i in sources))))
            if len(sources) * len(receivers) > len(best) * len(best_channels):
                best, best_channels = sources, receivers
    if len(best_channels) < 2:
        raise ValueError("No usable common measured rectangle")
    return np.asarray(best), np.asarray(best_channels)


def fresnel_audits(root, srap, rows, arrays):
    os.environ.setdefault("LA_WORKSPACE", str(root))
    model_path = srap / "legacy_la.py"
    legacy = _load_module("hc_cie_fresnel_history", model_path)
    for layout in ("contiguous", "interleaved"):
        start, cpu = time.perf_counter(), time.process_time()
        case = f"fresnel_twin_8GHz__{layout}"
        source = srap / f"data/public/{case}.npz"
        endpoint = srap / f"results/lazy_fresnel/{case}__dual.npz"
        with np.load(source, allow_pickle=False) as saved:
            meta = json.loads(str(saved["case_json"]))
            public = {k: saved[k] for k in ("E", "S", "Y", "centers", "k", "radius", "rx_ids", "valid_mask", "delta")}
        with np.load(endpoint, allow_pickle=False) as saved:
            material = saved["m_30_light"]
        receive_sets = [set(column[valid]) for column, valid in zip(public["rx_ids"], public["valid_mask"].T, strict=True)]
        sources, receivers = common_rectangle(receive_sets)
        y = np.empty((len(receivers), len(sources)), dtype=complex)
        for output, j in enumerate(sources):
            rows_by_id = {value: i for i, value in enumerate(public["rx_ids"][j]) if public["valid_mask"][i, j]}
            y[:, output] = public["Y"][[rows_by_id[i] for i in receivers], j]
        model = legacy.PaperModel(public, meta, legacy.Clock(1800))
        j, _ = model._solve_many(material, material[:, None] * public["E"][:, sources])
        selected_receive = public["S"][receivers]
        pred = selected_receive @ j
        radius = float(public["delta"])
        equation = material[:, None] * public["E"][:, sources]
        closure = np.linalg.norm(j - material[:, None] * model.green.apply(j) - equation) / np.linalg.norm(equation)
        common = dict(case=case, array=layout, dimension=2, model_size="60^2",
                      observation=str(source.relative_to(root)), operator=str(model_path.relative_to(root)),
                      retained_train_source_indices=sources.tolist(), original_rx_indices=receivers.tolist(),
                      train_channels_total=int(public["valid_mask"].sum()),
                      radius_scope="same 15% of all measured training channels, retained without rescaling")
        rows.append(dict(**common, anchor="air", endpoint=None,
                         **lower_bound(y, np.zeros_like(y), radius, components=1)))
        rows.append(dict(**common, anchor="Dual light_30", endpoint=str(endpoint.relative_to(root)),
                         original_equation_relative=float(closure),
                         original_counters=dict(model.clock.counts), wall_s=time.perf_counter() - start,
                         cpu_s=time.process_time() - cpu, **lower_bound(y, pred, radius, components=1)))
        arrays[f"{case}_Y"] = y
        arrays[f"{case}_Dual_recomputed"] = pred
        print(case, len(receivers), len(sources), rows[-2]["minimum_changed_cells"],
              rows[-1]["minimum_changed_cells"], flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, default=Path.cwd())
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    root = args.workspace.resolve()
    output = args.output or Path(__file__).resolve().parents[1] / "results/spectra"
    output.mkdir(parents=True, exist_ok=True)
    srap = root / "research/srap_cie_array_study_20261121"
    rows, arrays = [], {}
    clock, cpu = time.perf_counter(), time.process_time()
    vector_audits(root, srap, rows, arrays)
    fresnel_audits(root, srap, rows, arrays)
    np.savez_compressed(output / "physical_limits.npz", **arrays)
    with (output / "physical_limits.csv").open("w", newline="", encoding="utf8") as handle:
        names = ["case", "array", "anchor", "dimension", "minimum_changed_cells", "minimum_response_rank",
                 "receive_channels", "sources", "residual_over_radius", "radius", "wall_s", "cpu_s"]
        writer = csv.DictWriter(handle, fieldnames=names, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    (output / "physical_limits.json").write_text(json.dumps(dict(
        wall_s=time.perf_counter() - clock, cpu_s=time.process_time() - cpu,
        environment="native single-thread CPU, vector kernel reused in complex128",
        vector_padding="twice absolute port difference between tol=2e-5 and tol=1e-10; an empirical numerical margin",
        rows=rows), indent=2) + "\n")


if __name__ == "__main__":
    main()
