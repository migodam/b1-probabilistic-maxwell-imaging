"""Reacquire the unchanged lossless holdout on 72×72 cells."""

import json
from time import perf_counter, process_time

import numpy as np

from acquire import COARSE, FINE, QUADRATURE_SIDE, RESULTS, material, relative
from holdout import cases, forward_batch
from ndt import PLATE_PERMITTIVITY, PLATE_SIDE, host


FINE72 = 72
SUBPIXELS = FINE * QUADRATURE_SIDE // FINE72


def material_integral(images):
    step_mm = 1e3 * PLATE_SIDE / images.shape[-1]
    return (images.sum(axis=(1, 2)) * step_mm**2).tolist()


def main():
    wall, cpu = perf_counter(), process_time()
    definitions = cases()
    truth_fine = np.array([material(d["primitives"], FINE72, SUBPIXELS) for d in definitions])
    ratio = FINE72 // COARSE
    truth = truth_fine.reshape(len(definitions), COARSE, ratio, COARSE, ratio).mean(axis=(2, 4))

    # This is an offline data generator; its reference material is read only to
    # verify that the common geometry quadrature did not change the scoring truth.
    with np.load(RESULTS / "holdout.npz", allow_pickle=False) as reference:
        clean48 = reference["clean"]
        ref_truth = reference["truth"]
        first_integral = material_integral(reference["truth_fine"])
        kept = {key: reference[key] for key in reference.files if key not in {
            "clean", "clean_same", "background", "background_same", "truth", "truth_fine"
        }}
    if not np.array_equal(kept["case_ids"], [d["id"] for d in definitions]):
        raise ValueError("Reference and fine-grid batches differ in case order")
    truth_difference = float(np.max(np.abs(truth - ref_truth)))
    if truth_difference > 1e-12:
        raise ValueError(f"The common geometry quadrature changed 24-grid truth by {truth_difference:g}")

    clean = np.empty_like(clean48)
    same = np.empty_like(clean48)
    background = np.empty_like(clean48[0])
    background_same = np.empty_like(background)
    diag = []
    for i, f in enumerate(kept["frequencies_ghz"]):
        f_start = perf_counter()
        clean[:, i], same[:, i], background[i], background_same[i], stats = forward_batch(
            float(f), truth, truth_fine, PLATE_PERMITTIVITY, 1., host)
        stats["wall_seconds"] = perf_counter() - f_start
        diag.append(stats)
        print(f"{f:5.2f} GHz  {stats['fine_krylov_iterations']:5d} iterations  "
              f"residual {stats['max_fine_relative_state_residual']:.2e}  "
              f"{stats['wall_seconds']:.2f}s", flush=True)

    observations = dict(**kept, clean=clean, clean_same=same,
                        background=background, background_same=background_same)
    np.savez_compressed(RESULTS / "holdout72_observations.npz", **observations)
    np.savez_compressed(RESULTS / "holdout72.npz", **observations, truth=truth, truth_fine=truth_fine)

    def mismatch(other):
        return [None if not d["primitives"] else [relative(y0, y1) for y0, y1 in zip(a, b)]
                for d, a, b in zip(definitions, other, clean)]

    metadata = dict(
        case_definitions=definitions, frequency_ghz=kept["frequencies_ghz"].tolist(),
        coarse_cells_per_side=COARSE, fine_cells_per_side=FINE72,
        fine_cell_geometry_quadrature_samples=SUBPIXELS**2,
        geometry_global_samples_per_side=FINE72 * SUBPIXELS,
        common_scoring_material_24_max_absolute_difference=truth_difference,
        contrast_integral_mm2_on_48=first_integral, contrast_integral_mm2_on_72=material_integral(truth_fine),
        relative_cross_grid_48_to_72=mismatch(clean48), relative_cross_grid_24_to_72=mismatch(same),
        diagnostics=diag, wall_seconds=perf_counter() - wall, cpu_seconds=process_time() - cpu,
    )
    (RESULTS / "holdout72.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
