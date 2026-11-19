"""Use archive truth only inside the measurement simulator and the later scorer."""
from dataclasses import dataclass
import json
from time import perf_counter

import numpy as np

from data14 import ACTIONS, HALF_SPLIT_FIRST, R10, R12, ROOT, FrequencyTM, complex_reads, load_cases


@dataclass
class Generator:
    truth: np.ndarray
    alternatives: np.ndarray
    gain: complex
    noise: np.ndarray
    reused: dict


def archived_generator():
    generators = {}
    layouts = json.loads((R12 / "results/spectral_new_layouts.json").read_text())
    layout_index = {row["label"]: i for i, row in enumerate(layouts)}
    for i, row in enumerate(json.loads((R12 / "results/spectral_online.json").read_text())):
        case = f"spectral-{i:02d}"
        layout = layout_index[row["label"]]
        side, rep = row["truth_side"], row["rep"]
        pair = np.asarray(layouts[layout]["pair"])
        rng = np.random.default_rng(667029 + rep * 131 + layout * 13 + side)
        base = (rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))) / np.sqrt(2.)
        tail = (rng.normal(size=(2, 3)) + 1j * rng.normal(size=(2, 3))) / np.sqrt(2.)
        reused = {n: complex_reads(row["arms"][old]["extra_ri"])
                  for n, old in (("split_0.02", "fixed_wide_split"), ("repeat", "repeat"))}
        generators[case] = Generator(pair[side], pair[1 - side:2 - side],
                                     complex(*row["common_gain"]), np.vstack((base[1:], tail)), reused)
    for j, file in enumerate(("natural_inputs.npz", "confirm_inputs.npz")):
        with np.load(R10 / "results" / file) as table:
            for i, name in enumerate(sorted(table["case_ids"])):
                truth = table[f"{name}/truth"]
                alternative = table[f"{name}/clean_witnesses"].reshape(-1, 6)
                alternative = alternative[np.linalg.norm(alternative - truth, axis=1) > 1e-3]
                rng = np.random.default_rng(745007331 + j * 100 + i)
                noise = (rng.normal(size=(4, 3)) + 1j * rng.normal(size=(4, 3))) / np.sqrt(2.)
                generators[str(name)] = Generator(
                    truth, alternative, 1. + complex(table[f"{name}/gain_error"]), noise, {})
    return generators


def main():
    # A design from the old observations has been saved before this simulator reads truth.
    choices = json.loads((ROOT / "results/measure_choices.json").read_text())
    if {item["case_id"] for item in choices["cases"]} != {case.case_id for case in load_cases()}:
        raise ValueError("The saved choices must cover the acquired cases.")
    data = {}
    diagnostics = {}
    generators = archived_generator()
    for case in load_cases():
        generator = generators[case.case_id]
        model = FrequencyTM(case.geometry)
        t0 = perf_counter()
        unique_offsets = sorted({d for action in ACTIONS.values() for d in action.deltas})
        clean = {delta: model.response(generator.truth, delta)[0] for delta in unique_offsets}
        reused_max = 0.
        for action in ACTIONS.values():
            if action == HALF_SPLIT_FIRST:
                continue  # This action reuses an exact prefix of split_0.02, with no extra simulation.
            new = np.array([generator.gain * clean[delta] + case.sigma * eps / np.sqrt(e)
                            for delta, e, eps in zip(action.deltas, action.energies, generator.noise)])
            if action.name in generator.reused:
                value = generator.reused[action.name]
                reused_max = max(reused_max, float(np.max(abs(value - new))))
                new = value
            data[f"{case.case_id}/{action.name}"] = new.reshape(-1, case.old.size)
        diagnostics[case.case_id] = dict(
            solve_count=model.factorizations, wall_s=perf_counter() - t0,
            archived_reconstruction_max_abs=reused_max)
    np.savez_compressed(ROOT / "results/measure_new_reads.npz", **data)
    (ROOT / "results/measure_acquisition.json").write_text(
        json.dumps(diagnostics, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
