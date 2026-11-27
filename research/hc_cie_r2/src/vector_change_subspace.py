"""Two-sided rank-six localization, followed by independent full SRAP audits."""
import argparse
from itertools import combinations
import json
import os
from pathlib import Path
import time

import vector_change as vector
import numpy as np
from scipy.linalg import svd


ROOT = Path(__file__).resolve().parents[1]
CONFIG = json.loads((ROOT / "configs/vector_change_subspace.json").read_text())


def _basis_errors(o, u, observed, n):
    left, singular, vh = svd(observed, full_matrices=False, check_finite=False)
    limited = CONFIG["rank"]
    left_span, right_span = left[:, :limited], vh[:limited].conj().T
    qo = np.linalg.qr(o.T.reshape(n, 3, -1).transpose(0, 2, 1), mode="reduced")[0]
    qu = np.linalg.qr(u.reshape(n, 3, -1).conj().transpose(0, 2, 1), mode="reduced")[0]
    left_projected = np.einsum("ar,nrc->nac", left_span.conj().T, qo, optimize=True)
    right_projected = np.einsum("ar,nrc->nac", right_span.conj().T, qu, optimize=True)
    left_residual = qo - np.einsum("ra,nac->nrc", left_span, left_projected, optimize=True)
    right_residual = qu - np.einsum("ra,nac->nrc", right_span, right_projected, optimize=True)
    angles_left = np.sum(abs(left_residual)**2, axis=(1, 2))
    angles_right = np.sum(abs(right_residual)**2, axis=(1, 2))
    return angles_left, angles_right, singular


def search(physics, model, background, reference, measured, o, u, dictionary):
    lower, upper = -background, 4.0 - background
    left_angles, right_angles, singular = _basis_errors(o, u, measured, model.N)
    combined = left_angles + right_angles
    selected = np.argsort(combined, kind="stable")[:CONFIG["max_voxels"]]
    pairs = sorted((tuple(sorted(p)) for p in combinations(map(int, selected), 2)),
                   key=lambda p: (combined[p[0]] + combined[p[1]], p))
    shortlist = []
    cache = {}
    for p in pairs[:CONFIG["exact_local_finalists"]]:
        ids = list(p)
        local = vector.local_columns(ids)
        _, k = vector.get_cache(physics, model, background, cache, ids)
        o_pair, u_pair = o[:, local], u[local]
        born = np.linalg.lstsq(dictionary[ids].T, vector._real_vector(measured), rcond=None)[0]
        born = np.clip(born, lower[ids], upper[ids])
        closed, method_info = vector.closed_readout(measured, o_pair, u_pair, k, lower[ids], upper[ids])
        fitted, residual, cost = vector.bounded_refit(o_pair, u_pair, k, measured, born, closed, lower[ids], upper[ids])
        shortlist.append(dict(sites=ids, change=fitted.tolist(), residual=residual, readout=method_info, **cost))
    audited = []
    for pair in sorted(shortlist, key=lambda r: r["residual"])[:CONFIG["independent_full_audits"]]:
        estimate, loss, state_error = vector.independently_audit(physics, model, background, pair["sites"], pair["change"],
                                                                 reference, measured)
        pair.update(independent_loss=loss, independent_state_error=state_error)
        audited.append((loss, pair, estimate))
    _, winner, predicted = min(audited, key=lambda entry: entry[0])
    return winner, predicted, shortlist, left_angles, right_angles, singular, selected


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, default=Path.cwd())
    parser.add_argument("--output", type=Path, default=ROOT / "results/vector_change")
    args = parser.parse_args()
    source = args.workspace.resolve() / "research/srap_cie_array_study_20261121/vector/physics.py"
    output = args.output.resolve()
    os.environ.setdefault("SRAP_DEVICE", "cpu")
    physics = vector._load_module("hc_cie_two_sided_vector_physics", source)
    physics.CD = physics.torch.complex128
    clock, cpu = time.perf_counter(), time.process_time()
    models = {label: physics.Model(12, label) for label in vector.DESIGN["arrays"]}
    reports, setup = [], []
    arrays = {}
    for parent in vector.DESIGN["parents"]:
        for label, model in models.items():
            case = f"{label}_{parent:02d}"
            path = output / f"{case}.npz"
            started, c0 = time.perf_counter(), time.process_time()
            counters = model.count.copy()
            with np.load(path, allow_pickle=False) as saved:
                background = saved["background"]
                reference, current, error_ref = vector.forward(physics, model, background)
                o, u, dictionary, adj_error = vector.build_dictionary(physics, model, background, current)
                setup.append(dict(case=case, wall_s=time.perf_counter() - started, cpu_s=time.process_time() - c0,
                                  counters=dict(model.count - counters), state_error=error_ref, adjoint_error=adj_error))
                paid = []
                for mode in CONFIG["conditions"]:
                    started, c0 = time.perf_counter(), time.process_time()
                    counters = model.count.copy()
                    observation = saved[f"{mode}_observation"]
                    winner, predicted, shortlist, left_angles, right_angles, singular, selected = search(
                        physics, model, background, reference, observation, o, u, dictionary)
                    for ket, data in dict(left_sines_squared=left_angles, right_sines_squared=right_angles,
                                          differential_singular_values=singular, proposed=winner["sites"],
                                          proposed_change=winner["change"], prediction=predicted).items():
                        arrays[f"{case}_{mode}_{ket}"] = np.asarray(data)
                    paid.append(dict(case=case, array=label, condition=mode, winner=winner, shortlist=shortlist,
                                     voxel_screen=selected.tolist(), predicted=predicted, measured=observation,
                                     wall_s=time.perf_counter() - started, cpu_s=time.process_time() - c0,
                                     counters=dict(model.count - counters)))
                # Score all candidates after all predictions for this observation have been chosen.
                clean, truth, changes = saved["clean_differential"], saved["true_sites"], saved["true_change"]
                for info in paid:
                    estimate, measured = info.pop("predicted"), info.pop("measured")
                    radius = 0.0 if info["condition"] == "clean" else 0.1 * np.linalg.norm(clean)
                    winner = info["winner"]
                    info.update(vector.score_change(winner["sites"], winner["change"], truth, changes,
                                                    clean, estimate, measured, radius),
                                true_sites_in_screen=int(len(set(truth) & set(info["voxel_screen"]))))
                    reports.append(info)
                    print(case, info["condition"], info["support_overlap"],
                          round(info["material_rel_l2"], 5), flush=True)
    np.savez_compressed(output / "two_sided.npz", **arrays)
    (output / "two_sided.json").write_text(json.dumps(dict(
        rows=reports, setups=setup, config=CONFIG, wall_s=time.perf_counter() - clock,
        cpu_s=time.process_time() - cpu), indent=2) + "\n")


if __name__ == "__main__":
    main()
