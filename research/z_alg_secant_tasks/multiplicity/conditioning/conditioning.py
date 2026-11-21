"""Frozen symmetry-breaking conditioning control for the separated-root mechanism."""

from itertools import product
from pathlib import Path
import json
import sys
from time import perf_counter, process_time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from separated import FIXED, FullModel, l2, newton, root_gaps, write_json
from certify_separated import ExactModel, arb, certify_point, cross_proofs, ctx, complex_range

import numpy as np
from threadpoolctl import threadpool_limits


DISTANCES = (3, 10, 30, 100)
MATERIAL = np.array([2., 8.])
SWAPPED = MATERIAL[::-1].copy()


def component(rx_x):
    return {
        "geometry_decimal": {"centers": [["-0.18", "0"], ["0.18", "0"]],
                             "sources": [["0", "-0.85"]],
                             "receivers": [[rx_x, "0.95"]]},
        "disk_radius_decimal": "0.04",
        "known_material_contrast_decimal": ["0", "0"],
        "active_indices": [0, 1], "legal_material_range": [0.5, 30.0],
        "local_observation_binary64": [[0., 0.]],
    }


def inspect(model, active):
    fields, f, J, A, rhs, current = model.solve(active, differentiate=True)
    sv = np.linalg.svd(J, compute_uv=False)
    return fields[model.local], {
        "material": active.tolist(),
        "physical_scattered_port": [[float(x.real), float(x.imag)] for x in fields[model.local]],
        "distance_to_target_if_available": l2(f),
        "material_jacobian_singular_values": sv.tolist(),
        "material_jacobian_condition": float(sv[0] / sv[-1]),
        "state_sigma_min": float(np.linalg.svd(A, compute_uv=False)[-1]),
        "state_relative_residual": l2(A @ current - rhs) / l2(rhs),
    }


def make_exact(input_data, count, spacing, target=None):
    model = ExactModel(input_data, count, spacing)
    if target is None:
        y, _ = model.state([arb(2), arb(8)])
        target = y[0, 0]
    model.target = [target] * count
    return model, target


def run():
    begin, cpu = perf_counter(), process_time()
    old = json.loads(FIXED.read_text())
    old_model = FullModel(old, 1, 0)
    legacy = [inspect(old_model, np.array(old[f"component_seed_{branch}"]))[1] for branch in "AB"]
    sym_input = component("0")
    sym_f = FullModel(sym_input, 1, 0)
    y, sym_a = inspect(sym_f, MATERIAL)
    sym_f.target = y
    _, sym_a = inspect(sym_f, MATERIAL)
    _, sym_b = inspect(sym_f, SWAPPED)
    sym_proof, sym_exact = make_exact(sym_input, 1, 0)
    sym_cert_a, _, _ = certify_point(sym_proof, "A", MATERIAL, sym_f.lo, sym_f.hi)
    sym_cert_b, _, _ = certify_point(sym_proof, "B", SWAPPED, sym_f.lo, sym_f.hi)
    print("mirror symmetry, sigma_min A/B:",
          sym_a["material_jacobian_singular_values"][-1], sym_b["material_jacobian_singular_values"][-1],
          "historical A/B:", *[x["material_jacobian_singular_values"][-1] for x in legacy], flush=True)

    displaced = component("0.013")
    source = FullModel(displaced, 1, 0)
    y, _ = inspect(source, MATERIAL)
    source.target = y
    _, disp_a = inspect(source, MATERIAL)
    _, disp_swap = inspect(source, SWAPPED)
    partner, _, _, branch_b = newton(source, SWAPPED)
    _, disp_b = inspect(source, partner)
    displaced["local_observation_binary64"] = [[float(y[0].real), float(y[0].imag)]]
    target_model, exact_y = make_exact(displaced, 1, 0)
    certified_displaced, single_proof = False, []
    if branch_b["status"] == "numerical_root" and l2(partner - MATERIAL) > 1e-6:
        for key, candidate in (("A", MATERIAL), ("B", partner)):
            certificate, _, _ = certify_point(target_model, key, candidate, source.lo, source.hi)
            single_proof.append(certificate)
        certified_displaced = all(x["status"] == "certified_simple_real_root" for x in single_proof)
    print("displaced receiver", branch_b["status"], partner.tolist(),
          "sigma min", disp_a["material_jacobian_singular_values"][-1],
          disp_b["material_jacobian_singular_values"][-1], "both certified", certified_displaced, flush=True)

    result = {
        "scope": "fixed two-disk conditioning control; all distances in reference wavelengths",
        "displacement_decimal": "0.013", "target_definition": "exact full scattered F((2,8)) at displaced receiver in one component",
        "legacy_near_fold": legacy, "mirror_geometry": sym_input,
        "symmetric_pair": {"A": sym_a, "B": sym_b, "exact_target": complex_range(sym_exact),
                           "reflection_reason": "reflection about x=0 swaps disks and leaves both transducers fixed",
                           "root_certificates": [sym_cert_a, sym_cert_b],
                           "float_work": dict(sym_f.cost), "interval_solve_count": sym_proof.state_solve_count},
        "displaced_geometry": displaced, "displaced_pair": {
            "A": disp_a, "unadjusted_exchange": disp_swap, "B": disp_b, "exchange_attempt": branch_b,
            "exact_target": complex_range(exact_y), "root_certificates": single_proof,
            "both_regular_roots_certified": certified_displaced, "float_work": dict(source.cost),
            "interval_solve_count": target_model.state_solve_count,
        },
        "fixed_spacings": DISTANCES,
    }
    cases, proofs, arrays = [], [], {}
    if certified_displaced:
        seeds = {"A": MATERIAL, "B": partner}
        for count in (2, 3):
            proved_at_first = False
            for R in DISTANCES:
                model = FullModel(displaced, count, R)
                prefix = f"L{count}_R{R}"
                arrays[f"{prefix}__target_local"] = model.target
                arrays[f"{prefix}__local_port_selection"] = model.local
                for kind, points in model.xy.items():
                    arrays[f"{prefix}__{kind}"] = points
                accepted, attempts = [], []
                for labels in product("AB", repeat=count):
                    label = "".join(labels)
                    start = np.concatenate([seeds[x] for x in labels])
                    z, ports, jac, meta = newton(model, start)
                    meta["branch_seed"] = label
                    attempts.append(meta)
                    if meta["status"] == "numerical_root" and not any(
                            l2(z - root) <= 1e-6 for _, root, _ in accepted):
                        accepted.append((label, z, ports))
                    arrays[f"{prefix}__{label}__material"] = z
                    arrays[f"{prefix}__{label}__all_ports"] = ports
                    arrays[f"{prefix}__{label}__local_jacobian"] = jac
                data = {"id": prefix, **model.geometry, "distinct_numerical_root_count": len(accepted),
                        "attempts": attempts, "cost": dict(model.cost),
                        "physical_gaps": root_gaps(accepted, model.local)}
                cases.append(data)
                print(prefix, len(accepted), "/", 2**count, flush=True)
                if len(accepted) == 2**count and not proved_at_first:
                    t, c = perf_counter(), process_time()
                    exact, _ = make_exact(displaced, count, R, target=exact_y)
                    certificates, roots = [], []
                    for label, guess, _ in accepted:
                        certificate, box, fields = certify_point(exact, label, guess, model.lo, model.hi)
                        certificates.append(certificate)
                        if box is not None:
                            roots.append((label, box, fields))
                    pair_checks, solo = cross_proofs(exact, roots)
                    proofs.append({
                        "id": prefix, "number_certified": len(roots),
                        "certificates": certificates, "cross_port_gaps": pair_checks,
                        "individually_discriminating_cross_ports": solo,
                        "cross_ports_receiver_and_source_indices": exact.cross,
                        "all_root_boxes_disjoint": bool(pair_checks and all(
                            r["material_l2"][0] > 0 for r in pair_checks)),
                        "interval_state_count": exact.state_solve_count,
                        "wall_seconds": perf_counter()-t, "cpu_seconds": process_time()-c,
                    })
                    proved_at_first = True
    result.update({"coupled_cases": cases, "coupled_certificates_at_first_complete_spacing": proofs,
                   "precision_bits": ctx.prec, "wall_seconds": perf_counter()-begin,
                   "cpu_seconds": process_time()-cpu, "legacy_float_work": dict(old_model.cost)})
    np.savez_compressed(HERE/"materials_and_ports.npz", **arrays)
    write_json(HERE/"results.json", result)


if __name__ == "__main__":
    with threadpool_limits(limits=1, user_api="blas"):
        run()
