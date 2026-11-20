"""Enclose four frozen 24×24 TM witnesses using complex ball linear algebra."""

import argparse
import json
import math
from pathlib import Path
from time import perf_counter, process_time

import flint
from flint import acb, acb_mat, arb, ctx
import numpy as np

from designs import all_designs
from ndt import BACKGROUND, BOUNDS, C0, PLATE_SIDE, grid

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
N = 24
CELLS = N * N
DESIGNS = ("f2_t4_r8", "f4_t2_r8")
FREQUENCIES = (30., 2., 12., 20.)
I = acb(0, 1)


def read_inputs():
    modes = json.loads((RESULTS / "material_weak_modes.json").read_text())
    pairs = modes["cases"][:2]
    if {(p["hard_design"], p["easy_design"]) for p in pairs} != {DESIGNS, DESIGNS[::-1]}:
        raise ValueError("The first two records must be the two requested schedule orientations")
    ids = tuple(pairs[0]["pixel_ids"])
    if not all(tuple(pair["pixel_ids"]) == ids for pair in pairs):
        raise ValueError("The first two pairs must share the same 88-cell region")
    if len(ids) != 88 or len(set(ids)) != len(ids) or min(ids) < 0 or max(ids) >= CELLS:
        raise ValueError("Expected 88 distinct material indices in the 576-cell map")
    lower, upper = BOUNDS
    if not all(lower <= x <= upper for pair in pairs for member in ("minus", "plus") for x in pair[member]):
        raise ValueError("Frozen hypotheses exceed the material box")

    with np.load(RESULTS / "acquisition.npz", allow_pickle=False) as fields:
        sources = fields["sources_m"]
        receivers = fields["receivers_m"]
    sigma = json.loads((RESULTS / "material_witnesses.json").read_text())["noise_sigma"]
    return pairs, ids, sources, receivers, sigma


def exact_array(x):
    return [[arb(float(v)) for v in row] for row in x]


def hankel(x, order=0):
    arg = acb(x)
    return arg.bessel_j(order) + I * arg.bessel_y(order)


def arg_distance(x, y):
    return sum((a - b)**2 for a, b in zip(x, y, strict=True)).sqrt()


def host_blocks(freq, ids, tx_positions, rx_positions):
    start = perf_counter()
    k = 2 * arb.pi() * arb(freq) * 10**9 / arb(C0)
    width = float(PLATE_SIDE / N)
    ka = k * arb(width) / arb.pi().sqrt()
    constant = I * arb.pi() * ka / 2
    external = constant * acb(ka).bessel_j(1)
    diagonal = constant * hankel(ka, 1) - 1

    # Air.dense defines the displacement kernel by the separately stored axis offsets.
    # Each IEEE value is exact here. Norms and all special functions use balls.
    offsets = [arb(float(x)) for x in np.arange(N) * width]
    kernel = {}
    for row in range(N):
        for col in range(row + 1):
            value = diagonal if not row else external * hankel(
                k * (offsets[row]**2 + offsets[col]**2).sqrt()
            )
            kernel[row, col] = kernel[col, row] = value

    rows, cols = zip(*(divmod(i, N) for i in range(CELLS)), strict=True)
    g = [[kernel[abs(ri - rj), abs(ci - cj)] for rj, cj in zip(rows, cols, strict=True)]
         for ri, ci in zip(rows, cols, strict=True)]
    b = arb(BACKGROUND)
    operator = acb_mat([[int(i == j) - b * entry for j, entry in enumerate(row)]
                        for i, row in enumerate(g)])
    centers = exact_array(grid(N))
    n_tx = len(tx_positions)
    n_rx = len(rx_positions)
    tx = exact_array(tx_positions)
    rx = exact_array(rx_positions)
    source_scale = I * (arb(freq) / 16).sqrt() / 4
    rhs = acb_mat([
        [
            *[source_scale * hankel(k * arg_distance(x, position)) for position in tx],
            *[g[i][j] for j in ids],
            *[external * hankel(k * arg_distance(x, position)) for position in rx],
        ]
        for i, x in enumerate(centers)
    ])
    build_time = perf_counter() - start
    start = perf_counter()
    # Host matrix is exactly symmetric. Solving its action on S^T gives (S L^-1)^T.
    solved = operator.solve(rhs, algorithm="precond")
    solve_time = perf_counter() - start
    e = acb_mat([[solved[i, j] for j in range(n_tx)] for i in ids])
    k_host = acb_mat([[solved[i, n_tx + s] for s in range(len(ids))] for i in ids])
    s = acb_mat([[solved[i, n_tx + len(ids) + r] for i in ids] for r in range(n_rx)])
    return e, k_host, s, dict(build_wall_s=build_time, host_solve_wall_s=solve_time,
                              host_order=CELLS, host_rhs=n_tx + len(ids) + n_rx,
                              max_host_component_radius=component_radius(solved))


def component_radius(matrix):
    radius = max(max(z.real.rad().upper(), z.imag.rad().upper()) for z in matrix.entries())
    return math.nextafter(float(radius), math.inf)


def response(e, k_host, s, material):
    q = [arb(v) for v in material]
    operator = acb_mat([
        [int(i == j) - qi * k_host[i, j] for j in range(len(q))]
        for i, qi in enumerate(q)
    ])
    rhs = acb_mat([[qi * e[i, j] for j in range(e.ncols())] for i, qi in enumerate(q)])
    return s * operator.solve(rhs, algorithm="precond")


def outward(value):
    if not value.is_finite():
        raise ArithmeticError("A nonfinite ball cannot certify the result")
    lo = math.nextafter(float(value.lower()), -math.inf)
    hi = math.nextafter(float(value.upper()), math.inf)
    if not arb(lo) <= value.lower() or not arb(hi) >= value.upper():
        raise ArithmeticError("The binary64 endpoints did not enclose the Arb endpoints")
    return {"lower": lo, "upper": hi, "arb": str(value)}


def solve_frequency(f, pairs, ids, sources, receivers, arms):
    wall, cpu = perf_counter(), process_time()
    relevant = [d for d in arms.values() if f in d.freqs]
    tx_ids = sorted({t for d in relevant for t in d.tx})
    rx_ids = sorted({r for d in relevant for r in d.rx})
    e, k_host, s, record = host_blocks(f, ids, sources[tx_ids], receivers[rx_ids])
    t0 = perf_counter()
    differences = [response(e, k_host, s, p["plus"]) - response(e, k_host, s, p["minus"]) for p in pairs]
    selected = []
    for difference in differences:
        norms = {}
        for d in relevant:
            vals = [difference[rx_ids.index(r), tx_ids.index(t)] for r in d.rx for t in d.tx]
            norms[d.name] = sum(abs(v)**2 for v in vals)
        selected.append(norms)
    record.update(ghz=f, tx_ids=tx_ids, rx_ids=rx_ids, local_order=len(ids), local_solves=2 * len(pairs),
                  local_rhs=2 * len(pairs) * len(tx_ids), local_wall_s=perf_counter() - t0,
                  cpu_s=process_time() - cpu, wall_s=perf_counter() - wall,
                  max_signal_difference_component_radius=max(component_radius(d) for d in differences))
    return selected, record


def summarize(pairs, arms, frequencies, calculations, sigma):
    results = []
    used = set(frequencies)
    for pair, bands in zip(pairs, calculations, strict=True):
        rms = sum((arb(a) - arb(b))**2 for a, b in zip(pair["plus"], pair["minus"], strict=True)).sqrt() / N
        results.append(dict(hard_design=pair["hard_design"], easy_design=pair["easy_design"],
                            material_rms=outward(rms), designs={}))
        for name, design in arms.items():
            if not set(design.freqs) <= used:
                continue
            f_terms = [bands[f][name] for f in design.freqs]
            distance = (sum(f_terms) / design.shots).sqrt() / arb(sigma)
            if not distance > 0:
                raise ArithmeticError("Could not show the normalized separation is positive")
            probability = (distance / 2).erfc() / 2
            exact = results[-1]["designs"]
            exact[name] = dict(
                gap_over_sigma=outward(distance),
                equal_prior_error=outward(probability),
                rmse_lower_bound_value=outward(rms * probability),
                samples=design.samples, shots=design.shots, radiated_energy=1.,
                per_frequency_unweighted_squared_distance={
                    str(f): outward(term) for f, term in zip(design.freqs, f_terms, strict=True)
                },
            )
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--precision", type=int, default=112)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--probe", action="store_true", help="Enclose 30 GHz only and print its cost")
    args = parser.parse_args()
    ctx.prec = args.precision
    ctx.threads = args.threads
    wall, cpu = perf_counter(), process_time()
    pairs, ids, sources, receivers, sigma = read_inputs()
    designs = all_designs()
    arms = {name: designs[name] for name in DESIGNS}
    frequencies = FREQUENCIES[:1] if args.probe else FREQUENCIES
    calculation = [{} for _ in pairs]
    records = []
    for f in frequencies:
        norms, record = solve_frequency(f, pairs, ids, sources, receivers, arms)
        for target, new in zip(calculation, norms, strict=True):
            target[f] = new
        records.append(record)
        print(f"{f:g} GHz: build {record['build_wall_s']:.2f}s, "
              f"576-cell host solve {record['host_solve_wall_s']:.2f}s, "
              f"four 88-cell solves {record['local_wall_s']:.2f}s, "
              f"signal radius ≤ {record['max_signal_difference_component_radius']:.2e}", flush=True)
    if args.probe:
        return
    cases = summarize(pairs, arms, frequencies, calculation, sigma)
    payload = dict(
        cases=cases, precision_bits=ctx.prec, threads=ctx.threads, python_flint=flint.__version__,
        flint_version=flint.__FLINT_VERSION__, solve_algorithm="precond",
        frequencies=list(frequencies), host_cells=CELLS, selected_cells=len(ids),
        material_input="First two frozen records of results/material_weak_modes.json",
        noise_unit_shot_sigma=sigma, noise_unit_shot_sigma_binary64_hex=float(sigma).hex(),
        coordinates="Finite scalar TM disk-pulse implementation. Binary64 collocation centers, "
                    "displacement-offset axis, cell width, archived source and receiver positions and materials "
                    "are interpreted as exact. Pi, distances, Hankel kernels, and propagation "
                    "from these inputs are evaluated in Arb.",
        model_timing=records, wall_s=perf_counter() - wall, cpu_s=process_time() - cpu,
    )
    (RESULTS / "witness_ball.json").write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
    for case in cases:
        for design, value in case["designs"].items():
            print(f"{case['hard_design']} vs {case['easy_design']} · {design} · "
                  f"gap {value['gap_over_sigma']['arb']} · "
                  f"risk {value['equal_prior_error']['arb']} · "
                  f"risk bound {value['rmse_lower_bound_value']['arb']}", flush=True)


if __name__ == "__main__":
    main()
