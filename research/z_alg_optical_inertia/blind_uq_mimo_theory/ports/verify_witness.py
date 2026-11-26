"""Outward verification on exact stored FFT, low-rank, and real material arrays."""
import argparse
import csv
import json
from pathlib import Path
from time import perf_counter, process_time

import numpy as np
from flint import acb, acb_mat, arb, ctx

from model import BANDS, COUNTS, MAIN_STRENGTH
from sweep import tag

HERE = Path(__file__).resolve().parent


def as_acb(z):
    return acb(float(z.real), float(z.imag))


def ball(x):
    return acb_mat(x.shape[0], x.shape[1], [as_acb(z) for z in x.ravel()])


def frobenius(matrix):
    squared = sum((v.real*v.real+v.imag*v.imag for v in matrix.entries()), arb(0))
    return squared.nonnegative_part().sqrt()


def outward(x):
    return [float(np.nextafter(float(x.lower()), -np.inf)), float(np.nextafter(float(x.upper()), np.inf))]


def exact_kernel(spectrum):
    """Separable inverse DFT encloses the exact saved binary64 FFT coefficients."""
    rows = [acb.dft([as_acb(z) for z in row], inverse=True) for row in spectrum]
    cols = [acb.dft([row[j] for row in rows], inverse=True) for j in range(spectrum.shape[1])]
    return [cols[j][i] for i in range(spectrum.shape[0]) for j in range(spectrum.shape[1])]


def proof_forward(material, data, cache, band):
    start = perf_counter(), process_time()
    ids = np.flatnonzero(material)
    n = len(ids)
    m = [arb(float(v)) for v in material[ids]]
    if np.any(material < 0) or np.any(material > 4.5):
        raise ValueError("Material outside [0,4.5]")
    if "G" in data:
        full = ball(data["G"][np.ix_(ids, ids)])
        basis = "exact stored dense G includes network"
    else:
        if band not in cache:
            cache[band] = exact_kernel(data["spectrum"])
        height, width = data["pad"]
        y, x = np.unravel_index(ids, tuple(data["grid_shape"]))
        lags = (y[:, None]-y[None, :]) % height * width + (x[:, None]-x[None, :]) % width
        bare = acb_mat(n, n, [cache[band][int(k)] for k in lags.ravel()])
        full = bare + ball(data["low_rank_left"][ids]) * ball(data["low_rank_right"][:, ids])
        basis = "exact stored bare FFT spectrum plus exact low-rank network arrays"
    sys = acb_mat(n, n, [int(i == j)-m[i]*full[i, j] for i in range(n) for j in range(n)])
    E0 = ball(data["E0"][ids])
    rhs = acb_mat(n, E0.ncols(), [m[i]*E0[i, j] for i in range(n) for j in range(E0.ncols())])
    J = sys.solve(rhs)
    Y = ball(data["S"][:, ids]) * J
    return Y, dict(active_state_dimension=n, rhs=E0.ncols(), seconds_wall=perf_counter()-start[0],
                   seconds_CPU=process_time()-start[1], basis=basis)


def equal_power_dual_band(rows):
    joint = {}
    for q in COUNTS:
        entries = [rows.get(f"{band}_{q}") for band in BANDS]
        if any(entry is None for entry in entries):
            continue
        endpoints = list(zip(*(entry["separation_over_2r"] for entry in entries), strict=True))
        # Each signal distance is reduced by 1/sqrt(2) at half the source
        # power; the unchanged per-band hard radii combine as sqrt(2)*r.
        low, high = [(sum((arb(x)**2 for x in end), arb(0))).sqrt()/2 for end in endpoints]
        joint[q] = dict(separation_over_joint_two_r=[outward(low)[0],outward(high)[1]],
                        shared_joint_ball_impossible=bool(low>1))
    return joint


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bands", nargs="+", type=int, default=BANDS)
    p.add_argument("--counts", nargs="+", type=int, default=COUNTS)
    args = p.parse_args()
    ctx.prec, ctx.threads = 128, 1
    numbers = {}
    kernel = {}
    w = np.load(HERE/"results"/"mechanism_witness.npz")
    materials = [w["with_region"].ravel(), w["without_region"].ravel()]
    roi = w["roi"].ravel()
    sigma = arb(float(w["reference_sigma"]))
    for band in args.bands:
        for q in args.counts:
            op_id = "n32_f400_a08" if band == 400 and q == 8 else f"n32_f{band}_a{q:02}_c048"
            with np.load(HERE/"public"/"operators"/f"{op_id}.npz") as archive:
                data = dict(archive)
            audited = [proof_forward(m, data, kernel, band) for m in materials]
            radius = arb(3)*sigma*q/2
            identifier = tag(band, q, MAIN_STRENGTH, 0)
            y = ball(w[identifier+"_midpoint_Y"])
            sep = frobenius(audited[0][0]-audited[1][0])
            mid = [frobenius(forward - y)/radius for forward, _ in audited]
            numbers[f"{band}_{q}"] = dict(
                measurements_in_midpoint_ball=all(v <= 1 for v in mid),
                separation_over_2r=outward(sep/(2*radius)),
                distance_each_to_numeric_midpoint_over_r=[outward(v) for v in mid],
                shared_ball_impossible=bool(sep > 2*radius), solves=[audit for _, audit in audited])
            print(band, q, numbers[f"{band}_{q}"], flush=True)
    f = HERE/"results"/"mechanism_witness_arb.json"
    if f.exists():
        all_rows = json.loads(f.read_text())["conditions"]
        all_rows.update(numbers)
    else:
        all_rows = numbers
    f.write_text(json.dumps(dict(arithmetic_bits=ctx.prec, roi_cells=int(roi.sum()),
                roi_mean_difference=float(np.mean(materials[0][roi]-materials[1][roi])),
                conditions=all_rows, equal_power_dual_band=equal_power_dual_band(all_rows)), indent=2))


if __name__ == "__main__":
    main()
