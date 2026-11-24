"""Generate the prospectively frozen, four-material single-tile control."""

import argparse
import json
from pathlib import Path
import sys
import tempfile

import numpy as np
from flint import ctx
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
if not __package__:
    sys.path.insert(0, str(HERE.parents[3]))

from research.z_alg_xudong_cie.holdout.data import (
    APERTURES, ARCHIVE, FREQUENCY, K, NATIVE_SIDE, OPTICAL, WIDTH, WORKSPACE,
    elapsed, timers, write_json, write_npz,
)
from research.z_alg_xudong_cie.holdout.generation import padded_radius
from research.z_alg_xudong_cie.holdout.physics import NativeForward, centers
from research.z_alg_xudong_cie.holdout.verify_truth import (
    ExactFourierGeometry, StableActiveModel, verified_row,
)

DEFAULT_ROOT = WORKSPACE / ".artifacts/z_alg_xudong_cie/site_control"
PROTOCOL = HERE / "FROZEN_PROTOCOL.json"
OPERATOR_NAME = "geometry96.npz"
TILE_SIDE = 24
NUM_TILES = 16
OIDS = tuple(f"b3site_{i:04d}" for i in range(1, 5))


def plan(spec):
    for oid in OIDS:
        for aperture in spec["apertures"]:
            for noise in spec["relative_noise"]:
                case = f"{oid}__{aperture.lower()}__n{int(round(100*noise)):02d}"
                yield dict(case=case, object_id=oid, array=aperture, noise_relative=noise,
                           side=NATIVE_SIDE, tx_count=10, rx_count=20, frequency=FREQUENCY,
                           measurement_file=f"observations/{case}.npz")


def tile_material(site, family, seed, peak_seed, index):
    rng = np.random.default_rng(np.random.SeedSequence([seed, index]))
    height, width = rng.integers(15, 20, size=2)
    offset_y = int(rng.integers(1, TILE_SIDE - height))
    offset_x = int(rng.integers(1, TILE_SIDE - width))
    local_y = (np.arange(height) + .5) / height * 2 - 1
    local_x = (np.arange(width) + .5) / width * 2 - 1
    y, x = np.meshgrid(local_y, local_x, indexing="ij")
    theta = rng.uniform(-.18, .18)
    x, y = (x * np.cos(theta) - y * np.sin(theta), x * np.sin(theta) + y * np.cos(theta))
    x -= rng.uniform(-.08, .08)
    y -= rng.uniform(-.08, .08)
    if family == 0:
        support = x**2 + (y / .91)**2 < 1
        level = .2 + .8 * (1 - np.minimum(.95, x**2 + y**2))
    elif family == 1:
        support = (((x - .14) / 1.04)**2 + ((y + .09) / .92)**2 < 1) & (
                   ((x + .18) / 1.04)**2 + ((y - .03) / .95)**2 < 1)
        level = np.where(y - .3*x > .04, 1., .64)
    elif family == 2:
        support = (x**2 + (y/.95)**2 < 1) & ((x + .12)**2 + (y - .15)**2 > .46**2)
        support &= ~((x > .12) & (y > .24 + .17*x))
        level = np.where(x + y < .04, 1., .71)
    else:
        upright = (x > -.96) & (x < -.29) & (y > -.92) & (y < .94)
        base = (y > -.89) & (y < -.26) & (x > -.89) & (x < .94)
        bend = (x > -.47) & (x < .16) & (y > -.36) & (y < -.03)
        support = upright | base | bend
        level = np.where(x > -.14, 1., .73)
    row, col = divmod(int(site), 4)
    y0, x0 = row*TILE_SIDE + offset_y, col*TILE_SIDE + offset_x
    peak_rng = np.random.default_rng(np.random.SeedSequence([peak_seed, index]))
    peak = float(peak_rng.uniform(3.25, 4.))
    material = np.zeros((NATIVE_SIDE, NATIVE_SIDE), dtype=float)
    value = level[support]
    rescaled = 2.5 + (peak - 2.5) * (value - value.min()) / np.ptp(value)
    local = np.zeros((height, width))
    local[support] = rescaled
    material[y0:y0+height, x0:x0+width] = local
    yy, xx = np.nonzero(material)
    extent = tuple(int(np.ptp(k) + 1) for k in (yy, xx))
    tiles = np.unique(4*(yy//TILE_SIDE) + xx//TILE_SIDE)
    if tiles.tolist() != [site] or not all(11 <= k <= 19 for k in extent):
        raise ArithmeticError("Frozen control did not obey its single-tile extent contract.")
    if not np.all((rescaled >= 2.5) & (rescaled <= 4.0)):
        raise ArithmeticError("Frozen contrast floor or cap was violated.")
    return material.ravel(), extent, peak


def load_geometry():
    with np.load(ARCHIVE / "fine_geometry.npz", allow_pickle=False) as source:
        arrays = {key: source[key] for key in ("spectrum", *(f"{a}__{q}" for a in APERTURES for q in ("S","B")))}
    arrays.update(centers=centers(NATIVE_SIDE), k=K, radius=WIDTH/NATIVE_SIDE/np.sqrt(np.pi))
    return arrays


def build(root):
    started = timers()
    ctx.prec, ctx.threads = 128, 1
    protocol = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    if root.exists():
        raise FileExistsError(f"Site-control destination already exists: {root}")
    root.parent.mkdir(parents=True, exist_ok=True)
    staged = Path(tempfile.mkdtemp(prefix=root.name+".building.", dir=root.parent))
    public = staged / "public"
    private = staged / "offline"
    rows = list(plan(protocol))
    geom_work = timers()
    geom = load_geometry()
    write_npz(public / OPERATOR_NAME, **geom)
    exact = ExactFourierGeometry(public / OPERATOR_NAME)
    native = NativeForward(geom["spectrum"], NATIVE_SIDE)
    geometry_cost = elapsed(geom_work)
    sites = np.random.default_rng(protocol["site_seed"]).choice(NUM_TILES, len(OIDS), replace=False)
    families = np.random.default_rng(protocol["ordering_seed"]).permutation(len(OIDS))
    verdicts, key = [], {}
    for index, (oid, site, family) in enumerate(zip(OIDS, sites, families, strict=True)):
        t0 = timers()
        m, extent, peak = tile_material(site, family, protocol["shape_seed"],
                                        protocol["material_peak_seed"], index)
        y, numeric = native.stable_dense(m, exact.B, exact.S, rtol=2e-10)
        numeric.update(elapsed(t0))
        rp = timers()
        independent = StableActiveModel(exact, m)
        bound_inverse, inv_receipt = independent.inverse_bound()
        truth, interval_correction, interval_work = independent.enclose_field(bound_inverse, 128)
        proofs = []
        for ia, aperture in enumerate(APERTURES):
            rx, tx = exact.inds[aperture]
            clean = y[np.ix_(rx, tx)]
            rng = np.random.default_rng(np.random.SeedSequence([protocol["noise_seed"], index, ia]))
            direction = rng.normal(size=clean.shape) + 1j*rng.normal(size=clean.shape)
            direction /= np.linalg.norm(direction)
            for row in (row for row in rows if row["object_id"] == oid and row["array"] == aperture):
                raw = float(row["noise_relative"] * np.linalg.norm(clean))
                observed = clean + raw * direction
                delta = padded_radius(raw, observed)
                measured = dict(Y=observed, delta=delta, raw_radius=raw)
                write_npz(public / row["measurement_file"], **measured)
                proof = verified_row(row["case"], measured, (rx, tx), truth, interval_correction)
                proofs.append(proof)
        verdicts.extend(proofs)
        write_npz(private / "objects" / f"{oid}.npz", m_native=m)
        key[oid] = dict(
            site_zero_based_row_major=int(site), family=protocol["shape_families"][family],
            realized_support_extent_y_x=extent, peak=peak, native_simulation=numeric,
            independent_interval_work=dict(inverse=inv_receipt, primal_and_adjoint=interval_work,
                                           physical_LU=1, inverse_rhs=independent.n, **elapsed(rp)),
        )
    write_json(private / "truth_key.json", dict(objects=key, fit_proofs=verdicts))
    write_json(public / "manifest.json", dict(
        scope=protocol["scope"], native_operator=OPERATOR_NAME, cases=rows,
        independent_materials=len(OIDS), total_observations=len(rows),
    ))
    safe = dict(
        independent_materials=len(OIDS), public_observations=len(rows),
        truth_inside_padded_ball=sum(x["true_material_in_padded_ball"] for x in verdicts),
        minimum_interval_margin=min(x["verified_margin_lower"] for x in verdicts),
        generation_LU=len(OIDS), generation_forward_rhs=exact.B.shape[1] * len(OIDS),
        verification_LU=len(OIDS), verification_inverse_rhs=sum(
            value["independent_interval_work"]["inverse_rhs"] for value in key.values()),
        verification_primal_rhs=exact.B.shape[1] * len(OIDS),
        verification_adjoint_rhs=exact.S.shape[0] * len(OIDS),
        geometry_setup=geometry_cost, **elapsed(started),
    )
    write_json(private / "audit_aggregate.json", safe)
    staged.replace(root)
    print(f"TRAIN READY {root/'public'}; verified native reference fit for "
          f"{safe['truth_inside_padded_ball']}/{len(rows)} rows.", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    args = parser.parse_args()
    with threadpool_limits(1):
        build(args.root)
