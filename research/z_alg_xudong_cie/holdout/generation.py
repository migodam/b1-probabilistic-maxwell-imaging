"""Create the B3 synthetic fitting and offline assets after the explicit freeze."""

import argparse
from pathlib import Path
import sys
import tempfile

import numpy as np
from flint import arb, ctx
from scipy.ndimage import label
from threadpoolctl import threadpool_limits

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from research.z_alg_xudong_cie.holdout.data import (
    APERTURES, ARCHIVE, DEFAULT_ROOT, INVERSE_SIDE, K, NATIVE_SIDE, NOISE_TAGS, OBJECT_IDS,
    WIDTH, elapsed, manifest_path, observation_id, planned_cases, public_directory,
    public_geometry, require_freeze, timers, write_json, write_npz,
)
from research.z_alg_xudong_cie.holdout.physics import (
    NativeForward, centers, historical, independent_ports,
)

SRAP_SEED = 20261124031
NEW_SEED = 20261124032
AMPLITUDE_SEED = 20261124033
ORDER_SEED = 20261124034
NOISE_SEED = 20261124035
SRAP_FAMILIES = (0, 1, 2, 3, 4, 5, 0, 2)
HIDDEN_TX = 9. + 36. * np.arange(10)
HIDDEN_RX = 4.5 + 18. * np.arange(20)
PADDED_RELATIVE_ALLOWANCE = 1e-6


def padded_radius(raw, measurement):
    squared = sum((arb(float(x.real))**2 + arb(float(x.imag))**2
                   for x in measurement.ravel()), arb(0))
    radius = (arb(raw) + squared.sqrt() / 10**6).upper()
    return float(np.nextafter(float(radius), np.inf))


def target_description():
    """All material-family labels and this permutation stay outside the public manifest."""
    planned = [("srap", fam, i) for i, fam in enumerate(SRAP_FAMILIES)]
    planned += [("new", fam, i) for fam in range(2) for i in range(2)]
    permutation = np.random.default_rng(ORDER_SEED).permutation(len(planned))
    return {oid: planned[i] for oid, i in zip(OBJECT_IDS, permutation, strict=True)}


def new_shape(family, instance, xy):
    rng = np.random.default_rng(np.random.SeedSequence([NEW_SEED, family, instance]))
    angle = rng.uniform(-np.pi, np.pi)
    shift = rng.uniform(-.06, .06, 2)
    varied = rng.uniform(.9, 1.07, 2)
    rot = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    x, y = (((xy - shift) @ rot) / varied).T

    def ellipse(cx, cy, a, b):
        return ((x-cx)/a)**2 + ((y-cy)/b)**2 < 1

    if family == 0:
        dx, dy = rng.uniform(-.025, .025, 2)
        panel = (abs(x/.61)**4 + abs(y/.44)**4) < 1
        windows = (ellipse(-.30 + dx, -.07, .115, .12) |
                   ellipse(.01, .18 + dy, .106, .10) |
                   ellipse(.31 + dx, -.18 + dy, .108, .10))
        weight = .85 + .14 * (y >= .025)
        return weight * (panel & ~windows), "three-window-panel"
    offset = rng.uniform(-.025, .025)
    pair = ellipse(-.32, .10, .31, .32) | ellipse(.31, -.10, .29, .29)
    bridge = (abs(x) < .22) & (abs(y + .31*x - offset) < .092)
    void = (ellipse(-.34, .12 + offset, .128, .15) |
            ellipse(.33, -.12 + offset, .124, .132))
    weight = .80 + .19 * (x > 0)
    return weight * ((pair | bridge) & ~void), "bridged-perforated-pair"


def make_material(oid, item, indices, legacy):
    arm, family, instance = item
    xy = centers(NATIVE_SIDE)
    if arm == "srap":
        previous = legacy.MATERIAL_SEED
        try:
            legacy.MATERIAL_SEED = SRAP_SEED
            shape = legacy.random_material(xy, 20_000 + instance, family)
        finally:
            legacy.MATERIAL_SEED = previous
        family_name = legacy.FAMILY_NAMES[family]
    else:
        shape, family_name = new_shape(family, instance, xy)
        if label(shape.reshape(NATIVE_SIDE, NATIVE_SIDE) > 0)[1] != 1:
            raise ArithmeticError("The frozen connected-family construction did not stay connected.")
    m = np.zeros(NATIVE_SIDE**2)
    active = shape > 0
    values = shape[active]
    if values.size == 0:
        raise ValueError("Frozen material construction is empty.")
    rng = np.random.default_rng(np.random.SeedSequence([AMPLITUDE_SEED, indices[oid]]))
    peak = rng.uniform(3.05, 4.)
    span = float(np.ptp(values))
    m[active] = 2.5 + (peak - 2.5) * (values - values.min()) / span if span else peak
    if not np.all((m[active] >= 2.5) & (m[active] <= 4.0)):
        raise ArithmeticError("Active contrast violates the frozen [2.5, 4] range.")
    return m, family_name


def prepare_geometry(root, legacy):
    start = timers()
    with np.load(ARCHIVE / "fine_geometry.npz", allow_pickle=False) as saved:
        native = {key: saved[key] for key in ("spectrum", *(f"{ap}__{k}" for ap in APERTURES for k in ("S", "B")))}
    coarse = {}
    for ap in APERTURES:
        with np.load(ARCHIVE / f"{ap.lower()}.npz", allow_pickle=False) as saved:
            coarse.update({f"{ap}__{k}": saved[k] for k in ("S", "B")})
            fit_rx_degrees = np.rad2deg(np.arctan2(saved["rx_m"][:, 1], saved["rx_m"][:, 0])) % 360
            for held, overlap in ((HIDDEN_TX, saved["tx_plane_wave_deg"]), (HIDDEN_RX, fit_rx_degrees)):
                separation = abs((held[:, None] - np.asarray(overlap).ravel() + 180) % 360 - 180)
                if separation.min() < 1e-8:
                    raise ArithmeticError("Held angles overlap archived fitting angles.")
    for side, values in ((INVERSE_SIDE, coarse), (NATIVE_SIDE, native)):
        values.update(centers=centers(side), radius=WIDTH / side / np.sqrt(np.pi), k=K,
                      grid_shape=np.array([side, side]))
    coarse_green = legacy.EqualAreaGreen(coarse["centers"], (INVERSE_SIDE, INVERSE_SIDE),
                                         coarse["radius"], K)
    coarse["spectrum"] = coarse_green.spectrum
    write_npz(public_geometry(root, INVERSE_SIDE), **coarse)
    write_npz(public_geometry(root, NATIVE_SIDE), **native)
    s, b = independent_ports(NATIVE_SIDE, HIDDEN_TX, HIDDEN_RX)
    write_npz(root / "offline" / "hidden_geometry96.npz", S=s, B=b,
              tx_angles_deg=HIDDEN_TX, rx_angles_deg=HIDDEN_RX)
    native.update({"held__S": s, "held__B": b})
    named = (*APERTURES, "held")
    sensor, rx = np.unique(np.concatenate([native[f"{key}__S"] for key in named], axis=0),
                           axis=0, return_inverse=True)
    incident, tx = np.unique(np.concatenate([native[f"{key}__B"].T for key in named], axis=0),
                             axis=0, return_inverse=True)
    source_map = dict(zip(named, np.split(tx, np.arange(1, len(named))*10), strict=True))
    receiver_map = dict(zip(named, np.split(rx, np.arange(1, len(named))*20), strict=True))
    return native["spectrum"], sensor, incident.T, source_map, receiver_map, dict(
        aggregate_source_columns=incident.shape[0], aggregate_receiver_rows=sensor.shape[0],
        reused_native_fourier_spectrum=True, new_coarse_green_setups=1, new_hidden_geometry_setups=1,
        **elapsed(start))


def generation(root, freeze_record):
    freeze = require_freeze(freeze_record)
    ctx.prec = 120
    root = Path(root).resolve()
    if root.exists():
        raise FileExistsError(f"The B3 holdout destination exists: {root}")
    root.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f"{root.name}.building.", dir=root.parent))
    began = timers()
    history = historical()
    chosen = target_description()
    indices = {k: i for i, k in enumerate(OBJECT_IDS)}
    spec, sensor, incident, tx, rx, setup = prepare_geometry(staging, history)
    forward = NativeForward(spec, NATIVE_SIDE)
    records = {}
    for oid, item in chosen.items():
        started = timers()
        m, family = make_material(oid, item, indices, history)
        coarse = history.area_average(m, NATIVE_SIDE, INVERSE_SIDE)
        clean, work = forward.forward(m, incident, sensor, rtol=2e-10)
        parts = {ap: clean[np.ix_(rx[ap], tx[ap])] for ap in (*APERTURES, "held")}
        for ia, aperture in enumerate(APERTURES):
            target = parts[aperture]
            rng = np.random.default_rng(np.random.SeedSequence([NOISE_SEED, indices[oid], ia]))
            direction = rng.normal(size=target.shape) + 1j*rng.normal(size=target.shape)
            direction /= np.linalg.norm(direction)
            for noise in NOISE_TAGS:
                name = observation_id(oid, aperture, noise)
                raw = float(noise * np.linalg.norm(target))
                measured = target + raw * direction
                radius = padded_radius(raw, measured)
                write_npz(public_directory(staging) / "observations" / f"{name}.npz",
                          Y=measured, delta=radius, raw_radius=raw)
        write_npz(staging / "offline" / "objects" / f"{oid}.npz",
                  m_native=m, m_coarse=coarse, **{f"clean_{k}": v for k, v in parts.items()})
        records[oid] = dict(lineage=item, family=family, peak=float(np.max(m)),
                            active_native_cells=int(np.count_nonzero(m)), forward_cost=work,
                            **elapsed(started))
        print(oid, "native solve RHS", work["forward_rhs"], "seconds", work["wall_seconds"], flush=True)
    write_json(staging / "offline" / "generation_receipt.json", dict(
        scope="private generation receipt with no algorithm predictions used",
        freeze_record=str(freeze), seeds=dict(srap=SRAP_SEED, original=NEW_SEED, amplitude=AMPLITUDE_SEED,
                                             order=ORDER_SEED, noise=NOISE_SEED),
        geometry_cost=setup, objects=records, **elapsed(began),
    ))
    write_json(manifest_path(staging), dict(
        scope="blind fitting measurements; four observations per independent material",
        independent_materials=len(OBJECT_IDS), cases=list(planned_cases()),
        relative_generation_padding=PADDED_RELATIVE_ALLOWANCE,
    ))
    staging.replace(root)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--freeze-record", type=Path, required=True)
    args = parser.parse_args()
    with threadpool_limits(1):
        generation(args.root, args.freeze_record)
