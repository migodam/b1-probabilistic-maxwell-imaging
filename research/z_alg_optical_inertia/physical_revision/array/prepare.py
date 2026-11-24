"""Build public array inputs and audit all balls with offline generating materials."""

import json
from fractions import Fraction
from pathlib import Path
import platform
import sys
from time import perf_counter, process_time

import numpy as np
from flint import acb_mat, arb, ctx
from scipy.special import hankel1, jv
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
STUDY = HERE.parents[1]
sys.path.insert(0, str(STUDY / "synthetic"))
sys.path.insert(0, str(STUDY / "submission/curve"))
from nonvacuity import ExactFourierGeometry, StableActiveModel, norm, hi, lo
from core.certificate import ball
from search import save_json, times

SOURCES = (10, 16, 24, 32)
RECEIVERS = (20, 32)
NOISE_PCT = (1, 2)
CAPS = (.75, 1., 1.25, 1.5, 2.)
FREQUENCY = 400e6
LIGHT_SPEED = 3e8
WIDTH = 2.
RECEIVER_RADIUS = 3.75
NOISE_SEED = 20261124031
INHERITED_NOISE_SEED = 20261122092
GUARD = 1e-6
FINE = 96
PUBLIC = STUDY / "synthetic/public"
EVAL = STUDY / "synthetic/evaluation"
INHERITED = STUDY / "submission/curve/observations"


def load(path):
    return json.loads(path.read_text())


def directions(n):
    return [Fraction(i, n) for i in range(n)]


def xy(fractions):
    angle = np.deg2rad([float(f) * 360. for f in fractions])
    return np.column_stack((np.cos(angle), np.sin(angle)))


def candidate_arrays():
    freq = 2 * np.pi * FREQUENCY / LIGHT_SPEED
    with np.load(PUBLIC / "fine_geometry.npz", allow_pickle=False) as f:
        old_s, old_b = f["PW360_10x20__S"], f["PW360_10x20__B"]
    delta = WIDTH/FINE
    one = (np.arange(FINE) + .5) * delta - WIDTH/2
    x, y = np.meshgrid(one, one)
    centers = np.column_stack((x.ravel(), y.ravel()))
    with np.load(PUBLIC / "coarse_green.npz", allow_pickle=False) as coarse:
        coarse_centers = coarse["centers"]
    if not np.allclose(coarse_centers[:32, 1], coarse_centers[0, 1]):
        raise ValueError("Saved grid indexing does not increase x within a row.")
    generated_old = np.exp(1j*freq*(centers@xy(directions(10)).T))
    if not np.allclose(generated_old, old_b, atol=2.e-14, rtol=0):
        raise ValueError("New grid coordinates do not reproduce the original plane waves.")
    cell_radius = delta / np.sqrt(np.pi)
    sensor, incident = {}, {}
    old_rows = dict(zip(directions(20), old_s, strict=True))
    old_cols = dict(zip(directions(10), old_b.T, strict=True))
    for p in RECEIVERS:
        pos = RECEIVER_RADIUS * xy(directions(p))
        d = np.linalg.norm(pos[:, None] - centers[None, :], axis=2)
        rows = .5j*np.pi*freq*cell_radius*jv(1, freq*cell_radius)*hankel1(0, freq*d)
        if p == 20 and not np.allclose(rows, old_s, atol=2.e-14, rtol=0):
            raise ValueError("New formula does not reproduce the historical ring operator.")
        sensor[p] = np.vstack([old_rows.get(a, row) for a, row in zip(directions(p), rows)])
        if p != 20:
            sensor[p] *= np.sqrt(20/p)
    for s in SOURCES:
        columns = np.exp(1j*freq*(centers@xy(directions(s)).T)).T
        incident[s] = np.vstack([old_cols.get(a, row) for a, row in zip(directions(s), columns)]).T
        if s != 10:
            incident[s] *= np.sqrt(10/s)
    return sensor, incident


def shared_draw(index):
    nr = sorted({x for p in RECEIVERS for x in directions(p)})
    ns = sorted({x for s in SOURCES for x in directions(s)})
    rng = np.random.default_rng(np.random.SeedSequence([NOISE_SEED, index]))
    draws = rng.standard_normal((len(nr), len(ns))) + 1j*rng.standard_normal((len(nr), len(ns)))
    original = np.random.default_rng(np.random.SeedSequence([INHERITED_NOISE_SEED, index, 0, 1]))
    old_draw = original.standard_normal((20, 10)) + 1j*original.standard_normal((20, 10))
    rx = {v:i for i,v in enumerate(nr)}
    tx = {v:i for i,v in enumerate(ns)}
    draws[np.ix_([rx[x] for x in directions(20)], [tx[x] for x in directions(10)])] = old_draw
    return {(p,s): draws[np.ix_([rx[x] for x in directions(p)], [tx[x] for x in directions(s)])]
            for p in RECEIVERS for s in SOURCES}


def index_chunks(bank, axis):
    start = 0
    result = {}
    for label, value in bank.items():
        end = start + value.shape[axis]
        result[label] = np.arange(start, end)
        start = end
    return result


def main():
    ctx.prec, ctx.threads = 128, 1
    start, cpu = perf_counter(), process_time()
    public = HERE / "observations"
    public.mkdir(exist_ok=True)
    t, c = perf_counter(), process_time()
    sensor, incident = candidate_arrays()
    geo = ExactFourierGeometry(PUBLIC / "fine_geometry.npz")
    geo.S = np.vstack(list(sensor.values()))
    geo.B = np.column_stack(list(incident.values()))
    rx = index_chunks(sensor, 0)
    tx = index_chunks(incident, 1)
    np.savez_compressed(public / "arrays.npz", **{f"S{p}": a for p,a in sensor.items()},
                        **{f"B{s}": a for s,a in incident.items()})
    setup = times(t, c)
    old = {x["id"]: x for x in load(INHERITED / "manifest.json")["observations"]}
    old_audit = {x["id"]: x for x in load(INHERITED.parent / "observation_audit.json")["observations"]}
    rows, audit = [], []
    cost = {}
    for ii in range(6):
        object_id = f"optfresh_{ii+1:02d}"
        t, c = perf_counter(), process_time()
        with np.load(EVAL / f"{object_id}.npz", allow_pickle=False) as f:
            material = f["m_fine"]
            family = str(f["family"])
        model = StableActiveModel(geo, material)
        inv, inverse = model.inverse_bound()
        exact_y, forward_error, forward = model.enclose_field(inv, 128)
        direction = shared_draw(ii)
        peak = float(material.max())
        for p in RECEIVERS:
            for s in SOURCES:
                if (s,p) != (10,20):
                    yc_ball = acb_mat(p,s,[exact_y[int(i),int(j)] for i in rx[p] for j in tx[s]])
                    yc = np.array([complex(z) for z in yc_ball.entries()]).reshape(p,s)
                    cleaner_error = hi(norm(ball(yc)-yc_ball)+forward_error)
                    random = direction[p,s]
                    noise_2 = random * (.02*np.linalg.norm(yc) / np.linalg.norm(random))
                for pct in NOISE_PCT:
                    ident = f"{object_id}__pw360_{s}x{p}__pair{pct}pct"
                    if (s,p) == (10,20):
                        meta = old[ident]
                        check = old_audit[ident]
                        with np.load(INHERITED / meta["file"], allow_pickle=False) as f:
                            y = f["Y"]
                        with np.load(EVAL / f"{object_id}__pw360_10x20__n02.npz", allow_pickle=False) as f:
                            numerical_clean = f["Y_fine_clean"]
                        norm_clean = float(np.linalg.norm(numerical_clean))
                        if pct == 2:
                            drawn_original = direction[p,s] * (.02*norm_clean/np.linalg.norm(direction[p,s]))
                            roundoff = 2*np.finfo(float).eps*np.max(np.abs(y))
                            if np.max(np.abs(y-numerical_clean-drawn_original)) > roundoff:
                                raise ArithmeticError("Source draw does not match the original frozen 2% direction.")
                        radius, sample_radius = meta["padded_radius"], meta["sample_radius"]
                        witness = check["true_distance_upper"]
                        margin = check["admissibility_margin_lower"]
                        error = check["clean_forward_error_upper"]
                        source = str((INHERITED / meta["file"]).relative_to(STUDY))
                    else:
                        y = yc + (pct/2.)*noise_2
                        yb = ball(y)
                        sample_radius = hi(norm(yb-ball(yc)))
                        radius = hi(arb(sample_radius) + arb(1)*norm(yb)/10**6)
                        error = cleaner_error
                        distance = norm(yb-yc_ball)+forward_error
                        witness = hi(distance)
                        margin = lo(arb(radius)-distance)
                        norm_clean = float(np.linalg.norm(yc))
                        source = "new native exact-forward enclosure and coupled draw"
                    if not margin > 0 or not radius > witness:
                        raise ArithmeticError(f"Generating-material witness failed for {ident}.")
                    file = ident + ".npz"
                    np.savez_compressed(public / file, Y=y)
                    rows.append(dict(id=ident, object_id=object_id, file=file, sources=s, receivers=p,
                                     noise_pct=pct, sample_radius=sample_radius, padded_radius=radius,
                                     source=source))
                    audit.append(dict(id=ident, generating_material_peak=peak, generating_family=family,
                                      clean_frobenius_norm=norm_clean, clean_nominal_error_upper=error,
                                      feasible_residual_upper=witness, feasible_margin_lower=margin))
        cost[object_id] = dict(**times(t,c), inverse=inverse, forward=forward)
        print(object_id, model.n, round(cost[object_id]["wall_seconds"], 2), flush=True)
    save_json(public / "manifest.json", dict(observations=rows))
    save_json(HERE / "forward_audit.json", dict(
        scope="Exact represented native operator, capped only by the generating peak in each witness",
        host=platform.platform(), expected_geometry="optfresh full ring 400 MHz",
        original_operator=str((PUBLIC / "fine_geometry.npz").relative_to(STUDY)),
        forward_proof="synthetic/nonvacuity.py: active inverse enclosure, Arb primal/adjoint correction",
        precision=ctx.prec, receiver_radius_m=RECEIVER_RADIUS, light_speed_m_per_s=LIGHT_SPEED,
        frequency_hz=FREQUENCY, side_m=WIDTH, common_guard=GUARD,
        new_noise_seed=NOISE_SEED, original_noise_seed=INHERITED_NOISE_SEED,
        shared_setup=setup, objects=cost, rows=audit, **times(start,cpu)))


if __name__ == "__main__":
    with threadpool_limits(1):
        main()
