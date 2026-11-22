"""Create independent public fields and separately retained scoring labels."""

import json
import os
from pathlib import Path
import sys
from time import perf_counter, process_time

sys.dont_write_bytecode = True
for key in ("VECLIB_MAXIMUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[key] = "1"

import numpy as np
from flint import arb, ctx
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
PILOT = HERE.parent
sys.path.insert(0, str(PILOT / "src"))
from generate_srap_bands import (
    ExactActive,
    ExactFourierGeometry,
    FINE,
    ball,
    high,
    load_geometry,
    norm,
)


def material(index, *, seed=202611230571, high_count=24):
    rng = np.random.default_rng(np.random.SeedSequence([seed, index]))
    q = (np.arange(96) + 0.5) / 48 - 1
    x, y = np.meshgrid(q, q)
    x, y = x - rng.uniform(-0.23, 0.23), y - rng.uniform(-0.23, 0.23)
    angle = rng.uniform(-np.pi, np.pi)
    co, si = np.cos(angle), np.sin(angle)
    x, y = (
        (co * x + si * y) / rng.uniform(0.34, 0.62),
        (-si * x + co * y) / rng.uniform(0.31, 0.56),
    )
    phi = rng.uniform(-np.pi, np.pi)
    group = index % 4
    if group == 0:
        shape = abs(x) ** 4 + abs(y) ** 4 < 0.94
        for cx, cy in ((-0.35, -0.16), (0.19, 0.30), (0.46, -0.30)):
            cx, cy = np.array((cx, cy)) + rng.uniform(-0.09, 0.09, 2)
            shape &= ((x - cx) / rng.uniform(0.12, 0.25)) ** 2 + (
                (y - cy) / rng.uniform(0.14, 0.28)
            ) ** 2 > 1
        value = shape * (0.78 + 0.22 * np.cos(2 * x + y) ** 2)
    elif group == 1:
        a = abs(y + 0.27 - 0.28 * np.sin(3.1 * x + phi)) < rng.uniform(0.13, 0.27)
        b = abs(y - 0.32 - 0.15 * np.sin(4.2 * x - phi)) < rng.uniform(0.13, 0.23)
        value = (abs(x) < 1.06) * np.maximum(a, rng.uniform(0.65, 0.95) * b)
    elif group == 2:
        rr = x * x + y * y
        outer = rr < (0.96 + 0.07 * np.cos(3 * np.arctan2(y, x) + phi)) ** 2
        hx, hy = rng.uniform(-0.18, 0.18, 2)
        hole = ((x - hx) / 0.56) ** 2 + ((y - hy) / 0.68) ** 2 < 1
        core = ((x - hx - 0.15) / 0.25) ** 2 + ((y - hy + 0.15) / 0.33) ** 2 < 1
        value = np.maximum(rng.uniform(0.56, 0.86) * (outer & ~hole), core)
    else:
        value = np.zeros_like(x)
        for cx, cy in ((-0.50, -0.28), (0.33, -0.18), (-0.10, 0.41), (0.56, 0.45)):
            cx, cy = np.array([cx, cy]) + rng.uniform(-0.14, 0.14, 2)
            ellipse = ((x - cx) / rng.uniform(0.27, 0.52)) ** 2 + (
                (y - cy) / rng.uniform(0.22, 0.45)
            ) ** 2 < 1
            value = np.maximum(value, rng.uniform(0.59, 1) * ellipse)
    peak = rng.uniform(2.6, 4.2) if index < high_count else rng.uniform(0.25, 0.65)
    return (value * peak / value.max()).ravel(), group, peak


def native_geometry():
    _, sensor, incident = load_geometry()
    geometry = ExactFourierGeometry(FINE)
    geometry.S, geometry.B = sensor, incident
    return geometry


def observe(geometry, m, *, index, noise_seed, noise_pct=(1, 2, 10)):
    action = ExactActive(geometry, m)
    inverse, _ = action.inverse_bound()
    integral, remainder, _ = action.enclose_field(inverse, 128)
    clean = np.array([complex(v) for v in integral.entries()]).reshape(20, 10)
    rng = np.random.default_rng(np.random.SeedSequence([noise_seed, index]))
    noise = rng.standard_normal(clean.shape) + 1j * rng.standard_normal(clean.shape)
    noise *= np.linalg.norm(clean) / np.linalg.norm(noise)
    ys, rs, slacks = [], [], []
    for pct in noise_pct:
        y = clean + pct / 100 * noise
        nominal = norm(ball(y) - ball(clean)).upper()
        r = high(nominal + norm(ball(y)) / 10**6)
        actual = norm(ball(y) - integral).upper() + remainder
        slack = arb(r) - actual
        if not slack > 0:
            raise ArithmeticError("The certified field must lie inside its declared ball.")
        ys.append(y)
        rs.append(r)
        slacks.append(slack)
    return ys, rs, slacks


def main(archived8=False):
    if archived8:
        import data_followup
        import data2d

        family_order = np.random.default_rng(data_followup.ORDER_SEED).permutation(8)
        centers = data2d.history.grid(96, data2d.WIDTH)
    object_count = 8 if archived8 else 32
    ctx.prec, ctx.threads = 112, 1
    start, cpu = perf_counter(), process_time()
    geometry = native_geometry()
    y_all, r_all, audit, labels = [], [], [], []
    for i in range(object_count):
        t, tt = perf_counter(), process_time()
        if archived8:
            family, variant = divmod(int(family_order[i]), 4)
            m = data_followup.material(centers, i + 1, family, variant, 96)
            oid = f"Cx{i + 1:02d}"
            peak = float(np.max(m))
            # The private historical image is only an optional identity comparison
            # for the offline generator. The public material recipe determines m.
            reference = data2d.HERE / "data_followup/truth" / f"{oid}__pw180_10x20.npz"
            identity_error = None
            if reference.is_file():
                with np.load(reference, allow_pickle=False) as ar:
                    identity_error = float(np.max(
                        abs(ar["m"].ravel() - data2d.area_average(m, 96, 32).ravel())
                    ))
                if identity_error > 1e-12:
                    raise ValueError("Regenerated fine material differs from the archived material.")
        else:
            m, family, peak = material(i)
            oid = f"ip_port_{i + 1:02d}"
            identity_error = None
        ys, rs, slacks = observe(geometry, m, index=i, noise_seed=202611230617 if archived8 else 202611230583)
        y_all.append(ys)
        r_all.append(rs)
        audit.append(
            dict(
                oid=oid,
                radius_slack_lower=[float(slack.lower()) for slack in slacks],
                archived_coarse_identity_max_abs=identity_error,
                wall_s=perf_counter() - t,
                cpu_s=process_time() - tt,
            )
        )
        labels.append(
            dict(
                oid=oid,
                geometry=(
                    ("continuous_high_low", "thin_bridge_multibody")
                    if archived8
                    else (
                        "perforated",
                        "bent_pair",
                        "eccentric_shell",
                        "composite_islands",
                    )
                )[family],
                peak=float(peak),
                under_cap=bool(peak <= 0.75),
            )
        )
        print(oid, round(audit[-1]["wall_s"], 2), flush=True)
    data = (
        PILOT
        / "data"
        / ("followup_external8.npz" if archived8 else "followup_observations.npz")
    )
    np.savez_compressed(
        data,
        object_id=np.array([r["oid"] for r in audit]),
        Y=np.array(y_all),
        radius=np.array(r_all),
        noise_pct=[1, 2, 10],
    )
    (
        HERE / ("external8_offline_scores.json" if archived8 else "offline_scores.json")
    ).write_text(json.dumps(labels, indent=2) + "\n")
    (
        HERE / ("external8_generation.json" if archived8 else "generation.json")
    ).write_text(
        json.dumps(
            dict(
                objects=object_count,
                observations=3 * object_count,
                physical_grid=96,
                precision_bits=112,
                native_field_rhs=10 * object_count,
                wall_s=perf_counter() - start,
                cpu_s=process_time() - cpu,
                rows=audit,
            ),
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    with threadpool_limits(limits=1):
        import argparse

        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument("--archived8", action="store_true")
        main(archived8=parser.parse_args().archived8)
