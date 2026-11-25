"""Create saved observations separately from all localization inference."""
import os
import json

for var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[var] = "1"

import numpy as np
from threadpoolctl import threadpool_limits

from common import (FREEZE, HERE, STUDY, Geometry, clock, elapsed, geom_path, geometry32,
                    pattern, place, save_json)


def prepare():
    public = HERE / "public"
    offline = HERE / "offline"
    public.mkdir(exist_ok=True)
    offline.mkdir(exist_ok=True)
    acquisition = {"geometry": {}, "observations": {}}
    for side in (32, 96):
        start = clock()
        if side == 32:
            np.savez_compressed(geom_path(side), **geometry32())
        else:
            native = STUDY / FREEZE["aperture"]["archive"]
            with np.load(native, allow_pickle=False) as arrays:
                np.savez_compressed(geom_path(side), **dict(arrays))
        acquisition["geometry"][side] = elapsed(start)
    geometry = {side: Geometry(geom_path(side)) for side in (32, 96)}
    model = FREEZE["model"]
    truth_config=json.loads((HERE / "freeze.json").read_text())["truth"]
    roster=truth_config["objects"]
    for name, side, shape, center, amplitude, seed in roster:
        started = clock()
        g = geometry[side]
        offsets, profile = pattern(side, shape)
        cells, weights = place(side, center, offsets, profile)
        rng = np.random.default_rng(seed)
        h = truth_config["footprint_microtexture_halfwidth"]
        micro = rng.uniform(-h, h, size=len(cells))
        background = rng.uniform(0, model["background_upper"], size=g.n)
        noise = rng.normal(size=(g.S.shape[0], g.B.shape[1]))
        noise = noise + 1j * rng.normal(size=noise.shape)
        noise /= np.linalg.norm(noise)
        for scene, x in (("vacuum", np.zeros(g.n)), ("weak", background.copy())):
            x[cells] = weights * (amplitude + micro)
            clean, current, info = g.forward(x)
            assert info == 0
            vive = current - x[:, None] * (g.B + g.apply(current))
            born = g.S @ (x[:, None] * g.B)
            r0 = model["noise_fraction"] * np.linalg.norm(clean)
            Y = clean + r0 * noise
            r = r0 + model["noise_radius_arithmetic_slack_fraction"] * np.linalg.norm(clean)
            oid = f"ipcl-fov-f1__{name}__{scene}__7tx9rx__n03"
            np.savez_compressed(public / f"{oid}.npz", Y=Y, r=r)
            np.savez_compressed(offline / f"{oid}.npz", truth=x, clean=clean, born=born,
                                centre=center, cells=cells, profile=weights, physical_current=current)
            acquisition["observations"][oid] = dict(side=side, shape=shape, scene=scene,
                                                       noise_radius=float(r))
            print(oid, f"Born={np.linalg.norm(born-clean)/np.linalg.norm(clean):.3f}",
                  f"VIE={np.linalg.norm(vive):.2e}", flush=True)
        acquisition["observations"][f"ipcl-fov-f1__{name}__weak__7tx9rx__n03"]["paired_generation"] = elapsed(started)
    save_json(public / "acquisition.json", acquisition)


if __name__ == "__main__":
    with threadpool_limits(limits=1, user_api="blas"):
        prepare()
