"""Read only frozen noisy port fields and the corresponding published radii."""

from pathlib import Path
import json
import numpy as np

PAPER = Path(__file__).resolve().parents[2] / "z_alg_optical_inertia"


def observations(supplement=False):
    if supplement:
        pilot = Path(__file__).resolve().parents[1]
        with np.load(pilot / "data/srap_supplement.npz", allow_pickle=False) as fields:
            for i, oid in enumerate(fields["object_id"]):
                for j, noise in enumerate(fields["noise_pct"]):
                    yield (
                        str(oid),
                        int(noise),
                        float(fields["radius"][i, j]),
                        fields["Y"][i, j],
                        (0.75,),
                    )
        return
    base = PAPER / "submission/curve/observations"
    manifest = json.loads((base / "manifest.json").read_text())["observations"]
    for row in manifest:
        if row["noise_pct"] not in (1, 2, 10):
            continue
        with np.load(base / row["file"], allow_pickle=False) as f:
            y = f["Y"]
        yield (
            row["object_id"],
            int(row["noise_pct"]),
            float(row["padded_radius"]),
            y,
            (0.5, 0.75, 1.0),
        )
    base = PAPER / "validation/srap"
    meta = {
        v["id"]: v for v in json.loads((base / "online.json").read_text())["objects"]
    }
    with np.load(base / "public_observations.npz", allow_pickle=False) as f:
        for key, y in zip(f["ids"], f["Y"]):
            row = meta[str(key)]
            yield row["object_id"], 10, float(row["padded_radius"]), y, (0.5, 0.75)
