"""Score saved source-weight and restart experiments after all search arms have finished."""
import json
from pathlib import Path
import numpy as np
from high_score import ENDPOINTS, truth_for, material_scores

OUT = Path(__file__).resolve().parents[1] / "results"


def main():
    rows = []
    for path in sorted(OUT.glob("basin_*.npz")):
        label = path.stem.removeprefix("basin_").removeprefix("controls_")
        record = json.loads(path.with_suffix(".json").read_text())
        truth, support, amplitude, _ = truth_for(ENDPOINTS[label].case_id)
        with np.load(path) as data:
            for key in data.files:
                if key.startswith("beam_") or key.endswith("_temporary"):
                    continue
                metrics = material_scores(data[key], truth, support, amplitude)
                stats = {}
                if key == "initial":
                    stats = record.get("initial", record.get("geometry", {}))
                elif "branches" in record:
                    arm = record["branches"][key]
                    stats = {**arm, **arm["restored"], "RHS": arm["rhs"]}
                else:
                    name = key.removesuffix("_long")
                    arm = record["results"][name]
                    stats = dict(arm)
                    if key.endswith("_long"):
                        stats.update(arm["long"])
                        for field in ["LU", "RHS", "wall_s"]:
                            stats[field] += arm[field]
                rows.append(dict(file=path.stem, case=label, method=key, **metrics,
                                 **{k:stats[k] for k in ("phi","relative_residual","LU","RHS","wall_s") if k in stats}))
    (OUT / "basin_scores.json").write_text(json.dumps(rows, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
