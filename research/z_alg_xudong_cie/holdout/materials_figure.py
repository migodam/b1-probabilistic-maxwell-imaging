"""Draw sealed blind examples at one common contrast range."""

import argparse
import json
from math import ceil
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from research.z_alg_xudong_cie.holdout.data import DEFAULT_ROOT, NATIVE_SIDE, require_freeze
from research.z_alg_xudong_cie.holdout.physics import nearest_lift

EXAMPLES = ("b3m_0001__pw180_10x20__n02", "b3m_0002__pw180_10x20__n02")


def render(scores_path, output, freeze_record, data_root, examples):
    require_freeze(freeze_record)
    scored = json.loads(scores_path.read_text(encoding="utf-8"))["rows"]
    selected = [r for r in scored if r["status"] in ("scored", "material_only")
                and r["source"]["case"] in examples]
    arms = sorted({r["source"]["arm"] for r in selected})
    if not arms:
        raise ValueError("The fixed example IDs have no scored material maps.")
    plots = []
    for example in examples:
        maps = {}
        oid = example.split("__")[0]
        with np.load(data_root / "offline" / "objects" / f"{oid}.npz", allow_pickle=False) as source:
            maps["Reference"] = source["m_native"].reshape(NATIVE_SIDE, NATIVE_SIDE)
        for row in (r for r in selected if r["source"]["case"] == example):
            meta = row["source"]
            with np.load(row["source_resolved_path"], allow_pickle=False) as source:
                m = source[meta.get("array", "m")]
            if "index" in meta:
                m = m[meta["index"]]
            maps[meta["arm"]] = nearest_lift(m).reshape(NATIVE_SIDE, NATIVE_SIDE)
        plots.append(maps)
    panels = ["Reference", *arms]
    columns = min(6, len(panels))
    rows_per_case = ceil(len(panels) / columns)
    rows = rows_per_case * len(examples)
    size = (2.65 * columns, 2.55 * rows)
    fig, axs = plt.subplots(rows, columns, figsize=size,
                            squeeze=False, layout="constrained")
    scale_cap = 4.0
    saturated = any(np.max(m) > scale_cap for group in plots for m in group.values())
    for i, (mapping, example) in enumerate(zip(plots, examples, strict=True)):
        for j in range(rows_per_case * columns):
            row_index, col_index = divmod(j, columns)
            ax = axs[i * rows_per_case + row_index, col_index]
            if j >= len(panels):
                ax.axis("off")
                continue
            arm = panels[j]
            if arm in mapping:
                image = ax.imshow(mapping[arm], vmin=0, vmax=scale_cap, cmap="viridis",
                                  origin="lower", interpolation="nearest",
                                  extent=(-1, 1, -1, 1))
            else:
                ax.text(.5, .5, "Unresolved", transform=ax.transAxes, ha="center", va="center")
            ax.set_xticks([-1, 0, 1])
            ax.set_yticks([-1, 0, 1])
            label = arm.replace("_", " ")
            if arm in mapping and np.max(mapping[arm]) > scale_cap:
                label += f" (peak {np.max(mapping[arm]):.2f})"
            ax.set_title(label)
            if col_index == 0:
                ax.set_ylabel(f"{example.split('__')[0]}\ny (m)")
            if row_index == rows_per_case - 1:
                ax.set_xlabel("x (m)")
    fig.colorbar(image, ax=axs.ravel().tolist(), fraction=.022, label="Contrast χ",
                 extend="max" if saturated else "neither")
    fig.suptitle("Fixed blind examples · 180° aperture · 2% noise")
    for suffix in (".png", ".pdf"):
        fig.savefig(output.with_suffix(suffix), dpi=220, facecolor="white")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("scores", type=Path)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--freeze-record", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    render(args.scores, args.output, args.freeze_record, args.root, EXAMPLES)
