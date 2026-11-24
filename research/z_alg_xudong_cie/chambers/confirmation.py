"""Score a sealed confirmation prediction manifest without tuning or target refits."""

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from threadpoolctl import threadpool_limits

from score import load_material, score_material
from training import CASES, elapsed, timers


def score_confirmation(manifest, output, freeze_record):
    start = timers()
    rows = []
    entries = json.loads(manifest.read_text(encoding="utf-8"))["predictions"]
    for record in entries:
        name = record["case"]
        if name not in CASES or CASES[name].phase != "confirmation":
            continue
        path = Path(record["path"])
        if not path.is_absolute():
            path = manifest.parent / path
        row_start = timers()
        try:
            m = load_material(path, record.get("array", "m"), record.get("index"))
            scored = score_material(name, m, freeze_record=freeze_record)
            rows.append(dict(status="scored", source=record, source_resolved_path=str(path.resolve()), **scored))
        except (ArithmeticError, RuntimeError, OSError, ValueError, KeyError) as exc:
            rows.append(dict(status="unresolved", source=record, reason=str(exc), **elapsed(row_start)))
        output.write_text(json.dumps(dict(rows=rows, **elapsed(start)), indent=2, allow_nan=False) + "\n")
        print(name, record.get("arm", ""), rows[-1]["status"], flush=True)
    return rows


def plot_results(rows, path):
    good = { (r["case_id"], r["source"]["arm"]): r for r in rows if r["status"] == "scored"}
    cases = [name for name, case in CASES.items() if case.phase == "confirmation"]
    methods = list(dict.fromkeys(arm for _, arm in good))
    if not methods:
        return
    values = np.arange(len(methods))
    labels = ("Fresnel Int · 6 GHz", "Fresnel Ext · 6 GHz", "Manitoba E · 3.5 GHz")
    fig, axes = plt.subplots(1, len(cases), figsize=(9.5, 4.3), sharex=True, sharey=True,
                             layout="constrained")
    for ax, case_id, title in zip(axes, cases, labels, strict=True):
        for row_index, arm in enumerate(methods):
            result = good.get((case_id, arm))
            if result:
                train = result["train"]["relative_to_measured"]
                held = result["held_transmit"]["relative_to_measured"]
                ax.plot([train, held], [row_index, row_index], color="#88919a", linewidth=1)
                ax.scatter(train, row_index, color="#1a677b", marker="o", s=24)
                ax.scatter(held, row_index, color="#983b3c", marker="D", s=24)
        ax.set_title(title)
        ax.set_xlabel("Relative complex-field error")
        ax.set_yticks(values, methods)
        ax.grid(axis="x", alpha=.3)
    axes[0].invert_yaxis()
    axes[0].scatter([], [], color="#1a677b", marker="o", label="Training")
    axes[0].scatter([], [], color="#983b3c", marker="D", label="Held transmitter")
    axes[0].legend(frameon=False, loc="upper left")
    for suffix in (".png", ".pdf"):
        fig.savefig(path.with_suffix(suffix), dpi=220, facecolor="white")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    input_data = parser.add_mutually_exclusive_group(required=True)
    input_data.add_argument("--predictions", type=Path)
    input_data.add_argument("--scores", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--freeze-record", type=Path)
    parser.add_argument("--plot", type=Path)
    args = parser.parse_args()
    if args.scores:
        if args.plot is None:
            parser.error("--scores requires --plot.")
        records = json.loads(args.scores.read_text(encoding="utf-8"))["rows"]
    else:
        if args.output is None or args.freeze_record is None:
            parser.error("--predictions requires --output and --freeze-record.")
        with threadpool_limits(1):
            records = score_confirmation(args.predictions, args.output, args.freeze_record)
    if args.plot:
        plot_results(records, args.plot)
