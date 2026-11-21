"""Recover design cells from public, previously scored manifests. No truth or solver inputs opened."""
import argparse
import csv
import json
from collections import Counter
from pathlib import Path

from lookup import read_rows, find_input


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path)
    parser.add_argument("--repo", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    tables = {
        "r13_fresh48": "data_public/fresh48/MANIFEST.json",
        "r14_independent24": "data_public/independent24/MANIFEST.json",
    }
    rows = {r["dataset"]: r for r in read_rows("COHORTS.csv")}
    with args.output.open("w", newline="", encoding="utf8") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["dataset", "forward_side", "inverse_side", "rx_aperture",
                         "sources", "receivers", "relative_complex_RMS_noise", "material_count"])
        for key, relative in tables.items():
            path = find_input(rows[key], relative, args.workspace, args.repo)
            if path is None:
                raise FileNotFoundError(f"{key}/{relative}: provide --workspace or --repo")
            manifest = json.loads(path.read_text(encoding="utf8"))["case_files"]
            bins = Counter((r["forward_side"], r["inverse_side"], r["receiver_aperture"],
                            r["sources"], r["receivers"], r["relative_noise_RMS"]) for r in manifest)
            for condition, count in sorted(bins.items()):
                writer.writerow([key, *condition, count])


if __name__ == "__main__":
    main()
