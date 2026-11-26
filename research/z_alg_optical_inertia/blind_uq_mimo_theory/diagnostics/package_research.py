"""Create the three overlayable Blind UQ/MIMO evidence archives."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from tempfile import NamedTemporaryFile
from zipfile import ZIP_DEFLATED, ZIP_STORED, ZipFile

ROOT = Path(__file__).resolve().parents[1]
PREFIX = ROOT.name
ASSETS = {
    "study": "blind-uq-mimo-study.zip",
    "operators": "blind-uq-mimo-operators.zip",
    "scoring-truths": "blind-uq-mimo-scoring-truths.zip",
}
COHORT = {
    32: tuple(f"n32_c{i:02d}" for i in range(1, 9)),
    96: ("n96_c01", "n96_c02"),
}
STUDY_DIRS = (
    "bounds", "diagnostics", "figures", "imaging", "ports", "results",
    "theory", "vendor", "vendor_r8", "reviews/norm_before_correction",
)
PUBLIC_EXTENSIONS = {
    ".c", ".css", ".csv", ".html", ".js", ".json", ".jsonl", ".log",
    ".md", ".npy", ".npz", ".pdf", ".png", ".py", ".sh", ".svg", ".toml", ".txt",
}
ROOT_EXTENSIONS = {".css", ".html", ".js", ".md", ".txt"}
LOCAL_ROOT_FILES = {"ORIGINAL_REQUEST.txt", "TASK_RESUME.md"}
PUBLIC_REVIEW_FILES = (
    "RESPONSE_ZH.md", "neon_02.md", "neon_03.md", "neon_04.md",
    "RESPONSE_ZH.html", "neon_02.html", "neon_03.html",
    "index.html", "response.html", "review.html", "reader.html", "review.css",
    "review-02.html", "review-03.html", "review-04.html",
)
STORED_EXTENSIONS = {".npz", ".png", ".pdf"}
ACQUISITION_COST_PATH = Path("ports/results/reference_acquisition_costs.json")


def public_descendants(directory):
    for subdir, dirs, names in os.walk(directory, followlinks=False):
        folder = Path(subdir)
        relative = folder.relative_to(ROOT)
        omitted = {"__pycache__", "private"}
        if relative == Path("ports/public"):
            omitted.add("operators")
        dirs[:] = sorted(
            name for name in dirs
            if not name.startswith(".") and name not in omitted
            and not (folder / name).is_symlink()
        )
        for name in sorted(names):
            file = folder / name
            if name.startswith(".") or file.is_symlink() or file.suffix.lower() not in PUBLIC_EXTENSIONS:
                continue
            if "tomographic_dataset_" in name.casefold():
                continue
            if file.relative_to(ROOT) == ACQUISITION_COST_PATH:
                continue
            yield file.relative_to(ROOT)


def study_members():
    members = [
        f.relative_to(ROOT) for f in ROOT.iterdir()
        if f.is_file() and not f.is_symlink() and not f.name.startswith(".")
        and f.name not in LOCAL_ROOT_FILES and f.suffix.lower() in ROOT_EXTENSIONS
    ]
    for name in STUDY_DIRS:
        members.extend(public_descendants(ROOT / name))
    for name in PUBLIC_REVIEW_FILES:
        f = ROOT / "reviews" / name
        if f.is_file() and not f.is_symlink():
            members.append(f.relative_to(ROOT))
    return sorted(set(members))


def operator_members():
    directory = ROOT / "ports/public/operators"
    members = []
    for side in COHORT:
        members.extend(
            f.relative_to(ROOT) for f in directory.glob(f"n{side}_f*.npz")
            if f.is_file() and not f.is_symlink()
        )
    for side in COHORT:
        if Path(f"ports/public/operators/n{side}_f400_a08.npz") not in members:
            raise FileNotFoundError(f"Missing {side}² reference operator")
    return sorted(members)


def scoring_members():
    # Explicitly include exactly the scoring cohort and its generator. This
    # never walks the private source directory or follows its symlinks.
    filenames = [f"{case}.npz" for cases in COHORT.values() for case in cases]
    members = [Path("ports/private") / filename for filename in (*filenames, "generate.py")]
    if any(not (ROOT / f).is_file() or (ROOT / f).is_symlink() for f in members):
        raise FileNotFoundError("Missing a regular scoring-cohort source file")
    return sorted(members)


def acquisition_costs():
    # Transfer only the original numerical cost records into the public
    # evidence. Raw port-wave and truth files are selected separately.
    archived = ROOT / ACQUISITION_COST_PATH
    if archived.is_file() and not archived.is_symlink():
        return json.loads(archived.read_text())
    rows = []
    for cases in COHORT.values():
        for case in cases:
            private_record = ROOT / "ports/private" / f"{case}.json"
            source = json.loads(private_record.read_text())
            rows.append({k: source[k] for k in
                         ("case_id", "nominal_operator", "full_solve", "independent_bare_plus_loaded_solve")})
    return {
        "scope": "Original synthetic reference acquisition, one nominal VIE solve and independent loaded-network replay",
        "cases": rows,
    }


def packed_name(relative):
    return f"{PREFIX}/{relative.as_posix()}"


def build(part, output):
    members = {"study": study_members, "operators": operator_members,
               "scoring-truths": scoring_members}[part]()
    destination = output / ASSETS[part]
    records = [
        {"path": f.as_posix(), "bytes": (ROOT / f).stat().st_size}
        for f in members
    ]
    extra = {}
    if part == "study":
        extra[ACQUISITION_COST_PATH.as_posix()] = acquisition_costs()
    manifest = {
        "part": part,
        "relative_to": PREFIX,
        "archive": ASSETS[part],
        "entries": records,
        "generated_records": list(extra),
    }
    with NamedTemporaryFile(prefix=f".{part}-", suffix=".zip", dir=output, delete=False) as tmp:
        staging = Path(tmp.name)
    try:
        with ZipFile(staging, "w", allowZip64=True, compression=ZIP_DEFLATED, compresslevel=3) as archive:
            for f in members:
                compress = ZIP_STORED if f.suffix.lower() in STORED_EXTENSIONS else ZIP_DEFLATED
                archive.write(ROOT / f, packed_name(f), compress_type=compress)
            for name, value in extra.items():
                archive.writestr(packed_name(Path(name)), json.dumps(value, indent=2) + "\n")
            archive.writestr(packed_name(Path(f"archive-manifest-{part}.json")),
                             json.dumps(manifest, indent=2) + "\n")
        os.replace(staging, destination)
    finally:
        staging.unlink(missing_ok=True)
    print(f"{destination} ({len(members)} source files, {destination.stat().st_size:,} bytes)", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--part", action="append", choices=tuple(ASSETS), help="Only create the named part")
    parser.add_argument("--output-dir", type=Path, default=ROOT/"release")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for part in dict.fromkeys(args.part or ASSETS):
        build(part, args.output_dir)


if __name__ == "__main__":
    main()
