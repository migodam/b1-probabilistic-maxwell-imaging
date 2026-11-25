"""Pack and independently replay the fixed-operator article's scientific evidence."""
import argparse
from collections import Counter
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from zipfile import ZipFile, ZIP_DEFLATED


PAPER = Path(__file__).resolve().parents[1]
STUDY = PAPER.parent
REPOSITORY = PAPER.parents[2]
DEFAULT_OUTPUT = PAPER / "fixed_operator_reproduction.zip"
MANIFEST_FILE = "RESEARCH_ARCHIVE_MANIFEST.json"
REPRODUCTION_GUIDE = "REPRODUCE_FIXED_OPERATOR.md"
REQUIREMENTS_FILE = "requirements-fixed-operator.txt"
FILE_CAP_MIB = 95

# Keep imported paths in their original locations. The optional legacy box
# controls allow the inherited replayer's other declared input IDs as well.
STUDY_DEPENDENCIES = {
    "native_basis": (
        "synthetic/public/fine_geometry.npz",
        "physical_revision/array/observations/arrays.npz",
        "physical_revision/array/observations/manifest.json",
        "physical_revision/array/observations/optfresh_01__pw360_32x32__pair2pct.npz",
        "physical_revision/array/observations/optfresh_02__pw360_32x32__pair2pct.npz",
    ),
    "native_disk": (
        "ip_guided_v4/boxes/run.py",
        "ip_guided_v4/boxes/verify.py",
        "ip_guided_v4/boxes/results/followup_disk.json",
        "ip_guided_v4/boxes/results/followup_disk_witnesses.npz",
        "ip_guided_v4/boxes/results/disk_radius.json",
        "ip_unified/numerics/feasibility/strict_endpoints/control_disk075__pw360_32x32__pair2pct.npz",
        "ip_unified/numerics/feasibility/strict_endpoints/control_disk075__pw360_32x32__pair2pct.json",
    ),
    "native_ring": (
        "ip_final_v5/strong_control/run_strong_ring.py",
        "ip_final_v5/strong_control/strong_ring_arrays.npz",
        "ip_final_v5/strong_control/strong_ring_verified.json",
    ),
    "phase_A_witness": ("ip_guided_v4/pool/pool_common.py",),
    "scaling_and_transfer": (
        "validation/max_struct/summary.json",
        "validation/transfer/summary.json",
        "ip_guided_v4/routing/all_routing_cases.csv",
        "ip_guided_v4/routing/cie_diagnostics.csv",
        "ip_guided_v4/routing/early_exit_diagnostics.csv",
        "ip_guided_v4/routing/historical_l2.csv",
        "ip_guided_v4/routing/measured_port_cases.csv",
        "ip_guided_v4/routing/paired_arms.csv",
        "ip_guided_v4/routing/principal_28.csv",
    ),
    "legacy_box_inputs": (
        "synthetic/public/coarse_green.npz",
        "synthetic/public/pw360_10x20.npz",
        "ip_guided_v4/boxes/results/shiftbox32_01__pw360_10x20__n02.json",
        "ip_guided_v4/boxes/results/shiftbox32_01__pw360_10x20__n02.npz",
        "ip_guided_v4/boxes/results/shiftbox32_01__pw360_10x20__n02__offline.npz",
    ),
    "source_provenance": (
        "fixed_operator_20261126/Z_ALG_Fixed_Operator_Research.pdf",
        "fixed_operator_20261126/Z_ALG_Research_Package.zip",
        "ip_guided_v4/manuscript_ip.tex",
        "ip_guided_v4/theory.tex",
        "ip_guided_v4/results.tex",
        "ip_guided_v4/technical.tex",
        "ip_final_v5/manuscript_ip.tex",
        "ip_final_v5/theory.tex",
        "ip_final_v5/results.tex",
        "ip_final_v5/technical.tex",
    ),
}
REPOSITORY_DEPENDENCIES = (
    "research/z_alg_secant_tasks/highdim/public_vie_inputs.npz",
    "research/z_alg_secant_tasks/highdim/observation_radii.json",
    "research/z_alg_secant_tasks/highdim/results/r10_twod_0001.npz",
)
MANDATORY_LOCAL = (
    "src/replay_native.py", "src/replay_package.py", "protocol/selfcheck.py",
    "experiments/A/replay.py", "experiments/A/results/acquisition.npz",
    "experiments/A/results/DECISION_SEAL.json",
    "experiments/A/results/sealed_truth.npz",
    "experiments/A/psd_followup/positive_tail_summary.json",
    "experiments/A/psd_followup/positive_tail_upper_audit.json",
    "experiments/B/replay.py", "experiments/B/public/native_geometry.npz",
    "experiments/C/model.py", "experiments/D/error_audit.py",
    "experiments/native/plot.py", "results/native_replay.json",
)
EXCLUDED_PARTS = {
    "__pycache__", ".git", ".svn", ".hg", ".venv", "venv", ".cache",
    ".artifacts", ".pytest_cache", ".mypy_cache", ".ruff_cache", ".mplconfig",
    "node_modules", "private_prompts", "personal_prompts", "__MACOSX",
}
EXCLUDED_NAMES = {
    "ORIGINAL_REQUEST.txt", "prompt.md", "prompt.txt", "personal_prompt.md",
    "personal_prompt.txt", ".DS_Store", ".env", "fixed_operator_reproduction.zip",
}
EXCLUDED_SUFFIXES = {".pyc", ".pyo", ".part", ".tmp", ".aux", ".fdb_latexmk", ".fls"}
PINNED_REQUIREMENTS = (
    "numpy==2.5.5", "scipy==1.18.1", "python-flint==0.9.0",
    "threadpoolctl==3.8.0", "matplotlib==3.11.2", "cvxpy==1.9.3",
    "clarabel==0.11.1", "scs==3.3.2", "sympy==1.14.0",
)
README = """# Fixed-operator material sets — reproduction archive

Extract this archive into an empty directory and run commands from that directory.
Paths match their locations in the original repository. All code, arrays, failure
records and saved certificates needed for independent replay are inside the archive.
The tested runtime is Python 3.13 with the package versions listed in
`requirements-fixed-operator.txt`.

```sh
python -m venv .replay-venv
. .replay-venv/bin/activate
python -m pip install -r requirements-fixed-operator.txt
P=research/z_alg_optical_inertia/fixed_operator_paper
python -B "$P/src/package_research.py" --inventory
python -B "$P/src/package_research.py" --replay-a
python -B "$P/src/package_research.py" --replay-native
```

The Phase A command ignores previously written replay decisions. It recomputes the
strict certificates and witnesses for all 24 frozen pairs from the public acquisition
and saved frames and writes `results/archive_phase_a_replay.json` in the extracted
copy. To inspect a subset, use `--begin 1 --end 1`. Truth is not read by this command.

The native command recomputes the two stored 96×96 controls, including the original
field operator and one-sided Arb proof, and rewrites `results/native_replay.json` in
the extracted copy. It can take several minutes. Both commands preserve the ZIP.

Further entry points:

- `experiments/A/README.md` and `psd_followup/README.md` under Phase A explain the
  blinded evaluation and subsequent diagnostics.
- `experiments/B/replay.py` replays the material map results.
- `experiments/C` contains the active and supplemental source-direction evidence.
- `experiments/D/error_audit.py` checks the archived coarse/fine and laboratory
  summaries. Raw facility measurements remain with their original providers.
- `experiments/native/plot.py` redraws the native and historical scaling figures.
- `protocol/selfcheck.py` and `src/replay_package.py` replay the public API and
  supplied mechanism programs.
- `manuscript_ip.pdf`, `RESEARCH_REPORT_ZH.md`, and `NEXT_PRO_HANDOFF.md` explain
  the scope, outcomes, and limitations in the paper directory.

The archive includes the frozen synthetic truth for offline scoring and preserves
the archived failure records. The online replay does not read that truth. The
archive manifest lists the source paths and lengths so that the same known data may
be published as individual files. The combined release ZIP is a separate artifact.
"""


def excludes(path, destinations):
    name = path.name
    return (
        path in destinations
        or any(part in EXCLUDED_PARTS for part in path.relative_to(REPOSITORY).parts)
        or name in EXCLUDED_NAMES
        or name.startswith(("._", "fixed_operator_reproduction"))
        or path.suffix in EXCLUDED_SUFFIXES
        or name.endswith(".synctex.gz")
    )


def inventory(destinations=()):
    paths = {}
    found_missing = []
    for role, names in STUDY_DEPENDENCIES.items():
        for name in names:
            path = STUDY / name
            if path.is_file():
                paths[path] = role
            else:
                found_missing.append(path)
    for name in REPOSITORY_DEPENDENCIES:
        path = REPOSITORY / name
        if path.is_file():
            paths[path] = "legacy_box_inputs"
        else:
            found_missing.append(path)

    sample_paths = sorted((STUDY / "synthetic/evaluation").glob("*__*__n02.npz"))
    if len(sample_paths) != 12:
        raise ValueError(f"Expected 12 archived fine/coarse samples, found {len(sample_paths)}.")
    for path in sample_paths:
        paths[path] = "fine_coarse_samples"

    for name in MANDATORY_LOCAL:
        if not (PAPER / name).is_file():
            found_missing.append(PAPER / name)
    for pair in range(1, 25):
        for subdirectory, suffix in (("raw", ".json"), ("raw", "_arrays.npz"), ("replay", ".json")):
            path = PAPER / "experiments/A/results" / subdirectory / f"pair_{pair:02d}{suffix}"
            if not path.is_file():
                found_missing.append(path)
    if found_missing:
        raise FileNotFoundError(
            "Required portable reproduction files are missing:\n"
            + "\n".join(str(path.relative_to(REPOSITORY)) for path in found_missing)
        )

    for path in PAPER.rglob("*"):
        if path.is_file() and not excludes(path, destinations):
            paths[path] = "fixed_operator_paper"
    for path in paths:
        if path.resolve() != path:
            raise ValueError(f"Archive source must be a real repository-relative file: {path}")

    return [
        {"path": path.relative_to(REPOSITORY).as_posix(),
         "bytes": path.stat().st_size, "role": paths[path]}
        for path in sorted(paths)
    ]


def catalog(records, cap_mib):
    by_role = Counter(row["role"] for row in records)
    largest = max(records, key=lambda row: row["bytes"])
    over_cap = [row for row in records if row["bytes"] >= cap_mib * 1024**2]
    return {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "paths_are_repository_relative": True,
        "regular_files": len(records),
        "input_bytes": sum(row["bytes"] for row in records),
        "files_by_role": dict(sorted(by_role.items())),
        "largest_file": largest,
        "individual_file_cap_mib": cap_mib,
        "over_individual_file_cap": over_cap,
        "paths_for_individual_publication": records,
        "archive_is_separate_release_asset": True,
        "runtime_reference": dict(lineage="Python 3.13", requirements=PINNED_REQUIREMENTS),
    }


def json_text(value):
    return json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n"


def generated_files(contents):
    return {
        MANIFEST_FILE: json_text(contents),
        REPRODUCTION_GUIDE: README,
        REQUIREMENTS_FILE: "\n".join(PINNED_REQUIREMENTS) + "\n",
    }


def stage(destination, records, contents):
    destination.mkdir(parents=True, exist_ok=True)
    if any(destination.iterdir()):
        raise ValueError(f"The staging directory must be empty: {destination}")
    for row in records:
        relative = Path(row["path"])
        path = destination / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(REPOSITORY / relative, path)
    for name, text in generated_files(contents).items():
        (destination / name).write_text(text)


def archive(destination, records, contents):
    destination.parent.mkdir(parents=True, exist_ok=True)
    expected = {row["path"] for row in records} | generated_files(contents).keys()
    fd, staging_name = tempfile.mkstemp(prefix=".fixed-operator-", suffix=".part", dir=destination.parent)
    os.close(fd)
    try:
        with ZipFile(staging_name, "w", compression=ZIP_DEFLATED, compresslevel=4, allowZip64=True) as bundle:
            for row in records:
                source = REPOSITORY / row["path"]
                bundle.write(source, arcname=row["path"])
                if bundle.getinfo(row["path"]).file_size != row["bytes"]:
                    raise RuntimeError(f"Source size changed during packaging: {row['path']}")
            for name, text in generated_files(contents).items():
                bundle.writestr(name, text)
        with ZipFile(staging_name) as bundle:
            if set(bundle.namelist()) != expected or len(bundle.namelist()) != len(expected):
                raise RuntimeError("The archive membership differs from the inventory.")
            corrupt = bundle.testzip()
            if corrupt is not None:
                raise RuntimeError(f"The archive CRC read failed: {corrupt}")
        os.chmod(staging_name, 0o644)
        os.replace(staging_name, destination)
    finally:
        Path(staging_name).unlink(missing_ok=True)


def replay_phase_a(begin, end, output):
    directory = PAPER / "experiments/A"
    sys.path.insert(0, str(directory))
    from common import digest, load, PUBLIC, RESULTS, SEAL, write_json
    from replay import check_pair
    seal = json.loads(SEAL.read_text())
    if digest(PUBLIC) != seal["files"][PUBLIC.name]:
        raise RuntimeError("The packaged Phase A acquisition differs from its frozen receipt.")
    data = load(PUBLIC)
    decision_seal = json.loads((RESULTS / "DECISION_SEAL.json").read_text())
    records = []
    for index in range(begin - 1, end):
        for suffix in (".json", "_arrays.npz"):
            name = f"raw/pair_{index + 1:02d}{suffix}"
            if digest(RESULTS / name) != decision_seal["files"][name]:
                raise RuntimeError(f"The packaged decision material differs from its frozen receipt: {name}")
        record = check_pair(data, RESULTS / "raw", index, audit_inertia=True)
        records.append(record)
        write_json(output, {"phase": "A", "results": records, "saved_decisions_reused": False})
        print(f"Phase A pair {index + 1:02d}: {record['checked']}; pass {record['pass_all']}", flush=True)
        if not record["pass_all"]:
            raise RuntimeError(f"Phase A proof replay failed for pair {index + 1}.")


def replay_native():
    environment = dict(os.environ, PYTHONPATH="", PYTHONDONTWRITEBYTECODE="1")
    subprocess.run([sys.executable, "-B", str(PAPER / "src/replay_native.py")],
                   cwd=REPOSITORY, env=environment, check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--inventory", action="store_true", help="check input dependencies without writing a ZIP")
    mode.add_argument("--stage", type=Path, help="materialize the inventory in an empty directory")
    mode.add_argument("--replay-a", action="store_true", help="independent proof replay from extracted archive")
    mode.add_argument("--replay-native", action="store_true", help="replay the two extracted 96×96 controls")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--inventory-json", type=Path)
    parser.add_argument("--file-cap-mib", type=float, default=FILE_CAP_MIB)
    parser.add_argument("--begin", type=int, default=1)
    parser.add_argument("--end", type=int, default=24)
    parser.add_argument("--replay-output", type=Path, default=PAPER / "results/archive_phase_a_replay.json")
    args = parser.parse_args()
    if args.replay_a:
        if not 1 <= args.begin <= args.end <= 24:
            parser.error("Require 1 <= begin <= end <= 24.")
        replay_phase_a(args.begin, args.end, args.replay_output)
    elif args.replay_native:
        replay_native()
    else:
        ignored = {args.output.resolve()}
        if args.inventory_json:
            ignored.add(args.inventory_json.resolve())
        records = inventory(ignored)
        contents = catalog(records, args.file_cap_mib)
        if args.inventory_json:
            args.inventory_json.parent.mkdir(parents=True, exist_ok=True)
            args.inventory_json.write_text(json_text(contents))
        print(json_text({k: v for k, v in contents.items() if k != "paths_for_individual_publication"}), end="")
        if args.stage:
            stage(args.stage, records, contents)
        elif not args.inventory:
            archive(args.output, records, contents)
            print(f"Verified ZIP: {args.output} ({args.output.stat().st_size:,} bytes)", flush=True)


if __name__ == "__main__":
    main()
