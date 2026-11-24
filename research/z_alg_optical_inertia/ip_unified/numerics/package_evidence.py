#!/usr/bin/env python3
"""Package the unified numerical evidence and its public replay dependencies."""

import argparse
from collections import Counter
import importlib.util
from pathlib import Path
from shutil import copyfileobj
from zipfile import ZIP_DEFLATED, ZIP_STORED, ZipFile, ZipInfo


HERE = Path(__file__).resolve().parent
UNIFIED = HERE.parent
STUDY = UNIFIED.parent
ROOT = STUDY.parents[1]
ASSET = "Z_ALG_Unified_Evidence.zip"
OUTPUT = ROOT / ".artifacts/z_alg_unified/release" / ASSET
FIXED_ZIP_DATE = (2026, 11, 25, 0, 0, 0)

EVIDENCE_SUFFIXES = {
    ".csv", ".json", ".jsonl", ".log", ".md", ".npz", ".pdf", ".png", ".py", ".tsv", ".txt",
}
HISTORICAL_STUDY = {
    "synthetic": "generate.py FINE_FEASIBILITY_ZH.md",
    "submission/single": """
        local_moment.py local_moment.json local_moment.npz
        coherent_dual.py general_fourier_verify.py
    """,
    "ip_submission_v2": "ENERGY_EVIDENCE_ZH.md energy_port_control.py energy_port_control.json",
    "ip_submission_v2/energy_sdp": """
        REPORT_ZH.md COST_LEDGER.json probe.py search_rays.py optimize_mean.py verify.py
        alt025.npz alt025.json verified_alt025.json energy_grid.csv
    """,
    "ip_submission_v2/energy_bound": "README.md",
    "ip_submission_v2/energy_bound/source": """
        verify_physical_current.py search_physical_current.py verify_energy_bound.py
    """,
    "ip_submission_v2/energy_bound/receipts": """
        coarse_12x12.npz pixel_profile.npz physical_current.json physical_search_replay.json
        radiation_transfer.json energy_bound.json
    """,
    "ip_submission_v3": "experiments.tex",
    "submission/diagnostic": "audit.py run.py",
    "validation/max_struct": """
        REPRODUCE.md FREEZE.md frozen.json public_provenance.json public_inputs.npz
        offline_scalar_scores.json scan.py verify.py offline.py summarize.py
        summary.json verified.json numerical_grid.json paired_results.csv
    """,
}
EXTERNAL_CODE = {
    "research/z_alg_xudong_cie/minpeak_native": "search.py gn_recovery.py cie_seed.py",
    "research/srap_cie_array_study_20261121": "engine.py legacy_la.py data2d.py",
    "research/delegated/max_struct_r1/r13_worktree/research": """
        research_max_struct_r10/src/lowrank_denominator.py
        research_z_alg_r4_diversity/src/physics.py
    """,
}


def named(base, groups):
    return [base / folder / filename for folder, names in groups.items()
            for filename in names.split()]


def public_physical_revision():
    builder = STUDY / "physical_revision/package_release.py"
    spec = importlib.util.spec_from_file_location("physical_revision_package", builder)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    # Reuse the historical public allowlist with its R8 source fallback.
    return mod.members()


def evidence_files():
    return [path for base in (HERE, UNIFIED / "figures") for path in base.rglob("*")
            if path.is_file()
            and path.suffix in EVIDENCE_SUFFIXES
            and all(not p.startswith(".") and p != "__pycache__" and "neon" not in p.lower()
                    and "private" not in p.lower() for p in path.relative_to(base).parts)]


def source_paths():
    paths = public_physical_revision()
    local = [UNIFIED / "REPRODUCE.md", *evidence_files(),
             *named(STUDY, HISTORICAL_STUDY), *named(ROOT, EXTERNAL_CODE)]
    paths.update({path.relative_to(ROOT): path for path in local})
    return paths


def make_zip(output, paths):
    output = output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    temp = output.with_suffix(".zip.part")
    with ZipFile(temp, "w", compression=ZIP_DEFLATED, compresslevel=3) as archive:
        for name, source in sorted(paths.items(), key=lambda item: item[0].as_posix()):
            info = ZipInfo(name.as_posix(), FIXED_ZIP_DATE)
            info.compress_type = ZIP_STORED if source.suffix in {".npz", ".png", ".pdf"} else ZIP_DEFLATED
            info.file_size = source.stat().st_size
            info.external_attr = 0o100644 << 16
            with source.open("rb") as stream, archive.open(info, "w") as entry:
                copyfileobj(stream, entry, length=1024 * 1024)
    if temp.stat().st_size >= 2_000_000_000:
        raise ValueError("The evidence archive exceeds the intended release asset size.")
    temp.replace(output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    paths = source_paths()
    missing = [str(path) for path, source in paths.items() if not source.is_file()]
    if missing:
        parser.error("Missing evidence members:\n" + "\n".join(sorted(missing)))
    make_zip(args.output, paths)
    kinds = Counter(path.suffix for path in paths)
    print(f"{args.output.resolve()} · {len(paths)} files · {kinds['.npz']} array archives"
          f" · {args.output.stat().st_size:,} bytes")


if __name__ == "__main__":
    main()
