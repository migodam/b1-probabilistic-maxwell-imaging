"""Prepare a bounded reading ZIP, never the complete raw evidence archive.

Default: write an explicit plan only. Root must run --build --source-commit HEAD
after its final scientific commit. No prediction, model, physics or GT reads.
"""
from pathlib import Path
import argparse
import json
import os
import subprocess
import time
import zipfile

HERE = Path(__file__).resolve().parent
P = HERE.parent
REPO = P.parent.parent
DELIVERY = REPO.parents[2] / "research/delegated/r6_delivery_20261027"
PRUNE = {"raw", "inputs", "models", "checkpoints", "weights", "contexts",
         "predictions", "original_optimizer", "__pycache__", ".git", ".venv",
         "read", "VISUAL_QA_INITIAL", "jobs", "native96_only_contexts",
         "matched_anchors_TRAIN512", "matched_heads_TRAIN512", "matched_inputs_TRAIN512"}
SUFFIXES = {".md", ".py", ".json", ".csv", ".png", ".svg"}
LIMITS = {".md": 2_000_000, ".py": 2_000_000, ".json": 1_000_000,
          ".csv": 2_000_000, ".png": 5_000_000, ".svg": 2_000_000}
REQUIRED = ["RESEARCH_REPORT_ZH.md", "RESULTS_SUMMARY.json", "GO_NO_GO.md",
            "THEORY_AND_COUNTEREXAMPLES.md", "START_HERE.md", "REPRODUCE.md",
            "GPT_PRO_HANDOFF/COMPACT_NUMERICAL_MANIFEST.json",
            "GPT_PRO_HANDOFF/COMPACT_NUMERICAL_README_ZH.md",
            "GPT_PRO_HANDOFF/compact/DLOW_COMPACT_NUMERICAL.npz",
            "GPT_PRO_HANDOFF/compact/DHIGH_COMPACT_NUMERICAL.npz"]
RAW_REFERENCE_NAMES = ["STABLE_SOURCE_MANIFEST.json", "STABLE_DELIVERY_RECEIPT.json",
                       "STABLE_SIZE_SUMMARY.json", "STABLE_SIZE_COMBINED.json",
                       "STABLE_SHARED_ADDENDUM.json", "PUBLICATION_LAYOUT_TEMPLATE.json"]


def git(*args):
    return subprocess.check_output(["git", *args], cwd=REPO, text=True).strip()


def collect():
    files, excluded = {}, []
    def add(path, archive_name, required=False):
        if path.is_symlink():
            raise ValueError("No symlinks in reading package: " + str(path))
        if not path.is_file():
            if required:
                raise FileNotFoundError(path)
            return
        files[archive_name] = path
    # Prune before descent. In particular this never traverses raw, contexts,
    # model files, historical input archives or prediction directories.
    for folder, dirs, names in os.walk(P):
        dirs[:] = sorted(d for d in dirs if d not in PRUNE and not d.startswith(".venv"))
        for name in sorted(names):
            path = Path(folder) / name
            if name in {"READING_ZIP_PLAN.json", "READING_ZIP_RECEIPT.json"}:
                continue
            if path.suffix.lower() not in SUFFIXES:
                continue
            if path.stat().st_size > LIMITS[path.suffix.lower()]:
                excluded.append({"path": str(path.relative_to(REPO)), "reason": "above bounded reading size",
                                 "bytes": path.stat().st_size})
                continue
            add(path, str(path.relative_to(REPO)))
    for relative in REQUIRED:
        path = P / relative
        add(path, str(path.relative_to(REPO)), required=True)
    # Generic repository-root Track B sources are part of this actual R6 source
    # layout; archive them without renaming or overwriting historical sources.
    for path in sorted((REPO / "track_b").glob("*.py")):
        add(path, str(path.relative_to(REPO)))
    for name in ("index.html", "report.html", "appendices.html"):
        add(P / name, str((P / name).relative_to(REPO)))
    # Include existing readable appendix pages for the static HTML links. This
    # enumerates only the generated read/ directory, never the raw evidence tree.
    if (P / "read").is_dir():
        for folder, dirs, names in os.walk(P / "read"):
            dirs[:] = sorted(d for d in dirs if d != "__pycache__")
            for name in sorted(names):
                path = Path(folder) / name
                if path.suffix == ".html":
                    add(path, str(path.relative_to(REPO)))
    for name in RAW_REFERENCE_NAMES:
        add(DELIVERY / name, "RAW_RELEASE_REFERENCE/" + name)
    # Root may finalize Release references after this preparation. Include only
    # actual small top-level manifests/receipts that exist at build time; no
    # recursive inventory or provisional public URL is invented here.
    for pattern in ("*RELEASE*MANIFEST*.json", "*PUBLICATION*RECEIPT*.json",
                    "FINAL*MANIFEST*.json", "FINAL*RECEIPT*.json"):
        for path in sorted(DELIVERY.glob(pattern)):
            if path.stat().st_size <= 1_000_000:
                add(path, "RAW_RELEASE_REFERENCE/" + path.name)
    return files, excluded


def main(build=False, source_commit=None):
    wall = time.perf_counter()
    if not (HERE / "COMPACT_NUMERICAL_MANIFEST.json").exists():
        raise RuntimeError("First finish the separately verified numerical handoff")
    compact = json.loads((HERE / "COMPACT_NUMERICAL_MANIFEST.json").read_text())
    if compact["status"] != "OFFLINE_LOSSLESS_HANDOFF_VERIFIED":
        raise ValueError("Compact verification missing")
    head = git("rev-parse", "HEAD")
    if build:
        if not source_commit or git("rev-parse", source_commit) != head:
            raise RuntimeError("Root must explicitly supply its final scientific commit, equal to HEAD")
        # Ensure the key final reports and builder correspond to that commit.
        # Newly built ZIP/receipt and ignored compact data may remain untracked.
        check_paths = [str((P / r).relative_to(REPO)) for r in REQUIRED[:6]]
        check_paths += [str(Path(__file__).relative_to(REPO))]
        tracked = git("ls-files", "--", *check_paths).splitlines()
        if set(tracked) != set(check_paths):
            raise RuntimeError("Final reports and ZIP builder must be committed first")
        if git("diff", "HEAD", "--", *check_paths):
            raise RuntimeError("Final scientific files differ from the supplied commit")
    files, excluded = collect()
    total = sum(path.stat().st_size for path in files.values())
    if total > 150_000_000:
        raise RuntimeError("Reading package exceeds 150MB uncompressed guard; review the plan")
    reference = dict(scope="Reading ZIP is NOT the full raw evidence archive",
                     local_delivery_directory=str(DELIVERY),
                     release_source_manifest="RAW_RELEASE_REFERENCE/STABLE_SOURCE_MANIFEST.json",
                     complete_raw_restoration="Use actual Release manifests and the separately delivered lossless raw parts",
                     public_urls="Use only final root publication receipt; no draft/template URL is asserted live here",
                     exclusions="Large raw/current/field/J/H/optimizer traces, historical input ZIP, checkpoints and virtual environments",
                     GT_exposed=True, not_posterior_samples=True,
                     conditions_not_pooled=True, new_physics_RHS=0, new_training=0, new_inference=0)
    plan = dict(status="PLAN_ONLY_NO_ZIP" if not build else "BUILD_AUTHORIZED_BY_EXPLICIT_ROOT_COMMIT",
                science_commit=head, root_must_execute_after_final_scientific_commit=True,
                entries=[dict(path=name, bytes=path.stat().st_size) for name, path in sorted(files.items())],
                files=len(files), uncompressed_source_bytes=total, excluded_by_size=excluded,
                raw_reference=reference)
    (HERE / "READING_ZIP_PLAN.json").write_text(json.dumps(plan, ensure_ascii=False, indent=2) + "\n")
    if not build:
        print(json.dumps({"status": plan["status"], "files": len(files), "source_bytes": total}))
        return
    output = DELIVERY / "Research_P3LR6_Results_and_GPT_Pro_Handoff.zip"
    DELIVERY.mkdir(parents=True, exist_ok=True)
    if output.exists():
        raise FileExistsError(output)
    temporary = output.with_suffix(".tmp.zip")
    with zipfile.ZipFile(temporary, "x", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
        for name, path in sorted(files.items()):
            archive.write(path, name)
        archive.writestr("RAW_RELEASE_REFERENCE/READING_SCOPE.json",
                         json.dumps(reference, ensure_ascii=False, indent=2) + "\n")
        archive.writestr("READING_PACKAGE_MANIFEST.json",
                         json.dumps(plan, ensure_ascii=False, indent=2) + "\n")
    # Byte-for-byte unpack check for every included source; no extra scientific
    # computation or cryptographic checks are introduced by this ZIP builder.
    with zipfile.ZipFile(temporary) as archive:
        for name, path in files.items():
            if archive.read(name) != path.read_bytes():
                raise ValueError("Reading ZIP byte mismatch: " + name)
    os.replace(temporary, output)
    receipt = dict(status="READING_ZIP_BYTE_ROUNDTRIP_PASS", science_commit=head,
                   path=str(output), bytes=output.stat().st_size, source_files=len(files),
                   uncompressed_source_bytes=total, full_raw_evidence_included=False,
                   GT_exposed=True, not_posterior_samples=True,
                   new_physics_RHS=0, new_training=0, new_inference=0,
                   wall_seconds=time.perf_counter()-wall)
    (HERE / "READING_ZIP_RECEIPT.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", action="store_true", help="Root only, after the final scientific commit")
    parser.add_argument("--source-commit", help="Explicit final scientific commit, required only for --build")
    args = parser.parse_args()
    main(args.build, args.source_commit)
