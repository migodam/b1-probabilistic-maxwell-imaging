"""Archive R9 evidence losslessly; inventory paths/sizes without content hashes.

Run only after scientific writers have finished. Third-party reading caches and
Python caches are not research outputs and are deliberately not redistributed.
"""
from pathlib import Path
import argparse
import json
import time
import zipfile

ROOT = Path(__file__).resolve().parents[1]


def excluded(path):
    return bool(set(path.parts) & {"__pycache__", "primary_texts", "release_assets", ".git"}) or path.suffix == ".pyc"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-uncompressed-mb", type=int, default=96)
    args = parser.parse_args()
    started, cpu = time.time(), time.process_time()
    out = ROOT / "release_assets"
    out.mkdir(exist_ok=True)
    files = sorted(p for p in ROOT.rglob("*") if p.is_file() and not excluded(p.relative_to(ROOT)))
    # The generated inventory and final packaging receipt are kept in Git and
    # attached separately, avoiding a self-referential archive inventory.
    files = [p for p in files if p.name not in {"PUBLICATION_MANIFEST.json", "PACKAGING_RECEIPT.json", ".gitignore"}]
    groups, group, size = [], [], 0
    for path in files:
        n = path.stat().st_size
        if group and size + n > args.max_uncompressed_mb * 1024**2:
            groups.append(group)
            group, size = [], 0
        group.append(path)
        size += n
    if group:
        groups.append(group)
    manifest = {"format": 1, "scope": "All R9 scientific files, except third-party reading and Python caches", "archives": []}
    for i, group in enumerate(groups, 1):
        archive = out / f"max-struct-r9-evidence-{i:02d}.zip"
        entries = []
        with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as z:
            for path in group:
                relative = Path("research/research_max_struct_r9") / path.relative_to(ROOT)
                z.write(path, str(relative))
                entries.append({"path": str(relative), "bytes": path.stat().st_size})
        manifest["archives"].append({"name": archive.name, "bytes": archive.stat().st_size, "files": entries})
        print(archive.name, archive.stat().st_size, len(entries), flush=True)
    audit = ROOT / "audit"
    audit.mkdir(exist_ok=True)
    (audit / "PUBLICATION_MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n")
    receipt = {"CPU_seconds": time.process_time() - cpu, "wall_seconds": time.time() - started,
               "scientific_forward_calls": 0, "file_count": sum(len(x) for x in groups),
               "uncompressed_bytes": sum(p.stat().st_size for p in files),
               "archive_bytes": sum(x["bytes"] for x in manifest["archives"]),
               "third_party_fulltexts_redistributed": False}
    (audit / "PACKAGING_RECEIPT.json").write_text(json.dumps(receipt, indent=2) + "\n")
    # Keep compact outcomes, public inputs, source, reports and figures in Git.
    # Complete states and long histories remain in the lossless release above.
    ignore = ["# Full scientific evidence is retained in the same repository's R9 release.",
              "__pycache__/", "*.pyc", "*.zip", "/release_assets/", "/delegated/nonconvex_prior_art/primary_texts/"]
    for path in files:
        rel = path.relative_to(ROOT)
        raw_array = path.suffix in {".npz", ".npy"} and (rel.parts[0] not in {"data_public", "sealed_truth", "FIGURES"} or "generation_checkpoints" in rel.parts)
        long_record = path.suffix in {".json", ".jsonl", ".csv", ".log"} and path.stat().st_size > 350_000
        if raw_array or long_record or path.stat().st_size > 20_000_000:
            ignore.append("/" + str(rel))
    (ROOT / ".gitignore").write_text("\n".join(ignore) + "\n")
    print(json.dumps(receipt))


if __name__ == "__main__":
    main()
