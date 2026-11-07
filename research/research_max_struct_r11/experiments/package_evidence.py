"""Lossless R11 scientific archives, after all scientific writers finish.

Reuse the R10 release convention: paths and byte sizes, no content hashes.
Redundant result transports, Python caches and third-party reading copies are
excluded. Historical source/input/overlay archives are execution evidence and
remain included together with their launch/transport receipts.
"""
from pathlib import Path
import argparse
import json
import time
import zipfile

ROOT = Path(__file__).resolve().parents[1]


def excluded(path):
    historical_archive = (path.suffix == ".zip" and path.name != "RESULT.zip") or path.name in {
        "PUBLIC_INVERSE_INPUTS_V2.tgz", "independent48_generation_sources.tgz",
        "independent48_public_start.tgz", "inverse_only.tgz"}
    return (bool(set(path.parts) & {"__pycache__", "primary_texts", "release_assets", ".git"})
            or path.suffix == ".pyc"
            or (not historical_archive and (path.suffix in {".zip", ".tgz", ".tar"} or path.name.endswith(".tar.gz")))
            or path.name.startswith("undermind_search_")
            or path.name in {"vector_fulltext_extraction.json", "verified_metadata.json"}
            or path.name in {"PUBLICATION_MANIFEST.json", "PACKAGING_RECEIPT.json", ".gitignore"})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-uncompressed-mb", type=int, default=96)
    parser.add_argument("--completed-science", action="store_true", help="Root confirms all writers stopped and postseal scoring is complete")
    args = parser.parse_args()
    if not args.completed_science:
        parser.error("Package only after all writers finish and postseal scoring is complete")
    wall, cpu = time.perf_counter(), time.process_time()
    out = ROOT / "release_assets"
    out.mkdir(exist_ok=False)
    files = sorted(p for p in ROOT.rglob("*") if p.is_file() and not excluded(p.relative_to(ROOT)))
    groups, group, size = [], [], 0
    for path in files:
        n = path.stat().st_size
        if group and size+n > args.max_uncompressed_mb*1024**2:
            groups.append(group)
            group, size = [], 0
        group.append(path)
        size += n
    if group:
        groups.append(group)
    manifest = {"format": 1, "scope": "All R11 scientific files; transport duplicates, caches and third-party fulltexts excluded",
                "archives": [], "hash_verification": "NOT_PERFORMED_USER_INSTRUCTION"}
    for i, group in enumerate(groups, 1):
        archive = out / f"max-struct-r11-evidence-{i:02d}.zip"
        entries = []
        with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as z:
            for path in group:
                relative = Path("research/research_max_struct_r11") / path.relative_to(ROOT)
                z.write(path, str(relative))
                entries.append({"path": str(relative), "bytes": path.stat().st_size})
        manifest["archives"].append({"name": archive.name, "bytes": archive.stat().st_size, "files": entries})
        print(archive.name, archive.stat().st_size, len(entries), flush=True)
    audit = ROOT / "audit"
    audit.mkdir(exist_ok=True)
    (audit / "PUBLICATION_MANIFEST.json").write_text(json.dumps(manifest, indent=2)+"\n")
    receipt = {"CPU_seconds": time.process_time()-cpu, "wall_seconds": time.perf_counter()-wall,
               "scientific_forward_calls": 0, "file_count": len(files),
               "uncompressed_bytes": sum(p.stat().st_size for p in files),
               "archive_bytes": sum(x["bytes"] for x in manifest["archives"]),
               "third_party_fulltexts_redistributed": False, "monetary_cost": None,
               "monetary_cost_reason": "No hardware or tool billing rate supplied"}
    (audit / "PACKAGING_RECEIPT.json").write_text(json.dumps(receipt, indent=2)+"\n")
    ignore = ["# Complete scientific arrays/histories are retained in the same repository R11 release.",
              "__pycache__/", "*.pyc", "*.zip", "*.tgz", "*.tar", "*.tar.gz", "/release_assets/", "/audit/literature/primary_texts/", "/audit/literature/undermind_search_*.json", "/audit/literature/vector_fulltext_extraction.json", "/audit/literature/verified_metadata.json"]
    for path in files:
        rel = path.relative_to(ROOT)
        raw = path.suffix in {".npz", ".npy"} and rel.parts[0] not in {"data_public", "sealed_truth", "FIGURES"}
        long_record = path.suffix in {".json", ".jsonl", ".csv", ".log"} and path.stat().st_size > 350_000
        if raw or long_record or path.stat().st_size > 20_000_000:
            ignore.append("/"+str(rel))
    (ROOT / ".gitignore").write_text("\n".join(ignore)+"\n")
    print(json.dumps(receipt))


if __name__ == "__main__":
    main()
