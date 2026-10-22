"""Stream existing P3-L evidence into bounded, independently readable archives.

No action on import. No truth/NPZ/checkpoint interpretation, hash calculation,
physics, network, release upload or source mutation. See --help and the saved
AUDIT/evidence_packaging/README.md. Metadata-only planning is the default.
"""
from pathlib import Path
from contextlib import ExitStack
import argparse
import fcntl
import gzip
import itertools
import json
import os
import tarfile
import time

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "research/p3l_local_posterior_imaging"
DEST = OUT / "AUDIT/evidence_packaging"
MIB = 1024 * 1024
PART_LIMIT = 700 * MIB
PAYLOAD_LIMIT = 680 * MIB  # Reserve space for tar/PAX headers and gzip overhead.
FRAGMENT = 64 * MIB
VERSION = "P3L_STREAMING_EVIDENCE_V1"


def write(path, obj):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def excluded(path):
    relative = path.relative_to(ROOT)
    if any(p in ("__pycache__", ".git", ".venv", "node_modules") for p in relative.parts):
        return True
    if path.name == ".DS_Store" or path.suffix in (".lock", ".tmp", ".pyc"):
        return True
    if path.is_relative_to(DEST) and path.name != "README.md":
        return True
    return False


def walk(folder):
    for directory, folders, files in os.walk(folder, followlinks=False):
        base = Path(directory)
        folders[:] = sorted(name for name in folders if not excluded(base/name))
        for name in folders:
            if (base/name).is_symlink():
                raise RuntimeError(f"symlink directory refused: {(base/name).relative_to(ROOT)}")
        for name in sorted(files):
            path = base/name
            if excluded(path):
                continue
            if path.is_symlink() or not path.is_file():
                raise RuntimeError(f"only regular files allowed: {path.relative_to(ROOT)}")
            yield path


def inventory():
    # Campaign raw artifacts include every NPZ/checkpoint and complete JSONL.
    # Source snapshots are explicit project-owned trees, not external packages.
    for folder in (OUT, ROOT/"src", ROOT/"tools", ROOT/"configs"):
        if folder.exists():
            if folder.is_symlink():
                raise RuntimeError("source/campaign root symlink refused")
            for path in walk(folder):
                stat = path.stat()
                yield dict(path=path.relative_to(ROOT).as_posix(), size=stat.st_size,
                           mtime_ns=stat.st_mtime_ns, mode=stat.st_mode & 0o777)
    for name in ("pyproject.toml", "requirements.txt", "README.md", "AGENTS.md", "LICENSE"):
        path = ROOT/name
        if path.is_file():
            if path.is_symlink():
                raise RuntimeError("symlink source metadata refused")
            stat = path.stat()
            yield dict(path=name, size=stat.st_size, mtime_ns=stat.st_mtime_ns, mode=stat.st_mode & 0o777)


def records(path):
    with path.open() as stream:
        for line in stream:
            yield json.loads(line)


def stable_inventory():
    path = DEST/"INVENTORY.jsonl"
    if path.exists():
        for old, current in itertools.zip_longest(records(path), inventory()):
            if old != current:
                raise RuntimeError("inventory changed since registration; preserve existing parts and investigate before new package revision")
    else:
        temporary = path.with_suffix(".jsonl.tmp")
        with temporary.open("x") as stream:
            for row in inventory():
                stream.write(json.dumps(row, separators=(",", ":")) + "\n")
        temporary.replace(path)
    return path


def members(inventory_path):
    for row in records(inventory_path):
        # Large files are fragmented without interpreting bytes; the journal
        # remains complete and is reconstructed verbatim in path/offset order.
        count = max(1, (row["size"] + FRAGMENT - 1)//FRAGMENT) if row["size"] > PAYLOAD_LIMIT else 1
        for index in range(count):
            offset = index*FRAGMENT if count > 1 else 0
            size = min(FRAGMENT, row["size"]-offset) if count > 1 else row["size"]
            name = f'.p3l_fragments/{row["path"]}/{index:08d}' if count > 1 else row["path"]
            yield dict(**row, archive_member=name, offset=offset, member_size=size,
                       fragment_index=index, fragment_count=count)


def partitions(inventory_path):
    part, estimate = [], 0
    for member in members(inventory_path):
        # PAX path headers need bounded overhead even for long scientific paths.
        estimated = ((member["member_size"]+511)//512)*512 + 8192
        if part and estimate+estimated > PAYLOAD_LIMIT:
            yield part
            part, estimate = [], 0
        part.append(member); estimate += estimated
    if part:
        yield part


def hold_quiet_campaign(stack):
    closure_path = OUT/"CAMPAIGN_CLOSURE.json"
    if not closure_path.is_file():
        raise RuntimeError("parent CAMPAIGN_CLOSURE.json required before freezing source inventory")
    closure = json.loads(closure_path.read_text())
    if not all(isinstance(closure.get(key), str) and closure[key].strip()
               for key in ("scientific_verdict", "closed_at_utc", "probability_reference_status")):
        raise RuntimeError("parent campaign closure fields incomplete")
    snapshot = OUT/"RESOURCE_LEDGER.json"
    if snapshot.stat().st_size > 8*MIB:
        raise RuntimeError("unexpectedly large resource snapshot")
    resource = json.loads(snapshot.read_text())
    if resource.get("active_wall"):
        raise RuntimeError("active resource reservations remain; paid queue must close first")
    worker_path = OUT/"paid_worker.lock"
    if not worker_path.exists():
        raise RuntimeError("existing campaign paid-worker lock required")
    if worker_path.is_symlink():
        raise RuntimeError("worker lock symlink refused")
    worker = stack.enter_context(worker_path.open("r"))
    try:
        fcntl.flock(worker, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        raise RuntimeError("known paid worker remains live")
    for path in sorted(OUT.glob("[ABCDE]/jobs/*/job.lock")):
        if not path.exists():
            continue
        if path.is_symlink():
            raise RuntimeError("lock symlink refused")
        # Probe one at a time rather than keeping thousands of descriptors.
        with path.open("r") as stream:
            try:
                fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                raise RuntimeError(f"known live worker/job remains: {path.relative_to(ROOT)}")
    # Holding the global worker lock excludes new paid work during the stream.
    if json.loads(snapshot.read_text()).get("active_wall"):
        raise RuntimeError("resource reservations changed while acquiring locks")
    return resource


class SliceReader:
    def __init__(self, stream, remaining):
        self.stream, self.remaining = stream, remaining

    def read(self, size=-1):
        size = self.remaining if size < 0 else min(size, self.remaining)
        data = self.stream.read(size)
        self.remaining -= len(data)
        return data


def archive_part(path, entries):
    temporary = path.with_suffix(path.suffix+".tmp")
    if temporary.exists():
        temporary.unlink()  # Only the documented incomplete generated part.
    with temporary.open("xb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=6, mtime=0) as compressed:
            with tarfile.open(fileobj=compressed, mode="w|", format=tarfile.PAX_FORMAT) as archive:
                for row in entries:
                    source = ROOT/row["path"]
                    before = source.stat()
                    if before.st_size != row["size"] or before.st_mtime_ns != row["mtime_ns"]:
                        raise RuntimeError("source changed before streaming")
                    info = tarfile.TarInfo(row["archive_member"])
                    info.size = row["member_size"]; info.mode = row["mode"]
                    info.mtime = row["mtime_ns"]//1_000_000_000
                    with source.open("rb") as stream:
                        stream.seek(row["offset"])
                        archive.addfile(info, SliceReader(stream, row["member_size"]))
                    after = source.stat()
                    if after.st_size != before.st_size or after.st_mtime_ns != before.st_mtime_ns:
                        raise RuntimeError("source changed during streaming; partial part retained")
    if temporary.stat().st_size > PART_LIMIT:
        raise RuntimeError("compressed part exceeded 700MiB; preserve partial and revise partition size")
    temporary.replace(path)


def restore(destination, bundle=DEST):
    """Reconstruct into a new directory; never overwrite source evidence."""
    destination = Path(destination).resolve()
    if destination.exists():
        raise FileExistsError("restore destination must be new")
    bundle = Path(bundle).resolve()
    manifest = json.loads((bundle/"MANIFEST.json").read_text())
    if manifest.get("status") != "COMPLETE_SIZE_PATH_MANIFEST_NO_HASHES":
        raise RuntimeError("complete manifest required")
    inventory_path = bundle/"INVENTORY.jsonl"
    destination.mkdir(parents=True)
    for number, entries in enumerate(partitions(inventory_path), 1):
        expected = manifest["parts"][number-1]
        part = bundle/expected["file"]
        if Path(expected["file"]).name != expected["file"]:
            raise RuntimeError("unsafe part filename in manifest")
        if part.stat().st_size != expected["compressed_bytes"]:
            raise RuntimeError("archive part size differs from manifest")
        with tarfile.open(part, mode="r|gz") as archive:
            for row in entries:
                member = archive.next()
                if member is None or member.name != row["archive_member"] or member.size != row["member_size"] or not member.isfile():
                    raise RuntimeError("archive path/size/type differs from inventory")
                logical = Path(row["path"])
                if logical.is_absolute() or ".." in logical.parts:
                    raise RuntimeError("unsafe inventory path")
                target = destination/logical
                target.parent.mkdir(parents=True, exist_ok=True)
                mode = "xb" if row["offset"] == 0 else "ab"
                if row["offset"] and target.stat().st_size != row["offset"]:
                    raise RuntimeError("fragment offset differs from restored size")
                with target.open(mode) as output, archive.extractfile(member) as source:
                    while True:
                        data = source.read(MIB)
                        if not data:
                            break
                        output.write(data)
            if archive.next() is not None:
                raise RuntimeError("unexpected extra archive member")
    for row in records(inventory_path):
        path = destination/row["path"]
        if path.stat().st_size != row["size"]:
            raise RuntimeError("restored logical file size differs")
        os.chmod(path, row["mode"])
        os.utime(path, ns=(row["mtime_ns"], row["mtime_ns"]))
    print(json.dumps(dict(status="RESTORED_PATH_SIZE_CHECKED_NO_HASHES", destination=str(destination))))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--finalize", action="store_true", help="explicitly authorize archive creation after paid work closes")
    parser.add_argument("--restore", type=Path, metavar="NEW_DIRECTORY", help="restore completed archives into a new directory")
    parser.add_argument("--bundle-dir", type=Path, help="downloaded manifest/inventory/parts directory; restore only")
    args = parser.parse_args()
    if args.restore is not None:
        if args.finalize:
            raise ValueError("choose restore or finalize")
        restore(args.restore, args.bundle_dir or DEST); return
    if args.bundle_dir is not None:
        raise ValueError("--bundle-dir is only available with --restore")
    if not args.finalize:
        total = count = largest = 0
        for row in inventory():
            count += 1; total += row["size"]; largest = max(largest, row["size"])
        print(json.dumps(dict(status="METADATA_ONLY_PLAN_NO_ARCHIVES_CREATED", file_count=count,
            logical_bytes=total, largest_file_bytes=largest,
            minimum_payload_parts=(total+PAYLOAD_LIMIT-1)//PAYLOAD_LIMIT,
            additional_disk_upper_projection_bytes=total+max(1, (total+PAYLOAD_LIMIT-1)//PAYLOAD_LIMIT)*20*MIB,
            caveat="upper projection includes coarse archive overhead; actual gzip size depends on bytes; no content opened")))
        return
    DEST.mkdir(parents=True, exist_ok=True)
    with ExitStack() as stack:
        resource = hold_quiet_campaign(stack)
        inventory_path = stable_inventory()
        manifest_path = DEST/"MANIFEST.json"
        manifest = (json.loads(manifest_path.read_text()) if manifest_path.exists() else
                    dict(format=VERSION, status="IN_PROGRESS", inventory="INVENTORY.jsonl", parts=[],
                         source_root="worktree-relative paths", resource_snapshot_sequence=resource.get("sequence"),
                         hashes_used=False, complete_ledger_included=True,
                         exclusions="transient .lock/.tmp/.pyc/.DS_Store, __pycache__, .git/.venv/node_modules, generated packaging files except README; logs retained"))
        if manifest.get("format") != VERSION:
            raise RuntimeError("incompatible existing manifest")
        for index, entries in enumerate(partitions(inventory_path)):
            name = f"P3L_evidence_part{index+1:04d}.tar.gz"
            path = DEST/name
            logical_bytes = sum(e["member_size"] for e in entries)
            old = manifest["parts"][index] if index < len(manifest["parts"]) else None
            if old is not None:
                if old["file"] != name or old["logical_bytes"] != logical_bytes or not path.exists() or path.stat().st_size != old["compressed_bytes"]:
                    raise RuntimeError("completed part path/size differs; no silent replacement")
                continue
            if path.exists():
                raise RuntimeError("unregistered completed part exists; inspect before resuming")
            archive_part(path, entries)
            manifest["parts"].append(dict(file=name, compressed_bytes=path.stat().st_size,
                logical_bytes=logical_bytes, member_count=len(entries),
                first_member=entries[0]["archive_member"], last_member=entries[-1]["archive_member"]))
            write(manifest_path, manifest)
            print(json.dumps(manifest["parts"][-1]), flush=True)
        # Recheck every path/size/mtime and detect newly created source files.
        stable_inventory()
        manifest.update(status="COMPLETE_SIZE_PATH_MANIFEST_NO_HASHES",
                        completed_unix=time.time(), total_compressed_bytes=sum(p["compressed_bytes"] for p in manifest["parts"]),
                        total_logical_bytes=sum(p["logical_bytes"] for p in manifest["parts"]))
        write(manifest_path, manifest)
        print(json.dumps(dict(status=manifest["status"], parts=len(manifest["parts"]),
                              total_compressed_bytes=manifest["total_compressed_bytes"])))


if __name__ == "__main__":
    main()
