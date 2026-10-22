"""Stream frozen P3L-R2 evidence without importing campaign/physics modules.

Run only after other jobs stop. Under run_p3r2_accounted.py, pass its exact
--job as --accounted-job. Output must be a NEW directory outside the campaign.
Failed bundles and incomplete parts are retained; reruns need a new directory.
"""
from contextlib import ExitStack
from pathlib import Path
import argparse
import fcntl
import gzip
import itertools
import json
import os
import resource
import stat
import subprocess
import sys
import tarfile
import time

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'research/p3l_r2_risk_posterior_exploration'
MIB = 1024 * 1024
VERSION = 'P3LR2_STREAM_SPLIT_V1'
LEDGERS = {OUT / 'RESOURCE_LEDGER.jsonl', OUT / 'RESOURCE_LEDGER.json'}
RESTORE = (
    'Require receipt.json status COMPLETE and all listed parts with matching sizes. '
    'Concatenate evidence.tar.gz.partNNN in numeric order, then gzip-decompress and '
    'extract the resulting tar into a NEW empty directory (e.g. cat '
    'evidence.tar.gz.part[0-9][0-9][0-9]* | tar -xzf - -C NEW_DIRECTORY). '
    'These are byte splits of ONE gzip stream; individual parts are not archives. '
    'Symlinks are provenance-only in INVENTORY.json and are not restored. '
    'Compare restored regular-file sizes with INVENTORY.json. The final campaign '
    'RESOURCE_LEDGER.jsonl and RESOURCE_LEDGER.json must be supplied separately '
    'after the accounted packaging wrapper closes. No digest verification is used.'
)


def dump_new(path, obj):
    with path.open('x') as stream:
        json.dump(obj, stream, indent=2, allow_nan=False)
        stream.write('\n')


def exclusion(path):
    if path in LEDGERS:
        return 'final ledger supplied separately after packaging closes'
    if any(x in {'.git', '__pycache__', '.venv', 'node_modules'} for x in path.relative_to(ROOT).parts):
        return 'transient/cache directory'
    if path.name == '.DS_Store' or path.name.endswith(('~', '.lock', '.tmp', '.temp', '.pyc', '.partial')):
        return 'lock or temporary file'
    if path.name.startswith('.#') or path.name.endswith(('.tar', '.tar.gz', '.tgz', '.zip', '.7z')) or '.tar.gz.part' in path.name:
        return 'generated output archive'
    return None


def tree(path):
    """Do not resolve or descend into symlinks; never enumerate an old study."""
    reason = exclusion(path)
    info = path.lstat()
    row = dict(path=path.relative_to(ROOT).as_posix(), size=info.st_size,
               mtime_ns=info.st_mtime_ns, mode=stat.S_IMODE(info.st_mode))
    if stat.S_ISLNK(info.st_mode):
        yield dict(row, type='symlink', link_target=os.readlink(path), archived=False,
                   exclusion=reason or 'provenance only; target never followed')
    elif reason:
        yield dict(row, type='excluded', archived=False, exclusion=reason)
    elif stat.S_ISDIR(info.st_mode):
        yield dict(row, type='directory', archived=True)
        with os.scandir(path) as children:
            for child in children:
                yield from tree(Path(child.path))
    elif stat.S_ISREG(info.st_mode):
        yield dict(row, type='file', archived=True)
    else:
        raise RuntimeError(f'Unsupported special file: {row["path"]}')


def entries():
    for path in (OUT, ROOT/'src/p3r2', ROOT/'src/p3r2_mechanism', ROOT/'src/p1_material',
                 ROOT/'research/p1_material_target_20261018/inputs/geometry.json'):
        if not path.exists() or path.is_symlink():
            raise RuntimeError(f'Required real campaign/source/reference path missing: {path}')
        yield from tree(path)
    with os.scandir(ROOT/'tools') as children:
        for child in children:
            if 'p3r2' in child.name and child.name.endswith('.py'):
                yield from tree(Path(child.path))
    for name in ('pyproject.toml', 'requirements.txt', 'README.md', 'AGENTS.md', 'LICENSE'):
        path = ROOT/name
        if path.exists() or path.is_symlink():
            yield from tree(path)


def rows(path):
    with path.open() as stream:
        for line in stream:
            yield json.loads(line)


def latest_ledger():
    # Book._latest uses the last complete JSONL row, not a possibly stale snapshot.
    with (OUT/'RESOURCE_LEDGER.jsonl').open('rb') as stream:
        stream.seek(0, 2)
        size = stream.tell()
        stream.seek(max(0, size-8*MIB))
        lines = stream.read(8*MIB).splitlines()
    if not lines or not lines[-1].strip():
        raise RuntimeError('Missing authoritative ledger row')
    row = json.loads(lines[-1])
    if not isinstance(row.get('after', {}).get('active_wall'), dict):
        raise RuntimeError('Invalid authoritative ledger active_wall')
    return row


def hold_quiet(stack, accounted_job):
    # All existing Book workers lock budget.lock for reservations. Holding it
    # freezes admissions and ledger updates until packaging exits. No imports.
    for name in ('paid_worker.lock', 'budget.lock'):
        path = OUT/name
        if path.is_symlink():
            raise RuntimeError(f'Symlink lock refused: {name}')
        handle = stack.enter_context(path.open('r+'))
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError(f'Campaign lock busy: {name}')
    row = latest_ledger()
    expected = {accounted_job} if accounted_job else set()
    if set(row['after']['active_wall']) != expected:
        raise RuntimeError('Other active jobs, or accounted wrapper absent: ' +
                           repr(row['after']['active_wall']))
    if accounted_job:
        if row.get('job') != accounted_job or row.get('event') != 'JOB_START' or not row.get('analysis_only') or row.get('stage') != 'D':
            raise RuntimeError('Accounted job must be the current stage D analysis-only JOB_START')
    return row


class SplitWriter:
    """Bound compressed parts exactly, independent of input sizes/compression."""
    def __init__(self, destination, limit):
        self.destination, self.limit = destination, limit
        self.stream = None
        self.size = 0
        self.parts = []
        self.total = 0

    def write(self, data):
        view = memoryview(data)
        length = len(view)
        while view:
            if self.stream is None:
                number = len(self.parts)+1
                path = self.destination/f'evidence.tar.gz.part{number:03d}'
                self.stream = path.open('xb')
                self.parts.append(dict(file=path.name, size=0))
                self.size = 0
            count = min(len(view), self.limit-self.size)
            self.stream.write(view[:count])
            self.size += count
            self.total += count
            self.parts[-1]['size'] = self.size
            view = view[count:]
            if self.size == self.limit:
                self.flush_close()
        return length

    def flush(self):
        if self.stream is not None:
            self.stream.flush()

    def flush_close(self):
        if self.stream is not None:
            self.stream.flush()
            os.fsync(self.stream.fileno())
            self.stream.close()
            self.stream = None


def check_row(row):
    path = ROOT/row['path']
    current = path.lstat()
    if (current.st_size, current.st_mtime_ns, stat.S_IMODE(current.st_mode)) != (row['size'], row['mtime_ns'], row['mode']):
        raise RuntimeError(f'Source changed: {row["path"]}')
    if row['type'] == 'file' and not stat.S_ISREG(current.st_mode):
        raise RuntimeError(f'Source is no longer regular: {row["path"]}')
    return path


def inventory_json(destination, journal, metadata):
    with (destination/'INVENTORY.json').open('x') as stream:
        stream.write('{\n"metadata": ')
        json.dump(metadata, stream, allow_nan=False)
        stream.write(',\n"entries": [\n')
        first = True
        for row in rows(journal):
            if not first:
                stream.write(',\n')
            json.dump(row, stream, allow_nan=False)
            first = False
        stream.write('\n]}\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-directory', required=True, type=Path, help='NEW directory outside OUT')
    parser.add_argument('--part-bytes', type=int, default=1800*MIB, help='compressed part limit, strictly below 2GB')
    parser.add_argument('--accounted-job', help='exact active stage D analysis-only wrapper job name')
    args = parser.parse_args()
    destination = args.output_directory.expanduser().resolve()
    if destination == OUT.resolve() or OUT.resolve() in destination.parents:
        raise ValueError('Output directory must be outside OUT')
    if not 1 <= args.part_bytes < 2_000_000_000:
        raise ValueError('--part-bytes must be positive and below 2,000,000,000')
    if any(destination == root or root in destination.parents
           for root in (ROOT/'src/p3r2', ROOT/'src/p1_material', ROOT/'tools')):
        raise ValueError('Output directory must be outside archived source trees')
    if destination.exists():
        raise FileExistsError('Preserve prior evidence; choose a NEW output directory')
    start = time.perf_counter()
    writer = None
    receipt = dict(format=VERSION, status='FAILED', started_unix=time.time(),
                   output_directory=str(destination), part_bytes=args.part_bytes,
                   accounted_job=args.accounted_job, restore_instructions=RESTORE,
                   final_ledger_included=False, hashes_used=False)
    # Lock/eligibility failures do not create a misleading empty evidence bundle.
    with ExitStack() as stack:
        ledger = hold_quiet(stack, args.accounted_job)
        receipt['ledger_sequence_at_start'] = ledger['sequence']
        revision = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=ROOT,
                                  check=True, capture_output=True, text=True).stdout.strip()
        receipt['source_revision'] = revision
        destination.mkdir(parents=True, exist_ok=False)
        journal = destination/'INVENTORY_ENTRIES.jsonl'
        try:
            count = logical = 0
            with journal.open('x') as stream:
                for row in entries():
                    stream.write(json.dumps(row, separators=(',', ':'))+'\n')
                    count += 1
                    logical += row['size'] if row['type'] == 'file' else 0
            receipt.update(inventory_entries=count, logical_file_bytes=logical,
                           version_references=['src/p1_material/prior.py',
                               'src/p1_material/target.py',
                               'research/p1_material_target_20261018/inputs/geometry.json',
                               str((OUT/'configs/TARGET.json').relative_to(ROOT))],
                           source_snapshot='Archived working-tree bytes including uncommitted source; revision is provenance, not a clean-tree claim')
            writer = SplitWriter(destination, args.part_bytes)
            with gzip.GzipFile(filename='', fileobj=writer, mode='wb', compresslevel=1, mtime=0) as zipped:
                with tarfile.open(fileobj=zipped, mode='w|', format=tarfile.PAX_FORMAT, dereference=False) as archive:
                    archive.copybufsize = MIB
                    for row in rows(journal):
                        if not row['archived']:
                            continue
                        path = check_row(row)
                        info = tarfile.TarInfo(row['path'])
                        info.mode = row['mode']
                        info.mtime = row['mtime_ns']//1_000_000_000
                        if row['type'] == 'directory':
                            info.type = tarfile.DIRTYPE
                            archive.addfile(info)
                        else:
                            info.size = row['size']
                            # O_NOFOLLOW closes the lstat/open symlink race.
                            fd = os.open(path, os.O_RDONLY | os.O_NOFOLLOW)
                            with os.fdopen(fd, 'rb') as source:
                                opened = os.fstat(source.fileno())
                                if (opened.st_size, opened.st_mtime_ns) != (row['size'], row['mtime_ns']):
                                    raise RuntimeError(f'Source changed at open: {row["path"]}')
                                archive.addfile(info, source)
                            check_row(row)
                        # TarFile retains member metadata even in streaming mode.
                        archive.members.clear()
            writer.flush_close()
            # Metadata-only final comparison catches additions/removals as well
            # as changed file sizes/times. Traversal-order changes fail safely.
            for before, after in itertools.zip_longest(rows(journal), entries()):
                if before != after:
                    raise RuntimeError('Source inventory changed during packaging; preserve failed bundle')
            if latest_ledger()['sequence'] != ledger['sequence']:
                raise RuntimeError('Ledger changed despite held budget lock')
            receipt['status'] = 'COMPLETE'
        except BaseException as error:
            receipt['error'] = repr(error)
            raise
        finally:
            if writer is not None:
                writer.flush_close()
            receipt.update(parts=writer.parts if writer else [],
                           compressed_bytes=writer.total if writer else 0)
            try:
                if journal.exists():
                    inventory_json(destination, journal, receipt)
            except BaseException as error:
                receipt.update(status='FAILED', finalization_error=repr(error))
                raise
            finally:
                receipt.update(wall_seconds=time.perf_counter()-start,
                               completed_unix=time.time(),
                               peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * (1 if sys.platform == 'darwin' else 1024))
                dump_new(destination/'receipt.json', receipt)
    print(json.dumps(dict(status=receipt['status'], output_directory=str(destination),
                          parts=len(receipt['parts']), compressed_bytes=receipt['compressed_bytes'])))


if __name__ == '__main__':
    main()
