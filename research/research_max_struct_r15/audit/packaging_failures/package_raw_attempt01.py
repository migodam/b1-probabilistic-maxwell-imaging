"""Archive only terminal R15 main264/fork56 raw evidence and public sources.

No scientific imports, array loading, scoring, publication, hashes or deletion.
`plan` writes manifests only; `create` explicitly creates the three ZIP groups.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import re
import sys
import time
import zipfile

TARGET_BYTES = 600_000_000
HARD_BYTES = 650_000_000
CAMPAIGNS = (
    ('main264_raw', 'twofold_exposed24_v1', 'FROZEN_MANIFEST.json', 264),
    ('fork56_raw', 'fullstate_interventions_v1', 'FROZEN_TASKS.json', 56),
)
ROOT_DOCUMENTS = ('EXPERIMENT_CONTRACT.md', 'THEORY_DERIVATIONS.md', 'ORIGINAL_REQUEST_ZH.txt')
TEXT_SUFFIXES = {'.py', '.json', '.jsonl', '.csv', '.log', '.md', '.txt', '.yaml', '.yml',
                 '.toml', '.ini', '.cfg', '.sh', '.ps1', '.html', '.css', '.xml', '.tsv'}
SKIP_PARTS = {'.git', '__pycache__', '.venv', 'venv', 'node_modules', '.pytest_cache',
              'sealed_truth', 'private_operational', 'privateoperational', 'credentials', 'tokens'}
PRIVATE_NAMES = re.compile(r'credential|secret|token|password|passwd|private[_-]?operational', re.I)
SECRET_LITERAL = re.compile(
    rb'(?:gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{30,}|'
    rb'(?<![A-Za-z0-9])sk-[A-Za-z0-9_-]{30,}|Bearer\s+[A-Za-z0-9._~-]{24,}|'
    rb'(?:api[_-]?key|access[_-]?token|refresh[_-]?token|client[_-]?secret|password)'
    rb'\s*["\x27]?\s*[:=]\s*["\x27][A-Za-z0-9_./+=-]{16,})', re.I)


def now():
    return datetime.now(timezone.utc).isoformat()


def write_new(path, value):
    with path.open('x', encoding='utf-8') as handle:
        json.dump(value, handle, indent=2, ensure_ascii=False, allow_nan=False)
        handle.write('\n')


def safe_component(value):
    if not isinstance(value, str) or not value or value in ('.', '..') or any(c in value for c in '/\\:'):
        raise ValueError('Unsafe task/path component')
    return value


def exclusion(path, root):
    relative = path.relative_to(root)
    for component in relative.parts:
        lower = component.lower()
        if lower in SKIP_PARTS or PRIVATE_NAMES.search(component):
            return 'PRIVATE_OR_EXCLUDED_PATH_COMPONENT'
        if lower == '.env' or lower.startswith('.env.') or lower.endswith(('.pem', '.p12', '.pfx', '.key')):
            return 'CREDENTIAL_FILE_NAME'
    if path.name in ('.DS_Store', 'TASK_RESUME.md', 'CHAT_READBACK.json'):
        return 'OPERATIONAL_OR_GENERATED_FILE'
    if path.is_symlink():
        return 'SYMLINK_NOT_FOLLOWED'
    return None


def walk_files(directory, root, exclusions):
    if not directory.is_dir() or directory.is_symlink():
        raise FileNotFoundError('Required admitted directory missing or symlink: ' + str(directory))
    for base, dirs, names in os.walk(directory, followlinks=False):
        here = Path(base)
        keep = []
        for name in sorted(dirs):
            path = here / name
            reason = exclusion(path, root)
            if reason:
                exclusions.append(dict(path=path.relative_to(root).as_posix(), reason=reason, contents_opened=False))
            else:
                keep.append(name)
        dirs[:] = keep
        for name in sorted(names):
            path = here / name
            reason = exclusion(path, root)
            if reason:
                exclusions.append(dict(path=path.relative_to(root).as_posix(), reason=reason, contents_opened=False))
                continue
            if not path.is_file() or not path.resolve().is_relative_to(root):
                raise ValueError('Nonregular/outside-root source: ' + str(path))
            yield path


def scan_public_text(path):
    """Block potential literal leakage, without logging any matched content."""
    if path.suffix.lower() not in TEXT_SUFFIXES and path.name != '.gitignore':
        return
    tail = b''
    with path.open('rb') as handle:
        while block := handle.read(1024 * 1024):
            combined = tail + block
            if SECRET_LITERAL.search(combined):
                raise ValueError('Possible credential literal; archive blocked, review path: ' + str(path))
            tail = combined[-2048:]


def read_json(path, root):
    resolved = path.resolve()
    if not resolved.is_relative_to(root) or exclusion(path, root) or exclusion(resolved, root):
        raise ValueError('Refusing excluded JSON')
    return json.loads(path.read_text(encoding='utf-8'))


def terminal_group(root, name, campaign, frozen_name, expected_count, exclusions):
    queue = root / 'queues' / campaign
    result = root / 'RESULTS' / campaign
    terminal = read_json(queue / 'RECEIPT.json', root)
    frozen = read_json(queue / frozen_name, root)
    tasks = frozen['tasks']
    ids = [safe_component(t['task_id']) for t in tasks]
    jobs = terminal.get('jobs', [])
    if frozen.get('campaign') != campaign or len(ids) != expected_count or len(set(ids)) != expected_count:
        raise ValueError('Frozen task inventory differs from declared R15 group: ' + campaign)
    if (terminal.get('campaign') != campaign or terminal.get('requested_tasks') != expected_count or
            terminal.get('status') not in ('COMPLETED', 'COMPLETED_WITH_ERRORS') or
            terminal.get('running') or terminal.get('pending_tasks', 0) or terminal.get('pending_cases', 0) or
            any(lane.get('remaining_arms', 0) for lane in terminal.get('lanes', []))):
        raise RuntimeError('Require terminal queue with no pending/running work: ' + campaign)
    if len(jobs) != expected_count or {j['task_id'] for j in jobs} != set(ids):
        raise ValueError('Terminal receipt does not account for every frozen task: ' + campaign)
    missing_results = []
    for job in jobs:
        ident = safe_component(job['task_id'])
        if job.get('status') not in ('FINISHED', 'FAILED', 'LAUNCH_FAILED', 'PREFLIGHT_FAILED'):
            raise ValueError('Nonterminal/unknown job status: ' + ident)
        record = read_json(queue / 'TERMINAL_JOBS' / (ident + '.json'), root)
        if record.get('task_id') != ident or record.get('status') != job['status']:
            raise ValueError('Per-task terminal record conflicts with final queue receipt: ' + ident)
        folder = result / ident
        if not folder.is_dir() or not (folder / 'RECEIPT.json').is_file():
            if job['status'] == 'FINISHED':
                raise FileNotFoundError('FINISHED child lacks original result/receipt: ' + ident)
            missing_results.append(dict(task_id=ident, queue_status=job['status'],
                                        result_directory_exists=folder.is_dir(),
                                        reason='Failure/preflight/launch retained as terminal record and log; no result fabricated'))
        if not (queue / 'logs' / (ident + '.log')).is_file():
            raise FileNotFoundError('Original terminal task log missing: ' + ident)
    paths = [*walk_files(queue, root, exclusions)]
    if result.is_dir():
        paths.extend(walk_files(result, root, exclusions))
    elif any(j['status'] == 'FINISHED' for j in jobs):
        raise FileNotFoundError('Successful group RESULTS directory missing')
    metadata = dict(name=name, kind='FULL_TERMINAL_RAW', campaign=campaign,
                    frozen_file=f'queues/{campaign}/{frozen_name}', queue_status=terminal['status'],
                    expected_tasks=expected_count, terminal_jobs=len(jobs),
                    failed_jobs=sum(j['status'] != 'FINISHED' for j in jobs),
                    missing_failed_result_artifacts=missing_results,
                    completeness='Every admitted file under original RESULTS and queues, including all states/candidates/caches/receipts/failures/logs/status snapshots; no material-only export',
                    queue_cost_snapshot={k: terminal.get(k) for k in ('controller_CPU_seconds', 'finished_child_CPU_seconds',
                        'measured_finished_child_CPU_seconds', 'missing_external_CPU_tasks', 'queue_wall_seconds')})
    return metadata, paths


def inventory(root, paths):
    rows = []
    for path in sorted(set(paths)):
        if exclusion(path, root) or not path.resolve().is_relative_to(root):
            raise ValueError('Unsafe/excluded inventory source')
        scan_public_text(path)
        stat = path.stat()
        text = path.suffix.lower() in TEXT_SUFFIXES or path.name == '.gitignore'
        # A conservative raw-deflate upper bound plus local/central/Zip64/name overhead.
        bound = stat.st_size + (5 * math.ceil(stat.st_size / 16383) + 64 if text else 0)
        charge = bound + 1024 + 2 * len(path.relative_to(root).as_posix().encode('utf-8'))
        rows.append(dict(path=path.relative_to(root).as_posix(), bytes=stat.st_size, mtime_ns=stat.st_mtime_ns,
                         compression='ZIP_DEFLATED' if text else 'ZIP_STORED', compression_level=1 if text else None,
                         conservative_zip_bytes=charge, numerical_transform='NONE'))
    return rows


def partition(rows, group):
    parts, members, charged = [], [], 4096
    for index, row in enumerate(rows):
        charge = row['conservative_zip_bytes']
        if charge + 4096 > HARD_BYTES:
            raise ValueError('Single file cannot fit 650MB lossless volume; no implicit truncation/fragmentation: ' + row['path'])
        if members and charged + charge > TARGET_BYTES:
            parts.append(dict(archive=f'r15_{group}_part{len(parts)+1:03d}.zip', file_indices=members,
                              conservative_zip_bytes=charged, exceeds_target=charged > TARGET_BYTES))
            members, charged = [], 4096
        members.append(index)
        charged += charge
    if members:
        parts.append(dict(archive=f'r15_{group}_part{len(parts)+1:03d}.zip', file_indices=members,
                          conservative_zip_bytes=charged, exceeds_target=charged > TARGET_BYTES))
    return parts


def make_plan(root):
    groups, exclusions = [], []
    for specification in CAMPAIGNS:
        metadata, paths = terminal_group(root, *specification, exclusions)
        metadata['files'] = inventory(root, paths)
        groups.append(metadata)
    paths = []
    for name in ('src', 'runtime', 'configs', 'data_public'):
        paths.extend(walk_files(root / name, root, exclusions))
    # The only admitted audit tree is the original operator-only rank/basis cache.
    rank = root / 'audit/rank_all24'
    for domain in ('real', 'complex'):
        if not (rank / domain / 'RECEIPT.json').is_file() or not list((rank / domain).glob('operator_group_*/SPACES.npz')):
            raise FileNotFoundError('Original rank/basis cache incomplete: ' + domain)
        rank_receipt = read_json(rank / domain / 'RECEIPT.json', root)
        if rank_receipt.get('status') != 'COMPLETE' or rank_receipt.get('completed') != 24:
            raise ValueError('Expected original completed24 rank cache: ' + domain)
    paths.extend(walk_files(rank, root, exclusions))
    for name in ROOT_DOCUMENTS:
        path = root / name
        if not path.is_file():
            raise FileNotFoundError('Original contract/theory/request required: ' + name)
        paths.append(path)
    if (root / '.gitignore').is_file():
        paths.append(root / '.gitignore')
    groups.append(dict(name='source_public', kind='SOURCE_RUNTIME_CONFIG_PUBLIC_AND_RANK_BASIS',
                       files=inventory(root, paths),
                       completeness='src/runtime/configs/data_public plus original audit/rank_all24 real+complex basis arrays/metadata and explicit contract/theory/original-request documents',
                       excluded_scope='No other audit, analysis, release, sealed truth, credentials, private operations or root startup/source ZIPs'))
    for group in groups:
        group['parts'] = partition(group['files'], group['name'])
        group['file_count'] = len(group['files'])
        group['original_bytes'] = sum(row['bytes'] for row in group['files'])
    return dict(schema='R15_FULL_RAW_PACKAGE_V1', created_utc=now(), source_root=str(root),
                target_volume_bytes=TARGET_BYTES, hard_max_volume_bytes=HARD_BYTES, groups=groups,
                total_files=sum(g['file_count'] for g in groups), original_bytes=sum(g['original_bytes'] for g in groups),
                excluded_paths_not_opened=exclusions,
                excluded_entire_trees=['sealed_truth', 'analysis', 'audit except rank_all24', 'archives', 'release', 'private_operational'],
                root_startup_or_source_archives_admitted=False, hash_operations='NONE',
                numerical_arrays_loaded=0, scientific_imports=0, physics=0, optimization=0, GT_access='NONE',
                originals_retained=True, upload_or_publication='NOT_PERFORMED',
                compression='NPZ and other binary ZIP_STORED; text JSON/CSV/log/source lossless DEFLATE level1',
                oversized_file_policy='No fragmentation, deletion, or truncation: reject individual conservative ZIP bound >650MB; >600MB singleton allowed only <=650MB')


def write_plan_files(output, plan):
    write_new(output / 'PACKAGE_PLAN.json', plan)
    for group in plan['groups']:
        write_new(output / (group['name'] + '_MANIFEST.json'), group)
    with (output / 'FILE_INVENTORY.csv').open('x', newline='', encoding='utf-8') as handle:
        writer = csv.writer(handle)
        writer.writerow(['group', 'archive', 'original_path', 'bytes', 'mtime_ns', 'compression', 'compression_level'])
        for group in plan['groups']:
            for part in group['parts']:
                for index in part['file_indices']:
                    row = group['files'][index]
                    writer.writerow([group['name'], part['archive'], row['path'], row['bytes'], row['mtime_ns'], row['compression'], row['compression_level']])
    (output / 'RESTORE.md').write_text(
        '# R15 original raw files\n\n'
        'All ZIP members retain their R15-root-relative original paths. Extract every listed part into one NEW empty R15 directory. '
        'The main264_raw, fork56_raw and source_public groups are separate; do not substitute earlier remote-input/source ZIPs for them. '
        'Compare exact member paths, uncompressed sizes and compression with PACKAGE_PLAN.json / FILE_INVENTORY.csv. '
        'Normal ZIP CRC is part of lossless ZIP handling; no SHA is used. No arrays are converted or loaded. '
        'A PLANNED_NO_ARCHIVES receipt is not an archived release. A PACKAGE_ERROR means incomplete output; preserve it and choose a new destination for any later attempt. '
        'Originals are never removed. Analysis and old-endpoint GT-derived arrays are deliberately outside these packages.\n', encoding='utf-8')


def create_archives(root, output, plan, completed):
    for group in plan['groups']:
        for part in group['parts']:
            start_cpu, start_wall = time.process_time(), time.perf_counter()
            archive = output / part['archive']
            entries = [group['files'][index] for index in part['file_indices']]
            with zipfile.ZipFile(archive, 'x', compression=zipfile.ZIP_STORED, allowZip64=True) as bundle:
                for row in entries:
                    path = root / row['path']
                    if exclusion(path, root) or not path.resolve().is_relative_to(root):
                        raise ValueError('Excluded/unsafe original source')
                    before = path.stat()
                    if (before.st_size, before.st_mtime_ns) != (row['bytes'], row['mtime_ns']):
                        raise RuntimeError('Source changed since inventory: ' + row['path'])
                    bundle.write(path, row['path'], compress_type=getattr(zipfile, row['compression']),
                                 compresslevel=row['compression_level'])
                    after = path.stat()
                    if (after.st_size, after.st_mtime_ns) != (row['bytes'], row['mtime_ns']):
                        raise RuntimeError('Source changed while archiving: ' + row['path'])
            actual = archive.stat().st_size
            if actual > HARD_BYTES:
                raise RuntimeError('Actual archive exceeds hard 650MB limit; output not publishable: ' + archive.name)
            with zipfile.ZipFile(archive, 'r') as bundle:
                infos = bundle.infolist()
                if [item.filename for item in infos] != [row['path'] for row in entries]:
                    raise ValueError('ZIP directory differs from manifest')
                if any(item.file_size != row['bytes'] or item.compress_type != getattr(zipfile, row['compression'])
                       for item, row in zip(infos, entries)):
                    raise ValueError('ZIP sizes/compression differ from manifest')
            completed.append(dict(group=group['name'], archive=archive.name, actual_bytes=actual,
                                  files=len(entries), original_bytes=sum(row['bytes'] for row in entries),
                                  member_compressed_bytes=sum(item.compress_size for item in infos),
                                  CPU_seconds=time.process_time()-start_cpu, wall_seconds=time.perf_counter()-start_wall,
                                  status='WRITTEN_SIZE_AND_CENTRAL_DIRECTORY_CHECKED'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('plan', 'create'))
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path, help='New directory only; defaults to archives/r15_raw_plan_v1 or archives/r15_raw_v1')
    args = parser.parse_args()
    root = args.root.resolve()
    output = (args.output or root / 'archives' / ('r15_raw_v1' if args.command == 'create' else 'r15_raw_plan_v1')).resolve()
    if not root.is_dir():
        raise ValueError('Source root does not exist')
    if output.is_relative_to(root) and exclusion(output, root):
        raise ValueError('Invalid source/output path')
    if not output.is_relative_to(root / 'archives'):
        raise ValueError('Output must remain under this R15 root/archives; raw archives cannot be startup/source packages')
    if sys.platform == 'darwin' and not output.is_relative_to('/Volumes'):
        raise ValueError('Use the external /Volumes disk for large archives')
    output.mkdir(parents=True, exist_ok=False)
    start_cpu, start_wall = time.process_time(), time.perf_counter()
    completed, plan = [], None
    try:
        plan = make_plan(root)
        write_plan_files(output, plan)
        if args.command == 'create':
            create_archives(root, output, plan, completed)
            # Recheck complete source membership/stat inventory and terminal/frozen metadata.
            fresh = make_plan(root)
            for key in ('groups', 'excluded_paths_not_opened'):
                if fresh[key] != plan[key]:
                    raise RuntimeError('Terminal/source inventory changed during packaging; retain failed attempt')
        receipt = dict(status='ARCHIVED_ORIGINAL_BYTES_NOT_PUBLISHED' if args.command == 'create' else 'PLANNED_NO_ARCHIVES',
                       created_utc=now(), command=args.command, source_root=str(root), output=str(output),
                       groups=[g['name'] for g in plan['groups']], files=plan['total_files'], original_bytes=plan['original_bytes'],
                       archive_count=len(completed), actual_archive_bytes=sum(x['actual_bytes'] for x in completed),
                       parts=completed, CPU_seconds=time.process_time()-start_cpu,
                       whole_process_CPU_at_receipt_seconds=time.process_time(), wall_seconds=time.perf_counter()-start_wall,
                       failed_attempts_in_this_directory=0, originals_retained=True, deletion=0, hash_operations='NONE',
                       physics=0, optimization=0, GT_access='NONE', publication='NOT_PERFORMED',
                       verification='Source membership/size/mtime and terminal IDs, ZIP member paths/sizes/compression; no SHA or numerical transformation')
        write_new(output / 'PACKAGE_RECEIPT.json', receipt)
        print(json.dumps({k: receipt[k] for k in ('status', 'files', 'original_bytes', 'archive_count', 'actual_archive_bytes', 'CPU_seconds', 'wall_seconds')}))
        return 0
    except BaseException as exc:
        write_new(output / 'PACKAGE_ERROR.json', dict(status='FAILED_INCOMPLETE_ALL_ATTEMPTS_RETAINED', error=repr(exc),
                    completed_parts=completed, plan_written=(output / 'PACKAGE_PLAN.json').is_file(),
                    retained_zip_files=[dict(path=p.name, bytes=p.stat().st_size) for p in output.glob('*.zip')],
                    CPU_seconds=time.process_time()-start_cpu, whole_process_CPU_at_error_seconds=time.process_time(),
                    wall_seconds=time.perf_counter()-start_wall, originals_retained=True, hash_operations='NONE', GT_access='NONE'))
        raise


if __name__ == '__main__':
    raise SystemExit(main())
