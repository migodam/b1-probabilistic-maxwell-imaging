#!/usr/bin/env python3
"""Create a NEW delivery-only post-pack supplement after verified upload completion.

Uses tracked Git changes plus narrow saved-replay/delivery metadata scopes.
Never runs science, modifies evidence, follows symlinks, hashes, uploads, or
opens an existing ZIP. ZIP format CRC fields are incidental format metadata.
"""
import argparse
import datetime
import json
import os
import re
from pathlib import Path
import shutil
import stat
import subprocess
import sys
import time
import zipfile

ROOT = Path(__file__).resolve().parents[1]
REL = Path('research/p3l_r4_conditional_transport')
CAMPAIGN = ROOT / REL
PACK = CAMPAIGN / 'packaging'
ARCHIVED_SOURCE = '60066394392fe04fb6a5b30371ffa34c978ade9d'
COMPLETE = 'ALL_PARTS_UPLOADED_NAME_SIZE_VERIFIED'
PER_FILE = 16 * 1024 * 1024
TOTAL = 128 * 1024 * 1024
MAX_FILES = 20000
ARCHIVE_ENDINGS = ('.zip', '.tar', '.tar.gz', '.tgz', '.gz', '.bz2', '.xz', '.7z')


def git(*args):
    return subprocess.check_output(['git', '-C', str(ROOT), *args])


def excluded(path):
    return (any(v.endswith(('.writing', '.pending', '.lock')) for v in path.parts)
            or path.name.endswith(ARCHIVE_ENDINGS)
            or re.search(r'\.part[0-9]+$', path.name) is not None
            or any(v in ('entries', '__pycache__', '.git', '.venv', 'node_modules') for v in path.parts))


def dump(path, value):
    with path.open('x') as f:
        json.dump(value, f, indent=2, ensure_ascii=False)
        f.write('\n'); f.flush(); os.fsync(f.fileno())


def identity(path):
    s = path.lstat()
    if not stat.S_ISREG(s.st_mode):
        raise ValueError('nonregular/symlink source refused: ' + str(path))
    return {'bytes': s.st_size, 'mtime_ns': s.st_mtime_ns, 'mode': stat.S_IMODE(s.st_mode)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True,
                        help='NEW output directory strictly below campaign packaging; must not exist')
    parser.add_argument('--upload-receipt', type=Path,
                        default=PACK / 'delivery_upload_parallel/UPLOAD_RECEIPT.json')
    args = parser.parse_args()
    started = time.perf_counter()
    dest = args.out.expanduser().absolute()
    if dest.exists() or PACK.resolve() not in dest.resolve().parents:
        parser.error('--out must be a NEW directory below campaign packaging')
    # No mutation is performed until the final upload and pack gates pass.
    upload_path = args.upload_receipt.expanduser().resolve()
    upload = json.loads(upload_path.read_text())
    if upload.get('status') != COMPLETE or upload.get('error'):
        raise ValueError('final uploader receipt is not verified complete')
    if upload.get('tag') != 'research-p3lr4-20261025' or upload.get('repo') != 'migodam/b1-probabilistic-maxwell-imaging':
        raise ValueError('unexpected final uploader repository/tag')
    bundle = PACK / 'full_evidence_20261025_v2'
    receipt = json.loads((bundle / 'RECEIPT.json').read_text())
    dependency = json.loads((bundle / 'inventory/DEPENDENCIES.json').read_text())
    if receipt.get('status') != 'COMPLETE' or dependency.get('workspace_commit') != ARCHIVED_SOURCE:
        raise ValueError('original archive COMPLETE/source identity mismatch')
    expected = {v['file']: v['bytes'] for v in receipt['parts']}
    uploaded = upload.get('assets', {})
    if len(expected) != len(receipt['parts']) or set(uploaded) != set(expected):
        raise ValueError('final uploader asset map differs from archive parts')
    if any(uploaded[name].get('bytes') != size for name, size in expected.items()):
        raise ValueError('final uploader asset size differs from archive parts')
    current = git('rev-parse', 'HEAD').decode().strip()
    changed = git('diff', '--name-only', '-z', '--diff-filter=ACMRT', ARCHIVED_SOURCE,
                  '--', str(REL), 'tools').split(b'\0')
    selected = {}; skipped = []; deleted = []
    source_bytes = 0

    def add(path, origin):
        nonlocal source_bytes
        path = path.absolute()
        if path == dest or dest in path.parents:
            raise ValueError('output overlaps a selected metadata scope')
        rel = path.relative_to(ROOT)
        if excluded(rel):
            skipped.append({'path': str(rel), 'reason': 'archive/part/temp/lock/runtime exclusion'})
            return
        if path.is_symlink() or any(p.is_symlink() for p in path.parents if p != ROOT.parent):
            raise ValueError('symlink source/ancestor refused: ' + str(path))
        row = identity(path)
        if row['bytes'] > PER_FILE:
            raise ValueError('unexpected oversized metadata (>16MiB): ' + str(rel))
        if PACK in path.parents:
            arc = 'delivery_metadata/' + str(path.relative_to(PACK))
        else:
            arc = 'workspace/' + str(rel)
        if arc in selected:
            selected[arc]['origins'].append(origin)
            return
        source_bytes += row['bytes']
        if source_bytes > TOTAL or len(selected) >= MAX_FILES:
            raise ValueError('supplement exceeds 128MiB/20000-file budget')
        selected[arc] = dict(row, source=str(rel), archive_path=arc, origins=[origin])

    def narrow_tree(path, origin):
        if not path.exists():
            return
        if path.is_symlink():
            raise ValueError('scope symlink refused: ' + str(path))
        for base, dirs, names in os.walk(path, followlinks=False):
            for name in list(dirs):
                child = Path(base) / name
                if child.is_symlink():
                    raise ValueError('scope directory symlink refused: ' + str(child))
                if excluded(child.relative_to(ROOT)):
                    dirs.remove(name)
            for name in sorted(names):
                add(Path(base) / name, origin)

    for raw in changed:
        if not raw:
            continue
        rel = Path(os.fsdecode(raw))
        if rel.is_absolute() or '..' in rel.parts:
            raise ValueError('unsafe Git path')
        # Source/results metadata may be packaged, scientific model/input payloads may not.
        within = rel.relative_to(REL) if rel.is_relative_to(REL) else None
        if within is not None and within.parts[0] in ('models', 'inputs'):
            raise ValueError('unexpected changed scientific payload outside supplement scope: ' + str(rel))
        if rel.suffix in ('.npz', '.npy', '.pt', '.pth', '.bin', '.dat') and not rel.is_relative_to(REL / 'AUDIT/public_replay'):
            raise ValueError('unexpected scientific array/model changed outside saved public replay: ' + str(rel))
        path = ROOT / rel
        if not path.exists():
            deleted.append(str(rel)); continue
        add(path, 'tracked Git diff since archive source')
    deleted = list(dict.fromkeys(deleted + [os.fsdecode(p) for p in git(
        'diff', '--name-only', '-z', '--diff-filter=D', ARCHIVED_SOURCE,
        '--', str(REL), 'tools').split(b'\0') if p]))
    add(Path(__file__).resolve(), 'explicit post-pack delivery packager source')
    # These are saved replay bytes only: no numpy imports, decoding, scoring, or inference.
    narrow_tree(CAMPAIGN / 'AUDIT/public_replay', 'saved public replay outputs/docs, including small NPZs')
    narrow_tree(CAMPAIGN / 'AUDIT/delivery_layout', 'restore layout and small part manifests')
    for name in ('RECEIPT.json', 'DELIVERY_RECEIPT.json', 'CATEGORY_EVENTS.jsonl', 'PARTS_CLOSED.jsonl'):
        add(bundle / name, 'final fullpack metadata')
    for name in ('SUMMARY.json', 'DEPENDENCIES.json', 'EXCLUSIONS_AND_SYMLINKS.json'):
        add(bundle / 'inventory' / name, 'small fullpack inventory metadata; entries excluded')
    # Interrupted pre-archive and sequential-upload receipts are part of delivery
    # history too; retain them without copying abandoned tar streams/journals.
    interrupted = PACK / 'full_evidence_20261025_v1'
    for name in ('RECEIPT.json', 'DELIVERY_RECEIPT.json', 'CATEGORY_EVENTS.jsonl', 'PARTS_CLOSED.jsonl'):
        path = interrupted / name
        if path.is_file():
            add(path, 'preserved interrupted first-pack metadata')
    for name in ('delivery_upload', 'delivery_upload_v2', 'delivery_upload_parallel', 'abandoned_408_batch_recovery',
                 'abandoned_starter_recovery', 'safe_drain_resize', 'final_delivery_assets'):
        narrow_tree(PACK / name, 'bounded delivery logs/receipts/recovery metadata')
    for path in sorted(PACK.iterdir()):
        if path.is_file() and path.suffix in ('.log', '.jsonl', '.md'):
            add(path, 'top-level delivery log/accounting/documentation')
    if not selected:
        raise ValueError('empty supplement refused')
    upload_id = identity(upload_path)
    # Custom final receipt is also retained, even when outside default upload subtree.
    add(upload_path, 'explicit final uploader receipt')
    final_audit = upload_path.parent / 'runs' / upload.get('run_id', '') / 'REMOTE_ASSET_AUDIT.json'
    if not final_audit.is_file():
        raise ValueError('final remote name/size/state audit missing')
    audit = json.loads(final_audit.read_text())
    if audit.get('status') != 'NAME_SIZE_STATE_VERIFIED' or audit.get('expected') != expected:
        raise ValueError('final remote audit expected map/status mismatch')
    for name, size in expected.items():
        asset = audit.get('remote_assets', {}).get(name, {})
        if asset.get('bytes') != size or asset.get('state') != 'uploaded':
            raise ValueError('final remote audit asset name/size/state mismatch: ' + name)
    add(final_audit, 'final remote name/size/state audit')
    summary = {key: upload.get(key) for key in (
        'status', 'error', 'tag', 'repo', 'run_id', 'wall_seconds', 'aggregate_cost_seconds',
        'cumulative_wall_seconds', 'cumulative_aggregate_cost_seconds', 'max_parallel_upload_streams')}
    summary.update(parts=len(expected), part_bytes=sum(expected.values()))
    manifest = {'format': 'R4_POSTPACK_SUPPLEMENT_V1', 'archive_source_commit': ARCHIVED_SOURCE,
                'current_commit': current, 'current_worktree_may_be_dirty': True,
                'source_selection': 'actual tracked changed bytes plus explicit saved replay/delivery scopes',
                'created_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
                'scope': 'delivery supplement only; 0 RHS; no model/score changes or execution; original archive preserved',
                'research_RHS': 0, 'inference': 0, 'training': 0, 'hashes_used': False,
                'zip_format_CRC_only': True, 'per_file_limit_bytes': PER_FILE,
                'total_payload_limit_bytes': TOTAL, 'selected_source_bytes': source_bytes,
                'final_upload_summary': summary, 'deleted_tracked_paths_not_recreated': deleted,
                'excluded_candidates': skipped, 'files': [selected[k] for k in sorted(selected)]}
    read_first = '''# R4 post-pack delivery supplement\n\nThis ZIP is separate from the sealed 120-part evidence archive at source\n60066394392fe04fb6a5b30371ffa34c978ade9d. It does not replace or rewrite it.\n\nExtract this supplement into a NEW directory. workspace/ holds repository-relative\npost-pack changed files and saved public-replay outputs; delivery_metadata/ holds\npackaging-relative receipts, bounded logs, and final asset audit. Review MANIFEST.json\nfor exact file sizes, nanosecond mtimes, origins, source/current commits, deletions,\nexcluded paths, and final upload summary. Do not silently overlay sealed evidence.\n\nSaved public replay NPZs are copied as opaque saved bytes only. This packaging job\nperforms 0 RHS, 0 training, 0 inference, and no model changes, scoring, or new scientific\nqualification. Updated reports may describe already completed work. Dependencies and\nhistorical absolute-path caveats remain in the separate restore layout documentation.\n\nNo content hashes are generated or checked. ZIP-format CRC fields are produced by the\nstandard ZIP writer; these do not establish scientific reproducibility.\n'''
    generated = {'READ_FIRST.md': read_first.encode(),
                 'MANIFEST.json': (json.dumps(manifest, indent=2, ensure_ascii=False)+'\n').encode()}
    if any(len(v) > PER_FILE for v in generated.values()) or source_bytes + sum(map(len, generated.values())) > TOTAL:
        raise ValueError('manifest plus source payload exceeds supplement limits')
    # Only now create the NEW destination. Never replace an existing supplement.
    dest.mkdir(parents=False, exist_ok=False)
    zip_path = dest / 'Research_P3LR4_Postpack_Supplement.zip'
    pending_zip = dest / 'Research_P3LR4_Postpack_Supplement.zip.writing'
    result = {'status': 'FAILED', 'archive_source_commit': ARCHIVED_SOURCE, 'current_commit': current,
              'research_RHS': 0, 'inference': 0, 'training': 0, 'hashes_used': False}
    try:
        with zipfile.ZipFile(pending_zip, 'x', compression=zipfile.ZIP_DEFLATED, compresslevel=1) as z:
            for name, data in generated.items():
                z.writestr(name, data)
            for arc in sorted(selected):
                row = selected[arc]; path = ROOT / row['source']
                if identity(path) != {k: row[k] for k in ('bytes', 'mtime_ns', 'mode')}:
                    raise RuntimeError('source changed before copy: ' + row['source'])
                info = zipfile.ZipInfo.from_file(path, arcname=arc, strict_timestamps=False)
                info.compress_type = zipfile.ZIP_DEFLATED
                with path.open('rb') as src, z.open(info, 'w', force_zip64=True) as dst:
                    shutil.copyfileobj(src, dst, length=1024*1024)
                if identity(path) != {k: row[k] for k in ('bytes', 'mtime_ns', 'mode')}:
                    raise RuntimeError('source changed during copy: ' + row['source'])
        if identity(upload_path) != upload_id or git('rev-parse', 'HEAD').decode().strip() != current:
            raise RuntimeError('final uploader receipt or current commit changed during packaging')
        pending_zip.rename(zip_path)
        for name, data in generated.items():
            with (dest / name).open('xb') as f:
                f.write(data);f.flush();os.fsync(f.fileno())
        result.update(status='POSTPACK_SUPPLEMENT_CREATED', files=len(selected), source_bytes=source_bytes,
                      zip_bytes=zip_path.stat().st_size, final_upload_summary=summary)
    except BaseException as exc:
        result['error'] = repr(exc)
        raise
    finally:
        result['wall_seconds'] = time.perf_counter()-started
        dump(dest / 'PACKAGING_RECEIPT.json', result)
    print(json.dumps(result))


if __name__ == '__main__':
    main()
