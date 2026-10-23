#!/usr/bin/env python3
"""Post-CLOSED evidence inventory and independent ZIP assets, standard library only.

Never imports physics, arrays, Book, or inference. Never publishes. Data archives
are opaque byte streams copied only after CLOSED.json exists. --inventory-only
has the same closing gate and does not open any array/data payload.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import stat
import sys
import zipfile

MAX_ASSET_BYTES = 1_800_000_000  # decimal bytes, <=1.8 GB per Release asset
ZIP_RESERVE = 65_536
COPY_BLOCK = 2 * 1024 * 1024
PRUNED = {'.git', '.deps', '.venv', 'venv', '__pycache__', 'node_modules',
          '.mypy_cache', '.pytest_cache', '.ruff_cache', '.DS_Store',
          'release_evidence', 'dist', 'build'}
TEXT = {'.md', '.rst', '.txt', '.html', '.css', '.js', '.mjs', '.json', '.jsonl',
        '.csv', '.tsv', '.log', '.out', '.err', '.yaml', '.yml', '.toml', '.ini'}
ARRAY = {'.npz', '.npy', '.h5', '.hdf5', '.mat', '.nc', '.dat', '.bin', '.pkl',
         '.pickle', '.pt', '.pth', '.msh', '.vtk', '.vtu', '.xdmf', '.xmf'}
FIGURE = {'.png', '.svg', '.jpg', '.jpeg', '.webp', '.pdf', '.gif'}
CODE_DIRS = {'src', 'vendor', 'tools', 'tests'}
STAGE_DIRS = CODE_DIRS | {'configs', 'inputs', 'results', 'AUDIT', 'figures',
                         'tables', 'docs'}
HISTORY_DIRS = CODE_DIRS | {'configs', 'inputs', 'results', 'docs', 'figures',
                           'research'}
SOURCE_EXTENSIONS = {'.py', '.sh', '.c', '.cpp', '.cc', '.h', '.hpp', '.f', '.f90',
                     '.m', '.jl', '.r', '.R'}
DEPENDENCY_NAMES = {'pyproject.toml', 'uv.lock', 'poetry.lock', 'Pipfile',
                    'Pipfile.lock', 'environment.yml', 'environment.yaml'}


def closed_gate(root: Path) -> Path:
    marker = root / 'CLOSED.json'
    if marker.is_symlink() or not marker.is_file():
        raise RuntimeError(f'Packaging is disabled until a regular {marker} exists; '
                           'this tool never creates a closing marker.')
    return marker


def public_json(path: Path):
    """Only explicit public configuration/validation metadata, never data or GT."""
    if any(word in path.name.upper() for word in ('TRUTH', 'GENERATION_PLAN')):
        raise RuntimeError('Internal policy error: protected content is never parsed')
    if path.is_symlink():
        raise RuntimeError(f'Public metadata must not be a symlink: {path}')
    return json.loads(path.read_text(encoding='utf-8'))


def witness_subset(root: Path):
    """Use the frozen FIRST nomination, not a best-result or directory-order search."""
    nomination = root / 'configs/FIRST_JOINT_NOMINATION.json'
    cfg = public_json(nomination)
    name = cfg.get('id')
    if not isinstance(name, str) or Path(name).name != name or name in {'.', '..'}:
        raise ValueError('Invalid frozen first-joint nomination id')
    base = Path('results/B2/validation') / name
    gate = public_json(root / base / 'RESULTS.json')
    if gate.get('numerical_pass') is not True or gate.get('joint_physics_threshold') is not True:
        raise RuntimeError('The first nominated joint witness lacks the two recorded validation flags')
    if public_json(root / base / 'MANIFEST.json') != cfg:
        raise RuntimeError('First nomination and independent validation manifest differ')
    selected = []
    for kind in ('cal', 'joint8'):
        for world in (0, 1):
            endpoint = base / 'fields' / f'{kind}_world{world}_fine' / 'endpoint_0000'
            # Check names, presence and nonempty size; do not open NPZ headers/arrays.
            for filename in ('PARAMETERS.json', 'FREQUENCIES.json', 'RESPONSE.npz',
                             'COMPLETE.json', *(f'f{i}.npz' for i in range(5))):
                path = root / endpoint / filename
                if path.is_symlink() or not path.is_file() or path.stat().st_size == 0:
                    raise RuntimeError(f'Incomplete selected public full-field endpoint: {path}')
            selected.append(endpoint.as_posix())
    return dict(nomination=nomination.relative_to(root).as_posix(),
                validation=base.as_posix(), endpoints=selected,
                expected_frequency_archives=20,
                scope='first independently validated joint witness; cal/joint, '
                      'two worlds, fine grid, all five frequencies per endpoint; '
                      'whole field/material/mesh archives and all endpoint companions',
                payload_validation='filenames/presence/bytes only; NPZ members not inspected')


def opaque_data(rel: Path) -> bool:
    parts = rel.parts
    return len(parts) >= 4 and parts[0] == 'results' and parts[1] in {'P1A', 'P1B'} and parts[2] == 'data'


def protected_name(rel: Path) -> bool:
    return (any('GENERATION_PLAN' in part.upper() or part.upper() in {'OFFLINE_TRUTH', 'GT'} for part in rel.parts)
            or (rel.suffix.lower() in ARRAY and 'TRUTH' in rel.name.upper()))


def raw_field(rel: Path) -> bool:
    """Known raw-field producers, identified by names rather than a size cutoff.

    Generic ARRAYS.npz is deliberately *not* excluded: historical response,
    parameter and diagnostic aggregates use that name and must remain whole.
    """
    name = rel.name.lower()
    if rel.suffix.lower() not in ARRAY:
        return False
    return bool(re.fullmatch(r'f\d+\.(npz|npy)', name)
                or name.startswith('world_')
                or re.search(r'(^|[_-])(full[_-]?fields?|fields?|electric[_-]?fields?|'
                             r'scattered[_-]?fields?|total[_-]?fields?)([_\-.]|$)', name))


def classify(rel: Path, origin: str, subset: dict):
    parts = rel.parts
    name, suffix = rel.name, rel.suffix.lower()
    if name in PRUNED or suffix in {'.pyc', '.pyo', '.zip', '.tar', '.gz'} or (suffix == '.lock' and name not in DEPENDENCY_NAMES):
        return False, 'INFRASTRUCTURE_OR_PRIOR_PACKAGE', 'runtime/cache/lock/archive is not research evidence'
    # GT and plan content is never parsed. Only designated P1 data files may be
    # copied opaquely after closure, including the sealed offline evidence files.
    if protected_name(rel):
        if origin == 'stage' and opaque_data(rel) and suffix in ARRAY | TEXT:
            return True, 'POST_CLOSE_SEALED_DATA_OPAQUE_COPY', 'CLOSED-gated whole-file copy; never decode or inspect GT/plan'
        return False, 'PROTECTED_NON_DATA_PATH', 'no permission to disclose protected material outside designated P1 data files'
    if origin == 'stage' and any(rel.as_posix().startswith(p + '/') for p in subset['endpoints']):
        return True, 'PUBLIC_FULL_FIELD_SUBSET', 'entire selected endpoint, including field/material/mesh archives and companions'
    if rel.as_posix() == 'results/p0b/finer-endpoints/ARRAYS.npz':
        return False, 'LOCAL_MIXED_FULL_FIELD_ARCHIVE', 'historical mixed archive has raw *_field members; retain original locally and export nonfield members after closure'
    if suffix in ARRAY and re.search(r'(^|[_-])(mean|prediction|q|checkpoint|latent|iterate|diagnostic|error|difference|comparison)([_\-.]|$)', name.lower()):
        return True, 'ENDPOINT_MEAN_PREDICTION_Q_CHECKPOINT_COMPLETE', 'explicit endpoint/checkpoint evidence retained whole, regardless of size'
    if raw_field(rel):
        return False, 'LOCAL_FULL_FIELD', 'retained unchanged locally; public full-field subset is restricted to the first validated joint witness'
    top = parts[0]
    allowed_dirs = STAGE_DIRS if origin == 'stage' else HISTORY_DIRS
    if origin == 'historical_local':
        # Mutable external source never replaces the frozen replay source/vendor.
        allowed_dirs = {'results', 'figures', 'docs'}
    if len(parts) > 1 and top not in allowed_dirs:
        return False, 'OUTSIDE_ALLOWLIST', 'directory is outside the explicit evidence/source allowlist'
    if top in CODE_DIRS and origin != 'historical_local':
        return True, 'REPLAY_SOURCE_OR_VENDOR', 'whole source/vendor/tool/test file, excluding runtimes and bytecode'
    if top == 'configs' and origin != 'historical_local':
        return True, 'FROZEN_CONFIGURATION', 'configuration retained verbatim'
    if top == 'inputs' and origin != 'historical_local' and suffix in ARRAY | TEXT:
        return True, 'PUBLIC_INPUT', 'public fixed input retained whole'
    if origin == 'stage' and opaque_data(rel) and suffix in ARRAY | TEXT:
        return True, 'POST_CLOSE_DATA_OPAQUE_COPY', 'post-closure whole-file data evidence; no array/header parsing'
    if 'SOURCE_SNAPSHOT' in parts or 'SOURCE_CONTRACT' in parts:
        return True, 'RECORDED_SOURCE_SNAPSHOT', 'literal source contract/snapshot retained'
    if origin == 'stage' and len(parts) > 2 and parts[:3] in {
            ('results', 'P1A', 'rules'), ('results', 'P1B', 'full_bank')}:
        return True, 'PUBLIC_RULE_AGGREGATE_COMPLETE', 'every existing public aggregate and companion; no size cutoff'
    if origin == 'stage' and len(parts) > 2 and parts[1] == 'P1B' and parts[2] in {'inference', 'offline_scores', 'offline'}:
        return True, ('FROZEN_OFFLINE_SCORE_COMPLETE' if 'offline' in parts[2]
                      else 'INFERENCE_RAW_TRAJECTORY_AND_OUTPUT_COMPLETE'), 'complete method evidence subtree; raw trajectories retained regardless of size'
    if suffix in ARRAY:
        upper = name.upper()
        if any(t in upper for t in ('TRAJECTORY', 'CHAIN', 'PARTICLE', 'PROPOSAL')):
            role = 'CHAIN_PARTICLE_PROPOSAL_ARRAY_COMPLETE'
        elif any(t in upper for t in ('MEAN', 'PREDICTION', 'LATENT', 'ITERAT', 'CHECKPOINT', 'ACCEPTED', 'FINAL', 'CURRENT', 'PARAMETER')):
            role = 'ENDPOINT_AND_ITERATION_CHECKPOINT_COMPLETE'
        elif name.startswith('J_f'):
            role = 'JACOBIAN_DIAGNOSTIC_ARRAY_COMPLETE'
        else:
            role = 'RESPONSE_MATERIAL_DIAGNOSTIC_AGGREGATE_COMPLETE'
        return True, role, 'whole scientific array/archive; selected by role/extension, never by a small-byte threshold'
    if suffix in TEXT | FIGURE | SOURCE_EXTENSIONS or name in DEPENDENCY_NAMES:
        return True, 'REPORT_METADATA_LEDGER_LOG_OR_FIGURE', 'reports, metadata, ledgers, failures, dependency requirements and figures retained'
    if name.startswith(('requirements', 'LICENSE', 'NOTICE', 'COPYING')):
        return True, 'DEPENDENCY_REQUIREMENTS_OR_LICENSE', 'version requirements/licenses only; system runtime is excluded'
    return False, 'UNRECOGNIZED_FILE', 'not in explicit source/report/array/data allowlist; listed for owner review'


def snapshot(path: Path):
    s = path.lstat()
    return dict(bytes=s.st_size, mtime_ns=s.st_mtime_ns, device=s.st_dev,
                inode=s.st_ino, mode=stat.S_IMODE(s.st_mode))


def walk_inventory(base: Path, origin: str, prefix: str, subset: dict, output: Path):
    rows, sources = [], {}
    for directory, dirs, files in os.walk(base, followlinks=False):
        parent = Path(directory)
        for name in sorted(dirs[:]):
            child = parent / name
            if name in PRUNED or child.is_symlink() or child == output or output in child.parents:
                dirs.remove(name)
                rel = child.relative_to(base).as_posix()
                rows.append(dict(origin=origin, relative_path=rel, archive_path=f'{prefix}/{rel}',
                                 bytes=None, role='EXCLUDED_DIRECTORY', included=False,
                                 reason='runtime/cache/package-output/symlink subtree; not traversed',
                                 local_path=str(child), kind='directory'))
        for name in sorted(files):
            path = parent / name
            rel = path.relative_to(base)
            archive = f'{prefix}/{rel.as_posix()}'
            if path.is_symlink() or not stat.S_ISREG(path.lstat().st_mode):
                rows.append(dict(origin=origin, relative_path=rel.as_posix(), archive_path=archive,
                                 bytes=path.lstat().st_size, role='EXCLUDED_SPECIAL_FILE', included=False,
                                 reason='symlink or nonregular file is never followed/copied', local_path=str(path), kind='file'))
                continue
            included, role, reason = classify(rel, origin, subset)
            row = dict(origin=origin, relative_path=rel.as_posix(), archive_path=archive,
                       **snapshot(path), role=role, included=included, reason=reason,
                       local_path=str(path), kind='file')
            rows.append(row)
            if included:
                sources[archive] = (path, row)
    return rows, sources


def derive_historical_nonfield(output: Path, rows: list, sources: dict, inventory_only: bool):
    """Only the source-confirmed mixed historical archive, never P1 data/GT.

    Inspect ZIP member names/sizes, not NPY headers or array values. Full-field
    members are never decompressed. Copy other members as opaque NPY byte
    streams into a clearly marked derived NPZ, leaving its original untouched.
    """
    originals = [r for r in rows if r['role'] == 'LOCAL_MIXED_FULL_FIELD_ARCHIVE']
    for original in originals:
        path = Path(original['local_path'])
        check_frozen(path, original)
        archive = '_NONFIELD_DERIVED/' + original['archive_path']
        target = output / 'derived_nonfield' / original['archive_path']
        members = []
        with zipfile.ZipFile(path, 'r') as src:
            infos = src.infolist()
            if len({info.filename for info in infos}) != len(infos):
                raise RuntimeError('Duplicate historical NPZ member names: ' + str(path))
            selected = []
            for info in infos:
                name = PurePosixPath(info.filename)
                if name.is_absolute() or '..' in name.parts or not info.filename.endswith('.npy') or info.is_dir():
                    raise RuntimeError('Unexpected historical NPZ member: ' + info.filename)
                excluded = bool(re.search(r'_(field|dofs)\.npy$', info.filename))
                members.append(dict(member=info.filename, uncompressed_bytes=info.file_size,
                    compressed_bytes=info.compress_size, included_in_nonfield_derivative=not excluded,
                    reason='raw complex field/DOF member stays local' if excluded else 'response/material/diagnostic NPY member copied without decoding'))
                if not excluded:
                    selected.append(info)
            if not any(info.filename.endswith('_data.npy') for info in selected):
                raise RuntimeError('Known mixed archive has no expected response members: ' + str(path))
            if not inventory_only:
                target.parent.mkdir(parents=True, exist_ok=True)
                with zipfile.ZipFile(target, 'x', compression=zipfile.ZIP_STORED, allowZip64=True) as dest:
                    for info in selected:
                        zinfo = zipfile.ZipInfo(info.filename)
                        zinfo.compress_type = zipfile.ZIP_STORED
                        zinfo.file_size = info.file_size
                        with src.open(info, 'r') as inp, dest.open(zinfo, 'w', force_zip64=True) as out:
                            while True:
                                data = inp.read(COPY_BLOCK)
                                if not data:
                                    break
                                out.write(data)
        check_frozen(path, original)
        original['mixed_archive_members'] = members
        row = dict(origin='historical_nonfield_derivative', relative_path=original['archive_path'],
            archive_path=archive, bytes=None, role='HISTORICAL_RESPONSE_MATERIAL_NONFIELD_DERIVATIVE',
            included=not inventory_only, reason='planned member-only derivative; not created in inventory-only mode' if inventory_only else 'opaque NPY members excluding *_field.npy/*_dofs.npy; original mixed archive remains local',
            local_path=str(target), kind='file', derived_from=original['archive_path'],
            array_values_decoded=False, GT_or_plan_access=False)
        if not inventory_only:
            row.update(snapshot(target))
            sources[archive] = (target, row)
        rows.append(row)


def check_frozen(path: Path, row: dict):
    if path.is_symlink() or path.resolve() != path or not path.is_file():
        raise RuntimeError(f'Source changed type/path after inventory: {path}')
    now = snapshot(path)
    if any(now[key] != row[key] for key in ('bytes', 'mtime_ns', 'device', 'inode')):
        raise RuntimeError(f'Source changed after inventory (no mutable package accepted): {path}')


def payloads(sources: dict, limit: int):
    """Split opaque files only if needed; chunks preserve original bytes."""
    result, splits = [], []
    chunk_limit = limit - ZIP_RESERVE - 4096
    for archive, (path, row) in sorted(sources.items()):
        length = row['bytes']
        count = max(1, math.ceil(length / chunk_limit))
        segments = []
        for i in range(count):
            offset = i * chunk_limit
            size = min(chunk_limit, length - offset)
            member = archive if count == 1 else f'__chunks__/{archive}.part{i + 1:06d}'
            item = dict(path=path, row=row, member=member, offset=offset, length=size)
            result.append(item)
            segments.append(item)
        if count > 1:
            splits.append(dict(archive_path=archive, original_bytes=length, segments=segments))
    return result, splits


def plan_assets(items: list, limit: int, current=None):
    groups = [] if current is None else current
    for item in items:
        overhead = 2 * len(item['member'].encode('utf-8')) + 256
        needed = item['length'] + overhead
        if needed + ZIP_RESERVE > limit:
            raise RuntimeError(f'Entry cannot fit requested asset cap: {item["member"]}')
        if not groups or groups[-1]['estimate'] + needed > limit:
            groups.append(dict(name=f'w-evidence-part{len(groups) + 1:03d}.zip',
                               estimate=ZIP_RESERVE, entries=[]))
        item['asset'] = groups[-1]['name']
        groups[-1]['entries'].append(item)
        groups[-1]['estimate'] += needed
    return groups


def json_write(path: Path, value):
    with path.open('x', encoding='utf-8') as stream:
        json.dump(value, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write('\n')


def expected_rule_names(root: Path):
    cfg = public_json(root / 'configs/P1_PROBABILITY_PROTOCOL.json')
    expected = [f'results/P1A/rules/full_{n}.npz' for n in cfg['full_orders']]
    expected += [f'results/P1A/rules/cdf_{i}_{n}.npz'
                 for i in range(len(cfg['cdf_checkpoints']))
                 for n in cfg.get('cdf_evaluation_orders', cfg['reference_orders'])]
    return dict(expected=expected, missing_locally=[p for p in expected if not (root / p).is_file()],
                interpretation='missing is explicit, never represented as complete or manufactured; '
                               'all existing public P1A/P1B aggregates are included without a size threshold')


def emit_metadata(output: Path, rows: list, subset: dict, sources: dict,
                  splits: list, groups: list, rule_status: dict, inventory_only: bool, limit: int):
    local_fields = [r for r in rows if r['role'] in {'LOCAL_FULL_FIELD', 'LOCAL_MIXED_FULL_FIELD_ARCHIVE', 'PUBLIC_FULL_FIELD_SUBSET'}]
    inventory = dict(policy_version=1, status='POST_CLOSED_INVENTORY_ONLY' if inventory_only else 'POST_CLOSED_PACKAGE_SELECTION',
                     records=rows, selected_files=sum(r['included'] for r in rows),
                     selected_bytes=sum(r['bytes'] or 0 for r in rows if r['included']),
                     planned_derivatives_not_created=[r['archive_path'] for r in rows if r['origin']=='historical_nonfield_derivative' and not r['included']],
                     size_filter_used=False, array_payloads_decoded=False, scientific_gate_inferred=False,
                     hashes_used=False, public_field_subset=subset, public_rule_presence=rule_status)
    json_write(output / 'INVENTORY.json', inventory)
    columns = ['origin', 'relative_path', 'archive_path', 'bytes', 'role', 'included', 'reason', 'local_path', 'kind']
    with (output / 'INVENTORY.csv').open('x', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, extrasaction='ignore')
        writer.writeheader(); writer.writerows(rows)
    json_write(output / 'LOCAL_FULL_FIELDS_INVENTORY.json', dict(records=local_fields,
        full_fields_all_uploaded=False, local_originals_preserved=True, public_subset=subset,
        excluded_bytes=sum(r['bytes'] or 0 for r in local_fields if not r['included']),
        reason='Historical and other current raw complex fields remain in their original local paths; '
               'Release provides the explicit first-witness cal/joint fine-grid subset only. '
               'Complete response/parameter/chain/diagnostic aggregates are separately retained.'))
    reassembly = []
    for item in splits:
        reassembly.append(dict(archive_path=item['archive_path'], original_bytes=item['original_bytes'],
            chunks=[dict(asset=s['asset'], member=s['member'], offset=s['offset'], bytes=s['length']) for s in item['segments']]))
    json_write(output / 'REASSEMBLY.json', dict(bytewise_only=True, decoded_arrays=False, files=reassembly,
        restore='Extract every ZIP into one empty directory. For each listed file, concatenate its chunk members '
                'in listed offset order into archive_path (binary mode, exclusive output), then check original_bytes. '
                'Unsplit files already have their original relative paths. No content hashes are generated.'))
    json_write(output / 'PACKAGE_MANIFEST.json', dict(policy_version=1,
        one_existing_repository_release=True, publishing_performed=False,
        payload_files=len(sources), data_zip_count_planned=len(groups), max_asset_bytes=limit,
        format='independent ZIP_STORED archives; very large individual files split into opaque byte chunks',
        source_change_check='size, mtime_ns, device, inode before/after streaming; no cryptographic content claim',
        full_fields_all_uploaded=False, inventory_and_reassembly_are_release_assets=True,
        system_runtime_included=False, dependencies='version requirements/lock metadata and fixed vendor source only',
        closing_marker='research_w_p0r_to_p1/CLOSED.json is required; no scientific gate is inferred from it',
        new_physics=False, GT_and_generation_plan='whole-file streaming copy only from designated P1A/P1B data after CLOSED; no content inspection'))
    (output / 'RESTORE.md').write_text(
        '# W evidence release assets\n\nAll ZIP files belong to one Release in the existing repository. '
        'They are independent archives, not a split-ZIP format. Extract all into one empty directory. '
        'Consult REASSEMBLY.json for any single-file byte chunks, then concatenate those chunks in offset order '
        'using binary mode and verify the recorded byte count. Do not treat counts as a hash verification.\n\n'
        'The public full-field subset is limited to the first independently validated joint witness, '
        'cal/joint, both worlds, fine grid, five complete frequency archives and all endpoint companion files. '
        'Historical and remaining current raw fields remain local; LOCAL_FULL_FIELDS_INVENTORY.json gives '
        'their paths, bytes and omission reasons. No claim is made that all full fields are uploaded.\n\n'
        'Opaque offline data/plan files are post-close evidence, not online-inference inputs. '
        'Source and fixed historical vendor are included beside the stage; system Python, installed packages, '
        'virtual environments and downloaded runtimes are not. No tool in this package automatically runs physics.\n',
        encoding='utf-8')
    return [output / name for name in ('INVENTORY.json', 'INVENTORY.csv',
            'LOCAL_FULL_FIELDS_INVENTORY.json', 'REASSEMBLY.json', 'PACKAGE_MANIFEST.json', 'RESTORE.md')]


def copy_segment(zf: zipfile.ZipFile, item: dict):
    path, row = item['path'], item['row']
    check_frozen(path, row)
    flags = os.O_RDONLY | getattr(os, 'O_NOFOLLOW', 0)
    fd = os.open(path, flags)
    with os.fdopen(fd, 'rb') as stream:
        s = os.fstat(stream.fileno())
        if (s.st_size, s.st_mtime_ns, s.st_dev, s.st_ino) != (row['bytes'], row['mtime_ns'], row['device'], row['inode']):
            raise RuntimeError(f'Source replaced while opening: {path}')
        stream.seek(item['offset'])
        stamp = dt.datetime.fromtimestamp(row['mtime_ns'] / 1e9, dt.timezone.utc)
        year = min(2107, max(1980, stamp.year))
        info = zipfile.ZipInfo(item['member'], (year, stamp.month, min(stamp.day, 28), stamp.hour, stamp.minute, stamp.second))
        info.compress_type = zipfile.ZIP_STORED
        info.create_system = 3
        info.external_attr = (stat.S_IFREG | 0o644 | (row['mode'] & 0o111)) << 16
        info.file_size = item['length']
        remaining = item['length']
        with zf.open(info, 'w', force_zip64=True) as dest:
            while remaining:
                data = stream.read(min(COPY_BLOCK, remaining))
                if not data:
                    raise RuntimeError(f'Source truncated during byte copy: {path}')
                dest.write(data); remaining -= len(data)
        after = os.fstat(stream.fileno())
        if (after.st_size, after.st_mtime_ns, after.st_dev, after.st_ino) != (row['bytes'], row['mtime_ns'], row['device'], row['inode']):
            raise RuntimeError(f'Source changed during byte copy: {path}')
    check_frozen(path, row)


def package(args):
    root = args.root.resolve()
    marker = closed_gate(root)  # before scanning/writing anything, including inventory-only
    history = (args.historical_root or root.parent / 'research_w_continuation').resolve()
    old_local = (args.historical_local_root or root.parents[3] / 'Research_W/research-w-continuation').resolve()
    for path in (root / 'src', history / 'src', history / 'vendor', old_local / 'results'):
        if not path.is_dir():
            raise RuntimeError(f'Required fixed replay source/local historical evidence is missing: {path}')
    subset = witness_subset(root)
    rule_status = expected_rule_names(root)
    limit = args.max_asset_bytes
    if not 1_000_000 <= limit <= MAX_ASSET_BYTES:
        raise ValueError(f'--max-asset-bytes must be between 1000000 and {MAX_ASSET_BYTES}')
    stamp = dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    output = (args.output or root / 'release_evidence' / stamp).resolve()
    if output.exists():
        raise FileExistsError('Output must be a new directory; earlier inventories/assets are never overwritten')
    # Output may be inside stage only under its excluded release_evidence subtree.
    if output.is_relative_to(history) or output.is_relative_to(old_local):
        raise ValueError('Output must not be placed inside a historical source/evidence tree')
    if output.is_relative_to(root) and not output.is_relative_to(root / 'release_evidence'):
        raise ValueError('Output inside stage must be under release_evidence/')
    rows, sources = [], {}
    for base, origin, prefix in ((root, 'stage', 'research_w_p0r_to_p1'),
            (history, 'historical_frozen', 'research_w_continuation'),
            (old_local, 'historical_local', 'historical_local_evidence')):
        current, selected = walk_inventory(base, origin, prefix, subset, output)
        rows.extend(current); sources.update(selected)
    if not any(r['included'] and r['origin'] == 'historical_frozen' and r['relative_path'].startswith('vendor/') for r in rows):
        raise RuntimeError('Fixed historical vendor source was not selected')
    # Top-level project dependency requirements augment stage/historical requirements.
    for path in sorted(root.parents[1].iterdir()):
        if path.is_file() and not path.is_symlink() and (path.name in DEPENDENCY_NAMES or path.name.startswith('requirements')):
            archive = '_dependency_requirements/' + path.name
            row = dict(origin='project_dependency_requirements', relative_path=path.name, archive_path=archive,
                       **snapshot(path), role='DEPENDENCY_REQUIREMENTS_ONLY', included=True,
                       reason='version declaration, not installed runtime', local_path=str(path), kind='file')
            rows.append(row); sources[archive] = (path, row)
    closed_gate(root)
    rows.sort(key=lambda r: (r['origin'], r['relative_path']))
    output.mkdir(parents=True, exist_ok=False)
    derive_historical_nonfield(output, rows, sources, args.inventory_only)
    items, splits = payloads(sources, limit)
    groups = plan_assets(items, limit)
    metadata = emit_metadata(output, rows, subset, sources, splits, groups, rule_status, args.inventory_only, limit)
    # Metadata assets are also embedded; adding them cannot renumber existing data assets.
    meta_sources = {f'_PACKAGE_METADATA/{p.name}': (p, snapshot(p)) for p in metadata}
    meta_items, meta_splits = payloads(meta_sources, limit)
    if meta_splits:
        raise RuntimeError('Metadata exceeds one asset payload; choose a larger --max-asset-bytes')
    groups = plan_assets(meta_items, limit, groups)
    if args.inventory_only:
        json_write(output / 'ASSETS.json', dict(status='INVENTORY_ONLY_NO_ZIP_CREATED',
            archive_plan=[dict(name=g['name'], upper_estimated_bytes=g['estimate'], entries=len(g['entries'])) for g in groups],
            standalone_metadata_assets=[p.name for p in metadata], max_asset_bytes=limit,
            array_data_payloads_opened=0, published=False))
        print(json.dumps(dict(status='INVENTORY_ONLY', output=str(output), selected_files=len(sources),
                              planned_zip_assets=len(groups), zip_created=False, published=False)))
        return
    assets = []
    try:
        for group in groups:
            closed_gate(root)
            target = output / group['name']
            partial = target.with_suffix('.zip.partial')
            with zipfile.ZipFile(partial, 'x', compression=zipfile.ZIP_STORED, allowZip64=True) as zf:
                for item in group['entries']:
                    copy_segment(zf, item)
            size = partial.stat().st_size
            if size > limit:
                raise RuntimeError(f'ZIP exceeds asset cap ({size}>{limit}); partial is retained, no success manifest')
            partial.rename(target)
            assets.append(dict(name=target.name, bytes=size, entries=len(group['entries'])))
        # Detect mutation of early copied sources as well; no content hashes used.
        for path, row in sources.values():
            check_frozen(path, row)
        for row in rows:
            if row['role'] == 'LOCAL_MIXED_FULL_FIELD_ARCHIVE':
                check_frozen(Path(row['local_path']), row)
        closed_gate(root)
        json_write(output / 'ASSETS.json', dict(status='LOCAL_PACKAGE_COMPLETE_NOT_PUBLISHED',
            zip_assets=assets, standalone_metadata_assets=[p.name for p in metadata],
            max_asset_bytes=limit, one_existing_repository_release=True, published=False,
            no_all_full_fields_upload_claim=True))
    except BaseException as exc:
        json_write(output / 'PACKAGE_FAILURE.json', dict(status='PACKAGE_INCOMPLETE_DO_NOT_PUBLISH',
            error_type=type(exc).__name__, error=str(exc), already_written_zip_assets=assets,
            partial_assets_retained=True, source_evidence_modified=False))
        raise
    print(json.dumps(dict(status='LOCAL_PACKAGE_COMPLETE_NOT_PUBLISHED', output=str(output),
                          zip_assets=len(assets), payload_files=len(sources), published=False)))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1], help='stage root; must already have CLOSED.json')
    parser.add_argument('--historical-root', type=Path, help='fixed published research_w_continuation source/vendor tree')
    parser.add_argument('--historical-local-root', type=Path, help='original local Research_W/research-w-continuation evidence tree (metadata scan, opaque copies only)')
    parser.add_argument('--output', type=Path, help='new output directory; default ROOT/release_evidence/UTC timestamp')
    parser.add_argument('--inventory-only', action='store_true', help='still requires CLOSED.json; write metadata/ZIP plan only, never open array/data payloads')
    parser.add_argument('--max-asset-bytes', type=int, default=MAX_ASSET_BYTES, help='hard cap <=1800000000; smaller values permit portable tests after closure')
    args = parser.parse_args(argv)
    try:
        package(args)
    except Exception as exc:
        print(f'package_w_evidence: {type(exc).__name__}: {exc}', file=sys.stderr)
        return 2
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
