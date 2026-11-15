"""Finite no-GT R15 review export; deliberately NOT a complete raw archive.

Require both actual main264/fork56 queues terminal. Preserve all fork bytes and
selected full states; reduce main five-view NPZs to their real m/Y_prediction.
Only IO and saved-array equality checks: no solver, physics, GT or selection.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import importlib.util
import io
import json
from pathlib import Path
import sys
import time
import traceback
import zipfile

sys.dont_write_bytecode = True
MAIN = 'twofold_exposed24_v1'
FORK = 'fullstate_interventions_v1'
COMMON = ('twofold_hop', 'twofold_dr', 'twofold_ap', 'twofold_full_dr', 'twofold_reset')
FOCUS = tuple(f'r14_new_{i:04d}' for i in (9, 20, 21, 23))
VIEWS = ('initial', 'raw', 'corrected', 'primary', 'secondary')
NAME = 'R15_REVIEW_EXPORT_V1'
TEXT_SUFFIXES = {'.json', '.jsonl', '.csv', '.log', '.txt', '.md'}


def now():
    return datetime.now(timezone.utc).isoformat()


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def write_new(path, data):
    with path.open('x', encoding='utf-8') as handle:
        json.dump(data, handle, indent=2, ensure_ascii=False, allow_nan=False)
        handle.write('\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    args = parser.parse_args()
    start_cpu, start_wall = time.process_time(), time.perf_counter()
    root = args.root.resolve()
    output = root / 'archives'
    destinations = {key: output / (NAME + suffix) for key, suffix in (
        ('zip', '.zip'), ('manifest', '_MANIFEST.json'), ('receipt', '_RECEIPT.json'),
        ('error', '_ERROR.json'), ('csv', '_FILES.csv'), ('readme', '_README.md'))}
    if not root.is_dir() or output.is_symlink() or any(path.exists() for path in destinations.values()):
        raise FileExistsError('New review-export filenames required; never overwrite any existing attempt')
    if sys.platform == 'darwin' and not output.resolve().is_relative_to('/Volumes'):
        raise ValueError('Use the external /Volumes disk for review exports')
    raw = module(root / 'src/package_raw.py', 'r15_review_raw_boundaries')
    source_api = module(root / 'src/queue_interventions.py', 'r15_review_actual_v2_source_api')
    exclusions, queue_info, source_paths = [], {}, {}
    # Only terminal validation/stat inventory; no raw package/source scan is run.
    for specification in raw.CAMPAIGNS:
        metadata, paths = raw.terminal_group(root, *specification, exclusions)
        queue_info[metadata['campaign']] = metadata
        source_paths[metadata['campaign']] = paths
    frozen_main = json.loads((root / 'queues' / MAIN / 'FROZEN_MANIFEST.json').read_text())
    tasks = frozen_main['tasks']
    if (len(tasks) != 264 or len({t['task_id'] for t in tasks}) != 264 or
        {t['case_id'] for t in tasks} != {f'r14_new_{i:04d}' for i in range(24)}):
        raise ValueError('Actual frozen264 task inventory required')
    rows, expected, problems, selections = {}, [], [], []

    def register(path, role, transform='ORIGINAL_BYTES'):
        path = Path(path)
        if raw.exclusion(path, root) or not path.resolve().is_relative_to(root):
            raise ValueError('Excluded/outside-root review source')
        relative = path.relative_to(root).as_posix()
        if not path.is_file():
            expected.append(dict(path=relative, role=role, exists=False, source_bytes=None,
                                 status='MISSING_RETAINED', reason='Not produced in terminal original run; no rerun or substitute'))
            problems.append(dict(path=relative, role=role, issue='MISSING_ORIGINAL'))
            return
        st = path.stat()
        expected.append(dict(path=relative, role=role, exists=True, source_bytes=st.st_size))
        if relative not in rows:
            rows[relative] = dict(path=relative, source_filename=path.name, source_bytes=st.st_size,
                source_mtime_ns=st.st_mtime_ns, roles=[role], transform=transform,
                status='PLANNED', exported_bytes=None, original_fields=None, retained_fields=None,
                dropped_fields=None, field_array_equal=None)
        elif role not in rows[relative]['roles']:
            if rows[relative]['transform'] != transform:
                raise ValueError('Conflicting copy/transform role: ' + relative)
            rows[relative]['roles'].append(role)

    main_paths = source_paths[MAIN]
    for task in tasks:
        ident, case, method = task['task_id'], task['case_id'], task['method']
        folder = root / 'RESULTS' / MAIN / ident
        if ident != case + '__' + method or task['campaign'] != MAIN:
            raise ValueError('Unexpected frozen task identity')
        # Actual small JSON/JSONL/csv/log artifacts, not a fabricated uniform schema.
        for path in main_paths:
            if path.is_relative_to(folder) and path.suffix.lower() in TEXT_SUFFIXES:
                if path.name.startswith('NATIVE_HISTORY') and not (case in FOCUS and method in COMMON):
                    continue
                register(path, 'MAIN_ORIGINAL_TEXT')
        for name in ('RECEIPT.json', 'HISTORY.json'):
            if not (folder / name).is_file():
                register(folder / name, 'MAIN_REQUIRED_OR_MISSING_TEXT')
        register(folder / 'CANDIDATES.npz', 'MAIN_ORIGINAL_PAID_MATERIAL_POOL')
        for view in VIEWS:
            register(folder / (view + '.npz'), 'MAIN_FIVE_VIEW_M_AND_Y', 'RETAIN_REAL_M_AND_Y_PREDICTION_ONLY')
            if not (folder / (view + '.json')).is_file():
                register(folder / (view + '.json'), 'MAIN_VIEW_ORIGINAL_METADATA')
        if method in COMMON:
            for index in range(1, 5):
                for side in ('BEFORE', 'AFTER'):
                    register(folder / f'CYCLE_{index:03d}_{side}.npz', 'COMMON_HOP_FIRST4_FULL_CYCLE_STATE')
        if case in FOCUS and method in COMMON and not (folder / 'NATIVE_HISTORY.json').is_file():
            register(folder / 'NATIVE_HISTORY.json', 'FOCUS4_COMMON5_NATIVE_TRAJECTORY')
    # Every actual fork file is retained byte-for-byte, with its original path.
    for path in source_paths[FORK]:
        if path.is_relative_to(root / 'RESULTS' / FORK):
            register(path, 'FORK_ALL_ORIGINAL_FILES')
    for campaign in (MAIN, FORK):
        queue = root / 'queues' / campaign
        for path in source_paths[campaign]:
            if not path.is_relative_to(queue):
                continue
            relative = path.relative_to(queue)
            if path.name.startswith('STATUS_'):
                continue
            if (relative.parts[0] in ('TERMINAL_JOBS', 'logs') or path.name in
                ('RECEIPT.json', 'FROZEN_MANIFEST.json', 'FROZEN_TASKS.json', 'MAIN_TERMINAL_SNAPSHOT.json') or
                path.name.startswith(('RESOURCE', 'EVENTS'))):
                register(path, 'TERMINAL_QUEUE_RECORD_FROZEN_RESOURCE_EVENT_JOB_LOG')
    for case in FOCUS:
        selection = source_api.source_checkpoint(root, case)
        if selection.get('requested_iteration') != 512:
            raise ValueError('Actual source_checkpoint API is not the declared V2 preferred512 policy')
        selections.append(dict(case_id=case, **selection))
        folder = root / 'RESULTS' / MAIN / (case + '__twofold_hop')
        register(folder / 'FIXED_OPERATOR_CACHE.npz', 'FOCUS4_HOP_ORIGINAL_FIXED_DENOMINATOR_CACHE')
        if selection.get('path') is None:
            problems.append(dict(case_id=case, issue='NO_COMPLETE_V2_SOURCE_CHECKPOINT', details=selection))
            continue
        chosen = Path(selection['path'])
        if not chosen.resolve().is_relative_to(folder):
            raise ValueError('V2 source checkpoint outside original hop task')
        register(chosen, 'FOCUS4_HOP_V2_FULL_SOURCE_CHECKPOINT')
        # Preserve real relative dependencies used by the selected main/fork source.
        info = source_api.complete_checkpoint(chosen)
        if not info.get('complete'):
            raise ValueError('Selected V2 source no longer complete')
        for key in ('fixed_cache', 'background_cache'):
            if info.get(key):
                register(Path(info[key]), 'FOCUS4_HOP_SOURCE_RELATIVE_CACHE_DEPENDENCY')
    if len([x for x in expected if x['role'] == 'COMMON_HOP_FIRST4_FULL_CYCLE_STATE']) != 960:
        raise ValueError('Fixed24 x5 x4 BEFORE/AFTER slots must total960 files /480 pairs')
    output.mkdir(parents=True, exist_ok=True)
    exported, array_checks = [], []
    manifest = dict(schema='R15_REVIEW_EXPORT_V1', role='REDUCED_REVIEW_COPY_NOT_FULL_RAW_RELEASE',
        root=str(root), created_utc=now(), expected_main_tasks=264, expected_fork_tasks=56,
        cycle_pairs=480, expected_cycle_files=960, cycle_count_explanation='24 cases x5 methods x4 cycles x2 BEFORE/AFTER files',
        queue_info=queue_info, source_checkpoint_API=dict(path=str(root / 'src/queue_interventions.py'),
            bytes=(root / 'src/queue_interventions.py').stat().st_size,
            function='source_checkpoint(root,cid) / complete_checkpoint(path)',
            policy='V2 preferred complete512; latest complete numeric<=512; otherwise earliest complete numeric>512; metadata only'),
        source_checkpoint_selections=selections, expected_slots=expected, files=list(rows.values()),
        excluded_paths_not_opened=exclusions, problems=problems,
        main_view_transform=dict(retain=['m', 'Y_prediction'],
            omitted='j/Q and every other original field; original full state remains in package_raw release',
            equality='np.array_equal original versus exported m and Y_prediction, with exact shape and dtype'),
        omitted_scope=['Other main NATIVE_HISTORY except focus4/common5', 'Main checkpoints/final states except fixed cycle/source selections',
            'RUNNING queue STATUS snapshots', 'Public/setup inputs already present locally',
            'All audit/analysis/oldendpoint/truth/private trees', 'Root startup ZIPs and complete raw release'],
        raw_originals_unchanged=True, GT_access='NONE', physics=0, optimization=0, hash_operations='NONE',
        publication='NOT_PERFORMED', complete_raw='Separate package_raw archives retain original j/Q, all checkpoints/caches and full main native history')
    try:
        import numpy as np
        with zipfile.ZipFile(destinations['zip'], 'x', allowZip64=True) as bundle:
            for row in rows.values():
                path = root / row['path']
                st = path.stat()
                if (st.st_size, st.st_mtime_ns) != (row['source_bytes'], row['source_mtime_ns']):
                    raise RuntimeError('Source changed before review export: ' + row['path'])
                if row['transform'] == 'RETAIN_REAL_M_AND_Y_PREDICTION_ONLY':
                    try:
                        with np.load(path, allow_pickle=False) as original:
                            row['original_fields'] = list(original.files)
                            if not {'m', 'Y_prediction'} <= set(original.files):
                                raise ValueError('Actual view lacks required m/Y_prediction; no invented Y or m-only fallback')
                            values = {key: original[key].copy() for key in ('m', 'Y_prediction')}
                        if np.iscomplexobj(values['m']):
                            raise ValueError('Original material m is complex; no implicit real cast')
                        temporary = io.BytesIO()
                        np.savez_compressed(temporary, **values)
                        payload = temporary.getvalue()
                        bundle.writestr(row['path'], payload, compress_type=zipfile.ZIP_STORED)
                        row.update(status='WRITTEN_TRANSFORMED_PENDING_READBACK', exported_bytes=len(payload),
                            retained_fields=list(values), dropped_fields=[k for k in row['original_fields'] if k not in values])
                    except Exception as exc:
                        row.update(status='FAILED_VIEW_RETAINED_IN_MANIFEST', error=repr(exc))
                        problems.append(dict(path=row['path'], issue='VIEW_EXPORT_ERROR', error=repr(exc)))
                        continue
                else:
                    compression = zipfile.ZIP_DEFLATED if path.suffix.lower() in TEXT_SUFFIXES else zipfile.ZIP_STORED
                    bundle.write(path, row['path'], compress_type=compression,
                                 compresslevel=1 if compression == zipfile.ZIP_DEFLATED else None)
                    row.update(status='WRITTEN_ORIGINAL_BYTES', exported_bytes=row['source_bytes'])
                after = path.stat()
                if (after.st_size, after.st_mtime_ns) != (row['source_bytes'], row['source_mtime_ns']):
                    raise RuntimeError('Source changed while exporting: ' + row['path'])
                exported.append(row['path'])
        # Read back the actual emitted member, not merely the temporary buffer.
        with zipfile.ZipFile(destinations['zip'], 'r') as bundle:
            infos = bundle.infolist()
            if [item.filename for item in infos] != exported:
                raise ValueError('Export ZIP members differ from completed manifest rows')
            for info in infos:
                row = rows[info.filename]
                if info.file_size != row['exported_bytes']:
                    raise ValueError('Exported member size differs from manifest')
                if row['transform'] == 'RETAIN_REAL_M_AND_Y_PREDICTION_ONLY':
                    with np.load(root / row['path'], allow_pickle=False) as original, \
                         np.load(io.BytesIO(bundle.read(info)), allow_pickle=False) as review:
                        checks = {key: bool(np.array_equal(original[key], review[key]) and
                            original[key].shape == review[key].shape and original[key].dtype == review[key].dtype)
                            for key in ('m', 'Y_prediction')}
                        if set(review.files) != {'m', 'Y_prediction'} or not all(checks.values()):
                            raise ValueError('Actual exported m/Y array_equal/shape/dtype check failed')
                    row.update(status='WRITTEN_TRANSFORMED_ARRAY_EQUAL_CHECKED', field_array_equal=checks)
                    array_checks.append(dict(path=row['path'], **checks))
                row['zip_compressed_bytes'] = info.compress_size
                row['outer_zip_compression'] = 'ZIP_STORED' if info.compress_type == zipfile.ZIP_STORED else 'ZIP_DEFLATED_LEVEL1'
        # A terminal source cannot silently change while the review ZIP is made.
        for row in rows.values():
            st = (root / row['path']).stat()
            if (st.st_size, st.st_mtime_ns) != (row['source_bytes'], row['source_mtime_ns']):
                raise RuntimeError('Original file changed during export readback')
        manifest.update(status='COMPLETE_WITH_RETAINED_ISSUES' if problems else 'COMPLETE',
                        actual_archive_bytes=destinations['zip'].stat().st_size,
                        exported_file_count=len(exported), transformed_view_array_checks=array_checks)
        write_new(destinations['manifest'], manifest)
        with destinations['csv'].open('x', newline='', encoding='utf-8') as handle:
            writer = csv.DictWriter(handle, fieldnames=['path', 'source_filename', 'source_bytes', 'exported_bytes',
                'roles', 'transform', 'status', 'retained_fields', 'dropped_fields', 'field_array_equal'])
            writer.writeheader()
            for row in rows.values():
                writer.writerow({k: json.dumps(row.get(k), ensure_ascii=False) if isinstance(row.get(k), (dict, list)) else row.get(k)
                                 for k in writer.fieldnames})
        with destinations['readme'].open('x', encoding='utf-8') as handle:
            handle.write('# R15 finite review copy — not the complete raw release\n\n'
                'This ZIP preserves original R15-relative file paths. Extract into the local review tree with its existing public inputs, '
                'source/runtime and operator-rank basis caches. It does not itself include those setup/public files.\n\n'
                'Main264 candidate pools and original JSON/JSONL metadata are retained. Five-view NPZ files contain only the '
                'original real m and original Y_prediction; j/Q and other fields are deliberately omitted. Exact array_equal, shape '
                'and dtype checks compare the original arrays with the actual exported ZIP members. These are IO checks, not new science.\n\n'
                'All actual fork56 files are retained byte-for-byte. The common-hop main family retains first4 cycle BEFORE/AFTER '
                'full NPZs: 480 expected pairs /960 expected files. Missing states remain explicit, with no rerun. Four fixed hop '
                'tasks retain a preferred512 or actual V2 metadata-only fallback checkpoint and its original fixed/background cache '
                'dependencies. Other common-cycle caches may still require the separate raw release for solver restoration.\n\n'
                'Only focus0009/0020/0021/0023 x5 common-hop main methods retain NATIVE_HISTORY. Queue exports contain terminal/frozen '
                'resource/event/job/log records, excluding RUNNING STATUS snapshots. No audit/analysis/oldendpoint/truth tree is read.\n\n'
                'The MANIFEST and FILES CSV give original names/byte sizes, copied versus transformed roles, retained/dropped view '
                'fields, missing slots and readback checks. RECEIPT gives actual bytes and new IO/validation CPU/wall. Original inverse '
                'execution costs are reused provenance, never refunded. package_raw.py provides the complete original raw assets; '
                'do not relabel this review ZIP as a full raw archive or a startup/source package. No GT, physics, optimizer, SHA, '
                'publication or source deletion is performed.\n')
        receipt = dict(status=manifest['status'], role=manifest['role'], archive=str(destinations['zip']),
            actual_archive_bytes=manifest['actual_archive_bytes'], source_bytes_of_admitted_existing_files=sum(r['source_bytes'] for r in rows.values()),
            exported_member_bytes=sum(r['exported_bytes'] or 0 for r in rows.values()),
            exported_files=len(exported), copied_original_files=sum(r['status']=='WRITTEN_ORIGINAL_BYTES' for r in rows.values()),
            transformed_views=len(array_checks), requested_views=264*5,
            expected_cycle_pairs=480, expected_cycle_file_slots=960,
            actual_cycle_files=sum('COMMON_HOP_FIRST4_FULL_CYCLE_STATE' in r['roles'] and r['path'] in exported for r in rows.values()),
            missing_slots=sum(not x['exists'] for x in expected), issues=len(problems), failures=problems,
            CPU_seconds=time.process_time()-start_cpu, whole_process_CPU_at_receipt_seconds=time.process_time(),
            wall_seconds=time.perf_counter()-start_wall, original_raw_retained=True,
            complete_raw_release=False, physics=0, optimization=0, GT_access='NONE', hash_operations='NONE',
            original_execution_cost='Reused in queue_info; not recharged or refunded; this receipt pays only export IO/validation CPU',
            source_arrays_modified=False, publication='NOT_PERFORMED')
        write_new(destinations['receipt'], receipt)
        print(json.dumps({k: receipt[k] for k in ('status', 'actual_archive_bytes', 'exported_files', 'transformed_views', 'actual_cycle_files', 'issues', 'CPU_seconds', 'wall_seconds')}))
    except BaseException as exc:
        manifest.update(status='FAILED_INCOMPLETE_REVIEW_EXPORT', error=repr(exc), written_members=exported,
                        actual_archive_bytes=destinations['zip'].stat().st_size if destinations['zip'].is_file() else 0)
        if not destinations['manifest'].exists():
            write_new(destinations['manifest'], manifest)
        write_new(destinations['error'], dict(status='FAILED_RETAINED_ALL_PARTIAL_OUTPUTS', error=repr(exc),
            traceback=traceback.format_exc(), written_members=exported, CPU_seconds=time.process_time()-start_cpu,
            whole_process_CPU_at_error_seconds=time.process_time(), wall_seconds=time.perf_counter()-start_wall,
            original_raw_retained=True, physics=0, optimization=0, GT_access='NONE', hash_operations='NONE'))
        raise


if __name__ == '__main__':
    main()
