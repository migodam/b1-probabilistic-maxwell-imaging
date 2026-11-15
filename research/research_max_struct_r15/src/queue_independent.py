"""Finite foreground case-group queue; no retries, truth, daemon or other-job changes."""
from __future__ import annotations
import os
THREAD_KEYS = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
               'VECLIB_MAXIMUM_THREADS', 'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS')
for key in THREAD_KEYS:
    os.environ[key] = '1'
import argparse
from collections import OrderedDict, deque
from datetime import datetime, timezone
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import time
import traceback


def utc():
    return datetime.now(timezone.utc).isoformat()


def write_new(path, value):
    with path.open('x', encoding='utf-8') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')


def case_groups(data):
    groups, ids = OrderedDict(), set()
    campaign = data['campaign']
    if not campaign or Path(campaign).name != campaign or any(c in campaign for c in '/\\') or campaign in ('.', '..'):
        raise ValueError('Campaign must be a single safe directory name')
    for task in data['tasks']:
        ident, case = task['task_id'], task['case_id']
        if not ident or Path(ident).name != ident or any(c in ident for c in '/\\') or ident in ('.', '..'):
            raise ValueError('Task ID must be a single safe filename')
        if ident in ids:
            raise ValueError('Duplicate task ID: ' + ident)
        if task['campaign'] != campaign:
            raise ValueError('Task campaign differs from manifest campaign')
        ids.add(ident)
        groups.setdefault(case, []).append(task)
    if not groups:
        raise ValueError('Manifest contains no tasks')
    return [(case, deque(tasks)) for case, tasks in groups.items()]


def parse_cores(value):
    cores = [int(v.strip()) for v in value.split(',')]
    if not cores or len(set(cores)) != len(cores) or min(cores) < 0:
        raise ValueError('Cores must be distinct nonnegative logical CPU IDs')
    return cores


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--legacy-repo', type=Path, required=True)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--cores', default='0,2,4,6,8,10,12,14,16,18,20,22')
    parser.add_argument('--core-type-labels', default='', help='Optional JSON logicalCPU-to-type map from the host resource audit; unspecified types remain UNVERIFIED')
    args = parser.parse_args()
    root, legacy, manifest = args.root.resolve(), args.legacy_repo.resolve(), args.manifest.resolve()
    data = json.loads(manifest.read_text())
    pending = deque(case_groups(data))
    cores = parse_cores(args.cores)
    core_labels = json.loads(args.core_type_labels) if args.core_type_labels else {}
    core_types = {str(c): core_labels.get(str(c), 'UNVERIFIED_HOST_CORE_TYPE') for c in cores}
    runner = root / 'src/run_independent.py'
    cpu_source = legacy / 'research/research_max_struct_r11/experiments/twod_queue.py'
    if not runner.is_file() or not cpu_source.is_file():
        raise FileNotFoundError('Declared runner and inherited whole-child CPU helper are required')
    # Refuse historical output collisions before dispatching even the first task.
    for task in data['tasks']:
        path = root / 'RESULTS' / task['campaign'] / task['task_id']
        if path.exists():
            raise FileExistsError('Existing prediction output is never overwritten: ' + str(path))
    import psutil
    proc = psutil.Process()
    if not hasattr(proc, 'cpu_affinity'):
        raise RuntimeError('This host does not expose process CPU affinity; pinned campaign cannot run here')
    available = proc.cpu_affinity()
    if not set(cores).issubset(set(available)):
        raise ValueError('Selected logical CPUs are outside the controller allowed affinity')
    if os.name == 'nt':
        proc.nice(psutil.NORMAL_PRIORITY_CLASS)
    spec = importlib.util.spec_from_file_location('r14_inherited_external_cpu', cpu_source)
    old = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(old)
    q = root / 'queues' / data['campaign']
    q.mkdir(parents=True, exist_ok=False)
    (q / 'logs').mkdir()
    (q / 'TERMINAL_JOBS').mkdir()
    write_new(q / 'FROZEN_MANIFEST.json', data)
    # Children use this immutable manifest, not a subsequently edited source file.
    frozen_manifest = q / 'FROZEN_MANIFEST.json'
    others = []
    for other in psutil.process_iter(['pid', 'name']):
        try:
            if other.pid != os.getpid() and 'python' in (other.info['name'] or '').lower():
                others.append(dict(pid=other.pid, name=other.info['name'], affinity=other.cpu_affinity()))
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            continue
    write_new(q / 'RESOURCE_PROBE.json', dict(utc=utc(), selected_cores=cores,
        allowed_cores=available, core_type_by_logical=core_types, core_type_source='Host-supplied audit labels; unspecified types unverified', no_absolute_CPU_crosscohort_comparison=True, CPU_percent_by_logical=psutil.cpu_percent(interval=1, percpu=True),
        other_python=others, priority='NORMAL', threads_per_worker=1,
        exclusive_CPU_claim=False, other_processes_modified=False,
        external_cpu_helper=str(cpu_source)))
    env = os.environ.copy()
    env['PYTHONPATH'] = os.pathsep.join([str(legacy), str(legacy / 'vendor')])
    env['PYTHONUNBUFFERED'] = '1'
    lanes = {core: dict(case_id=None, tasks=deque()) for core in cores}
    running, jobs = {}, []
    start, controller_cpu_start = time.perf_counter(), time.process_time()
    sequence, last_signature = 0, None
    events = (q / 'EVENTS.jsonl').open('x', encoding='utf-8')

    def emit(row):
        row = dict(utc=utc(), **row)
        events.write(json.dumps(row, allow_nan=False) + '\n')
        events.flush()
        print(json.dumps(row), flush=True)

    def save(status):
        nonlocal sequence, last_signature
        signature = (status, len(pending), len(jobs),
                     tuple((p.pid, j['task']['task_id']) for p, j in running.items()),
                     tuple((c, l['case_id'], len(l['tasks'])) for c, l in lanes.items()))
        if signature == last_signature:
            return
        last_signature = signature
        values = [j['external'].get('child_process_CPU_seconds') for j in jobs]
        complete_cpu = all(v is not None for v in values)
        record = dict(status=status, campaign=data['campaign'], started_from_manifest=str(manifest),
            frozen_manifest=str(frozen_manifest), requested_tasks=len(data['tasks']),
            requested_cases=len(set(t['case_id'] for t in data['tasks'])),
            pending_cases=len(pending),
            pending_tasks=sum(len(ts) for _, ts in pending) + sum(len(l['tasks']) for l in lanes.values()),
            running=[dict(task_id=j['task']['task_id'], case_id=j['task']['case_id'],
                          pid=p.pid, cpu=j['cpu']) for p, j in running.items()],
            lanes=[dict(cpu=c, case_id=l['case_id'], remaining_arms=len(l['tasks'])) for c, l in lanes.items()],
            jobs=jobs, queue_wall_seconds=time.perf_counter() - start,
            controller_CPU_seconds=time.process_time() - controller_cpu_start,
            finished_child_CPU_seconds=sum(values) if complete_cpu else None,
            measured_finished_child_CPU_seconds=sum(v for v in values if v is not None),
            external_CPU_complete=complete_cpu, missing_external_CPU_tasks=[j['task_id'] for j in jobs
                if j['external'].get('child_process_CPU_seconds') is None],
            CPU_scope='Child external CPU plus controller CPU are disjoint; wall and child internal CPU never added',
            selected_cores=cores, core_type_by_logical=core_types, no_absolute_CPU_crosscohort_comparison=True, automatic_retry=False, truth_access='NONE', scoring='NOT_RUN',
            case_group_policy='All arms of a case sequential on one core; release lane only after case completion')
        sequence += 1
        path = q / ('RECEIPT.json' if status != 'RUNNING' else f'STATUS_{sequence:06d}.json')
        write_new(path, record)

    def terminal(row):
        write_new(q / 'TERMINAL_JOBS' / (row['task_id'] + '.json'), row)
        jobs.append(row)
        emit(dict(event='END', **row))

    while pending or running or any(l['tasks'] for l in lanes.values()):
        occupied = {j['cpu'] for j in running.values()}
        for core, lane in lanes.items():
            if core in occupied:
                continue
            if not lane['tasks']:
                lane['case_id'] = None
                if not pending:
                    continue
                lane['case_id'], lane['tasks'] = pending.popleft()
                emit(dict(event='CASE_ASSIGNED', case_id=lane['case_id'], cpu=core,
                          task_ids=[t['task_id'] for t in lane['tasks']]))
            task = lane['tasks'].popleft()
            argv = [str(runner), '--root', str(root), '--legacy-repo', str(legacy),
                    '--manifest', str(frozen_manifest), '--task-id', task['task_id'], '--cpu', str(core)]
            # Pin and set NORMAL before importing the runner/numerical libraries.
            boot = ('import os,sys,runpy,psutil;'
                    f'psutil.Process().cpu_affinity([{core}]);'
                    'psutil.Process().nice(psutil.NORMAL_PRIORITY_CLASS) if os.name=="nt" else None;'
                    f'sys.argv={argv!r};runpy.run_path({str(runner)!r},run_name="__main__")')
            log_path = q / 'logs' / (task['task_id'] + '.log')
            log = log_path.open('x', encoding='utf-8')
            try:
                process = subprocess.Popen([sys.executable, '-c', boot], env=env, cwd=root,
                                           stdout=log, stderr=subprocess.STDOUT)
            except Exception as exc:
                log.write(traceback.format_exc())
                log.close()
                terminal(dict(task_id=task['task_id'], case_id=task['case_id'],
                    method=task.get('method'), cpu=core, status='LAUNCH_FAILED', returncode=None,
                    external=dict(child_process_CPU_seconds=0., CPU_scope='No child was created'),
                    wall_seconds=0., child_status=None, error=repr(exc), log=str(log_path)))
                continue
            running[process] = dict(task=task, cpu=core, log=log, start=time.perf_counter(), log_path=log_path)
            emit(dict(event='START', task_id=task['task_id'], case_id=task['case_id'],
                      pid=process.pid, cpu=core))
        for process, run in list(running.items()):
            code = process.poll()
            if code is None:
                continue
            try:
                external = old.external_cpu(process)
            except Exception as exc:
                external = dict(child_process_CPU_seconds=None, CPU_scope='EXTERNAL_CPU_CAPTURE_FAILED',
                                error=repr(exc))
            run['log'].close()
            task = run['task']
            rpath = root / 'RESULTS' / task['campaign'] / task['task_id'] / 'RECEIPT.json'
            child, error = None, None
            try:
                if rpath.is_file():
                    child = json.loads(rpath.read_text())
            except Exception as exc:
                error = repr(exc)
            good = (code == 0 and child is not None and child.get('status') not in ('FAILED', None)
                    and (os.name != 'nt' or external.get('child_process_CPU_seconds') is not None))
            row = dict(task_id=task['task_id'], case_id=task['case_id'], method=task.get('method'),
                cpu=run['cpu'], pid=process.pid, status='FINISHED' if good else 'FAILED',
                returncode=code, external=external, wall_seconds=time.perf_counter() - run['start'],
                child_status=child.get('status') if child else None,
                child_internal_CPU_seconds=child.get('CPU_seconds') if child else None,
                receipt=str(rpath), receipt_read_error=error, log=str(run['log_path']))
            terminal(row)
            del running[process]
        save('RUNNING')
        if running:
            time.sleep(.5)
    final = 'COMPLETED' if all(j['status'] == 'FINISHED' for j in jobs) else 'COMPLETED_WITH_ERRORS'
    save(final)
    events.close()
    return 0 if final == 'COMPLETED' else 1


if __name__ == '__main__':
    raise SystemExit(main())
