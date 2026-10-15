"""Independent G accounting. No solver import and no historical balance reuse."""
from collections import Counter
from contextlib import contextmanager
from pathlib import Path
import fcntl
import json
import os
import time

ROOT = Path(__file__).resolve().parents[1]


class BudgetStop(RuntimeError):
    pass


def plain(value):
    if hasattr(value, 'tolist'):
        return value.tolist()
    if hasattr(value, 'item'):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(type(value).__name__)


def dump(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, default=plain, ensure_ascii=False,
                               allow_nan=False, indent=2) + '\n')


class GBudget:
    """Reserve states/RHS independently before kernels; failures stay charged.

    The A30 adapter supplies actual RHS counts, including adjoints and failed
    attempts. Spans are nested timing diagnostics, never summed as job wall.
    """
    def __init__(self, root=ROOT):
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)
        self.ledger = self.root / 'COST_LEDGER.jsonl'
        self.failures = self.root / 'FAILURE_LEDGER.jsonl'
        self.counts = Counter()
        self.stage = None
        self.current_job = None
        self.current_arm = None
        self.arm_limits = None
        self.started = None
        self.limits = dict(states=384, rhs=12288, wall=10800.)

    def rows(self):
        if not self.ledger.exists():
            return []
        return [json.loads(s) for s in self.ledger.read_text().splitlines()]

    def append(self, row):
        value = dict(row, utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()))
        value.setdefault('job', self.current_job)
        value.setdefault('stage', self.stage)
        value.setdefault('arm', self.current_arm)
        with self.ledger.open('a') as stream:
            stream.write(json.dumps(value, default=plain, allow_nan=False) + '\n')
            stream.flush()
            os.fsync(stream.fileno())

    def failure(self, reason, **metadata):
        value = dict(reason=reason, **metadata, job=self.current_job,
                     arm=self.current_arm, stage=self.stage)
        with self.failures.open('a') as stream:
            stream.write(json.dumps(value, default=plain, allow_nan=False) + '\n')
            stream.flush()
            os.fsync(stream.fileno())

    def totals(self, arm=None):
        rows = self.rows()
        if arm is not None:
            rows = [r for r in rows if r.get('arm') == arm]
        return {key: sum(r['amount'] for r in rows if r.get('event') == event)
                for key, event in [('states', 'STATE_RESERVED'), ('rhs', 'RHS_RESERVED')]}

    def check(self):
        prior = sum(r.get('wall_seconds', 0.) for r in self.rows()
                    if r.get('event') in ('JOB_END', 'ZERO_PHYSICS_ANALYSIS'))
        elapsed = time.perf_counter() - self.started if self.started is not None else 0.
        if prior + elapsed >= self.limits['wall']:
            raise BudgetStop('G_CONTROLLED_WALL_EXHAUSTED')

    def _reserve(self, key, event, amount, **metadata):
        self.check()
        if self.stage != 'E1' or self.current_job is None:
            raise BudgetStop('NEW_PHYSICS_FORBIDDEN_OUTSIDE_G_E1_JOB')
        if not isinstance(amount, int) or isinstance(amount, bool) or amount < 1:
            raise ValueError('POSITIVE_INTEGER_RESERVATION_REQUIRED')
        with (self.root / 'accounting.lock').open('a') as stream:
            fcntl.flock(stream, fcntl.LOCK_EX)
            try:
                global_count = self.totals()[key]
                if global_count + amount > self.limits[key]:
                    raise BudgetStop('G_GLOBAL_' + key.upper() + '_CAP')
                if self.current_arm is not None:
                    local_count = self.totals(self.current_arm)[key]
                    if local_count + amount > self.arm_limits[key]:
                        raise BudgetStop('G_ARM_' + key.upper() + '_CAP')
                self.append(dict(event=event, amount=amount, cumulative=global_count + amount,
                                 **metadata))
            finally:
                fcntl.flock(stream, fcntl.LOCK_UN)

    def reserve_state(self, material, *, resolution, purpose):
        self._reserve('states', 'STATE_RESERVED', 1, resolution=resolution,
                      purpose=purpose, material_shape=list(material.shape))

    def reserve_rhs(self, count, role):
        if isinstance(count, bool) or int(count) != count:
            raise ValueError('INTEGER_RHS_COUNT_REQUIRED')
        self._reserve('rhs', 'RHS_RESERVED', int(count), role=role)

    @contextmanager
    def span(self, label, **counts):
        self.check()
        self.counts.update(counts)
        wall, cpu = time.perf_counter(), time.process_time()
        status = 'FAIL'
        try:
            yield
            status = 'PASS'
        finally:
            self.append(dict(event='ACTION', label=label, status=status, counts=counts,
                             wall_seconds=time.perf_counter() - wall,
                             cpu_seconds=time.process_time() - cpu))

    @contextmanager
    def arm(self, round_id, case, method, states=12, rhs=288):
        if self.current_arm is not None:
            raise ValueError('NO_NESTED_ARM')
        self.current_arm = f'{round_id}/{case}/{method}'
        self.arm_limits = dict(states=min(12, int(states)), rhs=min(288, int(rhs)))
        try:
            yield self
        finally:
            self.current_arm = None
            self.arm_limits = None

    @contextmanager
    def job(self, name, stage='E1'):
        rows = self.rows()
        if self.current_job is not None or any(r.get('job') == name for r in rows):
            raise ValueError('JOB_ID_ALREADY_USED_OR_NESTED')
        begin = [r for r in rows if r.get('event') == 'JOB_BEGIN']
        end = [r for r in rows if r.get('event') == 'JOB_END']
        if len(begin) != len(end):
            raise BudgetStop('UNCLOSED_PREVIOUS_G_JOB_REQUIRES_EXPLICIT_REPAIR')
        if stage == 'E1':
            config = json.loads((self.root / 'E1_PROTOCOL_FROZEN.json').read_text())
            if config.get('execution_authorized') is not True or config.get('backend_ready') is not True:
                raise BudgetStop('FROZEN_IDENTICAL_PHYSICS_AND_LEGAL_BASE_REQUIRED')
            # This checks the registered experiment identity, not its scientific
            # correctness. The unchanged adapter then checks actual arrays and
            # the runner performs independently charged base validation.
            expected = dict(
                target_id='A30_24x24_TM_8src32rx_v1',
                A30_commit='717d3d26a290248bd2b9aeb55f3d03694b407b84',
                corrector_commit='f1b674eb4893e39126701ed48f3d3dcef3098c3c',
                corrector_id='ResearchA_E2_32realHaar_GN2_v1',
                packing='source-major-Re32-Im32', gain=1,
                whitening='sqrt(2)/original_sigma')
            backend = config.get('backend', {})
            if any(backend.get(k) != v for k, v in expected.items()):
                raise BudgetStop('G_INCOMPLETE_OR_DIFFERENT_FROZEN_BACKEND_METADATA')
            if abs(float(backend.get('data_threshold', -1))-23.785449004057032) > 1e-12:
                raise BudgetStop('G_DIFFERENT_DATA_THRESHOLD')
            cases = config.get('cases', [])
            if (not cases or len(cases) != len(set(cases)) or len(cases) > 4
                    or any(c not in {'p250069_weak_a1_snr20', 'p250073_weak_a1_snr20',
                                     'p250077_weak_a1_snr20'} for c in cases)):
                raise BudgetStop('G_REGISTERED_QUALIFIED_BASE_IDS_REQUIRED')
        self.current_job, self.stage = name, stage
        self.counts = Counter()
        self.started = time.perf_counter()
        cpu, initial = time.process_time(), self.totals()
        status, error = 'PASS', None
        try:
            try:
                os.nice(10)
            except OSError:
                pass
            self.append(dict(event='JOB_BEGIN'))
            self.check()
            yield self
        except BaseException as exc:
            status, error = 'FAIL', repr(exc)
            self.failure('JOB_FAILURE', error=error)
            raise
        finally:
            totals = self.totals()
            self.append(dict(event='JOB_END', status=status, error=error,
                             wall_seconds=time.perf_counter() - self.started,
                             cpu_seconds=time.process_time() - cpu,
                             new_states=totals['states'] - initial['states'],
                             new_rhs=totals['rhs'] - initial['rhs'],
                             aggregate_counts=dict(self.counts), GPU=0))
            self.started = None
            self.current_job, self.stage = None, None
            dump(self.root / 'RESOURCE_SUMMARY.json', self.receipt())

    def receipt(self):
        rows = self.rows()
        total = self.totals()
        jobs = [r for r in rows if r.get('event') in ('JOB_END', 'ZERO_PHYSICS_ANALYSIS')]
        action_counts = Counter()
        for row in jobs:
            action_counts.update(row.get('aggregate_counts', {}))
        return dict(new_full_state_attempts=total['states'], new_complex_full_solver_RHS=total['rhs'],
                    state_cap=384, RHS_cap=12288, controlled_wall_cap=10800,
                    controlled_wall_seconds=sum(r.get('wall_seconds', 0.) for r in jobs),
                    CPU_seconds=sum(r.get('cpu_seconds', 0.) for r in jobs),
                    failed_jobs=sum(r.get('status') == 'FAIL' for r in jobs),
                    action_counts=dict(action_counts),
                    GPU=0, T_states_borrowed=0, historical_ledgers_modified=False,
                    time_scope='Only JOB_END/unique zero-physics receipts; nested ACTION timings not summed')
