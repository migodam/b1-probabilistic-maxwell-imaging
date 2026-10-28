"""R7's independent, fail-closed paid physics adapter.

Coordinates passed to evaluate and returned Jacobians are Gaussianized x.
Importing this module loads the frozen prior, never a Maxwell state or solver.
"""
from __future__ import annotations
import os
for _name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
              'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS', 'BLIS_NUM_THREADS'):
    os.environ[_name] = '1'
import fcntl
import copy
import hashlib
import json
import math
import resource
import sys
import time
import uuid
from collections import Counter, OrderedDict
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
import numpy as np

REPO = Path(__file__).resolve().parents[3]
OUT = REPO / 'research/p3l_r7_material_solution_fibers'
R4 = REPO.parent / 'b1-p3l-r4-conditional-transport'
A30 = REPO.parents[1] / 'A30/a30-maxwell-material-manifold'
sys.path.insert(0, str(R4 / 'src'))
from p3r4.core import (Target as OriginalTarget, BudgetStop, pack, plain,
                       from_probit, to_probit, probit_derivative)
for _module, _expected in (('p3r4.core', R4/'src/p3r4/core.py'),
                            ('p1_material.prior', R4/'src/p1_material/prior.py')):
    if Path(sys.modules[_module].__file__).resolve() != _expected.resolve():
        raise ImportError('R7 nonfrozen source collision: '+_module)

CAMPAIGN = 'P3LR7-20261028'
CAPS = {'E0': 0, 'E1': 0, 'E2': 50000}
SIGMA = 0.00119211798514519


def _atomic(path, writer):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + '.' + uuid.uuid4().hex + '.tmp')
    try:
        with temp.open('wb') as stream:
            writer(stream)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temp, path)
        fd = os.open(path.parent, os.O_RDONLY)
        try: os.fsync(fd)
        finally: os.close(fd)
    finally:
        if temp.exists(): temp.unlink()


def write(path, obj):
    payload = (json.dumps(plain(obj), allow_nan=False, ensure_ascii=False,
                          indent=2) + '\n').encode()
    _atomic(path, lambda stream: stream.write(payload))


def save(path, **arrays):
    _atomic(path, lambda stream: np.savez_compressed(stream, **arrays))


def _empty():
    return dict(campaign_id=CAMPAIGN, sequence=0, RHS=0, forward_RHS=0,
                derivative_RHS=0, states=0, LU=0, failures=0,
                physics_wall_seconds=0., controlled_wall_seconds=0., CPU_seconds=0.,
                stages={s: dict(RHS=0, states=0, LU=0) for s in CAPS},
                active_actions={}, active_jobs={})


def _ledger(root, *, require_snapshot=True):
    """Validate every authoritative journal row; never silently repair damage."""
    root = Path(root)
    path, snapshot = root/'RESOURCE_LEDGER.jsonl', root/'RESOURCE_LEDGER_SNAPSHOT.json'
    previous = _empty()
    if path.exists():
        payload = path.read_bytes()
        if payload and not payload.endswith(b'\n'):
            raise RuntimeError('incomplete R7 journal: explicit reconciliation required')
        for raw in payload.splitlines():
            try: row = json.loads(raw)
            except Exception as exc: raise RuntimeError('invalid R7 journal') from exc
            after = row.get('after', {})
            if (row.get('campaign_id') != CAMPAIGN or after.get('campaign_id') != CAMPAIGN
                    or row.get('sequence') != previous['sequence'] + 1
                    or after.get('sequence') != row['sequence']):
                raise RuntimeError('foreign or discontinuous R7 journal')
            stage = row.get('stage')
            if stage not in CAPS: raise RuntimeError('invalid R7 ledger stage')
            for key in ('RHS', 'states', 'LU'):
                delta = row.get(key, 0)
                if type(delta) is not int or delta < 0 or after[key] != previous[key] + delta:
                    raise RuntimeError('invalid R7 reservation accounting')
                for s in CAPS:
                    if after['stages'][s][key] != previous['stages'][s][key] + (delta if s == stage else 0):
                        raise RuntimeError('invalid R7 stage accounting')
            if after['RHS'] > 50000 or any(after['stages'][s]['RHS'] > CAPS[s] for s in CAPS):
                raise RuntimeError('R7 historical cap violation')
            for key, kind in (('forward_RHS', 'forward'), ('derivative_RHS', 'derivative')):
                if after[key] != previous[key] + (row['RHS'] if row.get('kind') == kind else 0):
                    raise RuntimeError('invalid R7 RHS kind accounting')
            if after['RHS'] != after['forward_RHS'] + after['derivative_RHS']:
                raise RuntimeError('unclassified R7 RHS')
            previous = after
    if snapshot.exists():
        try: stored = json.loads(snapshot.read_text())
        except Exception as exc: raise RuntimeError('invalid R7 snapshot') from exc
        if stored != previous: raise RuntimeError('R7 journal/snapshot mismatch: reconcile explicitly')
    elif require_snapshot and previous['sequence']:
        raise RuntimeError('missing R7 snapshot: reconcile explicitly')
    return previous


def totals(root=OUT):
    return _ledger(root)


def reconcile_snapshot(root=OUT):
    """Rebuild only the R7 snapshot from its valid journal; preserve all charges.

    This does not clear interrupted jobs/actions or permit automatic replay.
    A damaged or foreign journal still raises; never import an R4/R6 snapshot.
    """
    root = Path(root)
    if root.resolve() != OUT.resolve():
        raise ValueError('reconciliation is restricted to the existing R7 campaign')
    root.mkdir(parents=True, exist_ok=True)
    with (root/'physics_worker.lock').open('a') as worker:
        fcntl.flock(worker, fcntl.LOCK_EX | fcntl.LOCK_NB)
        with (root/'ledger.lock').open('a') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            # Deliberately exclude an untrusted/mismatching snapshot from recovery.
            snapshot = root/'RESOURCE_LEDGER_SNAPSHOT.json'
            backup = None
            if snapshot.exists():
                backup = snapshot.with_name(snapshot.name + '.reconcile-' + uuid.uuid4().hex)
                os.replace(snapshot, backup)
            try:
                state = _ledger(root, require_snapshot=False)
                write(snapshot, state)
            except BaseException:
                if backup is not None: os.replace(backup, snapshot)
                raise
            return state


def reconcile_interrupted(job, reason, root=OUT):
    """Explicit dead-job audit, preserving RHS/LU and charging reserved wall.

    The worker PID must be dead. This adds an append-only receipt and cannot
    repair a damaged journal. Call reconcile_snapshot first after a torn snapshot.
    """
    root = Path(root)
    if root.resolve() != OUT.resolve() or not str(reason).strip():
        raise ValueError('existing R7 root and explicit audit reason required')
    with (root/'physics_worker.lock').open('a') as worker:
        fcntl.flock(worker, fcntl.LOCK_EX | fcntl.LOCK_NB)
        with (root/'ledger.lock').open('a') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            state = totals(root)
            meta = state['active_jobs'].get(job)
            if meta is None: raise ValueError('no active R7 job to reconcile')
            try: os.kill(int(meta['pid']), 0)
            except ProcessLookupError: pass
            else: raise RuntimeError('recorded worker PID still exists; fail closed')
            actions = [key for key, value in state['active_actions'].items() if value['job'] == job]
            for key in actions: del state['active_actions'][key]
            del state['active_jobs'][job]
            state['failures'] += max(1, len(actions))
            state['controlled_wall_seconds'] += float(meta['max_wall'])
            state['sequence'] += 1
            row = dict(campaign_id=CAMPAIGN, sequence=state['sequence'],
                       event='RECONCILE_INTERRUPTED', stage=meta['stage'], job=job,
                       utc=time.time(), RHS=0, states=0, LU=0, kind=None,
                       reason=str(reason), abandoned_actions=actions,
                       charged_reserved_wall_seconds=meta['max_wall'],
                       CPU_seconds_status='UNAVAILABLE', after=state)
            with (root/'RESOURCE_LEDGER.jsonl').open('ab') as stream:
                stream.write((json.dumps(row, allow_nan=False, separators=(',', ':'))+'\n').encode())
                stream.flush(); os.fsync(stream.fileno())
            write(root/'RESOURCE_LEDGER_SNAPSHOT.json', state)
            write(root/'raw/jobs'/job/'RECONCILED.json', row)
            return state


class Book:
    """Reserve before every LU/state/RHS action. Failures never refund."""
    def __init__(self, stage, job, parent=None, role='E2_PHYSICS', max_wall=7200,
                 *, root=OUT, test_mode=False):
        if stage not in CAPS: raise ValueError(stage)
        if stage != 'E2': raise BudgetStop('E0/E1 have zero physical RHS budget')
        self.root = Path(root).resolve()
        self.test_mode = bool(test_mode)
        if self.root != OUT.resolve() and not self.test_mode:
            raise ValueError('production ledger must be the R7 ledger')
        if self.test_mode and self.root == OUT.resolve():
            raise ValueError('fake tests must use a separate temporary ledger')
        self.root.mkdir(parents=True, exist_ok=True)
        self.workerlock = (self.root/'physics_worker.lock').open('a')
        try:
            fcntl.flock(self.workerlock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            self._authorize()
            if (self.root/'CLOSURE.json').exists(): raise BudgetStop('R7 closed')
            self.stage, self.job = stage, str(job)
            if not self.job or Path(self.job).name != self.job or self.job in ('.', '..'):
                raise ValueError('job must be one safe path component')
            self.max_wall = float(max_wall)
            if not math.isfinite(self.max_wall) or self.max_wall <= 0: raise ValueError('max_wall')
            self.context = dict(parent=parent, role=role, synthetic=self.test_mode)
            self.counts, self.credit, self.active = Counter(), 0, None
            self.finished = False
            self.started, self.cpu0 = time.perf_counter(), time.process_time()
            self.jobdir = self.root/'raw/jobs'/self.job
            if self.jobdir.exists(): raise RuntimeError('existing job: no automatic replay')
            state = totals(self.root)
            if state['active_actions'] or state['active_jobs']:
                raise RuntimeError('unreconciled interrupted R7 action/job; fail closed')
            self._ledger_state = state
            self._journal_identity = self._journal_stat()
            self.jobdir.mkdir(parents=True)
            self._event('JOB_BEGIN', pid=os.getpid(), max_wall=self.max_wall)
            write(self.jobdir/'START.json', dict(campaign_id=CAMPAIGN, stage=stage,
                  job=self.job, pid=os.getpid(), utc=time.time(), synthetic=self.test_mode))
        except BaseException:
            self.workerlock.close()
            raise

    def _authorize(self):
        auth = json.loads((self.root/'configs/AUTHORIZATION.json').read_text())
        freeze = json.loads((self.root/'configs/E2_FREEZE.json').read_text())
        if not isinstance(freeze, dict) or not freeze:
            raise BudgetStop('E2 must freeze before physical Book')
        if auth.get('campaign_id', CAMPAIGN) != CAMPAIGN or auth.get('explicit_user_authorization') is not True:
            raise BudgetStop('R7 explicit authorization missing')
        total = auth.get('total_complex_RHS')
        caps = auth.get('stage_caps', {})
        if (type(total) is not int or not 0 <= total <= 50000 or set(caps) != set(CAPS)
                or any(type(caps[s]) is not int or not 0 <= caps[s] <= CAPS[s] for s in CAPS)):
            raise BudgetStop('R7 authorization cannot expand frozen hard caps')
        self.auth = auth

    def _check(self):
        if self.finished: raise RuntimeError('closed R7 Book')
        if time.perf_counter() - self.started > self.max_wall:
            raise BudgetStop('R7 controlled wall cap before next action')
        # Credit-consuming callbacks also verify continuity, even when they
        # consume a reservation rather than append a new paid event.
        with (self.root/'ledger.lock').open('a') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            self._latest_live()

    def _journal_stat(self):
        path = self.root/'RESOURCE_LEDGER.jsonl'
        if not path.exists(): return None
        stat = path.stat()
        return (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns)

    def _latest_live(self):
        """Constant-sized continuity check under the ledger lock.

        Full replay occurred when this exclusive Book opened. Any outside
        journal mutation, replacement, truncation or append blocks this Book.
        The latest journal row and atomically written snapshot must both equal
        its last committed state. Public totals/reconciliation still replay all.
        """
        identity = self._journal_stat()
        if identity != self._journal_identity:
            raise RuntimeError('R7 journal changed outside the live Book; fail closed')
        state = self._ledger_state
        if state['sequence']:
            path = self.root/'RESOURCE_LEDGER.jsonl'
            with path.open('rb') as stream:
                size = identity[2]
                count = min(size, 65536)
                while True:
                    stream.seek(size-count)
                    tail = stream.read(count)
                    if not tail.endswith(b'\n'):
                        raise RuntimeError('incomplete R7 journal; fail closed')
                    start = tail.rfind(b'\n', 0, len(tail)-1)
                    if start >= 0 or count == size: break
                    count = min(size, count*2)
            try: row = json.loads(tail[start+1:])
            except Exception as exc: raise RuntimeError('invalid latest R7 journal row') from exc
            if (row.get('campaign_id') != CAMPAIGN or row.get('sequence') != state['sequence']
                    or row.get('after') != state):
                raise RuntimeError('R7 latest journal state mismatch; fail closed')
        snapshot = self.root/'RESOURCE_LEDGER_SNAPSHOT.json'
        if snapshot.exists():
            if json.loads(snapshot.read_text()) != state:
                raise RuntimeError('R7 live journal/snapshot mismatch; reconcile explicitly')
        elif state['sequence']:
            raise RuntimeError('missing R7 live snapshot; reconcile explicitly')
        return copy.deepcopy(state)

    def _event(self, event, *, rhs=0, states=0, lu=0, kind=None, **fields):
        with (self.root/'ledger.lock').open('a') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            t = self._latest_live()
            if any(type(v) is not int or v < 0 for v in (rhs, states, lu)):
                raise ValueError('R7 reservations are nonnegative integers; no refunds')
            if rhs and kind not in ('forward', 'derivative'): raise ValueError('RHS kind')
            v = t['stages'][self.stage]
            if t['RHS']+rhs > self.auth['total_complex_RHS'] or v['RHS']+rhs > self.auth['stage_caps'][self.stage]:
                raise BudgetStop('R7 50000 RHS hard cap')
            if event == 'JOB_BEGIN':
                limit = self.auth.get('controlled_CPU_wall_seconds')
                if limit is not None and t['controlled_wall_seconds']+fields['max_wall'] > float(limit):
                    raise BudgetStop('R7 controlled wall reservation cap')
                t['active_jobs'][self.job] = dict(stage=self.stage, **fields)
            elif event == 'PHYSICS_BEGIN':
                if t['active_actions']: raise RuntimeError('R7 physics action already live')
                t['active_actions'][fields['action']] = dict(job=self.job, stage=self.stage)
            elif event == 'PHYSICS_END':
                t['active_actions'].pop(fields['action'])
                t['physics_wall_seconds'] += fields['wall_seconds']
                t['failures'] += int(fields.get('error') is not None)
            elif event == 'JOB_END':
                t['active_jobs'].pop(self.job)
                t['controlled_wall_seconds'] += fields['wall_seconds']
                t['CPU_seconds'] += fields['CPU_seconds']
            for key, delta in (('RHS', rhs), ('states', states), ('LU', lu)):
                t[key] += delta
                v[key] += delta
            if kind: t[kind+'_RHS'] += rhs
            t['sequence'] += 1
            row = dict(campaign_id=CAMPAIGN, sequence=t['sequence'], event=event,
                       stage=self.stage, job=self.job, utc=time.time(), RHS=rhs,
                       states=states, LU=lu, kind=kind, context=self.context.copy(),
                       **fields, after=t)
            # Journal first, then snapshot. A crash between them blocks all work
            # until explicit snapshot reconciliation; it never loses a charge.
            with (self.root/'RESOURCE_LEDGER.jsonl').open('ab') as stream:
                stream.write((json.dumps(plain(row), allow_nan=False, separators=(',', ':'))+'\n').encode())
                stream.flush()
                os.fsync(stream.fileno())
            write(self.root/'RESOURCE_LEDGER_SNAPSHOT.json', t)
            self._ledger_state = t
            self._journal_identity = self._journal_stat()

    def reserve_state(self, material, *, resolution, purpose):
        self._check()
        if self.active is None: raise RuntimeError('state outside physical action')
        if resolution != 24: raise ValueError('only frozen TM24 is authorized')
        if self.credit: raise RuntimeError('unconsumed forward reservation')
        self._event('STATE_RESERVED', rhs=8, states=1, lu=1, kind='forward', purpose=purpose)
        self.credit += 8

    def reserve_rhs(self, count, role):
        self._check()
        if self.active is None: raise RuntimeError('solve outside physical action')
        if type(count) is not int or count <= 0: raise ValueError('RHS count')
        if role == 'full_forward':
            if self.credit != count: raise RuntimeError('missing exact forward reservation')
            self.credit -= count
            self.counts['actual_forward_RHS'] += count
        else:
            self._event('DERIVATIVE_RESERVED', rhs=count, kind='derivative', solver_role=role)
            self.counts['actual_derivative_RHS'] += count

    @contextmanager
    def span(self, label, **counts):
        start = time.perf_counter()
        self.counts.update(counts)
        try: yield
        finally: self.counts[label+'_seconds'] += time.perf_counter()-start

    @contextmanager
    def physical(self, label):
        self._check()
        if self.active is not None: raise RuntimeError('nested R7 physics')
        action = uuid.uuid4().hex
        self._event('PHYSICS_BEGIN', action=action, label=label)
        self.active = action
        start, cpu, error = time.perf_counter(), time.process_time(), None
        try: yield
        except BaseException as exc:
            error = repr(exc)
            raise
        finally:
            self._event('PHYSICS_END', action=action, label=label, error=error,
                        wall_seconds=time.perf_counter()-start, CPU_seconds=time.process_time()-cpu)
            self.active = None

    def close(self, status='COMPLETE', error=None):
        if self.finished: return
        if self.active is not None: raise RuntimeError('cannot close live physical action')
        receipt = dict(status=status, error=error, wall_seconds=time.perf_counter()-self.started,
                       CPU_seconds=time.process_time()-self.cpu0, actual_counters=dict(self.counts),
                       unused_prepaid_RHS=self.credit, synthetic=self.test_mode,
                       peak_RSS_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
                           * (1 if sys.platform == 'darwin' else 1024))
        self._event('JOB_END', **receipt)
        write(self.jobdir/'CLOSED.json', receipt)
        self.finished = True
        self.workerlock.close()

    def __enter__(self): return self
    def __exit__(self, kind, error, tb):
        self.close('FAILED' if error else 'COMPLETE', None if error is None else repr(error))


@dataclass
class EndpointState:
    """Saved endpoint arrays without an LU, matrix, or native state owner."""
    chi: np.ndarray
    current: np.ndarray
    field: np.ndarray
    data: np.ndarray
    residual: float


@dataclass
class Point:
    x: np.ndarray
    q: np.ndarray
    state: object
    key: str
    T: np.ndarray | None = None
    J: np.ndarray | None = None
    Jnative: np.ndarray | None = None


class Target(OriginalTarget):
    """Frozen P1 target, one live LU owner and at most two cached Points."""
    def __init__(self, book, *, backend_factory=None):
        if book.finished: raise RuntimeError('closed R7 Book')
        if book.test_mode != (backend_factory is not None):
            raise ValueError('fake backend required exactly for isolated test Books')
        super().__init__(book, keep_fields=True)
        if (self.sigma != SIGMA or self.cfg['gain'] != 1. or self.cfg['grid'] != 24
                or self.cfg['source_count'] != 8 or self.cfg['receiver_count'] != 32):
            raise ValueError('frozen original P1 target changed')
        self._backend_factory, self._live_point = backend_factory, None
        self.cached_points = OrderedDict()
        try:
            from threadpoolctl import threadpool_limits
            self._thread_limit = threadpool_limits(limits=1, user_api='blas')
        except ImportError:
            raise RuntimeError('threadpoolctl required to enforce one BLAS thread')

    def backend(self):
        if self.model is None:
            if self._backend_factory is not None:
                self.model = self._backend_factory(self.book)
            else:
                if os.environ.get('B1_P3L_A30_ROOT') not in (None, str(A30)):
                    raise ValueError('R7 forbids alternate A30 solver root')
                sys.path.insert(0, str(A30/'src'))
                from a30.physics import FrozenMaxwellAdapter
                self.model = FrozenMaxwellAdapter.from_background_npz(
                    R4/'research/p3l_r4_conditional_transport/inputs/background/weak.npz',
                    self.book, purpose_prefix='P3LR7_E2', device='cpu', backend='scipy')
                if Path(sys.modules['a30.physics'].__file__).resolve() != (A30/'src/a30/physics.py').resolve():
                    raise ImportError('R7 nonfrozen A30 module collision')
            if (self.model.n, self.model.q, self.model.s) != (576, 32, 8):
                raise ValueError('frozen TM24 dimensions changed')
            write(self.book.jobdir/'SOLVER_SOURCES.json', source_manifest())
        return self.model

    def _release_lu(self):
        if self.model is not None: self.model._transfers.clear()
        if self._live_point is not None:
            p = self._live_point
            s = p.state
            p.state = EndpointState(s.chi, s.current, s.field, s.data, s.residual)
            self._live_point = None

    def _persist(self, p, path=None):
        arrays = dict(x=p.x, q=p.q, material=np.asarray(p.state.chi, dtype=np.complex128),
                      chi=np.asarray(p.state.chi, dtype=np.complex128),
                      prediction=p.state.data, current=p.state.current, field=p.state.field,
                      residual=np.array(p.state.residual), backward_residual=np.array(p.state.residual))
        if p.J is not None: arrays.update(J=p.J, Jx=p.J, Jnative=p.Jnative)
        if p.T is not None: arrays['receiver_T'] = p.T
        save(path or self.book.jobdir/'endpoints'/f'{p.key}.npz', **arrays)

    def evaluate(self, x, purpose='target'):
        x = np.asarray(x, dtype=float)
        if x.shape != (66,) or not np.isfinite(x).all(): raise ValueError('finite Gaussianized x66 required')
        q = from_probit(x)
        chi, _ = self.prior.decode(q)
        self.evaluations += 1
        key = f'{self.book.job}_{self.evaluations:07d}'
        self.book.context.update(endpoint=key, purpose=purpose)
        save(self.book.jobdir/'attempts'/f'{key}.npz', x=x, q=q, material=np.asarray(chi, np.complex128))
        self._release_lu()
        try:
            with self.book.physical('full_forward_8_RHS'):
                state = self.backend().full_state(chi, purpose=purpose)
                if not np.isfinite(state.residual) or not 0 <= state.residual <= 1e-10:
                    raise RuntimeError('Maxwell backward residual contract')
            p = Point(x.copy(), q.copy(), state, key)
            self._live_point = p
            self._persist(p)
            self.cached_points[key] = p
            while len(self.cached_points) > 2: self.cached_points.popitem(last=False)
            return p
        except BaseException as exc:
            write(self.book.jobdir/'failures'/f'{key}.json', dict(endpoint=key, error=repr(exc), purpose=purpose))
            raise

    def evaluate_native(self, q, purpose='target'):
        return self.evaluate(to_probit(q), purpose)

    def material_tangent(self, x):
        x = np.asarray(x, dtype=float)
        if x.shape != (66,) or not np.isfinite(x).all(): raise ValueError('finite x66')
        return self.prior.derivative(from_probit(x)) * probit_derivative(x)[None, :]

    def white(self, data):
        data = np.asarray(data)
        if data.shape != (32, 8): raise ValueError('full receiver/source layout required')
        return pack(data)*np.sqrt(2)/self.sigma

    def jacobian(self, p):
        if p.J is None:
            if p.T is None:
                if p is not self._live_point:
                    raise RuntimeError('endpoint LU released; evaluate again or use persisted J/T')
                with self.book.physical('receiver_J_32_RHS'):
                    p.T = self.backend().receiver_transfer(p.state)
                # No backend cache may retain another LU owner.
                self.model._transfers.clear()
            C = np.einsum('rp,ps->rsp', p.T, p.state.field).reshape(256, 576)
            native = C @ self.prior.derivative(p.q)
            p.Jnative = np.column_stack([self.white(native[:, k].reshape(32, 8)) for k in range(66)])
            p.J = p.Jnative * probit_derivative(p.x)[None, :]
            self._persist(p)  # J is durable before returning to the optimizer.
        return p.J

    def jacobian_native(self, p):
        self.jacobian(p)
        return p.Jnative

    def solve_anchor(self, rhs, p, *, adjoint=False, role='jet'):
        """Paid live-LU solve, one RHS per complex column (eight per jet)."""
        if p is not self._live_point:
            raise RuntimeError('anchor LU released; paid reevaluation required')
        rhs = np.asarray(rhs, dtype=np.complex128)
        if rhs.ndim not in (1, 2) or rhs.shape[0] != 576 or not np.isfinite(rhs).all():
            raise ValueError('finite anchor RHS shape (576,) or (576,k) required')
        if role == 'full': raise ValueError('full role reserved for forward-credit consumption')
        with self.book.physical(str(role)+'_anchor_solve'):
            return self.backend().solve_anchor(rhs, p.state, adjoint=adjoint, role=str(role))

    def phi_gradient(self, p, y):
        return self.jacobian(p).T @ self.white(p.state.data-np.asarray(y))

    def objective(self, p, y):
        return self.phi(p, y) + .5*float(p.x @ p.x)

    def persist_cache(self, p, path): self._persist(p, path)

    def load_cache(self, path, owner=None):
        raise RuntimeError('R7 does not import historical live-LU caches; load saved endpoint/J arrays only')

    def close(self):
        self._release_lu()
        self.cached_points.clear()


def source_manifest():
    """Loaded source paths and BLAKE2b provenance; no solver is imported here."""
    paths = {'adapter': Path(__file__), 'original_core': R4/'src/p3r4/core.py',
             'original_prior': R4/'src/p1_material/prior.py',
             'original_target': R4/'research/p3l_r4_conditional_transport/configs/TARGET.json',
             'A30_adapter': A30/'src/a30/physics.py',
             'frozen_CPU_solver': A30/'vendor/frozen_a29/src/a24/physics.py',
             'frozen_contracts': A30/'vendor/frozen_a29/src/a24/contracts.py',
             'saved_background': R4/'research/p3l_r4_conditional_transport/inputs/background/weak.npz'}
    for name, module in tuple(sys.modules.items()):
        if name in ('a24.physics', 'a24.contracts', 'a30.physics') and getattr(module, '__file__', None):
            paths['loaded:'+name] = Path(module.__file__)
    return dict(campaign_id=CAMPAIGN, original_gain=1., original_sigma=SIGMA,
                forward_complex_RHS=8, first_receiver_J_complex_RHS=32,
                loaded_solver_modules=[name for name in ('a30.physics', 'a24.physics', 'a24.contracts') if name in sys.modules],
                files={key: dict(path=str(path.resolve()), bytes=path.stat().st_size,
                                blake2b=hashlib.blake2b(path.read_bytes()).hexdigest())
                       for key, path in paths.items()})
