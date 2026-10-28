"""Nonphysical bookkeeping/chain-rule fixtures. No A30/backend import or solve."""
import importlib.util
import json
import sys
import tempfile
import unittest
from unittest.mock import patch
from pathlib import Path
from types import SimpleNamespace
import numpy as np

REPO = Path(__file__).resolve().parents[4]
SOURCE = REPO/'research/p3l_r7_material_solution_fibers/src/physics.py'
spec = importlib.util.spec_from_file_location('r7_physics_fake_test', SOURCE)
physics = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = physics
spec.loader.exec_module(physics)


def authorize(root, cap=50000):
    (root/'configs').mkdir()
    physics.write(root/'configs/AUTHORIZATION.json', dict(campaign_id=physics.CAMPAIGN,
        explicit_user_authorization=True, total_complex_RHS=cap,
        stage_caps={'E0': 0, 'E1': 0, 'E2': cap}))
    physics.write(root/'configs/E2_FREEZE.json', {'fixture': 'NONPHYSICAL_ONLY'})


class FakeBackend:
    n, q, s = 576, 32, 8
    def __init__(self, book):
        self.book, self._transfers = book, {}
        self.T = (np.arange(32)[:, None] + 1)*(np.arange(576)[None, :] + 1)*1e-8*(1+.3j)
    def full_state(self, chi, *, purpose):
        self.book.reserve_state(chi, resolution=24, purpose=purpose)
        self.book.reserve_rhs(8, 'full_forward')
        field = np.tile(np.linspace(1, 2, 8), (576, 1)).astype(complex)
        current = chi[:, None]*field
        return SimpleNamespace(chi=chi, field=field, current=current, data=self.T@current,
                               residual=0., factor=('FAKE_ONLY',), matrix=None)
    def receiver_transfer(self, state):
        self.book.reserve_rhs(32, 'receiver_transfer_adjoint')
        return self.T.copy()
    def solve_anchor(self, rhs, state, *, adjoint=False, role):
        self.book.reserve_rhs(1 if rhs.ndim == 1 else rhs.shape[1],
                              role+('_adjoint' if adjoint else '_forward'))
        return rhs.copy()


class AdapterTests(unittest.TestCase):
    def book(self, root, job='test', cap=50000):
        authorize(root, cap)
        return physics.Book('E2', job, root=root, test_mode=True)

    def test_shapes_chain_rule_persistence_and_one_live_lu(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            with self.book(root) as book:
                target = physics.Target(book, backend_factory=FakeBackend)
                from threadpoolctl import threadpool_info
                self.assertTrue(all(item['num_threads'] == 1 for item in threadpool_info() if item['user_api'] == 'blas'))
                x = np.zeros(66); x[64:] = [.4, -.7]
                p = target.evaluate(x)
                J = target.jacobian(p)
                self.assertEqual(J.shape, (512, 66))
                self.assertEqual(physics.totals(root)['RHS'], 40)
                self.assertTrue(np.array_equal(J, target.jacobian_native(p)*physics.probit_derivative(x)))
                with np.load(book.jobdir/'endpoints'/f'{p.key}.npz') as z:
                    self.assertTrue({'x','q','material','prediction','current','field','residual','J','Jnative'} <= set(z.files))
                    self.assertEqual(z['material'].dtype, np.complex128)
                    self.assertTrue(np.array_equal(z['J'], J))
                direction = np.random.default_rng(21).normal(size=66)
                eps = 1e-6
                def fake_response(v):
                    chi = target.prior.decode(physics.from_probit(v))[0]
                    return target.model.T @ (chi[:, None]*p.state.field)
                finite = target.white((fake_response(x+eps*direction)-fake_response(x-eps*direction))/(2*eps))
                self.assertLess(np.linalg.norm(finite-J@direction)/np.linalg.norm(finite), 2e-7)
                dchi = (target.prior.decode(physics.from_probit(x+eps*direction))[0]
                        -target.prior.decode(physics.from_probit(x-eps*direction))[0])/(2*eps)
                self.assertLess(np.linalg.norm(dchi-target.material_tangent(x)@direction)/np.linalg.norm(dchi), 2e-7)
                target.solve_anchor(np.zeros((576, 8), complex), p, role='jet')
                self.assertEqual(physics.totals(root)['RHS'], 48)
                second = target.evaluate(x+.01)
                self.assertIsInstance(p.state, physics.EndpointState)
                self.assertFalse(hasattr(p.state, 'factor'))
                with self.assertRaises(RuntimeError): target.solve_anchor(np.zeros((576, 8)), p)
                third = target.evaluate(x+.02)
                self.assertEqual(len(target.cached_points), 2)
                self.assertFalse(hasattr(second.state, 'factor'))
                self.assertTrue(hasattr(third.state, 'factor'))
                target.close()
            self.assertFalse(physics.totals(root)['active_jobs'])
            self.assertEqual(physics.totals(root)['RHS'], 64)
            closed = json.loads((book.jobdir/'CLOSED.json').read_text())
            self.assertGreater(closed['peak_RSS_bytes'], 0)

    def test_cap_reserves_before_action_and_no_refunds(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            with self.book(root, cap=8) as book:
                with self.assertRaisesRegex(RuntimeError, 'intentional failure'):
                    with book.physical('failing_fake_state'):
                        book.reserve_state(None, resolution=24, purpose='fake')
                        raise RuntimeError('intentional failure before solver')
                self.assertEqual(physics.totals(root)['RHS'], 8)
                self.assertEqual(book.credit, 8)
                with book.physical('blocked_fake_derivative'):
                    with self.assertRaises(physics.BudgetStop): book.reserve_rhs(1, 'fake_derivative')
                self.assertEqual(physics.totals(root)['failures'], 1)

    def test_lock_and_fail_closed_snapshot(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            book = self.book(root)
            with self.assertRaises(BlockingIOError):
                physics.Book('E2', 'parallel', root=root, test_mode=True)
            book.close()
            snapshot = root/'RESOURCE_LEDGER_SNAPSHOT.json'
            data = json.loads(snapshot.read_text()); data['RHS'] = 0; data['sequence'] += 1
            physics.write(snapshot, data)
            with self.assertRaisesRegex(RuntimeError, 'mismatch'): physics.totals(root)

    def test_zero_stages_and_missing_freeze(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            authorize(root)
            for stage in ('E0', 'E1'):
                with self.assertRaises(physics.BudgetStop):
                    physics.Book(stage, 'zero', root=root, test_mode=True)
            (root/'configs/E2_FREEZE.json').unlink()
            with self.assertRaises(FileNotFoundError):
                physics.Book('E2', 'missing', root=root, test_mode=True)
            self.assertFalse((root/'RESOURCE_LEDGER.jsonl').exists())

    def test_interrupted_job_and_truncated_journal_block_replay(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            book = self.book(root)
            book.workerlock.close()  # Simulated abrupt process exit; no JOB_END.
            with self.assertRaisesRegex(RuntimeError, 'interrupted'):
                physics.Book('E2', 'retry', root=root, test_mode=True)
            with (root/'RESOURCE_LEDGER.jsonl').open('ab') as stream: stream.write(b'{')
            with self.assertRaisesRegex(RuntimeError, 'incomplete'): physics.totals(root)

    def test_explicit_reconciliation_retains_charge(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            book = self.book(root)
            # Simulate a dead worker with a reserved, interrupted state.
            book._event('PHYSICS_BEGIN', action='fake-dead-action', label='fake')
            book.active = 'fake-dead-action'
            book.reserve_state(None, resolution=24, purpose='interrupted fake')
            book.workerlock.close()
            (root/'RESOURCE_LEDGER_SNAPSHOT.json').unlink()
            with self.assertRaisesRegex(RuntimeError, 'missing'): physics.totals(root)
            with patch.object(physics, 'OUT', root):
                physics.reconcile_snapshot(root)
                with self.assertRaisesRegex(RuntimeError, 'PID still exists'):
                    physics.reconcile_interrupted('test', 'fake audit', root)
                with patch.object(physics.os, 'kill', side_effect=ProcessLookupError):
                    physics.reconcile_interrupted('test', 'fake audit, dead worker confirmed', root)
            recovered = physics.totals(root)
            self.assertEqual((recovered['RHS'], recovered['states'], recovered['LU']), (8, 1, 1))
            self.assertFalse(recovered['active_actions'])
            self.assertFalse(recovered['active_jobs'])
            self.assertEqual(recovered['controlled_wall_seconds'], 7200)

    def test_live_constant_size_checks_reject_external_changes(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            book = self.book(root)
            # A live event must not invoke full public journal replay.
            with patch.object(physics, 'totals', side_effect=AssertionError('unexpected full replay')):
                with book.physical('fake_constant_size'):
                    book.reserve_state(None, resolution=24, purpose='fake')
                    book.reserve_rhs(8, 'full_forward')
            snapshot = root/'RESOURCE_LEDGER_SNAPSHOT.json'
            original = snapshot.read_bytes()
            invalid = json.loads(original); invalid['RHS'] += 1
            physics.write(snapshot, invalid)
            with self.assertRaisesRegex(RuntimeError, 'snapshot mismatch'):
                book._event('FAKE_AUDIT')
            snapshot.write_bytes(original)
            with (root/'RESOURCE_LEDGER.jsonl').open('ab') as stream: stream.write(b'\n')
            with self.assertRaisesRegex(RuntimeError, 'changed outside'):
                book._event('FAKE_AUDIT')
            book.workerlock.close()


if __name__ == '__main__':
    before = 'a30.physics' in sys.modules
    result = unittest.main(exit=False).result
    if before or 'a30.physics' in sys.modules:
        raise AssertionError('fake audit must never import the Maxwell backend')
    sys.exit(not result.wasSuccessful())
