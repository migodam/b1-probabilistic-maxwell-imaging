"""Substantive geometry, accounting, and factor-invariant checks."""
import json
import multiprocessing
from pathlib import Path
import numpy as np
import unittest
import tempfile
from wslam.accounting import BudgetBook, BudgetExceeded, CAPS
from wslam.geometry import validate_geometry, enumerate_supports, disk_wall_gap
from wslam.factors import FactorRegistry, refresh_logweights

def _reserve_worker(path,count):
    book = BudgetBook(path)
    for _ in range(count):
        book.reserve('E0','parallel',rhs=1,factorizations=1)

class InfrastructureTests(unittest.TestCase):
    def setUp(self):
        self._temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self._temporary.cleanup)
        self.tmp_path = Path(self._temporary.name)


    def test_boundary_empty_ring_and_passivity(self):
        empty = validate_geometry(np.zeros((2,2)))
        assert empty['legal'] and empty['components']==0 and empty['holes']==0
        assert np.all(empty['epsinf']==1) and empty['area']==0
        ring = np.ones((3,3),dtype=bool)
        ring[1,1] = False
        v = validate_geometry(ring,epsinf=6,ell=.25)
        assert v['legal'] and v['holes']==1 and v['components']==1
        self.assertAlmostEqual(v['area'],2)
        self.assertAlmostEqual(validate_geometry(np.ones((2,2)))['area'],2.25)
        assert not validate_geometry(ring,epsinf=6.001)['legal']
        assert not validate_geometry(ring,ell=-.001)['legal']

    def test_disconnection_and_connected_checkerboard_are_distinct(self):
        assert not validate_geometry([[1,0],[0,1]])['legal']
        m = np.zeros((6,24),dtype=bool)
        # Connect the two diagonal cells around the other side; the local pinch
        # remains forbidden even when the global support is connected.
        m[1,1]=m[2,2]=True
        m[1,0]=m[2,0]=m[3,0]=m[3,1]=m[3,2]=True
        v = validate_geometry(m)
        assert v['components']==1
        assert 'checkerboard_diagonal_contact' in v['reasons']

    def test_exact_rectangle_disk_gap_and_target_bounds(self):
        full = np.ones((6,24),dtype=bool)
        self.assertAlmostEqual(disk_wall_gap(full,(0,1.6),.18),.17)
        assert validate_geometry(full,target_center=(0,1.6),target_radius=.18)['legal']
        assert not validate_geometry(full,target_center=(1.1,1.6),target_radius=.1)['legal']
        one = np.zeros((6,24),dtype=bool)
        one[0,0]=True
        # Closest point is the upper-right corner of the occupied rectangle.
        expected = np.hypot(.125,.125)-.08
        self.assertAlmostEqual(disk_wall_gap(one,(-1.25,.75),.08),expected)

    def test_conditioned_bernoulli_has_correct_relative_mass(self):
        rows = enumerate_supports(2)
        self.assertAlmostEqual(sum(p for _,p in rows),1)
        assert all(validate_geometry(s)['legal'] for s,_ in rows)
        p0 = next(p for s,p in rows if s.sum()==0)
        p1 = next(p for s,p in rows if s.sum()==1)
        self.assertAlmostEqual(p1/p0,.65/.35)
        assert len(rows)==14

    def test_atomic_reservation_and_failed_work_stays_charged(self):
        path = self.tmp_path/'ledger.json'
        ctx = multiprocessing.get_context('spawn')
        workers = [ctx.Process(target=_reserve_worker,args=(str(path),20)) for _ in range(3)]
        for p in workers:
            p.start()
        for p in workers:
            p.join(30)
            assert p.exitcode==0
        book = BudgetBook(path)
        assert book.snapshot('E1')['E0E1']['used']['rhs']==60
        book.event('E0','failed',kind='failure',exit_code=1)
        book.reserve('E1','fill',rhs=3940,factorizations=1140)
        with self.assertRaises(BudgetExceeded):
            book.reserve('E0','denied',rhs=1)
        assert book.snapshot('E0')['E0E1']['used']['rhs']==4000
        book.event('E0','cpu',process_cpu_seconds=7200)
        with self.assertRaises(BudgetExceeded):
            book.reserve('E0','even-zero-denied')
        assert len(json.loads(path.read_text())['events'])==63

    def test_factors_own_observations_and_replace_by_version(self):
        registry = FactorRegistry()
        registry.register('prior','prior')
        registry.register('obs1','likelihood',observation_id='observation-1')
        with self.assertRaises(ValueError):
            registry.register('duplicate','likelihood',observation_id='observation-1')
        with self.assertRaises(ValueError):
            registry.register('obs1','likelihood',observation_id='observation-1')
        registry.register('obs1','likelihood',version=2,observation_id='observation-1')
        assert len(registry)==2
        base = np.array([.2,.3,.5])
        old,new = np.array([.4,.2,.8]),np.array([.1,.9,.3])
        old_posterior = base*old
        old_posterior /= old_posterior.sum()
        refreshed = np.exp(refresh_logweights(np.log(old_posterior),np.log(old),np.log(new)))
        expected = base*new
        expected /= expected.sum()
        assert np.allclose(refreshed,expected,rtol=1e-13,atol=1e-15)

    def test_backend_diagnostics_never_doublecharge_and_scientific_metadata_serializes(self):
        path = self.tmp_path/'ledger.json'
        book = BudgetBook(path)
        book.reserve('E0E1','solve',rhs=3,factorizations=1)
        book.event('E0E1','backend_solve',rhs=np.int64(3),factorizations=np.int64(1),
                   scalar=np.float64(.5),array=np.array([1.,2.]),
                   complex_scalar=np.complex128(1+2j),complex_array=np.array([1+2j,3+4j]))
        assert book.totals('E0E1') == dict(rhs=3,factorizations=1,cpu_seconds=0.)
        event = json.loads(path.read_text())['events'][-1]
        assert event['observed_rhs']==3 and event['observed_factorizations']==1
        assert 'rhs' not in event and 'factorizations' not in event
        assert event['array']==[1.,2.] and event['scalar']==.5
        assert event['complex_scalar']==[1.,2.]
        assert event['complex_array']==[[1.,2.],[3.,4.]]


if __name__ == '__main__':
    unittest.main()
