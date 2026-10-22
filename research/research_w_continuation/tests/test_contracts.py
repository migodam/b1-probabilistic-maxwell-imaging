import tempfile,unittest,json
from pathlib import Path
import numpy as np
from wq.accounting import Book,SearchBook,BudgetStop
from wq.inverse_support import PixelGrid,project_parameters,pixel_material,occupied_area,realify_data,realify_jacobian

class Contracts(unittest.TestCase):
    def test_atomic_nonrefund(self):
        with tempfile.TemporaryDirectory() as d:
            b=Book(d);b.reserve('validation','all',rhs=3000)
            with self.assertRaises(BudgetStop):b.reserve('validation','overflow',rhs=1)
            self.assertEqual(b.totals()['rhs'],3000)
    def test_stage_isolation_and_phase(self):
        with tempfile.TemporaryDirectory() as d:
            q=Book(Path(d)/'q');b=SearchBook(Path(d)/'b');q.reserve('validation','q',rhs=100)
            b.reserve('search','cap',rhs=3200)
            with self.assertRaises(BudgetStop):b.reserve('search','extra',rhs=1)
            self.assertEqual(q.totals()['rhs'],100)
    def test_whitening(self):
        a=np.array([1+2j,3-4j]);self.assertAlmostEqual(np.linalg.norm(realify_data(a,noise_std=.01)),np.linalg.norm(a)/.01)
    def test_full_domain_legal_air(self):
        g=PixelGrid();rng=np.random.default_rng(1);p=project_parameters(rng.normal(size=g.n_parameters),g)
        self.assertLessEqual(occupied_area(p,g),3.)
        e=pixel_material(p,1.,g)(np.array([[0,2.],[0,-.3],[1.9,.2]]));self.assertEqual(e[0],1);self.assertEqual(e[1],1)
    def test_shared_conductivity(self):
        g=PixelGrid();p=project_parameters(np.r_[np.ones(g.n_pixels),np.full(g.n_pixels,.02)],g)
        pts=np.array([[-1.8,.1],[-1.7,.1]])
        a=pixel_material(p,.8,g)(pts);b=pixel_material(p,1.2,g)(pts)
        np.testing.assert_allclose(.8*a.imag,1.2*b.imag)
    def test_real_packing_matches_jacobian(self):
        rng=np.random.default_rng(7);j=rng.normal(size=(48,10))+1j*rng.normal(size=(48,10));v=rng.normal(size=10)
        np.testing.assert_allclose(realify_jacobian(j,noise_std=.01)@v,realify_data(j@v,noise_std=.01))

if __name__=='__main__':unittest.main()
