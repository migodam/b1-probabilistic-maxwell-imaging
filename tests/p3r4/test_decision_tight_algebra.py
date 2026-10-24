"""Diagonal closed-form fixtures: no inverse, LU, SVD or linear solve calls."""
import unittest
import numpy as np
from p3r4.decision_refresh_tight import inverse_norm_upper_bounds,neumann_constants,phi_interval_from_saved_residual

class TightAlgebra(unittest.TestCase):
 def test_exact_diagonal_norms_and_new_valid_window(self):
  n=9;B=np.eye(n,dtype=complex);b=inverse_norm_upper_bounds(B)
  self.assertEqual(b['frobenius'],3.)
  self.assertEqual(b['selected_upper'],1.)
  self.assertEqual(b['profile']['new_complex_RHS'],0)
  old=neumann_constants(0.,b['frobenius'],.4);new=neumann_constants(0.,b['selected_upper'],.4)
  self.assertIsNone(old['inverse_upper']);self.assertTrue(new['sufficient_condition_valid'])
  # L1=(1-.4/3)I, so its inverse norm is known by scalar arithmetic.
  self.assertGreaterEqual(new['inverse_upper'],1/(1-.4/3))

 def test_enclosure_without_any_solve(self):
  rng=np.random.default_rng(2026102511);n=9;sources=2;receivers=3;sigma=.2
  S=rng.normal(size=(receivers,n))+1j*rng.normal(size=(receivers,n));incident=rng.normal(size=(n,sources))+1j*rng.normal(size=(n,sources))
  for delta in (.01,.03,.08,.13333333333333333):
   chi=np.full(n,delta);exact_current=chi[:,None]*incident/(1-delta)
   approximate=exact_current+.001*(rng.normal(size=(n,sources))+1j*rng.normal(size=(n,sources)))
   observation=rng.normal(size=(receivers,sources))*.1
   rho=chi[:,None]*incident-(1-delta)*approximate
   phi=float(np.linalg.norm(S@exact_current-observation)**2/sigma**2)
   constants=neumann_constants(0.,1.,np.linalg.norm(chi))
   bound=phi_interval_from_saved_residual(np.sqrt(2)/sigma*np.linalg.norm(S@approximate-observation),np.linalg.norm(rho),sigma,np.linalg.norm(S),constants)
   self.assertTrue(bound['available']);self.assertLessEqual(bound['phi_lower'],phi);self.assertGreaterEqual(bound['phi_upper'],phi)

if __name__=='__main__':unittest.main()
