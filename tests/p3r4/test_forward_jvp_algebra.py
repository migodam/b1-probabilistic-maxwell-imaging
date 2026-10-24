"""Packing and derivative/unit tests only; no neural optimizer or physics."""
import unittest
import numpy as np
from p3r4.forward_jvp_surrogate import (whiten,unwhiten,fixed_directions,
 input_tangent_matrix,parameter_count)

class ForwardJVPAlgebra(unittest.TestCase):
 def test_reversible_complex_packing_and_phi_units(self):
  rng=np.random.default_rng(2026102518);a=rng.normal(size=(32,8))+1j*rng.normal(size=(32,8));y=rng.normal(size=(32,8))+1j*rng.normal(size=(32,8));sigma=.00119211798514519
  w=whiten(a,sigma);np.testing.assert_allclose(unwhiten(w,sigma),a,atol=1e-15)
  self.assertAlmostEqual(.5*np.linalg.norm(whiten(a-y,sigma))**2,np.vdot(a-y,a-y).real/sigma**2,delta=1e-6)
 def test_fixed_unit_raw_coordinate_directions(self):
  for d in (64,66):
   V=fixed_directions(d);np.testing.assert_allclose((V*V).sum(axis=0),1.)
   np.testing.assert_array_equal(V,fixed_directions(d))
   self.assertLess(parameter_count(d),800000)
 def test_normalization_chain_rule_with_fixed_context(self):
  # Literal linear tensors, not a neural model/optimization.
  import torch
  torch.set_num_threads(1);rng=np.random.default_rng(3);d=2;m=3;J0=rng.normal(size=(m,d));A=input_tangent_matrix(J0);n=A.shape[0]
  W=rng.normal(size=(m,n));input_scale=rng.uniform(.5,2,size=n);output_scale=rng.uniform(.5,2,size=m);h=rng.normal(size=d);V=rng.normal(size=d)
  tA=torch.tensor(A,dtype=torch.float64);tW=torch.tensor(W,dtype=torch.float64);sx=torch.tensor(input_scale,dtype=torch.float64);sy=torch.tensor(output_scale,dtype=torch.float64)
  function=lambda delta:sy*(tW@((tA@delta)/sx))
  _,jvp=torch.autograd.functional.jvp(function,torch.tensor(h,dtype=torch.float64),torch.tensor(V,dtype=torch.float64))
  expected=output_scale*(W@((A@V)/input_scale));np.testing.assert_allclose(jvp.numpy(),expected,atol=1e-12)

if __name__=='__main__':unittest.main()
