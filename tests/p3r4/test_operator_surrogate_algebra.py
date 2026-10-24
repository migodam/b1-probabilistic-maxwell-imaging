import unittest
import numpy as np
from scipy.special import ndtr
from p3r4.operator_surrogate import (compose_from_transfer,one_secant,gauss_newton_action,
 coordinate_derivative,ordered_evolution,ranker_parameter_count,FEATURE_NAMES)
from p3r4.core import pack

class OperatorAlgebra(unittest.TestCase):
 def test_exact_cached_composition_packing(self):
  rng=np.random.default_rng(2026102514);T=rng.normal(size=(4,9))+1j*rng.normal(size=(4,9));E=rng.normal(size=(9,3))+1j*rng.normal(size=(9,3));D=rng.normal(size=(9,5))+1j*rng.normal(size=(9,5));sigma=.2
  J=compose_from_transfer(T,E,D,sigma)
  expected=np.column_stack([pack(T@(E*D[:,k,None]))*np.sqrt(2)/sigma for k in range(5)])
  np.testing.assert_allclose(J,expected,rtol=1e-13,atol=1e-12)

 def test_gaussian_prior_gradient_is_in_action(self):
  J=np.diag([3.,4.]);r=np.array([.2,-.4]);x=np.array([.7,-.2]);D=np.eye(2,dtype=complex)
  delta,material=gauss_newton_action(J,r,x,D)
  np.testing.assert_allclose(delta,-(np.diag(J)*r+x)/(1+np.diag(J)**2))
  np.testing.assert_allclose(material,np.r_[delta,np.zeros(2)]/np.sqrt(2))

 def test_probit_property_chain_and_secant(self):
  rng=np.random.default_rng(7);A=rng.normal(size=(9,66));prior=type('LinearPrior',(),{'derivative':lambda self,q:A})()
  x=rng.normal(size=66);q=x.copy();q[64:]=ndtr(x[64:]);D=coordinate_derivative(prior,q,x,66)
  for k in (2,64,65):
   h=1e-5;a=x.copy();b=x.copy();a[k]+=h;b[k]-=h;a[64:]=ndtr(a[64:]);b[64:]=ndtr(b[64:])
   np.testing.assert_allclose(D[:,k],A@(a-b)/(2*h),rtol=1e-8,atol=1e-9)
  J=rng.normal(size=(17,5));h=rng.normal(size=5);df=rng.normal(size=17);updated,u,v=one_secant(J,h,df)
  np.testing.assert_allclose(updated@h,df,atol=1e-14)

 def test_telescope_and_cross_terms_and_capacity(self):
  rng=np.random.default_rng(11);J0=rng.normal(size=(17,5));d=rng.normal(size=(17,5));f=-.7*d;p=rng.normal(size=(17,5))*.1
  result=ordered_evolution(J0,J0+d,J0+d+f,J0+d+f+p)
  self.assertLess(result['telescoping_error'],1e-13);self.assertLess(result['squared_norm_identity_error'],1e-12)
  self.assertLess(result['cross_inner_products'][-1],0)
  self.assertLessEqual(ranker_parameter_count(),800000);self.assertEqual(len(FEATURE_NAMES),1176)

if __name__=='__main__':unittest.main()
