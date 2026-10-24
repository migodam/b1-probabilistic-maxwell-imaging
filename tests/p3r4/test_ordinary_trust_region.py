"""Pure TR KKT and whitened/probit gradient algebra, no Maxwell."""
import unittest
import numpy as np
from scipy.special import ndtr
from p3r4.core import pack,probit_derivative
from p3r4.ordinary_trust_region import exact_psd_trust_region

class OrdinaryTrustRegion(unittest.TestCase):
 def test_interior_and_singular_boundary_KKT(self):
  H=np.diag([2.,8.]);g=np.array([-.1,.2]);s,m=exact_psd_trust_region(H,g,1.)
  np.testing.assert_allclose(s,-np.linalg.solve(H,g),atol=1e-13);self.assertFalse(m['boundary'])
  H=np.diag([0.,4.]);g=np.array([1.,0.]);s,m=exact_psd_trust_region(H,g,.5)
  np.testing.assert_allclose(s,[-.5,0.],atol=1e-12);self.assertLess(m['KKT_residual'],1e-12)
 def test_general_boundary_psd_sufficient_KKT(self):
  rng=np.random.default_rng(20261025201);A=rng.normal(size=(9,5));H=A.T@A;g=rng.normal(size=5)
  s,m=exact_psd_trust_region(H,g,.08)
  self.assertLessEqual(np.linalg.norm(s),.08*(1+1e-12));self.assertGreater(m['multiplier'],0)
  self.assertLess(np.linalg.norm(g+H@s+m['multiplier']*s),1e-11)
  self.assertGreaterEqual(np.linalg.eigvalsh(H+m['multiplier']*np.eye(5)).min(),-1e-12)
 def test_whitened_Phi_and_full66_probit_chain(self):
  rng=np.random.default_rng(8);A=rng.normal(size=(32,8,2))+1j*rng.normal(size=(32,8,2));y=rng.normal(size=(32,8))+1j*rng.normal(size=(32,8));x=rng.normal(size=66);sigma=2.
  F=lambda z:A[:,:,0]*z[0]+A[:,:,1]*ndtr(z[64])
  white=lambda z:pack(z)*np.sqrt(2)/sigma
  J=np.zeros((512,66));J[:,0]=white(A[:,:,0]);J[:,64]=white(A[:,:,1])*probit_derivative(x)[64]
  for lam in (0,1):
   objective=lambda z:.5*np.dot(white(F(z)-y),white(F(z)-y))+.5*lam*np.dot(z,z)
   grad=J.T@white(F(x)-y)+lam*x
   for k in (0,5,64):
    h=np.zeros(66);h[k]=1e-6;fd=(objective(x+h)-objective(x-h))/(2e-6)
    self.assertAlmostEqual(fd,grad[k],delta=1e-7)
if __name__=='__main__':unittest.main()
