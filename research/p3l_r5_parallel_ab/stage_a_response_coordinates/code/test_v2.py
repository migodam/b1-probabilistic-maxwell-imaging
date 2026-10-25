import unittest
import numpy as np
from run_v2 import Conditional,nested_original,response_bounded,trunc_normal_nodes,boxline
class Toy:
    def full(self,z,purpose):
        z=np.asarray(z);f=np.zeros(512);f[:2]=np.array([.3,.1])*z
        return dict(z=z.copy(),f=f,chi=np.full(576,z[0]+1j*z[1]),key='toy')
    def jac(self,z,purpose,V=None):
        J=np.zeros((512,2));J[:2]=np.diag([.3,.1]);return J if V is None else J@V
class Tests(unittest.TestCase):
    def test_linear_original_response(self):
        toy=Toy();sig=.1;J=toy.jac(np.zeros(2),'')/sig;m=dict(z=np.zeros(2),J=J,C=np.linalg.inv(np.eye(2)+J.T@J),r=toy.full(np.zeros(2),''));c=Conditional(toy,m,np.zeros(512),sig)
        o,_,_=nested_original(c,24,9);r,_,dg=response_bounded(c,24,9)
        self.assertLess(abs(np.expm1(o['logZ']-r['logZ'])),.005)
        self.assertLess(np.linalg.norm(o['mean']-r['mean']),1e-8)
        self.assertEqual(dg['inverse_failures'],0)
    def test_stable_positive_gaussian_tail(self):
        s,w,lm=trunc_normal_nodes(8.,9.,7);self.assertTrue(np.all(np.asarray(s)>=8));self.assertTrue(np.all(np.asarray(s)<=9));self.assertAlmostEqual(sum(w),1)
        self.assertTrue(np.isfinite(lm))
    def test_box_partition(self):
        V=np.array([.6,.8]);W=np.array([-.8,.6]);lo,hi=boxline(np.zeros(2),V,W,1.)
        for b in [lo,hi]:self.assertLessEqual(np.abs(V+W*b).max(),4+1e-12)
if __name__=='__main__':unittest.main()
