"""Substantive small-field consistency and Euler self-defect tests, no physics."""
import unittest
import numpy as np
from p3r4.dynamic_geometry import heun_corrector


def circuit(x,X,Y,h,method='Heun'):
    z=np.asarray(x,float).copy()
    for F,sign in ((X,1),(Y,1),(X,-1),(Y,-1)):
        v=F(z);predictor=z+sign*h*v
        z=predictor if method=='Euler' else heun_corrector(z,v,F(predictor),sign*h)
    return z


class HeunAlgebra(unittest.TestCase):
    def test_commuting_fields_remove_order_h2_self_defect(self):
        x=np.array([.7,-.3]);X=lambda z:np.array([z[0],0.]);Y=lambda z:np.array([0.,z[1]])
        normalized=[]
        for h in (.1,.05,.025,.0125):
            np.testing.assert_allclose((circuit(x,X,Y,h,'Euler')-x)/h**2,-x,atol=2e-12)
            normalized.append(np.linalg.norm((circuit(x,X,Y,h)-x)/h**2))
        self.assertLess(normalized[-1],normalized[0]/30)

    def test_noncommuting_linear_fields_recover_signed_bracket(self):
        A=np.diag([1.,0.]);B=np.array([[0.,1.],[0.,0.]])
        X=lambda z:A@z;Y=lambda z:B@z;x=np.array([.7,-.3]);bracket=(B@A-A@B)@x
        errors=[]
        for h in (.1,.05,.025,.0125):
            errors.append(np.linalg.norm((circuit(x,X,Y,h)-x)/h**2-bracket))
        self.assertTrue(all(a>b for a,b in zip(errors[:-1],errors[1:])))
        self.assertLess(errors[-1],errors[0]/5)
        euler=(circuit(x,X,Y,.0125,'Euler')-x)/.0125**2
        self.assertGreater(np.linalg.norm(euler-bracket),.5)
        reverse=(circuit(x,Y,X,.0125)-x)/.0125**2
        self.assertLess(np.linalg.norm(reverse+bracket),.01)

if __name__=='__main__':unittest.main()
