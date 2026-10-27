import unittest
from decimal import Decimal
import numpy as np
from bgga.decimal_bounds import normalized_intervals,independent_reference

class DecimalBounds(unittest.TestCase):
    def test_all_known_encloses_independent_high_precision(self):
        p=np.array([.1,.2,.7]);ll=np.array([-32.1,-33.3,-1000.])
        r=normalized_intervals(p,ll,np.ones(3,bool))
        for l,u,q in zip(r['lower'],r['upper'],independent_reference(p,ll)):
            self.assertLessEqual(Decimal.from_float(l),q);self.assertLessEqual(q,Decimal.from_float(u))
        self.assertEqual(r['unknown_mass_upper'],0.)
    def test_unknown_positive_prior_not_deleted(self):
        p=np.array([.1,.2,.7]);ll=np.array([-32.1,-33.3,-1000.])
        r=normalized_intervals(p,ll,np.array([True,False,True]))
        ref=independent_reference(p,ll)
        self.assertGreater(r['unknown_mass_upper'],.99)
        for l,u,q in zip(r['lower'],r['upper'],ref):
            self.assertLessEqual(Decimal.from_float(l),q);self.assertLessEqual(q,Decimal.from_float(u))
