import unittest
import numpy as np
from p3r4.best_query_candidate import ordinary_best_query,material_best_query

class BestQueryCandidate(unittest.TestCase):
 def test_rejected_but_forward_valid_query_is_allowed(self):
  q=[dict(UID='anchor',Phi=10.,chi=np.array([0j]),prior_legal=True,forward_valid=True,callback=True),
     dict(UID='paid',Phi=5.,chi=np.array([1+0j]),prior_legal=True,forward_valid=True,callback=False),
     dict(UID='invalid',Phi=1.,chi=np.array([1+0j]),prior_legal=True,forward_valid=False,callback=True)]
  self.assertEqual(ordinary_best_query(q)['UID'],'paid')
 def test_material_selector_keeps_original_noiseball_and_tie_rule(self):
  q=[dict(UID='anchor',Phi=10.,chi=np.array([0j]),prior_legal=True,forward_valid=True),
     dict(UID='b',Phi=12.,chi=np.array([1+0j]),prior_legal=True,forward_valid=True),
     dict(UID='a',Phi=12.,chi=np.array([1+0j]),prior_legal=True,forward_valid=True),
     dict(UID='outside',Phi=300.,chi=np.array([1+0j]),prior_legal=True,forward_valid=True)]
  self.assertEqual(material_best_query(q,np.array([1+0j]))['UID'],'a')
if __name__=='__main__':unittest.main()
