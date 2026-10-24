"""Prospective budget/roster/release guard; zero evaluator calls."""
import unittest,sys,tempfile,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'tools'))
import run_r4_interval_chain_boundary as driver

class BoundaryPreparation(unittest.TestCase):
 def test_fixed_grid_and_immutable_cache_ownership(self):
  p=driver.protocol()
  self.assertEqual(len(p['cases']),36);self.assertEqual(sum(driver.CAPS.values()),300000)
  self.assertEqual(len(set(c['parent'] for c in p['cases'])),4)
  self.assertEqual(p['transitions_per_arm'],9216);self.assertEqual(p['original_R0_setup_RHS_perapproxarm'],2304)
  self.assertIn('one immutable original',p['original_R0_lifetime'])
  self.assertEqual([(c['beta'],c['seed'],c['parent_index']) for c in driver.CASES],sorted((c['beta'],c['seed'],c['parent_index']) for c in driver.CASES))
 def test_false_release_stops_before_protocol_or_physics(self):
  with tempfile.TemporaryDirectory() as folder:
   path=Path(folder)/'release.json';path.write_text(json.dumps(dict(released=False,experiment=driver.EXPERIMENT,campaign_id='P3LR4-20261025')))
   with self.assertRaises(PermissionError):driver.check_release(path)

if __name__=='__main__':unittest.main()
