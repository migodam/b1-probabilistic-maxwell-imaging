"""Zero-physics prospective wrapper contract and pre-Book release guards."""
import unittest,tempfile,sys,json
from pathlib import Path
from unittest.mock import patch
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'tools'))
import run_r4_interval_chain_locked as driver

class LockedFourArmPreparation(unittest.TestCase):
 def test_frozen_distinct_stage_roster_and_caps(self):
  p=driver.protocol()
  self.assertEqual(p['ledger_stage'],'IV');self.assertEqual(p['seed'],20261025171)
  self.assertEqual(p['parents'],tuple(f'R4L1_{i:05d}' for i in range(8)))
  self.assertEqual(sum(p['arm_complex_RHS_caps'].values()),300000)
  self.assertEqual(driver.RESULT,driver.OUT/'results/IV/interval_chain_locked_v1')
 def test_missing_locked_gate_prevents_launch_without_book(self):
  with tempfile.TemporaryDirectory() as folder:
   with patch.object(driver,'PREPARED',Path(folder)/'missing_prepared'),patch.object(driver,'LOCKED',Path(folder)/'missing_locked'):
    with patch.object(driver.engine,'run',side_effect=AssertionError('must remain inert')) as launch:
     status=driver.preflight()
     self.assertEqual(status['status'],'NOT_READY_MISSING_PREREQUISITES');self.assertEqual(status['new_RHS'],0)
     launch.assert_not_called()

if __name__=='__main__':unittest.main()
