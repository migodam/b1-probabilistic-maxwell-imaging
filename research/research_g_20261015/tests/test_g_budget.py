"""Accounting fixtures only: ledger reservations never run a physical solver."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
from g_budget import GBudget, BudgetStop
ROOT=Path(__file__).resolve().parents[1]
class Field:
    shape=(24,24)
class BudgetContract(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(); self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.config=json.loads((ROOT/'E1_PROTOCOL_FROZEN.json').read_text())
        self.write_config(); self.book=GBudget(self.root)
        self.nice=patch('g_budget.os.nice'); self.nice.start(); self.addCleanup(self.nice.stop)
    def write_config(self):
        (self.root/'E1_PROTOCOL_FROZEN.json').write_text(json.dumps(self.config))
    def state(self):self.book.reserve_state(Field(),resolution=24,purpose='fixture-only')
    def test_e0_rejects_all_new_physics(self):
        with self.book.job('e0',stage='E0'):
            with self.assertRaises(BudgetStop):self.state()
            with self.assertRaises(BudgetStop):self.book.reserve_rhs(1,'fixture')
        self.assertEqual(self.book.totals(),{'states':0,'rhs':0})
    def test_e1_requires_complete_frozen_metadata(self):
        original=copy.deepcopy(self.config)
        invalid=[{}, {'execution_authorized':True,'backend_ready':True}]
        for key in original['backend']:
            if key in {'target_id','A30_commit','corrector_commit','corrector_id','packing','gain','whitening','data_threshold'}:
                item=copy.deepcopy(original); del item['backend'][key]; invalid.append(item)
        for key in ['execution_authorized','backend_ready','cases']:
            item=copy.deepcopy(original);del item[key];invalid.append(item)
        for index,config in enumerate(invalid):
            with self.subTest(index=index):
                self.config=config;self.write_config()
                with self.assertRaises((BudgetStop,ValueError)):
                    with self.book.job('bad'+str(index)):pass
        self.config=original; self.write_config()
        with self.book.job('complete'):pass
    def test_global_state_cap_independent_from_rhs(self):
        self.book.append({'event':'STATE_RESERVED','amount':383})
        with self.book.job('state-cap'):
            self.state()
            with self.assertRaisesRegex(BudgetStop,'GLOBAL_STATES_CAP'):self.state()
            self.book.reserve_rhs(1,'independent')
        self.assertEqual(self.book.totals(),{'states':384,'rhs':1})
    def test_global_rhs_cap_independent_from_states(self):
        self.book.append({'event':'RHS_RESERVED','amount':12287})
        with self.book.job('rhs-cap'):
            self.book.reserve_rhs(1,'last')
            with self.assertRaisesRegex(BudgetStop,'GLOBAL_RHS_CAP'):self.book.reserve_rhs(1,'overflow')
            self.state()
        self.assertEqual(self.book.totals(),{'states':1,'rhs':12288})
    def test_arm_hard_caps_cannot_be_raised(self):
        arm='first/case/method'
        self.book.append({'event':'STATE_RESERVED','amount':11,'arm':arm})
        self.book.append({'event':'RHS_RESERVED','amount':287,'arm':arm})
        with self.book.job('arm'):
            with self.book.arm('first','case','method',states=999,rhs=999):
                self.state();self.book.reserve_rhs(1,'last')
                with self.assertRaisesRegex(BudgetStop,'ARM_STATES_CAP'):self.state()
                with self.assertRaisesRegex(BudgetStop,'ARM_RHS_CAP'):self.book.reserve_rhs(1,'overflow')
    def test_failed_reserved_attempt_is_not_refunded(self):
        with self.assertRaisesRegex(RuntimeError,'synthetic failure'):
            with self.book.job('failed'):
                self.state();self.book.reserve_rhs(8,'fixture');raise RuntimeError('synthetic failure')
        self.assertEqual(self.book.totals(),{'states':1,'rhs':8})
        self.assertEqual(self.book.receipt()['failed_jobs'],1)
    def test_wall_cap_10800(self):
        self.book.append({'event':'ZERO_PHYSICS_ANALYSIS','wall_seconds':10800,'cpu_seconds':0})
        with self.assertRaisesRegex(BudgetStop,'WALL_EXHAUSTED'):
            with self.book.job('wall'):self.fail('must not execute')
        self.assertEqual(self.book.totals(),{'states':0,'rhs':0})
    def test_reopening_does_not_reset_ledger(self):
        with self.book.job('first'):
            self.state();self.book.reserve_rhs(8,'fixture')
        reopened=GBudget(self.root)
        self.assertEqual(reopened.totals(),{'states':1,'rhs':8})
        with reopened.job('second'):reopened.reserve_rhs(1,'fixture')
        self.assertEqual(reopened.totals(),{'states':1,'rhs':9})
    def test_bool_rhs_is_rejected(self):
        with self.book.job('bool'):
            for value in [True,False]:
                with self.assertRaises(ValueError):self.book.reserve_rhs(value,'fixture')
        self.assertEqual(self.book.totals()['rhs'],0)
if __name__=='__main__':unittest.main()
