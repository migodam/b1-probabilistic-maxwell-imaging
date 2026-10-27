import json,tempfile,unittest
from pathlib import Path
from bgga.accounting import BudgetBook,ResourceAuthorizationRequired,BudgetExceeded

class Accounting(unittest.TestCase):
    def test_saved_diagnostic_matrix_serialization(self):
        import numpy as np
        from bgga.accounting import _default
        text=json.dumps(dict(matrix=np.eye(2),scalar=np.int32(2)),default=_default)
        self.assertEqual(json.loads(text),dict(matrix=[[1.,0.],[0.,1.]],scalar=2))
    def test_no_implicit_resource_grant(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);(p/'AUTHORIZATION.json').write_text(json.dumps(dict(new_physics_authorized=False,authorized_caps={})))
            b=BudgetBook(p)
            with self.assertRaises(ResourceAuthorizationRequired):b.reserve('G0','test',rhs=1)
            self.assertFalse(b.path.exists())
    def test_synthetic_grant_is_not_replenished(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);(p/'AUTHORIZATION.json').write_text(json.dumps(dict(new_physics_authorized=True,
                authorized_caps=dict(G0=dict(rhs=2,full_fe=1,cpu_seconds=10)))))
            b=BudgetBook(p);b.reserve('G0','synthetic',rhs=2,factorizations=1)
            with self.assertRaises(BudgetExceeded):BudgetBook(p).reserve('G0','again',rhs=1)
            b.event('G0','kernel',rhs=999,factorizations=999)
            self.assertEqual(b.snapshot('G0')['rhs'],2)
    def test_aggregate_cpu_only_once(self):
        with tempfile.TemporaryDirectory() as d:
            b=BudgetBook(Path(d));b.event('PREFLIGHT','kernel',kernel_cpu_seconds=.2)
            b.event('PREFLIGHT','job_completion',process_cpu_seconds=.4)
            self.assertEqual(b.snapshot('PREFLIGHT')['process_cpu_seconds'],.4)

if __name__=='__main__':unittest.main()
