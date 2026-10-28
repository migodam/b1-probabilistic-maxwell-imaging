"""Review the preserved runner snapshot without executing physics or GT reads."""
import ast
import importlib.util
import json
import sys
import tempfile
from collections import Counter
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import numpy as np

P=Path(__file__).resolve().parents[2]
sys.dont_write_bytecode=True
sys.path.insert(0,str(P/'src'))
source=P/'audit/prelocked/source_snapshot/run_e2.py'
spec=importlib.util.spec_from_file_location('r7_prelocked_snapshot',source)
runner=importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)
tree=ast.parse(source.read_text())
functions={f.name:f for f in tree.body if isinstance(f,ast.FunctionDef)}
readers={}
for name,f in functions.items():
    readers[name]=[ast.unparse(c.args[0]) for c in ast.walk(f)
        if isinstance(c,ast.Call) and isinstance(c.func,ast.Attribute)
        and isinstance(c.func.value,ast.Name) and c.func.value.id=='np'
        and c.func.attr=='load' and c.args]

class FakeTarget:
    def __init__(self):
        self.book=SimpleNamespace(counts=Counter())
        self._live_point=None;self.evaluations=0
    def evaluate(self,x,purpose):
        self.evaluations+=1
        self.book.counts['actual_forward_RHS']+=8 # FAKE bookkeeping only
        p=SimpleNamespace(x=x.copy(),J=None,T=None,key=f'FAKE_{self.evaluations}',
                          state=SimpleNamespace(data=np.zeros((32,8),complex)))
        self._live_point=p
        return p
    def objective(self,p,y):return .5*float(p.x@p.x)
    def phi(self,p,y):return 0.
    def white(self,z):return np.r_[z.real.T,z.imag.T].ravel()
    def jacobian(self,p):p.J=np.zeros((512,66));return p.J
    def persist_cache(self,p,path):pass

def fake_optimizer(fun,x,jac,**kwargs):
    for k in range(64):fun(np.full(66,k/100))
    return SimpleNamespace(status=0,success=False,message='FAKE cap fixture',nfev=64,njev=0,optimality=1.)

with tempfile.TemporaryDirectory() as temp:
    target=FakeTarget()
    with patch.object(runner,'least_squares',side_effect=fake_optimizer):
        runner.direct_anchor(target,np.zeros((32,8),complex),'FAKE_PARENT',Path(temp),cap=64)
    status=json.loads((Path(temp)/'ANCHOR_STATUS.json').read_text())
    anchor=dict(frozen_cap=64,reported_fresh_state_attempts=status['fresh_state_attempts'],
                actual_fake_evaluate_calls=target.evaluations,
                finding='restoration evaluations bypass fresh cap' if target.evaluations>64 else 'cap respected')

# The saved receipt is measurement/derivative QA, never offline truth.
probe=json.loads((P/'audit/prelocked/source_snapshot/PROBE_RESULT.json').read_text())
cost=dict(base_forward=8,three_jet_solves=3*8,source_adjoint=8,
          FD_3directions_3steps_plusminus=3*3*2*8,finite_jet_endpoints=2*8)
assert sum(cost.values())==200==probe['actual_cost']['RHS']
assert probe['actual_cost']['forward_RHS']==168
assert probe['actual_cost']['derivative_RHS']==32
report=dict(status='HOLD_PENDING_PRELOCKED_FIXES',new_Maxwell_RHS=0,GT_access=False,
    source=str(source),np_load_paths_by_function=readers,anchor_cap_fixture=anchor,
    first_probe_200_RHS_decomposition=cost,
    first_probe_finite_jet_errors=[t['relative_error'] for t in probe['finite_jet_trials']],
    loaded_solver_modules=[name for name in ('a30.physics','a24.physics') if name in sys.modules])
assert not report['loaded_solver_modules']
(P/'audit/prelocked/RESULT.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
