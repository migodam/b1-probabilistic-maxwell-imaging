"""Twelve predeclared index quantiles, independent field grid, same discrete prior."""
from pathlib import Path
import json,gc
import numpy as np
from wslam.accounting import BudgetBook,json_default
from wslam.forward import Evaluator
from wslam.world import World
root=Path(__file__).resolve().parents[1];source=root/'results/e2e3/finite_bank/online'
out=root/'results/e2e3/bank_validation';out.mkdir(parents=True,exist_ok=False)
grid=np.load(source/'FINITE_GRID.npz',allow_pickle=False)
bank=np.load(source/'FINITE_BANK.npz',allow_pickle=False)
if int(bank['completed'])!=2950:raise RuntimeError('Incomplete prior bank')
ids=np.linspace(0,2949,12).astype(int);book=BudgetBook();rows=[]
for i in ids:
    world=World(grid['support'][i].reshape(2,2),float(grid['epsinf'][i]),float(grid['ell'][i]),grid['pose'][i])
    solver=Evaluator(book,'E2E3',h=.17,n=384,grid_refinement=3,tx=bank['tx'],rx=bank['rx'])
    ev=solver.evaluate(world)
    error=float(np.linalg.norm(ev.data-bank['data'][i])/.01)
    np.savez_compressed(out/f'ROW_{i:04d}.npz',primary=bank['data'][i],independent=ev.data,
                        state_index=i,support=world.support,pose=world.pose)
    rows.append(dict(index=int(i),white_difference=error,residual=ev.state.residual,world=world.metadata()))
    del ev;solver.clear();gc.collect()
    (out/'RESULTS.json').write_text(json.dumps(dict(rows=rows,status='IN_PROGRESS'),indent=2,default=json_default))
result=dict(decision='PASS' if all(r['white_difference']<=.1 for r in rows) else 'FAIL',
    frozen_numerical_tolerance=.1,max_white_difference=max(r['white_difference'] for r in rows),rows=rows,
    status='QUALIFIED_FIXED_GRID_AUDIT_NOT_UNIFORM_ERROR_BOUND',continuous_posterior_certificate=False,
    all_prior_states_independent_recomputed=False,stage_cost_snapshot=book.totals('E2E3'))
(out/'RESULTS.json').write_text(json.dumps(result,indent=2,default=json_default));print(json.dumps(result,default=json_default))
