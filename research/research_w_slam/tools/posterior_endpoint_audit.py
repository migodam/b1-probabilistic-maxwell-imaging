"""Pre-data all-case final MAP refinement; no truth and no posterior retuning."""
from pathlib import Path
import gc,json
import numpy as np
from wslam.accounting import BudgetBook,json_default
from wslam.forward import Evaluator
from wslam.world import World

ROOT=Path(__file__).resolve().parents[1]

def main():
    protocol=json.loads((ROOT/'POSTERIOR_ENDPOINT_GRID_AUDIT.json').read_text())
    source=ROOT/'results/e2e3/finite_bank/online'
    with np.load(ROOT/'results/e2e3/inference/POSTERIORS.npz',allow_pickle=False) as z:
        if not bool(z['complete']):raise RuntimeError('Complete online reference required')
        final=z['weights'][:,:,-1].copy()
    with np.load(source/'FINITE_GRID.npz',allow_pickle=False) as z:grid={k:z[k].copy() for k in z.files}
    with np.load(source/'FINITE_BANK.npz',allow_pickle=False) as z:
        bank=z['data'].copy();tx=z['tx'].copy();rx=z['rx'].copy()
    out=ROOT/'results/e2e3/posterior_endpoint_audit';out.mkdir(parents=True,exist_ok=False)
    book=BudgetBook();rows=[];paid={}
    for parent,repeat in np.ndindex(final.shape[:2]):
        i=int(np.argmax(final[parent,repeat]))
        if i not in paid:
            world=World(grid['support'][i].reshape(2,2),float(grid['epsinf'][i]),float(grid['ell'][i]),grid['pose'][i])
            solver=Evaluator(book,'E2E3',h=.17,n=384,grid_refinement=3,tx=tx,rx=rx)
            ev=solver.evaluate(world)
            paid[i]=dict(state_index=i,white_difference=float(np.linalg.norm(ev.data-bank[i])/.01),residual=ev.state.residual)
            np.savez_compressed(out/f'ROW_{i:04d}.npz',primary=bank[i],independent=ev.data,state_index=i,support=world.support,pose=world.pose)
            del ev;solver.clear();gc.collect()
        rows.append(dict(parent=parent,repeat=repeat,**paid[i]))
        result=dict(status='RUNNING',protocol=protocol,rows=rows,unique_paid_states=len(paid),truth_access=False,posterior_changed=False)
        (out/'RESULTS.json').write_text(json.dumps(result,indent=2,default=json_default)+'\n')
    result.update(status='COMPLETE',decision='PASS' if all(r['white_difference']<=.1 for r in rows) else 'FAIL',
        max_white_difference=max(r['white_difference'] for r in rows),uniform_error_bound=False,posterior_perturbation_certificate=False)
    (out/'RESULTS.json').write_text(json.dumps(result,indent=2,default=json_default)+'\n')
    print(json.dumps(result,default=json_default))

if __name__=='__main__':main()
