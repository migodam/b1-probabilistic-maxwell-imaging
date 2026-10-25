"""Supplementary physical source/receiver swap, same mean units; no fitted gain."""
from pathlib import Path
import json,gc
import numpy as np
from wslam.accounting import BudgetBook,json_default
from wslam.forward import Evaluator
from wslam.world import World
from wslam.physics import acquisition
root=Path(__file__).resolve().parents[1];out=root/'results/e0/reciprocity'
out.mkdir(parents=True,exist_ok=False)
old=np.load(root/'results/e1/rank_audit/JACOBIANS.npz',allow_pickle=False)
tx,rx=acquisition(outer=True);book=BudgetBook();rows=[]
world=World(np.ones((2,2)),np.array([[3.,4.],[3.,4.]]),.08,(.2,2.))
for si,ri in [(0,8),(1,3),(2,12),(3,5)]:
    values=[]
    for n,refinement,key in [(256,2,'primary_data'),(384,3,'independent_data')]:
        solver=Evaluator(book,'E0E1',tx=rx[[ri]],rx=tx[[si]],n=n,grid_refinement=refinement)
        ev=solver.evaluate(world);mean=complex(ev.data[0,0]);original=complex(old[key][si,ri])
        values.append(mean)
        rows.append(dict(source_index=si,receiver_index=ri,n=n,field_refinement=refinement,
            forward_mean=original,swapped_mean=mean,reciprocity_absolute_white=float(abs(mean-original)/.01),
            residual=ev.state.residual,unit_amplitude=True,complex_scale_fitted=False,
            interpretation='complex transpose reciprocity, not conjugation'))
        del ev;solver.clear();gc.collect()
    rows[-1]['swapped_two_grid_white']=float(abs(values[0]-values[1])/.01)
result=dict(rows=rows,decision='PASS' if all(r['reciprocity_absolute_white']<=.1 for r in rows) else 'FAIL',
    tolerance='existing0.1 numerical output uncertainty; no new physical/statistical gate',
    new_source_receiver_pairs=8,continuous_PDE_error_certificate=False)
(out/'RESULTS.json').write_text(json.dumps(result,indent=2,default=json_default));print(json.dumps(result,default=json_default))
