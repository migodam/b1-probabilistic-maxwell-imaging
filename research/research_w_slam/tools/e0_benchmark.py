"""First real fixed-contract full-wave cost and interface measurement."""
from pathlib import Path
import json,time,resource
import numpy as np
from wslam.accounting import BudgetBook
from wslam.world import World
from wslam.forward import Evaluator,volume_extraction
from wslam.physics import paired_data_decomposition,paired_feedback_residual
root=Path(__file__).resolve().parents[1]
if not (root/'EXPERIMENT_CONTRACT.json').exists():raise RuntimeError('Freeze first')
out=root/'results/e0/benchmark';out.mkdir(parents=True,exist_ok=False)
book=BudgetBook();solver=Evaluator(book,'E0E1',h=.22,n=128,grid_refinement=1)
w=World(np.ones((2,2),bool),np.array([[3.,4.],[3.,4.]]),.08,(.2,2.))
t=time.perf_counter();wall=solver.evaluate(w,include_target=False);full=solver.evaluate(w,jacobian=True,kind='left_right')
total=full.data;dtarget=total-wall.data
numeric=paired_feedback_residual(full.backend,wall.state,full.state)
np.savez_compressed(out/'ARRAYS.npz',wall=wall.data,total=total,target_induced=dtarget,J=full.jacobian,points=full.state.material_points,epsilon=full.state.epsilon,field=full.state.field,mesh_p=full.backend.mesh.p,mesh_t=full.backend.mesh.t)
row=dict(ndofs=full.backend.ndofs,wall_seconds=time.perf_counter()-t,counts=book.totals('E0E1'),maximum_residual=max(wall.state.residual,full.state.residual),
         point_volume_white=float(np.linalg.norm(volume_extraction(full)-full.data)/.01),decomposition_relative=float(np.linalg.norm(total-wall.data-dtarget)/np.linalg.norm(total)),
         feedback_equation=numeric,world=w.metadata(),jacobian_column_names=full.names,
         use='cost/interface benchmark only; not imaging or posterior performance')
(out/'RESULTS.json').write_text(json.dumps(row,indent=2,default=lambda a:a.item() if isinstance(a,np.generic) else np.asarray(a).tolist()))
print(json.dumps(row,default=lambda a:a.item() if isinstance(a,np.generic) else np.asarray(a).tolist()),flush=True)
