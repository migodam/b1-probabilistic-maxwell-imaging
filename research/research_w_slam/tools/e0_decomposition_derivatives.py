"""Paid full-wave derivatives of wall-only and wall-dependent target response."""
from pathlib import Path
from dataclasses import replace
import json,gc
import numpy as np
from wslam.accounting import BudgetBook,json_default
from wslam.forward import Evaluator
from wslam.world import World
from wslam.physics import pack,packj,acquisition
root=Path(__file__).resolve().parents[1]
out=root/'results/e0/decomposition_derivatives';out.mkdir(parents=True,exist_ok=False)
book=BudgetBook();tx,rx=acquisition()
base=World(np.ones((2,2),bool),np.array([[3.,4.],[3.,4.]]),.08,(.2,2.))
full=Evaluator(book,'E0E1');wall=Evaluator(book,'E0E1')
f=full.evaluate(base,jacobian=True,kind='left_right')
w=wall.evaluate(base,jacobian=True,kind='left_right',include_target=False)
jfull=f.jacobian.copy();jwall=w.jacobian.copy();df=f.data-w.data
np.savez_compressed(out/'DECOMPOSITION.npz',wall=w.data,full=f.data,target_induced=df,
    full_J=jfull,wall_J=jwall,target_induced_J=jfull-jwall,names=np.asarray(f.names))
del f,w;rows=[]
for label,column,h in [('wall_left',0,.001),('pose_x',2,.002),('intrinsic_contrast',4,.001)]:
    ff=[];ww=[]
    for sign in (1.,-1.):
        if column==0:
            eps=np.asarray(base.epsinf).copy();eps[:,0]+=sign*h;v=replace(base,epsinf=eps)
        elif column==2:v=replace(base,pose=np.asarray(base.pose)+[sign*h,0])
        else:v=replace(base,contrast=base.contrast+sign*h)
        a=full.evaluate(v);b=wall.evaluate(v,include_target=False)
        ff.append(a.data.copy());ww.append(b.data.copy())
        np.savez_compressed(out/f'{label}_{sign:+}.npz',full=a.data,wall=b.data,delta=a.data-b.data)
    fd_f=(ff[0]-ff[1]).ravel()/(2*h);fd_w=(ww[0]-ww[1]).ravel()/(2*h)
    fd_delta=fd_f-fd_w;exact_delta=jfull[:,column]-jwall[:,column]
    rows.append(dict(parameter=label,step=h,
        full_derivative_relative_error=float(np.linalg.norm(fd_f-jfull[:,column])/max(np.linalg.norm(jfull[:,column]),1e-300)),
        wall_derivative_norm=float(np.linalg.norm(jwall[:,column])),
        target_induced_wall_dependence_norm=float(np.linalg.norm(exact_delta)),
        delta_derivative_relative_error=float(np.linalg.norm(fd_delta-exact_delta)/max(np.linalg.norm(exact_delta),1e-300)),
        target_induced_FD_white=float(np.linalg.norm(fd_delta)/.01),
        discrete_wall_pose_FD_aliasing_white=float(np.linalg.norm(fd_w)/.01) if column==2 else None))
    gc.collect()
# Vacuum is the same declared source/receiver convention, not fitted scaling.
vac=full.evaluate(replace(base,support=np.zeros((2,2),bool)),include_target=False)
vac_error=float(np.linalg.norm(vac.data-vac.state.direct)/.01) if hasattr(vac.state,'direct') else float(np.linalg.norm(vac.data-vac.backend._direct)/.01)
result=dict(rows=rows,vacuum_total_minus_analytic_white=vac_error,
    source_normalization='unit iH0^(1)/4, no fitted complex multiplier',
    decomposition_scope='full wave wall-only plus target-induced response, whose wall derivative is nonzero',
    state_cost_snapshot=book.totals('E0E1'))
(out/'RESULTS.json').write_text(json.dumps(result,indent=2,default=json_default,allow_nan=False))
print(json.dumps(result,default=json_default))
