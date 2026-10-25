"""Paid independent-grid Jacobian uncertainty, followed by parent mechanism review."""
from pathlib import Path
import json
import numpy as np
from wslam.accounting import BudgetBook,json_default
from wslam.forward import Evaluator
from wslam.world import World
from wslam.physics import acquisition,packj
from wslam.information import generalized_retention,svd_audit,dynamic_information,schur
from scipy.linalg import block_diag

root=Path(__file__).resolve().parents[1]
out=root/'results/e1/rank_audit';out.mkdir(parents=True,exist_ok=False)
book=BudgetBook();tx,rx=acquisition(outer=True)
base=World(np.ones((2,2),bool),np.array([[3.,4.],[3.,4.]]),.08,(.2,2.))
fine=Evaluator(book,'E0E1',h=.17,n=384,grid_refinement=3,tx=tx,rx=rx)
ev=fine.evaluate(base,1.,jacobian=True,kind='left_right')
old=np.load(root/'results/e1/dynamic/pose0_f1.npz')
j0=old['J'].reshape(4,16,7);j1=ev.jacobian.reshape(4,16,7)
np.savez_compressed(out/'JACOBIANS.npz',primary=j0,independent=j1,primary_data=old['data'],independent_data=ev.data)
rows=[]
for label,sources,receivers in [('single_port',[0],[8]),('two_port',[0,1],[8]),('full_aperture',[0,1],list(range(16)))]:
    a0=packj(j0[np.ix_(sources,receivers,range(2))].reshape(-1,2))
    a1=packj(j1[np.ix_(sources,receivers,range(2))].reshape(-1,2))
    b0=packj(j0[np.ix_(sources,receivers,[2,3])].reshape(-1,2))
    b1=packj(j1[np.ix_(sources,receivers,[2,3])].reshape(-1,2))
    ea=float(np.linalg.norm(a1-a0,2));eb=float(np.linalg.norm(b1-b0,2))
    # Measured grid discrepancy is an empirical uncertainty, not a rigorous bound.
    audit0=generalized_retention(a0,b0,atol=max(ea,eb))
    audit1=generalized_retention(a1,b1,atol=max(ea,eb))
    rows.append(dict(case=label,wall_operator_grid_difference=ea,pose_operator_grid_difference=eb,
        empirical_threshold=max(ea,eb),primary_A_singular=np.linalg.svd(a0,compute_uv=False),
        primary_B_singular=np.linalg.svd(b0,compute_uv=False),
        retained_A_rank=audit0['wall_rank_audit']['rank'],retained_B_rank=audit0['nuisance_rank_audit']['rank'],
        primary_retention=audit0['eigenvalues'],independent_retention=audit1['eigenvalues'],
        interpretation='empirical rank sensitivity, no continuum error bound'))
intrinsic=[]
aa=[];bb=[]
for s in (0,1):
    aa.append(packj(j0[s,[8],:2]))
    bb.append(packj(j0[s,[8],2:7]))
a=np.vstack(aa);b=np.vstack(bb);bind=block_diag(*bb)
d=np.zeros((2,10));d[:,:2]=-np.eye(2);d[:,5:7]=np.eye(2)
q=.03**2*np.eye(2)
ds=dynamic_information(a,bind,d,q)
rs=generalized_retention(a,b)
intrinsic.append(dict(case='minimal_shared_unknown_intrinsic',
    names=['pose_x','pose_y','target_contrast','target_loss','target_radius'],
    B_rank=svd_audit(b)['rank'],rho_shared=rs['eigenvalues'],
    dynamic_D2_formula=ds['formula'],motion_information=ds['motion_information'],
    independent_pose_and_intrinsic_data_information=ds['data_only_information'],
    null_motion_rank=ds['null_motion_rank_audit']['rank'],
    interpretation='free intrinsic material can undo the known-intrinsic two-port information recovery; motion prior on pose alone need not help'))
np.savez_compressed(out/'MINIMAL_INTRINSIC_D2.npz',A=a,B_shared=b,B_independent=bind,D=d,Q=q,
                    information=ds['information'],null_motion_projector=ds['null_motion_projector'])
profile=json.loads((root/'results/e1/dynamic/RESULTS.json').read_text())['rows']
finite=[r for r in profile if r['case']=='finite_compensation_profile' and abs(r['wall_step'])==.01]
result=dict(status='NUMERICALLY_VERIFIED_LOCAL_CONTROLLED_MECHANISM',rows=rows,intrinsic_negative_control=intrinsic,finite_profiles=finite,
    residual=float(ev.state.residual),new_cost=book.totals('E0E1'),
    general_unknown_wall_reconstruction=False,global_nonidentifiability_proved=False,
    empirical_grid_uncertainty_only=True)
(out/'RESULTS.json').write_text(json.dumps(result,indent=2,default=json_default,allow_nan=False))
print(json.dumps(result,default=json_default,allow_nan=False))
