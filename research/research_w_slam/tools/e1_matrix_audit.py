"""D1/D2/D3 evaluated on paid actual Maxwell Jacobians; no new PDE."""
from pathlib import Path
import json
import numpy as np
from wslam.information import dynamic_information,schur,conditional_innovation,svd_audit
from wslam.accounting import json_default
from wslam.physics import pack,packj
root=Path(__file__).resolve().parents[1];source=root/'results/e1'
out=source/'matrix_audit';out.mkdir(parents=True,exist_ok=False)
rows=[]
for label,path in [('D1',source/'dynamic/full_joint_base_right.npz'),
                   ('D2',source/'rank_audit/MINIMAL_INTRINSIC_D2.npz')]:
    z=np.load(path);a=z['A'];b=z['B_independent'];d=z['D'];q=z['Q']
    actual=dynamic_information(a,b,d,q)
    ordinary=schur(a,b,d.T@np.linalg.solve(q,d))
    error=float(np.linalg.norm(actual['information']-ordinary['information'])/max(1.,np.linalg.norm(ordinary['information'])))
    audit=svd_audit(b);v=audit['parameter_basis'];s=audit['singular_values'][:audit['rank']]
    # Deliberately incorrect comparator, only to exhibit why D2 is needed.
    gdag=(v/s**2)@v.T;c=audit['pseudoinverse']@a;r=a-b@c
    naive=r.T@r+(d@c).T@np.linalg.solve(q+d@gdag@d.T,d@c)
    np.savez_compressed(out/(label+'.npz'),A=a,B=b,D=d,Q=q,
        corrected=actual['information'],augmented_schur=ordinary['information'],naive_incorrect=naive,
        DC=actual['DC'],projected_conflict=actual['projected_conflict'])
    rows.append(dict(case=label,formula=actual['formula'],B_rank=audit['rank'],B_columns=b.shape[1],
        normalized_corrected_vs_augmented_error=error,
        naive_incorrect_vs_augmented_error=float(np.linalg.norm(naive-ordinary['information'])),
        absolute_information_eigenvalues=np.linalg.eigvalsh(actual['information']),
        scope='matrix identity on actual Maxwell Jacobian, not finite PDE uniqueness'))
z=np.load(source/'dynamic/pose0_f1.npz');j=z['J'].reshape(4,16,7)
j0=packj(j[0,:,:4]);j1=packj(j[1,:,:4])
prior=np.eye(4)*100.**-2
h0=prior+j0.T@j0;h1=h0+j1.T@j1
d3=conditional_innovation(h0,2,j1[:,:2],j1[:,2:])
s1=h1[:2,:2]-h1[:2,2:]@np.linalg.solve(h1[2:,2:],h1[2:,:2])
rng=np.random.default_rng(2026102512);noise0=rng.normal(size=j0.shape[0]);noise1=rng.normal(size=j1.shape[0])
cov0=np.linalg.solve(h0,np.eye(4));cov1=np.linalg.solve(h1,np.eye(4))
mu0=cov0@(j0.T@noise0);mu1=cov1@(j0.T@noise0+j1.T@noise1)
innovation=noise1-j1@mu0;innovation_cov=np.eye(len(noise1))+j1@cov0@j1.T
np.savez_compressed(out/'D3.npz',old_H=h0,new_H=h1,statistical_prior=prior,
    old_covariance=cov0,new_covariance=cov1,old_mean=mu0,new_mean=mu1,
    old_cross_covariance=cov0[:2,2:],new_cross_covariance=cov1[:2,2:],
    conditional_wall_jacobian=d3['conditional_jacobian'],innovation=innovation,
    innovation_covariance=innovation_cov,information_increment=d3['increment'])
rows.append(dict(case='D3',normalized_schur_increment_error=float(np.linalg.norm(d3['new_information']-s1)/max(1.,np.linalg.norm(s1))),
    increment_eigenvalues=np.linalg.eigvalsh(d3['increment']),
    statistical_prior_standard_deviation=100.,LM_as_prior=False,
    innovation_mahalanobis_squared=float(innovation@np.linalg.solve(innovation_cov,innovation)),
    old_wall_variance=np.diag(cov0)[:2],new_wall_variance=np.diag(cov1)[:2],
    scope='proper linear Gaussian audit using actual Maxwell Jacobians and independent simulated noise; not nonlinear posterior'))
raw=json.loads((source/'dynamic/RESULTS.json').read_text())['rows']
finite=[]
for row in raw:
    if row['case']!='finite_compensation_profile' or row['mode']!='independent_twoTX':continue
    poses=np.asarray(row['endpoint_pose']).reshape(2,2)
    difference=poses[1]-poses[0]
    finite.append(dict(step=row['wall_step'],data_white=row['endpoint_white'],
        actual_pose_difference=difference,motion_incompatibility_squared=float(difference@difference/.03**2),
        scope='cost of saved data-only compensating endpoint under zero-mean motion prior; not a new motion-regularized optimum'))
result=dict(rows=rows,finite_motion_conflicts=finite,new_pde_calls=0,
    labels='D1/D2/D3 matrix verifications and finite endpoint prior incompatibility are separate evidence')
(out/'RESULTS.json').write_text(json.dumps(result,indent=2,default=json_default,allow_nan=False))
print(json.dumps(result,default=json_default))
