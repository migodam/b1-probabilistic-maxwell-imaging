import json
from pathlib import Path
import numpy as np
from lxz.local2d import MaterialPrior
ROOT=Path('/Volumes/migodam\'s-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-finite-feedback')
OUT=ROOT/'research/finite_feedback_fiber_diagnostic_20261016'
R=OUT/'results'
def whiten(z,sigma):
    return np.concatenate((z.real.T,z.imag.T),axis=1).ravel()*np.sqrt(2)/sigma
def geom(C,chart,u,sigma):
    D=chart.analytic_derivative(u);DM=np.r_[D.real,D.imag]/24
    _,s,Vh=np.linalg.svd(DM,full_matrices=False);rank=int(np.sum(s>s[0]*1e-7));E=Vh[:rank].T/s[:rank]
    vj=C.reshape(256,576)@D
    A=np.column_stack([whiten(vj[:,j].reshape(32,8),sigma) for j in range(D.shape[1])])
    _,sr,V=np.linalg.svd(A@E,full_matrices=False);v=E@V[-1];v/=np.linalg.norm(v)
    if v[np.argmax(abs(v))]<0:v=-v
    return dict(D=D,A=A,v=v,rank=rank,sr=sr)
def angle(x,y):
    c=float(np.real(np.vdot(x,y))/(np.linalg.norm(x)*np.linalg.norm(y)))
    return dict(cosine=c,degrees=float(np.degrees(np.arccos(np.clip(c,-1,1)))))
def correct(A,v,r,cap):
    Q=np.eye(len(v))-np.outer(v,v)/(v@v)
    d=-Q@np.linalg.pinv(A@Q,rcond=1e-8)@r
    return d*min(1,cap/max(np.linalg.norm(d),1e-300))
report=dict(scope='saved array and decoder algebra replay only; no Maxwell import, factorization, solve, or ledger access',new_physics_states=0,new_physics_RHS=0,anchors=[])
for key in ['250033_weak','250034_ring']:
    bg=key.split('_')[-1]
    meta=json.loads((OUT/'inputs'/f'{key}.json').read_text());chart=MaterialPrior(meta['generator'],bg)
    a=np.load(OUT/'inputs'/f'{key}.npz');sigma=float(a['sigma'])
    base=np.load(R/'endpoints'/f'{key}_base.npz');x0=base['chi'];pred0=base['prediction']
    saved=np.load(R/f'{key}_CHART.npz');u0=chart.coordinates();decoded=chart.decode(u0)[0]
    C=np.load(R/f'{key}_JACOBIAN.npz')['C'];g0=geom(C,chart,u0,sigma)
    rec=dict(anchor=key,decode_anchor_max_abs=float(np.max(abs(decoded-x0))),decode_anchor_RMS=float(np.linalg.norm(decoded-x0)/24),saved_u0_max_abs=float(np.max(abs(u0-saved['u0']))),initial_v_max_abs=float(np.max(abs(g0['v']-saved['v']))),steps=[])
    eps=.0005;v0=g0['v'];xp=chart.decode(u0+eps*v0)[0];xm=chart.decode(u0-eps*v0)[0]
    d2base=(xp-2*x0+xm)/eps**2;d2decoded=(xp-2*decoded+xm)/eps**2
    rec['curvature_center_stencil_relative_difference']=float(np.linalg.norm(d2base-d2decoded)/max(np.linalg.norm(d2decoded),1e-300))
    rec['curvature_center_material_RMS_contamination']=float(np.linalg.norm(d2base-d2decoded)/24)
    u=u0.copy();prev=v0.copy();prevD=g0['D'];g=g0
    for step in range(1,4):
        if step>1:
            C=np.load(R/f'{key}_updated_{step}_JACOBIAN.npz')['C'];g=geom(C,chart,u,sigma)
        vraw=g['v'].copy()
        wrongdot=float(np.real(np.vdot(g['D']@vraw,g0['D']@prev)))
        truedot=float(np.real(np.vdot(g['D']@vraw,prevD@prev)))
        wrongsign=-1 if wrongdot<0 else 1;truesign=-1 if truedot<0 else 1
        v=wrongsign*vraw
        ep=np.load(R/'endpoints'/f'{key}_pc_{step}_predict.npz');cp=np.load(R/'endpoints'/f'{key}_pc_{step}_correct.npz')
        up=u+.25*v;dx=correct(g['A'],v,whiten(ep['prediction']-pred0,sigma),.125);uc=up+dx
        xpred=chart.decode(up)[0];xcorr=chart.decode(uc)[0]
        rec['steps'].append(dict(step=step,frozen_source_sign=wrongsign,true_transport_sign=truesign,sign_choices_identical=wrongsign==truesign,frozen_source_orientation_dot=wrongdot,true_transport_orientation_dot=truedot,true_material_direction_angle=angle(g['D']@v,prevD@prev),initial_D_proxy_angle=angle(g['D']@v,g0['D']@prev),prior_coordinate_direction_angle=angle(v,prev),predictor_material_max_abs_replay_error=float(np.max(abs(xpred-ep['chi']))),corrector_material_max_abs_replay_error=float(np.max(abs(xcorr-cp['chi']))),corrector_delta_norm=float(np.linalg.norm(dx)),u_corrected=uc.tolist()))
        prevD=g['D'];prev=v;u=uc
    rec['actual_result_unchanged_by_orientation_fix']=all(s['sign_choices_identical'] for s in rec['steps'])
    rec['max_endpoint_replay_abs_error']=max(max(s['predictor_material_max_abs_replay_error'],s['corrector_material_max_abs_replay_error']) for s in rec['steps'])
    report['anchors'].append(rec)
report['all_sign_choices_identical']=all(a['actual_result_unchanged_by_orientation_fix'] for a in report['anchors'])
report['maximum_endpoint_replay_abs_error']=max(a['max_endpoint_replay_abs_error'] for a in report['anchors'])
Path('/tmp/fff_review_numeric.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
