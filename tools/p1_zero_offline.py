"""Isolated offline pair evaluator. Never constructs online proposals."""
from pathlib import Path
import sys,json,time
import numpy as np
from scipy.stats import norm
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p1_material.prior import MaterialPrior
SRC=ROOT/'research/p1_material_target_20261018';OUT=ROOT/'research/p1_zero_physics_audit_20261018'


def pack(z):return np.concatenate((z.real.T,z.imag.T),axis=1).ravel()
def plain(v):
    if isinstance(v,np.ndarray):return v.tolist()
    if isinstance(v,np.generic):return v.item()
    if isinstance(v,dict):return {k:plain(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [plain(x) for x in v]
    return v


def main():
    if not (OUT/'ONLINE_AUDIT_CLOSED.json').exists():raise RuntimeError('freeze online directions before GT access')
    if (OUT/'OFFLINE_EVALUATOR_CLOSED.json').exists():raise RuntimeError('offline evaluation already sealed')
    start=time.perf_counter();cfg=json.loads((SRC/'configs/ONLINE_TARGET.json').read_text());P=MaterialPrior(cfg,json.loads((SRC/'inputs/geometry.json').read_text())['points']);sigma=cfg['noise']['sigma']
    with np.load(OUT/'A_ONLINE_GEOMETRY.npz') as z:q=z['q'].copy();x=z['chi'].copy();pred=z['prediction'].copy();Q=z['Q'].copy();K=z['K'].copy();C=z['C'].copy();A=z['A'].copy()
    with np.load(SRC/'results/physical/endpoints/state_0000.npz') as z:qg=z['q'].copy();xg=z['chi'].copy();pg=z['prediction'].copy()
    with np.load(SRC/'inputs/online/P1_obs_A.npz') as z:y=z['observation'].copy()
    h=xg-x;hm=np.r_[h.real,h.imag]/24;delta=pack(pg-pred)*np.sqrt(2)/sigma;d=float(np.linalg.norm(delta));ra=pack(pred-y)*np.sqrt(2)/sigma;rg=pack(pg-y)*np.sqrt(2)/sigma
    matlin=pack((C@h).reshape(32,8))*np.sqrt(2)/sigma;latlin=A@(qg-q);U,s,Vh=np.linalg.svd(K,full_matrices=False);coord=Vh@(Q.T@hm);normh=float(np.linalg.norm(hm));ortho=hm-Q@(Q.T@hm)
    decomposition=[]
    for threshold in [.01,.001,.0001]:
        strong=s>s[0]*threshold;decomposition.append(dict(relative_data_threshold=threshold,strong_count=int(sum(strong)),weak_count=int(sum(~strong)),strong_energy_of_full_h=float(np.dot(coord[strong],coord[strong])/normh**2),weak_energy_of_full_h=float(np.dot(coord[~strong],coord[~strong])/normh**2),outside_local_prior_tangent_energy=float(np.dot(ortho,ortho)/normh**2)))
    taskweak=s*.02<=1
    def cosine(a,b):return float(np.dot(a,b)/(np.linalg.norm(a)*np.linalg.norm(b))) if np.linalg.norm(a)*np.linalg.norm(b)>0 else None
    phiA=float(.5*ra@ra);phiG=float(.5*rg@rg)
    result=dict(analysis_type='OFFLINE_ORACLE_PAIR_DIAGNOSTIC_NOT_DISCOVERY',target_id=cfg['target_id'],source_online_directions='A_ONLINE_GEOMETRY.npz frozen first',material_RMS=normh,material_complex_MAE=float(np.mean(abs(h))),Re_MAE=float(np.mean(abs(h.real))),Im_MAE=float(np.mean(abs(h.imag))),whitened_prediction_distance=d,observed_Phi_anchor=phiA,observed_Phi_GT=phiG,observed_log_likelihood_ratio_anchor_over_GT=phiG-phiA,fixed_pair_equal_prior_optimal_classification_error=float(norm.cdf(-d/2)),fixed_pair_error_semantics='prospective independent noisy repetition with these two now-fixed means only; not observed-data p-value, not composite-model power or posterior calibration',data_fitted_anchor=True,selection_adjusted_significance='NOT_AVAILABLE',both_noise_ball_compatible=bool(np.linalg.norm(ra)<=23.785449004057032 and np.linalg.norm(rg)<=23.785449004057032),noise_ball_compatibility_not_indistinguishability=True,finite_prior_secant_legal=P.legal(q) and P.legal(qg),latent_line_prior_legal='yes: Gaussian x convex [0,1]^2; no intermediate Maxwell evaluations',local_prior_tangent_energy=float(1-np.dot(ortho,ortho)/normh**2),strong_weak_decomposition=decomposition,task_snr_weak=dict(material_RMS=.02,white_response_threshold=1.,weak_count=int(sum(taskweak)),weak_energy_of_full_h=float(np.dot(coord[taskweak],coord[taskweak])/normh**2),outside_local_prior_tangent_energy=float(np.dot(ortho,ortho)/normh**2)),finite_response_terms=dict(actual_white_norm=d,full_material_J_h_norm=float(np.linalg.norm(matlin)),latent_J_delta_q_norm=float(np.linalg.norm(latlin)),finite_Maxwell_remainder_norm=float(np.linalg.norm(delta-matlin)),finite_decoder_secant_deviation_response_norm=float(np.linalg.norm(matlin-latlin)),actual_material_linear_cosine=cosine(delta,matlin),first_order_Maxwell_remainder_cosine=cosine(matlin,delta-matlin)),latent_point_logprior_anchor=P.logprior(q),latent_point_logprior_GT=P.logprior(qg),posterior_event_mass=None,posterior_reference='NOT_AVAILABLE',GT_used_online=False,new_solver_calls=0,new_RHS=0)
    (OUT/'OFFLINE_GT_PAIR_EVALUATION.json').write_text(json.dumps(plain(result),ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    np.savez_compressed(OUT/'OFFLINE_GT_PAIR_ARRAYS.npz',GT_material=xg,anchor_material=x,material_difference=h,GT_prediction=pg,anchor_prediction=pred,white_prediction_difference=delta,material_linear_prediction=matlin,latent_linear_prediction=latlin,prior_tangent_projected_material_difference=Q@(Q.T@hm),orthogonal_material_difference=ortho)
    receipt=dict(event='ISOLATED_OFFLINE_EVALUATOR_COMPLETED',wall_seconds=time.perf_counter()-start,new_Maxwell_states=0,new_complex_RHS=0,new_LU=0,GT_used_only_here=True,online_outputs_not_modified=True)
    (OUT/'OFFLINE_EVALUATOR_CLOSED.json').write_text(json.dumps(receipt,indent=2)+'\n')
    with (OUT/'RESOURCE_LEDGER.jsonl').open('a') as f:f.write(json.dumps(receipt)+'\n')
    print(json.dumps(plain(result),ensure_ascii=False))


if __name__=='__main__':main()
