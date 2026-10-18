"""Cached linear algebra and material decoder only. No GT or physical solver."""
from pathlib import Path
import sys,json,time
import numpy as np
from scipy.optimize import lsq_linear
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p1_material.prior import MaterialPrior,structural_change
from p1_material.proposals import frozen_latent_step
SRC=ROOT/'research/p1_material_target_20261018';OUT=ROOT/'research/p1_zero_physics_audit_20261018'
READS=[]

def readjson(path):READS.append(str(path.relative_to(ROOT)));return json.loads(path.read_text())
def arrays(path,keys):
    if 'offline_materials' in path.parts:raise RuntimeError('OFFLINE_INPUT_FORBIDDEN')
    READS.append(str(path.relative_to(ROOT)))
    with np.load(path,allow_pickle=False) as z:return {k:z[k].copy() for k in keys}
def plain(v):
    if isinstance(v,np.ndarray):return v.tolist()
    if isinstance(v,np.generic):return v.item()
    if isinstance(v,dict):return {k:plain(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [plain(x) for x in v]
    return v
def write(name,value):
    path=OUT/name;path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(plain(value),ensure_ascii=False,indent=2,allow_nan=False)+'\n')
def pack(z):return np.concatenate((z.real.T,z.imag.T),axis=1).ravel()
def matpack(z):return np.r_[z.real,z.imag]/24


def load_online():
    cfg=readjson(SRC/'configs/ONLINE_TARGET.json');P=MaterialPrior(cfg,readjson(SRC/'inputs/geometry.json')['points']);states=[]
    rows=[json.loads(x) for x in (SRC/'PHYSICS_LEDGER.jsonl').read_text().splitlines()];READS.append(str((SRC/'PHYSICS_LEDGER.jsonl').relative_to(ROOT)))
    for row in rows:
        if row['event']=='ENDPOINT_PERSISTED' and row['purpose'].startswith(('P1_obs_A','P1_obs_B')):
            state=arrays(ROOT/row['array'],['q','chi','prediction']);state['path']=row['array'];state['purpose']=row['purpose'];states.append(state)
    jac=[]
    for path in sorted((SRC/'results/physical/jacobians').glob('*.npz')):
        z=arrays(path,['q','A','C','D']);z['path']=str(path.relative_to(ROOT));jac.append(z)
    def J(q):return next(x for x in jac if np.array_equal(x['q'],q))
    return cfg,P,states,J


def main():
    if (OUT/'ONLINE_AUDIT_CLOSED.json').exists():raise RuntimeError('online audit sealed; do not overwrite')
    start=time.perf_counter();protocol=readjson(OUT/'configs/AUDIT_PROTOCOL.json');cfg,P,states,J=load_online();sigma=cfg['noise']['sigma']
    aqual=readjson(SRC/'results/physical/P1_obs_A_qualification.json');q=np.array(aqual['q']);anchor=next(x for x in states if x['purpose'].startswith('P1_obs_A') and np.array_equal(x['q'],q));j=J(q);A,D,C=j['A'],j['D'],j['C'];M=np.r_[D.real,D.imag]/24
    with np.load(SRC/'inputs/online/P1_obs_A.npz') as z:y=z['observation'].copy()
    READS.append(str((SRC/'inputs/online/P1_obs_A.npz').relative_to(ROOT)));r=pack(anchor['prediction']-y)*np.sqrt(2)/sigma
    U,s,Vh=np.linalg.svd(M,full_matrices=False);rank=sum(s>s[0]*protocol['primary_material_SVD_threshold']);T=Vh[:rank].T/s[:rank];Q=U[:,:rank];K=A@T
    Uk,sk,Vkh=np.linalg.svd(K,full_matrices=False);_,sa,Vha=np.linalg.svd(A,full_matrices=False)
    Hprior=np.diag(np.r_[np.ones(64),np.zeros(2)]);B=np.r_[A,np.eye(66)[:64]];_,sp,Vph=np.linalg.svd(B@T,full_matrices=False)
    check=np.column_stack([pack((C@D[:,k]).reshape(32,8))*np.sqrt(2)/sigma for k in range(66)])
    _,aux=P.decode(q);sat=aux['fraction'][P.mask];stats=dict(cells_in_ROI=int(sum(P.mask)),exact_host=int(sum(sat==0)),exact_pure=int(sum(abs(sat-1)<1e-14)),transition=int(sum((sat>0)&(sat<1-1e-14))))
    rawweak=[]
    for i in range(10):
        v=Vha[-1-i];rawweak.append(dict(order=i,unit_latent_material_RMS=float(np.linalg.norm(M@v)),unit_latent_white_response=float(np.linalg.norm(A@v)),white_response_per_material_RMS=float(np.linalg.norm(A@v)/max(np.linalg.norm(M@v),1e-300))))
    sensitivity=[]
    for tol in protocol['material_SVD_relative_thresholds']:
        keep=s>s[0]*tol;TT=Vh[keep].T/s[keep];ss=np.linalg.svd(A@TT,compute_uv=False);sensitivity.append(dict(relative_threshold=tol,material_rank=int(sum(keep)),weakest_response_per_unit_material=float(ss[-1]),strongest=float(ss[0])))
    directions={};trials=[];family={"raw_latent_weak":Vha[-3:][::-1].T,"material_metric_weak":T@Vkh[-3:][::-1].T,"prior_aware_material_weak":T@Vph[-3:][::-1].T}
    for name,vectors in family.items():
        for k in range(3):
            v=vectors[:,k].copy();v/=np.linalg.norm(M@v)
            # Current loss coordinate is at its upper bound: choose inward half-ray.
            if v[65]>0:v=-v
            key=f'{name}_{k}';directions[key]=v
            for distance in protocol['finite_decoder_RMS_diagnostics']:
                new,matched=frozen_latent_step(P,q,v,distance);x,_=P.decode(new);dq=new-q;delta=x-anchor['chi'];material_linear=pack((C@delta).reshape(32,8))*np.sqrt(2)/sigma;coordinate_linear=A@dq
                trials.append(dict(direction=key,target_RMS=distance,actual_RMS=P.material_rms(x,anchor['chi']),radius_matched=bool(matched),latent_Gaussian_displacement=float(np.linalg.norm(dq[:64])),latent_prior_log_change=P.logprior(new)-P.logprior(q),predicted_DeltaPhi_material_linear=float(r@material_linear+.5*np.dot(material_linear,material_linear)),predicted_DeltaPhi_coordinate_linear=float(r@coordinate_linear+.5*np.dot(coordinate_linear,coordinate_linear)),material_linear_white_drift=float(np.linalg.norm(material_linear)),coordinate_linear_white_drift=float(np.linalg.norm(coordinate_linear)),decoder_tangent_error_RMS=float(np.linalg.norm(matpack(delta)-M@dq)),structural_change=structural_change(P,q,new),prior_legal=P.legal(new),Maxwell_validated=False,posterior_mass=None))
                directions[key+f'_s{distance}_q']=new;directions[key+f'_s{distance}_chi']=x
    g=A.T@r+Hprior@q;pg=g.copy();pg[65]=max(g[65],0) if q[65]>=1 else g[65]
    metrics=dict(target_id=cfg['target_id'],source_anchor=anchor['path'],material_singular_values=s,latent_response_singular_values=sa,material_normalized_response_singular_values=sk,prior_aware_material_singular_values=sp,material_condition_number=float(s[0]/s[-1]),latent_response_condition_number=float(sa[0]/sa[-1]),material_normalized_response_condition_number=float(sk[0]/sk[-1]),rank_threshold_sensitivity=sensitivity,decoder_plateaus=stats,weak_latent_directions=rawweak,linear_task_drift_at_0p02=sk*.02,proper_prior_local_drift_at_0p02=sp*.02,gradient_likelihood=A.T@r,gradient_posterior=g,projected_gradient_inf=float(max(abs(pg))),loss_upper_boundary_gradient=float(g[65]),cached_chain_rule_relative_error=float(np.linalg.norm(check-A)/np.linalg.norm(A)),no_exact_fiber_certificate=True,finite_Maxwell_region='NOT_EVALUATED',new_Maxwell_states=0,new_complex_RHS=0)
    write('A_ONLINE_GEOMETRY.json',metrics);write('A_DECODER_ONLY_TRIALS.json',trials)
    np.savez_compressed(OUT/'A_ONLINE_GEOMETRY.npz',q=q,chi=anchor['chi'],prediction=anchor['prediction'],residual=r,M=M,A=A,C=C,Q=Q,T=T,K=K,material_response_singular_values=sk,material_response_right_vectors=Vkh,**directions)
    # Reconstruct the old B trajectory using only measured endpoints and saved Js.
    with np.load(SRC/'inputs/online/P1_obs_B.npz') as z:yb=z['observation'].copy()
    READS.append(str((SRC/'inputs/online/P1_obs_B.npz').relative_to(ROOT)));bs=[x for x in states if x['purpose'].startswith('P1_obs_B')];brows=[]
    for i,st in enumerate(bs):
        qb=st['q'];Ab=J(qb)['A'];rb=pack(st['prediction']-yb)*np.sqrt(2)/sigma;prior=np.r_[qb[:64],0.,0.];grad=Ab.T@rb+prior;lam=max(1e-6,.001*float(np.median(np.diag(Ab.T@Ab)+np.diag(Hprior))))
        MM=np.r_[Ab,Hprior,np.sqrt(lam)*np.eye(66)];rhs=np.r_[-rb,-prior,np.zeros(66)];lo=np.r_[np.full(64,-np.inf),-qb[64:]];hi=np.r_[np.full(64,np.inf),1-qb[64:]]
        fit=lsq_linear(MM,rhs,bounds=(lo,hi),tol=1e-10,max_iter=1000);raw=fit.x;alpha=min(1,.5/max(np.linalg.norm(raw[:64]),1e-30),.2/max(np.linalg.norm(raw[64:],np.inf),1e-30));obj=.5*rb@rb+.5*qb[:64]@qb[:64]
        row=dict(index=i,source=st['path'],Phi=float(.5*rb@rb),Gaussian_prior_penalty=float(.5*qb[:64]@qb[:64]),projected_gradient_inf=float(max(abs(grad))),unscaled_GN_Gaussian_step=float(np.linalg.norm(raw[:64])),uniform_postsolve_scale=float(alpha),damping=lam,raw_GN_linear_predicted_objective=float(.5*np.linalg.norm(rb+Ab@raw)**2+.5*np.linalg.norm(qb[:64]+raw[:64])**2),active_lower=np.flatnonzero(qb[64:]<=1e-10),active_upper=np.flatnonzero(qb[64:]>=1-1e-10))
        if i+1<len(bs):
            new=bs[i+1];d=new['q']-qb;rn=pack(new['prediction']-yb)*np.sqrt(2)/sigma;objn=.5*rn@rn+.5*np.linalg.norm(new['q'][:64])**2;pred=.5*np.linalg.norm(rb+Ab@d)**2+.5*np.linalg.norm(qb[:64]+d[:64])**2
            row.update(actual_objective_reduction=float(obj-objn),linear_predicted_reduction=float(obj-pred),trust_agreement_ratio=float((obj-objn)/(obj-pred)),actual_Gaussian_step=float(np.linalg.norm(d[:64])),full_secant_relative_linear_error=float(np.linalg.norm((rn-rb)-Ab@d)/np.linalg.norm(rn-rb)),material_RMS_step=P.material_rms(new['chi'],st['chi']))
        brows.append(row)
    write('B_OPTIMIZER_CACHED_AUDIT.json',dict(rows=brows,ground_truth_used=False,no_new_physics=True,unexecuted_GN_predictions_not_full_endpoints=True))
    forbidden=[x for x in sys.modules if x=='a30' or x.startswith(('a30.','a24.','a25.'))]
    write('ONLINE_INPUT_MANIFEST.json',dict(read_paths=sorted(set(READS)),offline_materials_reads=0,dataset_truth_endpoint_reads=0,backend_modules_loaded=forbidden,full_solver_calls=0))
    receipt=dict(event='ONLINE_CACHED_AUDIT_COMPLETED',wall_seconds=time.perf_counter()-start,new_Maxwell_states=0,new_complex_RHS=0,new_LU=0,decoder_trial_endpoints=len(trials),stored_Jacobians_reused=24,GT_used=False)
    write('ONLINE_AUDIT_CLOSED.json',receipt)
    with (OUT/'RESOURCE_LEDGER.jsonl').open('a') as f:f.write(json.dumps(receipt)+'\n')
    print(json.dumps(plain(dict(material_condition=metrics['material_condition_number'],data_condition=metrics['material_normalized_response_condition_number'],weakest_material_data=sk[-5:],weakest_prior_aware=sp[-5:],B_trust_ratios=[x.get('trust_agreement_ratio') for x in brows],receipt=receipt)),ensure_ascii=False))


if __name__=='__main__':main()
