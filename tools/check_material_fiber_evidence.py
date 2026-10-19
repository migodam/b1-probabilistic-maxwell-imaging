"""Zero-solve checks of the paid arrays, pushforward metric and MH algebra."""
from pathlib import Path
from collections import Counter
import sys,json,time
import numpy as np
from scipy.special import logsumexp
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from material_fiber.core import OUT,write,save,mpack
from p1_material.prior import MaterialPrior,topology,structural_change

def main():
    start=time.perf_counter();cfg=json.loads((OUT/'configs/TARGET.json').read_text());P=MaterialPrior(cfg,json.loads((ROOT/'research/p1_material_target_20261018/inputs/geometry.json').read_text())['points']);anchors=json.loads((OUT/'A/ANCHORS.json').read_text());key=anchors['P2_train_000']
    with np.load(OUT/f'cache/{key}.npz') as a:q=a['q'];Trecv=a['receiver_T'];field=a['field']
    D=P.derivative(q);M=np.r_[D.real,D.imag]/24;C=np.einsum('qp,ps->qsp',Trecv,field).reshape(256,576)
    from material_fiber.core import pack
    J=np.column_stack([pack((C@D[:,k]).reshape(32,8))*np.sqrt(2)/cfg['noise']['sigma'] for k in range(66)])
    H=np.diag(np.r_[np.ones(64),np.zeros(2)])
    def mat(M,J,H):
        Q,s,Vh=np.linalg.svd(M,full_matrices=False);T=Vh.T/s;K=J@T;A=K.T@K+T.T@H@T;W=np.linalg.solve(np.eye(66)+.02**2*A,np.eye(66));return Q@W@Q.T,W,A
    old,W,A=mat(M,J,H);change=np.diag(np.linspace(.5,2,66));new,_,_=mat(M@change,J@change,change.T@H@change)
    coordinate_error=float(np.linalg.norm(old-new)/np.linalg.norm(old))
    B=A+.01*np.eye(66);WB=np.linalg.solve(np.eye(66)+.02**2*B,np.eye(66));identity=np.linalg.norm(WB-W+.02**2*WB@(B-A)@W)
    rng=np.random.default_rng(2026101980);mh=[]
    for _ in range(100):
        a=P.draw(rng)[0];b=P.local_prior_reversible(a,rng,beta=.15);rho=np.sqrt(1-.15**2)
        lf=-np.sum((b[:64]-rho*a[:64])**2)/(2*.15**2);lr=-np.sum((a[:64]-rho*b[:64])**2)/(2*.15**2)
        mh.append(float(P.logprior(b)-P.logprior(a)+lr-lf))
    def reflected_logpdf(y,x,sd=.04):
        k=np.arange(-8,9);z=np.r_[y-x+2*k,y+x+2*k];return logsumexp(-.5*(z/sd)**2)-np.log(sd*np.sqrt(2*np.pi))
    prop_error=max(abs(reflected_logpdf(a,b)-reflected_logpdf(b,a)) for a,b in [(0.,.04),(1.,.98),(.2,.8),(.6,.4)])
    with np.load(OUT/'inputs/background/weak.npz') as a:G=a['G'];S=a['S'];incident=a['incident']
    maxima=Counter();n=0
    for stage in 'ABCDE':
        for path in (OUT/stage/'endpoints').glob('*.npz'):
            with np.load(path) as a:chi=a['chi'];current=a['current'];ef=a['field'];pred=a['prediction'];qq=a['q'];back=float(a['backward_residual'])
            if not P.legal(qq) or not np.array_equal(P.decode(qq)[0],chi):raise RuntimeError('material support audit failed')
            maxima['shared_material_error']=max(maxima['shared_material_error'],float(np.linalg.norm(current-chi[:,None]*ef)/max(1e-30,np.linalg.norm(current))))
            maxima['field_error']=max(maxima['field_error'],float(np.linalg.norm(ef-incident-G@current)/max(1e-30,np.linalg.norm(ef))))
            maxima['data_error']=max(maxima['data_error'],float(np.linalg.norm(pred-S@current)/max(1e-30,np.linalg.norm(pred))))
            maxima['backward_residual']=max(maxima['backward_residual'],back);n+=1
    package_checks=[]
    for path in (OUT/'MaterialHypothesisPackage').glob('*.json'):
        p=json.loads(path.read_text());hs=p['hypotheses'];ps=[]
        for h in hs:
            if h['posterior_mass'] is not None or not h['probability_status'].endswith('CANDIDATE_ONLY'):raise RuntimeError('unsupported probability claim')
            with np.load(OUT/'MaterialHypothesisPackage'/h['material_array']) as a:ps.append(a['chi'].copy())
        minimum=min([P.material_rms(a,b) for i,a in enumerate(ps) for b in ps[:i]]+[1e9]);package_checks.append(dict(parent=p['observation_id'],count=len(hs),minimum_material_distance=None if minimum==1e9 else minimum,probability_mass_null=True))
    stability=[]
    for path in (OUT/'MaterialHypothesisPackage').glob('*.json'):
        p=json.loads(path.read_text());hs=p['hypotheses']
        if not hs:continue
        with np.load(OUT/'MaterialHypothesisPackage'/hs[0]['material_array']) as a:q0=a['q'].copy()
        for h in hs[1:]:
            with np.load(OUT/'MaterialHypothesisPackage'/h['material_array']) as a:qq=a['q'].copy();frac=a['fraction'].copy()
            changes={str(th):topology(P.decode(q0)[1]['fraction'],threshold=th)['event']!=topology(frac,threshold=th)['event'] for th in [.25,.5,.75]}
            st=structural_change(P,q0,qq);stability.append(dict(parent=p['observation_id'],key=h['endpoint_key'],structure_pass=st['task_structure_pass'],topology_change_by_threshold=changes,structure_and_middle_topology=bool(st['task_structure_pass'] and changes['0.5'])))
    write(OUT/'TOPOLOGY_STABILITY.json',stability)
    result=dict(material_coordinate_invariance_relative_error=coordinate_error,soft_filter_resolvent_identity_absolute_error=float(identity),pCN_prior_cancellation_max=max(abs(x) for x in mh),reflected_property_symmetry_max=prop_error,real_physical_endpoints_audited=n,physical_equations=dict(maxima),candidate_package_checks=package_checks,zero_new_states=0,zero_new_RHS=0,zero_new_LU=0,wall_seconds=time.perf_counter()-start,meaning='cached algebra and full-array consistency; not new solver performance or posterior calibration')
    write(OUT/'SCIENTIFIC_NUMERICAL_CHECKS.json',result)
    if coordinate_error>1e-5 or max(maxima.values())>1e-10 or max(abs(x) for x in mh)>1e-10:raise RuntimeError('scientific evidence consistency failed')
    print(json.dumps(result))
if __name__=='__main__':main()
