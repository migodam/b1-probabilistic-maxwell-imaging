"""Synthetic and paid-cache A17-X audit. No physical backend or LU application.

Only the synthetic() function solves linear systems, all dimension 18 or less.
cached() reads selected already-paid array keys and performs matrix products.
"""
from pathlib import Path
import json, os, sys, time
import numpy as np
from scipy.linalg import null_space
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from a17x.finite import q_operator,galerkin,exchange,tail_bound,source_chart,row_fit
SAVED=ROOT/'research/a17x_f3_competition_audit'
OUT=Path(os.environ.get('B1_A17X_AUDIT_OUT',str(SAVED)))
legacy=ROOT.parent/'b1-research-f3/research/research_f3_quotient_exploration'
F3=Path(os.environ.get('B1_F3_ROOT',str(legacy if legacy.exists() else ROOT/'research/research_f3_quotient_exploration')))
THIN=SAVED/'evidence/reused_f3'

def plain(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    if isinstance(x,dict):return {k:plain(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [plain(v) for v in x]
    return x
def dump(path,x):path.write_text(json.dumps(plain(x),ensure_ascii=False,indent=2,allow_nan=False)+'\n')
def rel(a,b):return float(np.linalg.norm(a-b)/max(np.linalg.norm(a),np.linalg.norm(b),1e-30))

def synthetic():
    rng=np.random.default_rng(20261018)
    def cm(shape):return rng.normal(size=shape)+1j*rng.normal(size=shape)
    n,s,p,m=18,7,2,2
    g=cm((n,n));g*=.4/np.linalg.norm(g,2);x=.4+.2j+.02*cm((n,));l=np.eye(n)-x[:,None]*g
    u=np.eye(n)[:,:s];sensor=cm((m,n))/n;b=cm((s,p));P=np.linalg.solve(l,u);H=u.conj().T@g@P;O=sensor@P
    # Ordinary goal-oriented primal/dual basis; exact setup used only in synthetic controls.
    primal=P[:,:3];dual=np.linalg.solve(l.conj().T,sensor.conj().T)
    V=np.linalg.qr(np.c_[primal,dual])[0]
    candidate=(np.eye(n)-V@V.conj().T)@cm((n,3));C=np.linalg.qr(candidate)[0]
    pv,dp,bar,gamma,nu=exchange(l,V,C,u);pn=galerkin(l,np.c_[V,C])@u
    d=.04*cm((s,));D=u@np.diag(d)@u.conj().T
    full=sensor@np.linalg.solve(l-D@g,u@np.diag(d)@b)
    exact=O@q_operator(d,H)@b
    old=sensor@pv@q_operator(d,u.conj().T@g@pv)@b
    goal=sensor@pn@q_operator(d,u.conj().T@g@pn)@b
    xp=pv+dp;a17=sensor@xp@q_operator(d,u.conj().T@g@xp)@b
    result=dict(seed=20261018,scope='SYNTHETIC_ONLY',matrix_dimension=n,exact_port_vs_full=rel(exact,full),same_basis_A17_vs_direct_goal_ROM=rel(a17,goal),same_basis_propagator_identity=rel(xp,pn),finite_errors=dict(old=float(np.linalg.norm(old-full)),goal=float(np.linalg.norm(goal-full)),a17=float(np.linalg.norm(a17-full))))
    # Receiver darkness does not ensure multi-source feedback compatibility.
    bh=cm((s,2));hh=.08*cm((s,s));oh=cm((2,s))
    obs=[];mis=[]
    for i in range(s):
        tt=np.eye(s)[:,i,None]*bh;ht=hh@tt;ga=np.sum(ht*bh.conj(),axis=1)/np.sum(abs(bh)**2,axis=1)
        obs.append((oh@tt).ravel());mis.append((ht-ga[:,None]*bh).ravel())
    obs=np.column_stack(obs);joint=np.r_[obs,np.column_stack(mis)]
    result['generic_two_source_kernel']=[dict(rcond=tol,receiver_dark_complex_dimension=null_space(obs,rcond=tol).shape[1],joint_complex_dimension=null_space(joint,rcond=tol).shape[1]) for tol in [1e-7,1e-10,1e-12]]
    # Per-component fits are not a legal isotropic 3-vector material.
    ee=np.array([[1,2],[2,1],[1,1]],complex);tt=np.array([.1,.2,.3])[:,None]*ee
    _,rr,defined=row_fit(tt,ee,block_size=3)
    result['vector_block_counterexample']=dict(componentwise_fit_error=0.,shared_isotropic_block_residual=float(np.linalg.norm(rr)/np.linalg.norm(tt)),source_rank=int(np.linalg.matrix_rank(ee)))
    try:source_chart(np.zeros(2),np.zeros(2),np.eye(2))
    except ValueError as e:result['zero_field_chart_rejected']=str(e)
    # Sufficient perturbation bound can be unavailable while full feedback is invertible.
    qr=q_operator(np.array([2.]),np.zeros((1,1)))
    result['unavailable_bound_not_singularity']=dict(bound=tail_bound(qr,np.ones((1,1)),np.ones((1,1)),1.,.1),full_feedback_min_singular_value=1.,full_feedback_inverse=-1.)
    result['checks_pass']=bool(result['exact_port_vs_full']<1e-12 and result['same_basis_A17_vs_direct_goal_ROM']<1e-12 and all(r['joint_complex_dimension']==0 and r['receiver_dark_complex_dimension']>0 for r in result['generic_two_source_kernel']) and result['vector_block_counterexample']['shared_isotropic_block_residual']>.01)
    dump(OUT/'checks/INDEPENDENT_ALGEBRA.json',result)
    return result

def cached():
    records=json.loads((F3/'results/RECORDS.json').read_text());eligible=json.loads((F3/'results/ELIGIBLE_ENDPOINTS.json').read_text())['unique_material_endpoints']
    rows=[r for r in records if r['phase']=='F2_C1']
    for e in eligible:
        if not any(r['array']==e['array'] for r in rows):rows.append(next(r for r in records if r['array']==e['array']))
    contexts={};provenance=[];results=[];thin=OUT/'evidence/reused_f3';thin.mkdir(exist_ok=True)
    for bg in ['ring','weak']:
        case=f'p250073_{bg}_a1_snr20';summary=json.loads((F3/f'results/F15/{case}_FINAL.json').read_text());q=np.array(summary['q'])
        if (THIN/f'{case}_anchor.npz').exists():
            with np.load(THIN/f'{case}_anchor.npz') as z:a={k:z[k] for k in z.files}
            with np.load(THIN/f'{case}_jacobian.npz') as z:jac={k:z[k] for k in z.files}
            with np.load(THIN/f'{case}_background.npz') as z:background={k:z[k] for k in z.files}
        else:
            for path in (F3/'results/endpoints').glob(case+'_F15_*.npz'):
                with np.load(path) as z:
                    if np.array_equal(z['q'],q):a={k:z[k] for k in ['q','chi','current','field','prediction']};break
            for path in (F3/'results/jacobians').glob(case+'*.npz'):
                with np.load(path) as z:
                    if np.array_equal(z['q'],q):jac={k:z[k] for k in ['q','T','C','D']};break
            with np.load(F3/f'inputs/background/{bg}.npz') as z:background={k:z[k] for k in ['G','S']}
        G=background['G'];S=background['S'];L=np.eye(576)-a['chi'][:,None]*G
        contexts[case]=(a,jac,G,S,L)
        for label,content in [('anchor',a),('jacobian',jac),('background',background)]:np.savez_compressed(thin/f'{case}_{label}.npz',**content)
        provenance.append(dict(case=case,anchor='F3 matched F15 q',full_physics_prepaid=True,new_LU_or_RHS=0,keys_not_read=['LU','pivots'],exact_receiver_propagator_cached=True,full_R0U_or_H_cached=False))
    for index,row in enumerate(rows):
        a,jac,G,S,L=contexts[row['case']]
        path=THIN/f'endpoint_{index:02}.npz'
        if not path.exists():path=F3/Path(row['array']).relative_to('research/research_f3_quotient_exploration')
        with np.load(path) as z:end={k:z[k] for k in ['q','chi','current','field','prediction']}
        h=end['chi']-a['chi'];delta_y=end['prediction']-a['prediction'];k=end['current']-a['current'];T=L@k;GK=G@k
        exact_relation=rel(T,h[:,None]*(a['field']+GK))
        order=np.argsort(abs(h))[::-1];amp=abs(h[order]);jh=jac['C']@h;compression=[]
        for rank in [1,7,16,32,64,128,256,576]:
            omitted=h.copy();omitted[order[:rank]]=0;operr=float(np.max(abs(omitted)))
            compression.append(dict(ports=rank,deltaX_spectral_error=operr,relative_operator_error=operr/max(abs(h)),cached_linear_receiver_error=float(np.linalg.norm(jac['C']@omitted)),relative_cached_linear_receiver_error=float(np.linalg.norm(jac['C']@omitted)/max(np.linalg.norm(jh),1e-300)),finite_receiver_error=None,P0_exact_membership_of_truncation='NOT_ESTABLISHED; full original material retained'))
        paths=[]
        for t in [.25,.5,.75]:
            E=a['field']+t*GK;d,res,defined=row_fit(t*T,E)
            paths.append(dict(t=t,multi_source_row_fit_relative_residual=float(np.linalg.norm(res)/max(np.linalg.norm(t*T),1e-300)),fixed_phase_relative_defect=float(np.linalg.norm(np.imag(d/((1+.2j)/abs(1+.2j))))/max(np.linalg.norm(d),1e-300)),undefined_field_rows=int(sum(~defined)),synthetic_material_ratio_only=True,full_state_validated=False,P0_exact='NOT_ESTABLISHED'))
        rr=dict(index=index,case=row['case'],method=row['method'],array=row['array'],delta_Phi=row['delta_Phi'],eligible=any(e['array']==row['array'] for e in eligible),material_RMS=float(np.linalg.norm(h)/24),floating_nonzero_deltaX_diagonal=int(sum(abs(h)>0)),numerical_rank_diagonal={str(tol):int(sum(abs(h)>tol*max(abs(h)))) for tol in [1e-14,1e-10,1e-6]},ports_for_relative_operator_error={str(tol):int(sum(abs(h)>tol*max(abs(h)))) for tol in [.1,.01,.001]},exact_paid_state_relation_relative_error=exact_relation,cached_receiver_current_relation_error=rel(S@k,delta_y),linear_vs_finite_relative_error=rel((jac['C']@h).reshape(32,8),delta_y),compression=compression,current_straight_interpolation=paths)
        results.append(rr);np.savez_compressed(thin/f'endpoint_{index:02}.npz',**end)
    dump(OUT/'results/CACHED_PORT_AND_SOURCE_AUDIT.json',dict(scope='PAID_F3_CACHE_ALGEBRA_NOT_NEW_PHYSICS',fixed_commit='a323ef223731870b390f875c55f041aa8cc855b4',parent_count=1,backgrounds=2,actual_P1_target=None,rows=results,provenance=provenance,physical_calls=0,new_RHS=0,new_LU=0,GT_read=False,finite_truncation_certificate=None))
    return results

def main():
    for sub in ['checks','results','evidence']:(OUT/sub).mkdir(parents=True,exist_ok=True)
    if (OUT/'AUDIT_CLOSED.json').exists():raise RuntimeError('audit already closed; preserve results')
    t=time.perf_counter();a=synthetic();s=time.perf_counter();rows=cached();elapsed=time.perf_counter()-t
    original=json.loads((SAVED/'checks/algebra_results.json').read_text());receipt=json.loads((SAVED/'checks/REPLAY_RECEIPT.json').read_text())
    record=dict(scope='ZERO_NEW_PHYSICS',synthetic_extended_seconds=s-t,cached_algebra_and_IO_seconds=time.perf_counter()-s,total_extended_wall_seconds=elapsed,package_replay_wall_seconds=receipt['wall_seconds'],new_full_states=0,new_complex_RHS=0,new_Maxwell_LU=0,NN_training=0,Mie=0,new3D=0,package_checks=original,independent_checks=a,cached_audit='results/CACHED_PORT_AND_SOURCE_AUDIT.json',cached_rows=len(rows),physical_competition='NOT_RUN',P1_status='UNKNOWN_NOT_READY_IN_INSPECTED_ASSETS')
    dump(OUT/'ZERO_PHYSICS_ALGEBRA_CHECKS.json',record)
    with (OUT/'RESOURCE_LEDGER.jsonl').open('a') as f:
        f.write(json.dumps(plain(dict(event='AUDIT_COMPLETED',**record)))+'\n')
    dump(OUT/'AUDIT_CLOSED.json',dict(new_physics=0,new_RHS=0,closed=True,reason='P1_target_and_separate_new_physics_authorization_missing'))
    print(json.dumps(dict(checks_pass=a['checks_pass'],cached_rows=len(rows),wall_seconds=elapsed,states=0,RHS=0)))

if __name__=='__main__':main()
