"""Fixed legal material evaluation: exact ports, ordinary ROM/Krylov and exchange."""
import json,time
import numpy as np
from scipy.linalg import lu_factor,lu_solve
from .core import OUT,write,save
from a17x.finite import finite_state_exchange

METHODS=['full','exact_T','goal_ROM','Krylov','A17_exchange']

def orth(a):
    if a.shape[1]==0:return a
    u,s,_=np.linalg.svd(a,full_matrices=False);return u[:,s>s[0]*1e-10]

def small_solve(book,A,b,role):
    with book.span(role,small_LU=1,small_RHS=b.shape[1],small_dimension=A.shape[0]):
        if np.linalg.cond(A)>1e12:raise np.linalg.LinAlgError('small core ill-conditioned')
        return lu_solve(lu_factor(A),b)

def galerkin(book,L,V,b):
    return V@small_solve(book,V.conj().T@L@V,V.conj().T@b,'projected_core')

def approximate(t,anchor,q,method,V):
    model=t.backend();chi,_=t.prior.decode(q);delta=chi-anchor.state.chi;L=np.eye(576)-chi[:,None]*model.G;b=delta[:,None]*anchor.state.field
    k=galerkin(t.book,L,V,b);same_space=[];actions=0
    for iteration in range(2):
        residual=b-L@k
        # Identical pool/reuse rights; ordinary goal residual and signed exchange utility differ.
        pool=[residual,L@residual,L.conj().T@(model.S.conj().T@(model.S@residual))]
        t.book.counts['ROM_L_columns']+=24;t.book.counts['ROM_adjoint_L_columns']+=8;t.book.counts['ROM_receiver_columns']+=8
        incoming=[]
        for candidate in pool:
            C=orth(candidate-V@(V.conj().T@candidate))
            if C.shape[1]:incoming.append(C[:,:min(8,C.shape[1])])
        if not incoming:break
        scored=[]
        for C in incoming:
            newV=orth(np.c_[V,C]);kn=galerkin(t.book,L,newV,b);actions+=C.shape[1]
            if method=='A17_exchange':
                pv,inc,_,_,_=finite_state_exchange(L,V,C,b);exchange=pv+inc
                same_space.append(float(np.linalg.norm(exchange-kn)/max(1e-30,np.linalg.norm(kn))))
                prediction=anchor.state.data+model.S@kn;y=t.observation(t.parent)
                score=float(np.vdot(prediction-y,prediction-y).real/t.sigma**2)
            elif method=='goal_ROM':score=float(np.linalg.norm(anchor.T@(b-L@kn)))
            else:score=float(np.linalg.norm(b-L@kn))
            scored.append((score,newV,kn))
        _,V,k=min(scored,key=lambda a:a[0])
        if V.shape[1]>=32:break
    return anchor.state.data+model.S@k,V,dict(rank=V.shape[1],pool_columns_scored=actions,same_space_relative_errors=same_space,finite_state_residual=float(np.linalg.norm(b-L@k)/max(1e-30,np.linalg.norm(b))))

def stage_E(t):
    parents=json.loads((OUT/'configs/PARENTS.json').read_text());anchors=json.loads((OUT/'A/ANCHORS.json').read_text());br=[json.loads(x) for x in (OUT/'B/candidates.jsonl').read_text().splitlines()];panel=[]
    for r in parents:
        if r['split']=='heldout':continue
        p=t.load(anchors[r['parent']]);y=t.observation(r['parent'])
        if not t.qualified(p,y):continue
        chosen=[]
        for row in br:
            if row['parent']!=r['parent']:continue
            path=OUT/f'B/endpoints/{row["key"]}.npz'
            with np.load(path) as a:q=a['q'].copy();pred=a['prediction'].copy();chi=a['chi'].copy()
            if any(t.prior.material_rms(chi,x['chi'])<1e-6 for x in chosen):continue
            chosen.append(dict(key=row['key'],q=q,pred=pred,chi=chi))
            if len(chosen)==4:break
        if len(chosen)==4:panel.append((r['parent'],p,chosen))
        if len(panel)==4:break
    write(OUT/'E/FROZEN_QUERY_PANEL.json',dict(selection='first four qualified train/validation anchors; first four chronological material-distinct B queries, including failures; no truth or accuracy selection',parents=[dict(parent=p,anchor=a.key,queries=[q['key'] for q in qq]) for p,a,qq in panel],query_count=sum(len(x[2]) for x in panel)))
    rows=[];details=[]
    for parent,paid,queries in panel:
        t.parent=parent;model=t.backend();setups={};indices=np.flatnonzero(t.prior.mask);U=np.eye(576,dtype=complex)[:,indices]
        for method in METHODS:
            if method=='full':setups[method]=dict(anchor=None,warm={});continue
            start=time.perf_counter();a=t.evaluate(paid.q,f'E_{parent}_{method}_cold_anchor')
            setup=dict(anchor=a,warm={})
            if method=='exact_T':
                with t.book.physical('exact_port_setup'):
                    P=model.solve_anchor(U,a.state,role='exact_T_ports');setup.update(P=P,H=model.G[indices]@P,O=model.S@P)
            if method in ['goal_ROM','A17_exchange']:
                t.jacobian(a);setup['V']=orth(np.c_[a.state.current,a.T.conj().T[:,:8]])[:,:16]
            elif method=='Krylov':setup['V']=orth(a.state.current)[:,:8]
            setup['cold_setup_seconds']=time.perf_counter()-start;setups[method]=setup
        for repetition in range(3):
            order=METHODS[repetition:]+METHODS[:repetition]
            for query in queries:
                for method in order:
                    setup=setups[method];q=query['q'];chi=query['chi'];before=t.book.used['E'].copy();start=time.perf_counter();detail={};error=None
                    try:
                        with t.book.physical('fixed_material_'+method):
                            if method=='full':
                                if query['key'] not in setup['warm']:
                                    p=t.evaluate(q,f'E_full_{parent}_{query["key"]}_{repetition}');setup['warm'][query['key']]=p;pred=p.state.data
                                else:
                                    p=setup['warm'][query['key']];t.book.reserve_state(chi,resolution=24,purpose='E_full_factor_reuse');j=model.solve_anchor(chi[:,None]*model.incident,p.state,role='full_factor_reuse');pred=model.S@j
                                    t.book.counts['warm_full_G_columns']+=8;t.book.counts['warm_full_S_columns']+=8
                            elif method=='exact_T':
                                t.book.reserve_state(chi,resolution=24,purpose='E_exact_T_evaluation');a=setup['anchor'];d=(chi-a.state.chi)[indices]
                                Z=small_solve(t.book,np.eye(len(d))-d[:,None]*setup['H'],d[:,None]*a.state.field[indices],'exact_port_core');pred=a.state.data+setup['O']@Z
                            else:
                                t.book.reserve_state(chi,resolution=24,purpose='E_'+method+'_evaluation');pred,V,detail=approximate(t,setup['anchor'],q,method,setup['V']);setup['V']=V
                    except np.linalg.LinAlgError as e:error=repr(e);pred=None
                    approximation_seconds=time.perf_counter()-start;y=t.observation(parent);ref=query['pred'];white=float(np.linalg.norm(t.white(pred-ref))) if pred is not None else None
                    phi_error=float(abs(np.vdot(pred-y,pred-y).real-np.vdot(ref-y,ref-y).real)/t.sigma**2) if pred is not None else None
                    fallback=pred is None or white>.05 or phi_error>.05
                    if fallback:
                        f=t.evaluate(q,f'E_fallback_{parent}_{query["key"]}_{method}_{repetition}');delivered=f.state.data
                    else:delivered=pred
                    row=dict(parent=parent,query=query['key'],method=method,repetition=repetition,cache='cold_sequence' if repetition==0 else 'warm_factors_and_own_basis_no_prediction_memoization',approximation_seconds=approximation_seconds,total_seconds=time.perf_counter()-start,white_prediction_error=white,absolute_Phi_error=phi_error,full_fallback=fallback,error=error,delivered_white_error=float(np.linalg.norm(t.white(delivered-ref))),states=t.book.used['E']['STATE']-before['STATE'],RHS=t.book.used['E']['RHS']-before['RHS'],**detail)
                    save(OUT/f'E/evaluations/{parent}_{query["key"]}_{method}_{repetition}.npz',q=q,chi=chi,prediction=delivered,approximation=np.zeros((0,0),complex) if pred is None else pred,reference_prediction=ref)
                    rows.append(row)
        details.append(dict(parent=parent,setup_seconds={m:s.get('cold_setup_seconds',0) for m,s in setups.items()},exact_port_count=len(indices),latent_dimension=66,port_count_not_latent_dimension=True))
        write(OUT/'E/RESULTS.json',dict(rows=rows,setups=details,same_space_identity_is_control=True,certified_rejection=False));print(json.dumps(dict(stage='E',parent=parent,consumed=dict(t.book.used['E']))),flush=True)
    write(OUT/'E/SUMMARY.json',dict(rows=rows,setups=details,probability_claims=False,accuracy_controls=dict(white_error=.05,DeltaPhi_error=.05),new_material_discovery=False))
