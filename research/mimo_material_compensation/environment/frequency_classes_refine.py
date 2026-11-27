"""Tighten the 32² finite-box estimate by retaining its signed quadratic term."""
from frequency_classes import (HERE, LAYOUT, AuditedState, VerifiedWave, bounds, check_direction,
                               choose_direction, down, elapsed, explicit_difference, fro, material,
                               now, np, opnorm, rd, released_operator, up, write_json)
from scipy.linalg import eigh


def q_lower(state,u,e,lam,sign):
    m=lam[:64].reshape(8,8)+1j*lam[64:].reshape(8,8)
    g=state.k_matrix
    un,en=fro(u),fro(e)
    ug=u.T@m.conj()
    mtx=ug@e.T
    uv_error=up(state.u_error*(en+state.e_error)+un*state.e_error+
                rd(16)*un*fro(m)*en)
    elementwise_bound=up(float(np.max(np.abs(g)))*(1.+rd(4)))
    quadratic=np.real(g*mtx)*sign
    q_error=up(state.k_matrix_error*(fro(mtx)+uv_error)+elementwise_bound*uv_error+
               rd(6)*elementwise_bound*fro(mtx))
    symmetric=(quadratic+quadratic.T)*.5
    q_error=up(q_error+rd(4)*fro(quadratic))
    values,v=eigh(symmetric,check_finite=False,driver="evd")
    # An approximate orthogonal expansion still yields a one-sided spectral
    # bound when the factorization and orthogonality residuals are included.
    basis=up(np.sqrt(1.+fro(v.T@v-np.eye(len(v)))+rd(len(v))*(fro(v)**2+np.sqrt(len(v)))))
    factor_error=up(fro(symmetric-(v*values)@v.T)+rd(len(v)+4)*
                    (fro(symmetric)+fro(v)**2*fro(values)))
    neg=max(0.,-float(values[0]))
    return up(neg*basis*basis+factor_error+q_error),dict(numerical_min_eigenvalue=float(values[0]),
                 symmetric_quadratic_error_upper=q_error,eigenbasis_factorization_error_upper=factor_error,
                 verified_basis_norm_upper=basis)


def main():
    start=now()
    t,name,r,obs,ix,roi,x,y=material(32)
    _,op=released_operator(32,460)
    w=VerifiedWave(op)
    states=[AuditedState(w,op["S"],op["E0"],m) for m in (x,y)]
    jacs=[s.jacobian(LAYOUT) for s in states]
    d,ed=explicit_difference(states[0].p,states[1].p)
    ed=up(ed+sum(s.response_error for s in states))
    initial=d/np.linalg.norm(d)
    rows=[]
    arrays={}
    for index,beta in enumerate((.015,.02,.025,.03)):
        proposal,_=choose_direction(d,[s[0] for s in jacs],x,y,beta)
        ranks=[]
        for name,l in (("empirical",initial),("optimized",proposal)):
            old,lam,proj,terms=check_direction(d,ed,l,states,jacs,beta)
            refined=[]
            for sign,st,(_,err,u,e),term in zip((1.,-1.),states,jacs,terms):
                q,log=q_lower(st,u,e,lam,sign)
                lo,hi=bounds(st.x,beta)
                delta=fro(np.maximum(-lo,hi))
                second=up(q*delta**2)
                third=up(term["nonlinear"]*term["k_times_max_edit"])
                subtract=up(min(term["nonlinear"],second+third))
                refined.append(dict(signed_quadratic_lower_penalty=second,
                                     third_and_higher_upper=third, refinement_deduction=subtract,
                                     **term,**log))
            cut=up(sum(s["safe_support"]+s["derivative_term_error"]+s["refinement_deduction"] for s in refined))
            lower=max(0.,down(proj-cut-rd(6)*(abs(proj)+cut)))
            ranks.append((lower,dict(direction=name, lower_over_two_r=max(0.,down(lower/(2*r))),
                                     first_order_remainder_only_lower_over_two_r=max(0.,down(old/(2*r))),
                                     projection_lower=proj,endpoint_terms=refined),lam))
        _,best,chosen=max(ranks,key=lambda q:q[0])
        key=f"lambda_{index:02d}"
        arrays[key]=chosen
        best.update(beta=beta,lambda_key=key,status="verified_excluded" if best["lower_over_two_r"]>1 else "unresolved")
        rows.append(best)
        print(beta,best["lower_over_two_r"],best["first_order_remainder_only_lower_over_two_r"],
              [(v["signed_quadratic_lower_penalty"],v["third_and_higher_upper"]) for v in best["endpoint_terms"]],flush=True)
    write_json(HERE/"frequency_classes_32_refine.json",dict(observation_id=obs,rows=rows,
                                                           dense_eigenfactorizations=16,**elapsed(start)))
    with (HERE/"frequency_classes_32_refine.npz").open("wb") as f:
        np.savez_compressed(f,**arrays)


if __name__=="__main__":
    main()
