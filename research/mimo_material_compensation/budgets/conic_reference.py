"""Restricted two-endpoint SOCP benchmark with separately checked directions."""
from bound import (CFG, HERE, Response, audited_state, block, directional, down, elapsed, explicit_difference,
                   fro, load_case, now, np, rd, real, support, up, write_json)
from budget_proof import VerifiedWave
from scipy.optimize import minimize


def scaled_residual(w,j,x,radius):
    y=(w+j@x)/radius
    return y,j/radius


def main(upper_mean=0.0415):
    start=now()
    side,base,exc,rec,roi,c,r,observation_id,observed=load_case("joint32")
    rec=rec[:16]
    wave=VerifiedWave(dict(base,S=rec,E0=exc))
    physical=Response(side,base,exc,rec)
    # Restricted benchmark: 12 predeclared evenly spaced exterior indices.
    # Both endpoints can differ on those same coordinates independently.
    ext_ids=np.flatnonzero(~roi)
    picks=ext_ids[np.linspace(0,len(ext_ids)-1,12,dtype=int)]
    act=np.zeros(len(c),bool)
    act[picks]=True
    beta=0.005
    a,y=(c,c+upper_mean*roi),observed["A"]
    proofs=[audited_state(wave,physical,rec,exc,x,derivatives=True) for x in a]
    ap=[block(w.p,"A") for w in proofs]
    bp=[block(w.p,"B") for w in proofs]
    old=[explicit_difference(p,y) for p in ap]
    b_delta,dround=explicit_difference(bp[1],bp[0])
    b_err=up(sum(s.response_error for s in proofs)+dround)
    jacs=[{k:s.jacobian(CFG["layout"][k]) for k in ("A","B")} for s in proofs]
    lo=[np.maximum(-beta,-x[picks]) for x in a]
    hi=[np.minimum(beta,np.nextafter(4.5-x[picks],np.inf)) for x in a]
    k=[s.k_upper(act,act,dense=True) for s in proofs]
    remainders=[]
    j_errors=[]
    for i,s in enumerate(proofs):
        rem={}
        er={}
        for lay in ("A","B"):
            j,j_error,u,e=jacs[i][lay]
            dummy,val=directional(np.zeros(128),np.zeros(128),j[:,picks],j_error,u[:,picks],e[picks],
                                  s.u_error,s.e_error,k[i],beta,"max",lo[i],hi[i],0.)
            rem[lay]=dummy["norm_control_nonlinear"]
            er[lay]=up(dummy["derivative_term_error"]+rem[lay])
        remainders.append(rem)
        j_errors.append(er)
    aug=[up(r+proofs[i].response_error+old[i][1]+j_errors[i]["A"]) for i in range(2)]
    A=[jacs[i]["A"][0][:,picks] for i in range(2)]
    B=[jacs[i]["B"][0][:,picks] for i in range(2)]
    n=len(picks)
    J=np.column_stack((-B[0],B[1]))
    D=[np.column_stack((A[0],np.zeros_like(A[0]))),np.column_stack((np.zeros_like(A[1]),A[1]))]
    dA=[o[0] for o in old]
    def cone(v,m,w,stepsize,limit):
        rv,scaled=scaled_residual(w,m,v[:2*n],stepsize)
        rn=np.linalg.norm(rv)
        gradient=rv@scaled/max(rn,1e-30)
        return ((v[-1] if limit is None else limit)-rn,
                np.r_[-gradient, 1. if limit is None else 0.])
    constraints=[]
    for mat,d,limit in [(J,b_delta,None)]+list(zip(D,dA,np.array(aug)/r)):
        constraints.append(dict(type="ineq",
            fun=lambda x,m=mat,w=d,h=limit:cone(x,m,w,r,h)[0],
            jac=lambda x,m=mat,w=d,h=limit:cone(x,m,w,r,h)[1]))
    init=np.r_[np.zeros(2*n),np.linalg.norm(b_delta)/r]
    opt=minimize(lambda x:x[-1],init,method="SLSQP",jac=lambda x:np.r_[np.zeros(2*n),1.],
                 constraints=constraints,bounds=list(zip(np.r_[lo[0],lo[1],0.],np.r_[hi[0],hi[1],np.inf])),
                 options=dict(ftol=1e-13,maxiter=500))
    direction=b_delta+J@opt.x[:2*n]
    direction=direction/(up(fro(direction))*(1.+rd(4)))
    av=[d+M@opt.x[:2*n] for d,M in zip(dA,D)]
    directionsA=[p/(up(fro(p))*(1.+rd(4))) for p in av]
    def dual(weights):
        alphas=[w*p for w,p in zip(weights,directionsA)]
        parts=[]
        center_term=float(direction@b_delta)
        err=up(rd(len(direction))*fro(direction)*fro(b_delta))
        for i,alpha in enumerate(alphas):
            oldterm=float(alpha@dA[i])
            center_term+=oldterm
            err=up(err+rd(len(alpha))*fro(alpha)*fro(dA[i]))
            qi=(-1 if i==0 else 1)*(B[i].T@direction) + A[i].T@alpha
            rounding=up(rd(2*len(direction)+8)*(fro(B[i])*fro(direction)+fro(A[i])*fro(alpha)))
            supp,_,_,_=support(-qi,lo[i],hi[i],beta,"max")
            boxdelta=fro(np.maximum(-lo[i],hi[i]))
            parts.append(up(aug[i]*fro(alpha)+supp+boxdelta*rounding))
        linear=max(0.,down(center_term-err*(1.+rd(4))-sum(parts)))
        lower=max(0.,down(linear-b_err-sum(q["B"] for q in j_errors)))
        return linear,lower,parts,alphas
    ad=minimize(lambda w:-dual(w)[0]/r,np.zeros(2),method="Nelder-Mead",bounds=[(0,20),(0,20)],
                options=dict(xatol=1e-7,maxiter=400))
    candidates=[(dual(w),w) for w in (np.zeros(2),ad.x)]
    line,cert,parts,alphas=max(candidates,key=lambda x:x[0][0])[0]
    xlower=dual(np.zeros(2))[0]
    rho=np.linalg.norm((b_delta+J@opt.x[:2*n])/r)
    datum=dict(observation_id=observation_id,r=r,beta=beta,roi_means=[0.,upper_mean],
               scope="restricted local benchmark on 12 exterior pixels at each endpoint, fixed nominal physics",
               optimizer="SLSQP evaluates norm-cone slack and gradients; status is numerical",
               sites=picks.tolist(),old_enlarged_radius_over_r=(np.array(aug)/r).tolist(),
               old_numerical_residual_over_r=[float(np.linalg.norm(v)/r) for v in av],
               epigraph_slack=[float(fun["fun"](opt.x)) for fun in constraints],
               primal_success=bool(opt.success),iterations=int(opt.nit),status=opt.message,
               numerical_linear_primal_over_r=float(rho),linear_dual_over_r=line/r,
               no_old_multiplier_linear_lower_over_r=xlower/r,
               certified_nonlinear_lower_over_r=cert/r,
               endpoint_remainders=remainders,endpoint_jacobian_remainder_errors=j_errors,
               optimized_normalization_upper=fro(direction),dual_costs=parts,**elapsed(start))
    write_json(HERE/"two_endpoint.json",datum)
    with (HERE/"two_endpoint_raw.npz").open("wb") as f:
        np.savez_compressed(f,anchor0=a[0],anchor1=a[1],site_indices=picks,
                            x=opt.x,linear_A=np.stack(A),linear_B=np.stack(B),
                            A_delta0=old[0][0],A_delta1=old[1][0],B_delta=b_delta,
                            lambda_B=direction,alpha_A0=alphas[0],alpha_A1=alphas[1])
    print(datum,flush=True)


if __name__=="__main__":
    import sys
    main(float(sys.argv[1])) if len(sys.argv)>1 else main()
