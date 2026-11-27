"""Directional response bounds and feasible exterior changes for the registered task."""
from quick import CFG, block, load_case, Response, HERE, write_json, elapsed, now, np
from budget_proof import VerifiedState, VerifiedNativeState, VerifiedWave, up, down, fro, norm_lower, rd
from scipy.optimize import minimize
from scipy.ndimage import label


def real(a):
    return np.concatenate((a.real.ravel(),a.imag.ravel()))


def legal(contrast, beta, exterior, kind):
    ext=contrast[exterior]
    if kind=="max":
        return up(float(np.max(np.abs(ext))))<=beta if beta else not np.any(ext)
    if kind=="l2":
        return fro(ext)<=down(beta*np.sqrt(len(ext))) if beta else not np.any(ext)
    if kind=="l1":
        return up(np.sum(np.abs(ext))*(1.+rd(len(ext))))<=down(beta*len(ext)) if beta else not np.any(ext)
    raise ValueError(kind)


def support(q, lower, upper, beta, kind):
    """Conservative box/ball support, plus an optimizer's support subgradient."""
    cap=np.where(q>=0,upper,-lower)
    n=len(cap)
    if not beta:
        return 0.,np.zeros_like(q),0.,0.
    if kind=="max":
        a=np.minimum(cap,beta)
        e=np.copysign(a,q)
        return up(np.dot(np.abs(q),a)*(1.+rd(n))),e,fro(a),min(beta,float(np.max(cap)))
    if kind=="l1":
        order=np.argsort(-np.abs(q),kind="stable")
        v=cap[order]
        total=up(beta*n)
        cum=np.cumsum(v)
        cum_prior=np.r_[0.,cum[:-1]]
        # Round upper support allocations up. Numerically exact subgradient
        # proposals use the same greedy order with ordinary binary64 arithmetic.
        safe_prior=cum_prior/(1.+rd(n))
        eabs=np.minimum(v,np.maximum(0.,(total-safe_prior)*(1.+rd(3))))
        e=np.zeros_like(q)
        e[order]=np.copysign(np.minimum(v,np.maximum(0.,beta*n-cum_prior)),q[order])
        d=max(float(np.max(cap)),max(float(np.max(-lower)),float(np.max(upper))))
        max_bound=min(d,total)
        delta2=min(fro(np.maximum(-lower,upper)),up(np.sqrt(max_bound*total)))
        return up(np.dot(eabs,np.abs(q[order]))*(1.+rd(n))),e,delta2,max_bound
    d=up(beta*np.sqrt(n))
    # Intersection support bounded by relaxing either the ball or its box.
    direction=np.where(cap>0,q,0.)
    ndir=fro(direction)
    if not np.any(direction):
        e=np.zeros_like(q)
    else:
        e=direction*(d/ndir)
    box_upper=up(np.dot(np.abs(q),cap)*(1.+rd(n)))
    if up(d*ndir)<=box_upper:
        val=up(d*ndir)
    else:
        val=box_upper
        e=np.copysign(cap,q)
    return val,e,min(fro(np.maximum(-lower,upper)),d),min(d,max(float(np.max(-lower)),float(np.max(upper))))


def optimize_direction(d, j, lower, upper, beta, kind):
    if norm_lower(d)==0.:
        return np.zeros_like(d),np.zeros_like(lower)
    initial=d/np.linalg.norm(d)
    scale=np.linalg.norm(d)
    # Primal least-squares solution provides a dual direction. For L2/L1 the
    # empirical direction is used and still rechecked against a safe support.
    if kind!="max" or beta==0:
        return initial,np.zeros_like(lower)
    lo=np.maximum(lower,-beta)
    hi=np.minimum(upper,beta)
    def objective(e):
        v=(d+j@e)/scale
        return .5*np.dot(v,v), (j.T@v)/scale
    result=minimize(objective,np.zeros(len(lo)),jac=True,method="L-BFGS-B",bounds=list(zip(lo,hi)),
                    options=dict(maxiter=35,maxls=15,ftol=1e-10,gtol=1e-9))
    e=np.clip(result.x,lo,hi)
    v=d+j@e
    direction=v/np.linalg.norm(v) if norm_lower(v)>0 else initial
    return direction,e


def directional(d, lam, j, err_j, u, e, err_u, err_e, k, beta, kind, lo, hi, err_d):
    """Enclose the response of every box- and budget-feasible exterior edit."""
    scale=up(max(1.,fro(lam)))
    lam=lam/(scale*(1.+rd(4)))
    if fro(lam)>1.:
        raise ArithmeticError("Response direction exceeds unit norm")
    q=-(j.T@lam)
    qerr=up(err_j+rd(j.shape[0])*fro(j)*fro(lam))
    s,_,e2,dmax=support(q,lo,hi,beta,kind)
    # Bound errors for all edits in the original (possibly asymmetric) box.
    if kind=="max" and beta:
        e2=fro(np.minimum(np.maximum(-lo,hi),beta))
        dmax=min(beta,float(np.max(np.maximum(-lo,hi))))
    eu=min(up(dmax*(fro(u)+err_u)),up(e2*(peak(u,axis=0)+err_u)))
    ee=min(up(dmax*(fro(e)+err_e)),up(e2*(peak(e,axis=1)+err_e)))
    gain=up(k*dmax) if dmax else 0.
    # Use Lambda to weight one side before bounding the nonlinear term.
    matrix=lam[:d.size//2].reshape(u.shape[0],e.shape[1])+1j*lam[d.size//2:].reshape(u.shape[0],e.shape[1])
    lu=matrix.conj().T@u
    le=e@matrix.conj().T
    lue=up(err_u+rd(u.shape[0])*fro(matrix)*fro(u))
    lee=up(err_e+rd(e.shape[1])*fro(matrix)*fro(e))
    elu=min(up(dmax*(fro(lu)+lue)),up(e2*(peak(lu,axis=0)+lue)))
    ele=min(up(dmax*(fro(le)+lee)),up(e2*(peak(le,axis=1)+lee)))
    remainder=up(min(elu*ee,eu*ele)*k/(1.-gain)) if gain<1 else None
    uniform=up(eu*ee*k/(1.-gain)) if gain<1 else None
    p=down(float(np.dot(lam,d))-rd(len(d))*fro(lam)*fro(d)-err_d)
    val=max(0.,down(p-s-e2*qerr-remainder)) if remainder is not None else 0.
    # Norm-subtraction control discards direction and asymmetric box geometry.
    unifj=min(fro(j)*e2,
              dmax*np.sum(np.linalg.norm(j,axis=0))*(1.+rd(j.shape[1]+j.shape[0])))
    old=max(0.,down(norm_lower(d)-err_d-unifj-e2*err_j-uniform)) if uniform is not None else 0.
    return dict(lower=val, norm_control_lower=old, projection_lower=p, safe_support=s,
                derivative_term_error=up(e2*qerr), nonlinear=remainder,
                norm_control_nonlinear=uniform, k_times_max_edit=gain),lam


def peak(a,axis):
    return up(float(np.max(np.linalg.norm(a,axis=axis)))*(1.+rd(2*a.shape[axis])))


def optimize_witness(physical, anchor, target, layout, e, ext, lo, hi, beta, kind, scale, maxiter=12):
    if beta==0 or kind!="max":
        return anchor.copy()
    low=np.maximum(lo,-beta)
    high=np.minimum(hi,beta)
    initial=np.clip(e,low,high)
    def objective(z):
        x=anchor.copy()
        x[ext]+=z
        st=physical.full(x,derivatives=True)
        d=block(st["p"],layout)-target
        dual=CFG["layout"][layout]
        u,en=st["u"][dual["rx"]][:,ext],st["e"][ext][:,dual["tx"]]
        grad=np.einsum("ij,ik,kj->k",d.conj(),u,en,optimize=True).real
        return .5*(np.linalg.norm(d)/scale)**2,grad/scale**2
    opt=minimize(objective,initial,jac=True,method="L-BFGS-B",bounds=list(zip(low,high)),
                 options=dict(maxiter=maxiter,maxls=14,ftol=1e-10,gtol=1e-8))
    x=anchor.copy()
    # Stay strictly within the represented beta despite endpoint addition rounding.
    xmax=np.clip(opt.x,low,high)
    xmax=np.clip(xmax, -max(0.,beta-5e-15),max(0.,beta-5e-15))
    x[ext]+=xmax
    np.clip(x,0.,4.5,out=x)
    return x


def audited_state(wave, physical, rec, exc, x, derivatives=False):
    source_cost=physical.model.cost
    active=int(np.count_nonzero(x))
    if len(x)==32**2 or active<=1100:
        state=VerifiedState(wave,rec,exc,x,full_k=len(x)==32**2)
        if active:
            source_cost["verified_active_lu"]+=1
            source_cost["verified_inverse_rhs"]+=active
            source_cost["verified_point_rhs"]+=exc.shape[1]
        return state
    state=VerifiedNativeState(wave,rec,exc,x,physical.full(x,derivatives=derivatives))
    details=state.log["inverse"]
    source_cost["verified_native_inverse_calls"]+=1
    source_cost["verified_native_fourier_columns"]+=details.get("rank",0)
    source_cost["verified_native_state_rhs"]+=exc.shape[1]
    if derivatives:
        source_cost["verified_native_adjoint_rhs"]+=len(rec)
    return state


def elapsed_cost(physical):
    return {key: value for key,value in physical.model.cost.items()}


def explicit_difference(a,b):
    return real(a-b),up(rd(3)*(fro(a)+fro(b)))


def strict_norm_difference(a,b,error):
    d,err=explicit_difference(a,b)
    return up(fro(d)+err+error)


def task_amplitudes(key):
    # First pass was frozen before either layout was evaluated. This grid
    # refines only the previously found threshold neighborhoods.
    if key=="joint32":
        return [0.03,0.035,0.04,0.045,0.05,0.06,0.125,0.20,0.21,0.22,0.23,0.24,0.25,0.375,0.5,1,2]
    if key=="air32":
        return [0.25,0.50,0.54,0.55,0.555,0.56,0.57,0.625,1,2]
    if key=="air96":
        return [0.25,0.50,0.55,0.56,0.57,0.58,0.59,0.625,1,2]
    if key=="native96":
        return [0.04,0.05,0.125]
    raise ValueError(key)


def evaluate(key, focused=False):
    start=now()
    side,base,exc,rec_all,roi,c,r,obs_id,observed=load_case(key)
    rec=rec_all[:16]
    wave=VerifiedWave(dict(base,S=rec,E0=exc))
    physical=Response(side,base,exc,rec)
    ref=audited_state(wave,physical,rec,exc,c)
    ext=~roi
    rows=[]
    anchors=[]
    arrays=dict(reference=c,roi=roi,p_reference=ref.p,j_reference=ref.j,radius=r)
    for tag,y in observed.items():
        arrays["observed_"+tag]=y
    budgets=CFG["budgets"] if side==32 else CFG["budgets"][:-2]
    amplitudes=task_amplitudes(key)
    if focused:
        amplitudes=[0.252,0.255,0.26,0.27,0.28]
        budgets=[0.,0.0001,0.0002,0.0003,0.0004,0.0005,0.0007,0.001,0.0012,0.0015,0.002,0.0025,0.003]
    for tid,t in enumerate(amplitudes):
        work=now()
        anchor=c+t*roi
        state=audited_state(wave,physical,rec,exc,anchor,derivatives=True)
        k=state.k_upper(ext,ext,dense=side==32)
        prefix=f"t{tid:02d}_"
        arrays[prefix+"anchor"]=anchor
        arrays[prefix+"j_anchor"]=state.j
        arrays[prefix+"p_anchor"]=state.p
        arrays[prefix+"u_anchor"]=state.u
        arrays[prefix+"e_anchor"]=state.e
        lo=-anchor[ext]
        hi=np.nextafter(4.5-anchor[ext],np.inf)
        hi[anchor[ext]==4.5]=0.
        estimates={}
        anchors.append(dict(t=t,k_exterior_upper=k,**state.log))
        print(key,t,"anchor",state.log,"k_EE",k,flush=True)
        for layout in ("A","B"):
            jac,ej,u,en=state.jacobian(CFG["layout"][layout])
            jp=jac[:,ext]
            uu,ee=u[:,ext],en[ext]
            for contract in (("prospective",) if focused else ("prospective","actual")):
                if contract=="actual" and layout not in observed:
                    continue
                rix=2*r if contract=="prospective" else r
                target=block(ref.p,layout) if contract=="prospective" else observed[layout]
                ref_error=ref.response_error if contract=="prospective" else 0.
                d,err=explicit_difference(block(state.p,layout),target)
                ed=up(err+state.response_error+ref_error)
                point_hi=strict_norm_difference(block(state.p,layout),target,up(state.response_error+ref_error))
                for kind in (("max",) if focused else ("max","l1","l2")):
                    for bi,beta in enumerate(budgets):
                        lam,z=optimize_direction(d,jp,lo,hi,beta,kind)
                        report,chosen=directional(d,lam,jp,ej,uu,ee,state.u_error,state.e_error,k,beta,kind,lo,hi,ed)
                        raw_lam=d/np.linalg.norm(d)
                        raw,raw_chosen=directional(d,raw_lam,jp,ej,uu,ee,state.u_error,state.e_error,k,beta,kind,lo,hi,ed)
                        if raw["lower"]>report["lower"]:
                            report,chosen=raw,raw_chosen
                        receipt=dict(t=t,beta=beta,kind=kind,layout=layout,contract=contract,
                              anchor_point_upper_over_threshold=up(point_hi/rix),
                              lower_over_threshold=max(0.,down(report["lower"]/rix)),
                              empirical_direction_lower_over_threshold=max(0.,down(raw["lower"]/rix)),
                              norm_control_lower_over_threshold=max(0.,down(report["norm_control_lower"]/rix)),
                              k_exterior_upper=k,t_index=tid,**report)
                        ident=f"{prefix}{layout}_{contract}_{kind}_b{bi:02d}"
                        receipt["dual_key"]=ident+"_lambda"
                        arrays[receipt["dual_key"]]=chosen
                        estimates[layout,contract,kind,bi]=dict(receipt=receipt,e=z,target=target,
                                                               forward_error=ref_error,threshold=rix)
        # Attempt nonlinear primal proposals only around an unhideability
        # threshold and use the zero edit immediately if already feasible.
        for (layout,contract,kind,bi),item in estimates.items():
            receipt=item["receipt"]
            beta=receipt["beta"]
            point_hi=receipt["anchor_point_upper_over_threshold"]
            denied=receipt["lower_over_threshold"]>1
            hideable=point_hi<1
            hi_now=point_hi
            witness=anchor
            if (not denied and not hideable and kind=="max" and beta>0 and
                 ((key=="joint32" and t<=0.5 and beta<=0.03) or
                  (key in ("air32","air96") and t<=0.625 and beta<=0.003))):
                candidate=optimize_witness(physical,anchor,item["target"],layout,item["e"],ext,
                                           lo,hi,beta,kind,item["threshold"],maxiter=8)
                wave_candidate=physical.full(candidate)
                candidate_ratio=np.linalg.norm(block(wave_candidate["p"],layout)-item["target"])/item["threshold"]
                receipt["nonlinear_candidate_numeric_over_threshold"]=float(candidate_ratio)
                if candidate_ratio<0.9999 and legal(candidate-anchor,beta,ext,kind):
                    candidate_state=audited_state(wave,physical,rec,exc,candidate)
                    hi_now=up(strict_norm_difference(block(candidate_state.p,layout),item["target"],
                         up(candidate_state.response_error+item["forward_error"]))/item["threshold"])
                    if hi_now<1:
                        hideable=True
                        witness=candidate
                        suffix=f"{prefix}{layout}_{contract}_{kind}_b{bi:02d}"
                        receipt["x_key"]=suffix+"_x"
                        arrays[receipt["x_key"]]=candidate
                        arrays[suffix+"_j"]=candidate_state.j
                        arrays[suffix+"_p"]=candidate_state.p
            if denied and hideable:
                raise ArithmeticError("Primal/dual disagreement")
            delta=(witness-anchor)[ext]
            area=(2/side)**2
            receipt.update(witness_upper_over_threshold=hi_now,
                           exterior_norm=dict(max=float(np.max(np.abs(delta))),
                                              weighted_l1=float(area*np.sum(np.abs(delta))),
                                              weighted_l2=float(np.sqrt(area)*np.linalg.norm(delta))),
                           status="hideable" if hideable else "unhideable" if denied else "unresolved")
            rows.append(receipt)
        print(key,t,"s",elapsed(work),"counts",elapsed_cost(physical),flush=True)
    stem=("focus_" if focused else "")+key
    with (HERE/f"raw_{stem}.npz").open("wb") as f:
        np.savez_compressed(f,**arrays)
    report=dict(key=key,side=side,cells=int(roi.sum()),area_m2=float(roi.sum())*(2/side)**2,
                components=int(label(roi.reshape(side,side))[1]),r=r,obs_id=obs_id,reference_audit=ref.log,
                actual_reference_upper={tag:strict_norm_difference(block(ref.p,tag),y,ref.response_error)/r
                                        for tag,y in observed.items()},
                anchor_audits=anchors,records=rows,**elapsed(start),cost=elapsed_cost(physical))
    write_json(HERE/f"budget_{stem}.json",report)
    positives=[(z["t"],z["beta"],z["kind"],z["contract"]) for z in rows
               if z["layout"]=="B" and z["status"]=="unhideable" and
               any(w["layout"]=="A" and w["status"]=="hideable" and
                   tuple(w[f] for f in ("t","beta","kind","contract"))==
                   tuple(z[f] for f in ("t","beta","kind","contract")) for w in rows)]
    print(key,"B positive",positives,"cost",report["wall_s"],flush=True)


if __name__=="__main__":
    import sys
    args=sys.argv[1:]
    if args==["--focused"]:
        evaluate("joint32",focused=True)
    else:
        for key in args or list(CFG["reference_tasks"]):
            evaluate(key)
