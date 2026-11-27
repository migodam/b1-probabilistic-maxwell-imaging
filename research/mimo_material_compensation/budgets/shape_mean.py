"""Single-endpoint support bounds for all legal ROI shapes at specified means."""
from bound import (CFG, HERE, Response, audited_state, block, down, elapsed, explicit_difference, fro, load_case,
                   now, np, optimize_direction, rd, real, support, up, write_json)
from budget_proof import VerifiedWave


def best_task_multiplier(q, roi_lo, roi_hi, t_minus_anchor, inequality):
    """Use the multiplier of the mean constraint to tighten the box support."""
    n=len(q)
    theta_candidates=np.r_[0.,n*q]
    if inequality:
        theta_candidates=theta_candidates[theta_candidates>=0]
    candidates=[]
    for theta in theta_candidates:
        direction=theta/n-q
        a=np.maximum(direction*roi_hi,direction*roi_lo)
        err=rd(n+5)*(np.sum(np.abs(a))+ (abs(theta)/n+np.sum(abs(q))) *
                     np.max(np.maximum(np.abs(roi_lo),np.abs(roi_hi))))
        lower=down(float(theta*t_minus_anchor)-up(np.sum(a)+err))
        candidates.append((lower,float(theta)))
    return max(candidates)


def verify_shape(key):
    start=now()
    side,base,ex,rec,roi,c,r,obs_id,obs=load_case(key)
    if side!=32:
        raise ValueError("Shape-free experiment uses the finite six-cell query")
    rec=rec[:16]
    ext=~roi
    wave=VerifiedWave(dict(base,S=rec,E0=ex))
    physical=Response(side,base,ex,rec)
    ref=audited_state(wave,physical,rec,ex,c)
    rows=[]
    arrays=dict(reference=c,roi=roi,radius=r,p_reference=ref.p)
    tasks=[(t,t,"fixed_mean") for t in (0.20,0.21,0.22,0.23,0.24,0.25,0.50,1.,2.,4.0)]
    tasks += [(t,4.5,"at_least_mean") for t in (0.25,2.,4.,4.25,4.4,4.45,4.48,4.49,4.5)]
    cache={}
    for idx,(t,a,kind) in enumerate(tasks):
        if a not in cache:
            center=c+a*roi
            state=audited_state(wave,physical,rec,ex,center,derivatives=True)
            kparts=[state.k_upper(ra,rb,dense=True) for ra in (roi,ext) for rb in (roi,ext)]
            cache[a]=state,kparts
        state,kparts=cache[a]
        n=int(roi.sum())
        if kind=="fixed_mean":
            roi_lo=np.full(n,max(-a,down((n-1)*(a-4.5))))
            roi_hi=np.full(n,min(4.5-a,up((n-1)*a)))
        else:
            roi_lo=np.full(n,max(-a, down(n*(t-a))))
            roi_hi=np.zeros(n)
        ma_roi=float(np.max(np.maximum(-roi_lo,roi_hi)))
        nprefix=f"task{idx:02d}_"
        arrays[nprefix+"anchor"]=state.x
        arrays[nprefix+"p_anchor"]=state.p
        for lay in ("A","B"):
            j,ej,u,e=state.jacobian(CFG["layout"][lay])
            constant,dround=explicit_difference(block(state.p,lay),block(ref.p,lay))
            err_c=up(dround+state.response_error+ref.response_error)
            lo=-c[ext]
            hi=np.nextafter(4.5-c[ext],np.inf)
            hi[c[ext]==4.5]=0.
            for bi,beta in enumerate((0.,0.0001,0.0003,0.001,0.003,0.01)):
                p,_=optimize_direction(constant,j[:,ext],lo,hi,beta,"max")
                p=p/(up(max(1.,fro(p)))*(1.+rd(4)))
                if fro(p)>1.:
                    raise ArithmeticError("Shape-bound direction exceeds unit norm")
                q=j.T@p
                task_linear,theta=best_task_multiplier(q[roi],roi_lo,roi_hi,t-a,kind=="at_least_mean")
                sup,_,_,_=support(-q[ext],lo,hi,beta,"max")
                de=fro(np.minimum(np.maximum(-lo,hi),beta)) if beta else 0.
                dmax_ext=min(beta,float(np.max(np.maximum(-lo,hi))))
                d2 = up(np.hypot(fro(np.maximum(-roi_lo,roi_hi)),de))
                qerr=up(ej+rd(len(p))*fro(j)*fro(p))
                scales=np.r_[np.sqrt(ma_roi),np.sqrt(dmax_ext)]
                scaled_k=np.array(kparts).reshape(2,2)*np.outer(scales,scales)
                contraction=fro(scaled_k)
                weights=np.where(roi,scales[0],scales[1])
                matrix=p[:64].reshape(8,8)+1j*p[64:].reshape(8,8)
                wu=u*weights
                we=e*weights[:,None]
                e_u=up(max(scales)*(state.u_error+rd(4)*fro(u)))
                e_e=up(max(scales)*(state.e_error+rd(4)*fro(e)))
                elu=up(fro(matrix.conj().T@wu)+e_u+rd(8)*fro(wu)*fro(matrix))
                ele=up(fro(we@matrix.conj().T)+e_e+rd(8)*fro(we)*fro(matrix))
                rem=up(min(elu*(fro(we)+e_e),ele*(fro(wu)+e_u))*contraction/(1.-contraction)) if contraction<1. else None
                projection=down(np.dot(p,constant)-err_c-rd(len(p))*fro(p)*fro(constant))
                lower=(max(0.,down(projection+task_linear-sup-qerr*d2-rem)) if rem is not None else 0.)
                witness=cache[t][0] if t in cache else state
                flat_constant,flat_round=explicit_difference(block(witness.p,lay),block(ref.p,lay))
                flat_upper=up((fro(flat_constant)+flat_round+witness.response_error+ref.response_error)/(2*r))
                result=dict(layout=lay,task=kind,t=t,anchor_mean=a,beta=beta,theta=theta,
                       projection_lower=projection,roi_linear_lower=task_linear,external_support=sup,
                       rounding_derivative_penalty=up(qerr*d2),nonlinear_upper=rem,
                       contraction_upper=contraction,part_k_upper=kparts,
                       lower_over_two_r=max(0.,down(lower/(2*r))),flat_witness_over_two_r=flat_upper,
                       status="unhideable" if lower>2*r else "hideable" if flat_upper is not None and flat_upper<1 else "unresolved")
                tag=f"{nprefix}{lay}_b{bi:02d}_lambda"
                arrays[tag]=p
                result["lambda_key"]=tag
                rows.append(result)
    with (HERE/f"shape_raw_{key}.npz").open("wb") as f:
        np.savez_compressed(f,**arrays)
    write_json(HERE/f"shape_{key}.json",dict(key=key,obs_id=obs_id,radius=r,records=rows,**elapsed(start)))
    print(key,[(a["layout"],a["t"],a["beta"],a["lower_over_two_r"]) for a in rows if a["status"]=="unhideable"],flush=True)


if __name__=="__main__":
    for case in ("joint32","air32"):
        verify_shape(case)
