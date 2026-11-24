"""Center-aware matrix SOS. Each row lives in a two-dimensional affine span."""
import argparse
from time import perf_counter, process_time

import cvxpy as cp
from flint import arb, ctx
import numpy as np
from scipy.linalg import eigh
from threadpoolctl import threadpool_limits

from replay import CASES, HERE, ball, gershgorin, herm, high, load_json, low, save_json, star, times
from outward import f_norm


def directions_at_capacity(a,b,frame_limit):
    # Fixed directions are proposed by two aggregate quadratic forms followed
    # by exact per-pixel scalar quadratic minimization for each direction.
    _,rotation = eigh(herm((a-b/2).conj().T@(a-b/2)))
    n = a.shape[1]
    candidates = np.column_stack((np.eye(n),rotation))
    av,bv = a@candidates,b@candidates
    xx = np.clip(np.real(av.conj()*bv)/np.maximum(abs(bv)**2,np.finfo(float).tiny),0,1)
    per_row = np.abs(av-xx*bv)**2
    max_scores = per_row.sum(axis=0)*frame_limit
    return candidates,max_scores


def add_poly_row(weight, basis):
    # An affine row r(x)=(1,x) R, with x in [0,1]. Embed a 2-by-2
    # symmetric polynomial lower matrix through the exact affine coefficients.
    lo = cp.Variable((2,2), symmetric=True)
    sos = cp.Variable((4,4), symmetric=True)
    rem = cp.Variable((2,2), symmetric=True)
    q0 = np.array([[1.,0.],[0.,0.]])
    q1 = np.array([[0.,1.],[1.,0.]])
    q2 = np.array([[0.,0.],[0.,1.]])
    constraints = [sos >> 0, rem >> 0, weight*q0-lo == sos[:2,:2],
                   weight*q1 == sos[:2,2:]+sos[2:,:2]+rem,
                   weight*q2 == sos[2:,2:]-rem]
    return cp.conj(basis.T)@lo@basis, constraints, (lo,sos,rem)


def outward_psd_lower(a):
    ab = ball(a)
    _,v = eigh(herm(a))
    vb = ball(v)
    q = star(vb)*ab*vb
    gram = star(vb)*vb
    g0 = gershgorin(gram,upper=False)
    g1 = gershgorin(gram,upper=True)
    lo = gershgorin(q,upper=False)
    if not g0>0:
        raise ArithmeticError("PSD proposal basis is singular.")
    return (lo/(g0 if lo < 0 else g1)).lower()


def certify_poly(ws,lo,sos,rem):
    q0 = np.array([[1.,0.],[0.,0.]])
    q1 = np.array([[0.,1.],[1.,0.]])
    q2 = np.array([[0.,0.],[0.,1.]])
    w = arb(float(ws))
    low_s = min(arb(0),outward_psd_lower(sos))
    low_r = min(arb(0),outward_psd_lower(rem))
    e0 = w*ball(q0)-ball(lo)-ball(sos[:2,:2])
    e1 = w*ball(q1)-ball(sos[:2,2:])-ball(sos[2:,:2])-ball(rem)
    e2 = w*ball(q2)-ball(sos[2:,2:])+ball(rem)
    error = (f_norm(e0)+f_norm(e1)+f_norm(e2)-2*low_s-low_r/4).upper()
    return error


def frame_cap(gamma, gram_shift):
    # Congruence to a represented full eigenbasis keeps the certified
    # scalar cap close to 1/lambda_max rather than a row-sum estimate.
    _,v=eigh(herm(gamma))
    vb=ball(v)
    g0=gershgorin(star(vb)*vb,upper=False)
    if not g0>0:
        raise ArithmeticError("Endpoint Gram eigenbasis is singular.")
    upper = gershgorin(star(vb)*ball(gamma)*vb,upper=True)/g0+arb(gram_shift)
    limit=low(arb(.9999)/upper)
    return limit,high(upper)


def main():
    t,c=perf_counter(),process_time()
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("case",choices=CASES)
    p.add_argument("--object",type=int,required=True)
    p.add_argument("--rows-factor",type=int,choices=(1,2,4),default=4)
    args=p.parse_args()
    dat=load_json(HERE/f"{args.case}.json")
    rec=dat["rows"][args.object-1]
    arrays=np.load(HERE/f"{args.case}.npz",allow_pickle=False)
    planes=[arrays[x] for x in ("center","center_field","hs")]
    j,e,h=(a[args.object-1] for a in planes)
    s=h.shape[1]
    modes_path=HERE/f"{args.case}_conditioned.npz"
    modes=np.load(modes_path,allow_pickle=False) if modes_path.exists() else arrays
    n=modes["negative_vec"]/np.sqrt(-modes["negative_eig"])
    nu=n.shape[1]
    coords=args.rows_factor*nu
    pixels=modes["pixel_index"][:coords]
    ii=np.r_[np.arange(coords),len(modes["pixel_index"])+np.arange(coords)]
    gamma=modes["endpoint_gram"][np.ix_(ii,ii)]
    strict_path=HERE/f"{args.case}_strict.npz"
    strict_meta_path=HERE/f"{args.case}_strict.json"
    strict=None
    if strict_path.exists():
        strict=load_json(strict_meta_path)
        if "endpoint_gram" in strict:
            gamma_strict=np.load(strict_path,allow_pickle=False)["endpoint_gram_midpoint"][np.ix_(ii,ii)]
            wlim,gram_upper=frame_cap(gamma_strict,strict["endpoint_gram"]["upper_scalar_shift"])
        else:
            strict=None
    if strict is None:
        wlim,gram_upper=frame_cap(gamma,0.)
    scale,cap=dat["profile"]["scale"],dat["profile"]["cap"]
    a,b=n[pixels],cap*modes["negative_gn"][pixels]
    candidate_dirs,scores=directions_at_capacity(a,b,wlim)
    best=int(np.argmax(scores))
    vals,v=eigh(h)
    v,vals=v[:,::-1],vals[::-1]
    ell=max(nu+1,int(np.count_nonzero(vals>0)))
    q=v[:,:ell]
    hs=herm(q.conj().T@h@q)
    rj, re = np.sqrt(scale)*j[pixels],np.sqrt(scale)*cap*e[pixels]
    result=dict(case=args.case,id=rec["id"],pixels=pixels.tolist(),frame_max_scalar_weight=wlim,
                endpoint_gram_upper=gram_upper,outward_frame_available=strict is not None,
                single_direction_capacity_candidates=scores.tolist(),max_single_candidate_capacity=float(scores[best]),
                baseline_ratio=rec["baseline_ratio"],optimized=[])
    archive=dict(directions=candidate_dirs,pixels=pixels,candidate_capacity=scores)
    for dim in ([0,1] if scores[best]>1.001 else [0]):
        trial_start=perf_counter(),process_time()
        ww=cp.Variable(coords,nonneg=True)
        qdim=dim+ell
        tvals=cp.Variable(ell)
        weights=ww<=wlim
        if dim:
            direction=candidate_dirs[:,best:best+1]
            da,db=a@direction,-b@direction
        else:
            da,db=np.zeros((coords,0)),np.zeros((coords,0))
        # Rescale sources to improve conditioning. The result reported below
        # divides the source variables and lower matrix back to physical units.
        rescale=10.
        r0=np.column_stack((da,rescale*rj@q))
        r1=np.column_stack((db,-rescale*re@q))
        symbols=[]
        lmi_full=np.zeros((qdim,qdim),dtype=complex)
        constraints=[weights,tvals[:-1]>=tvals[1:]]
        for i in range(coords):
            form, on, symbol=add_poly_row(ww[i],np.stack((r0[i],r1[i])))
            lmi_full = lmi_full+form
            constraints.extend(on)
            symbols.append(symbol)
        correction=cp.bmat([[np.eye(dim), np.zeros((dim,ell))],
                     [np.zeros((ell,dim)), cp.diag(tvals)-rescale**2*hs]]) if dim else (
                     cp.diag(tvals)-rescale**2*hs)
        constraints.append(lmi_full-correction >> 0)
        problem=cp.Problem(cp.Maximize(cp.sum(tvals[nu-dim:])),constraints)
        status=None
        try:
            problem.solve(solver="SCS",max_iters=16000,eps_abs=2e-7,eps_rel=2e-7)
            status=problem.status
        except cp.SolverError as err:
            result["optimized"].append(dict(c=dim,status=str(err),**times(*trial_start)))
        if status is None or ww.value is None:
            continue
        radius=rec["padded_radius"]
        target=max(0.,problem.value/rescale**2)
        ws=np.clip(ww.value,0,wlim)
        lo,sos,rem=(np.stack([entry[i].value for entry in symbols]) for i in range(3))
        entry=dict(c=dim,status=status,objective_squared=target,ratio=float(np.sqrt(target)/radius),
                   source_projection_dimension=ell,positive_weights=int(np.count_nonzero(ws>1.e-8)),
                   negative_center_dimension_optimized=dim,**times(*trial_start))
        if dim==0 and strict is not None:
            proof=(strict.get("outward_centers") or [r["outward_baseline"] for r in dat["rows"]])[args.object-1]
            strict_l=[]
            t0=perf_counter(),process_time()
            # Use the full source space for the final certificate. For c=0
            # define MY directly as minus the aggregate certified source form.
            herror=arb(proof["source_form_saved_to_trial_F_error_upper"])+arb(proof["absolute_gap_defect_correction_upper"])
            jerror=arb(proof["center_error_frobenius_upper"])+arb(proof["center_saved_to_trial_F_error_upper"])
            eerror=arb(strict["projector"]["green_norm_upper"])*arb(proof["center_error_frobenius_upper"])+arb(
                proof["total_center_field_saved_to_trial_F_error_upper"])
            mu=arb(scale).sqrt()
            rb=np.column_stack((j[pixels],e[pixels]))
            hb=ball(h)
            hcenter_err=arb(0)
            for i in range(coords):
                err=certify_poly(ws[i],lo[i],sos[i],rem[i])
                li=ball(lo[i])
                li[0,0]-=err
                li[1,1]-=err
                # Using represented J,E inside a ball with explicitly
                # bounded complete-center discrepancies.
                r=ball(np.stack((rb[i,:s],rb[i,s:])))*mu
                for k in range(s):
                    r[1,k]*=-arb(cap)
                hb+=star(r)*li*r
                rerr=mu*(jerror+arb(cap)*eerror)
                hcenter_err+=f_norm(li)*(2*f_norm(r)*rerr+rerr**2)
                strict_l.append(np.array([[complex(li[k,l]) for l in range(2)] for k in range(2)]))
            y=np.array([[complex(hb[k,l]) for l in range(s)] for k in range(s)])
            val,vec=eigh(herm(y))
            vb=ball(vec)
            shifted=star(vb)*hb*vb
            for i,sig in enumerate(val):
                shifted[i,i]-=arb(float(sig))
            e0=f_norm(shifted)
            vg=star(vb)*vb
            qlo,qhi=gershgorin(vg,upper=False),gershgorin(vg,upper=True)
            if not qlo>0:
                raise ArithmeticError("Final source trial eigenbasis is singular.")
            eigs_lower=[]
            for value in val[:s-nu]:
                dlo=(arb(float(value))-e0).lower()
                elo=(dlo/(qhi if dlo>0 else qlo)-herror-hcenter_err).lower()
                eigs_lower.append(low(elo))
            d2=sum((arb(max(0.,lam)) for lam in eigs_lower),arb(0))
            new_d=low(d2.sqrt()) if d2>0 else 0.
            entry.update(strict_distance_lower=new_d,strict_ratio_lower=low(arb(new_d)/arb(radius)),
                         strict_gain_over_baseline_upper=low(arb(new_d)-arb(
                             rec["outward_baseline"]["distance_upper"])),
                         outward_center_error_in_source_form=high(hcenter_err),
                         outward_check=times(*t0))
            archive["c0_polynomial_lower_certified"]=np.array(strict_l)
        archive.update({f"c{dim}_weight":ws,f"c{dim}_lo":lo,f"c{dim}_sos":sos,f"c{dim}_rem":rem})
        result["optimized"].append(entry)
        print(args.case,dim,entry,flush=True)
    result.update(total=times(t,c))
    out=HERE/f"{args.case}_center_sdp_{args.object}_{args.rows_factor}"
    save_json(out.with_suffix(".json"),result)
    np.savez_compressed(out.with_suffix(".npz"),**archive)


if __name__=="__main__":
    ctx.prec,ctx.threads=112,2
    with threadpool_limits(1):
        main()
