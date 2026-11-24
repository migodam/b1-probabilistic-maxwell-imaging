"""Outward replay of endpoint Grams and a support-cardinality obstruction."""
import argparse
from time import perf_counter, process_time

from flint import acb, acb_mat, arb, ctx
import numpy as np
from threadpoolctl import threadpool_limits

from replay import (CASES, HERE, ball, configure, frobenius_squared, gershgorin,
                    high, load_json, low, mid, outward_sources, save_json, star, times)


def f_norm(a):
    return frobenius_squared(a).sqrt().upper()


def sq(a):
    return a.real*a.real+a.imag*a.imag


def rowscale(a, scale):
    return acb_mat(a.nrows(), a.ncols(), [a[i,j]*scale[i]
                    for i in range(a.nrows()) for j in range(a.ncols())])


def top_sum(values, amount):
    return sum((arb(float(r)) for r in sorted(values, reverse=True)[:amount]), arb(0))


def safe_projector(native, state, vec, eigen):
    start = perf_counter(), process_time()
    vb = ball(vec)
    ng = vb.ncols()
    diagonal = [arb(float(x)) for x in eigen]
    av = native.exact_action(vb)
    vlam = acb_mat(vb.nrows(), ng, [vb[i,j]*diagonal[j]
                    for i in range(vb.nrows()) for j in range(ng)])
    residual = av-vlam
    gplus = arb(state["positive_gap_lower"])
    positive_leak = f_norm(acb_mat(vb.nrows(), ng, [
        residual[i,j]/(gplus-diagonal[j]) for i in range(vb.nrows()) for j in range(ng)]))
    gram = star(vb)*vb
    gmin, gmax = gershgorin(gram, upper=False), gershgorin(gram, upper=True)
    ident = acb_mat(ng,ng,[1 if i==j else 0 for i in range(ng) for j in range(ng)])
    gram_error = f_norm(gram-ident)
    umin = (gmin-positive_leak**2).lower()
    if not umin > 0:
        raise ArithmeticError("Projected negative columns do not provably span the exact negative space.")
    symbol = [arb(1)+value.real for row in native.exact_symbol for value in row]
    max_negative = max(arb(0), -min(x.lower() for x in symbol))
    max_green = max(abs(a).upper() for row in native.g.spec_ball for a in row)
    # U = Pi_- V. On ran Pi_-, ||I-U U*|| <= ||V*V-I||+||Pi_+ V||².
    # Expand (-K)_+ - U diag(-lambda) U*, then compare the U and V outer products.
    correction = (max_negative*(gram_error+positive_leak**2)
                  + f_norm(residual)*gmax.sqrt()
                  + 2*gmax.sqrt()*positive_leak*max(-v for v in diagonal)).upper()
    diagonal_free = sum(symbol, arb(0))/len(symbol)
    mu = arb(native.w.scale)
    weight_upper = []
    for i in range(vb.nrows()):
        norms = [sq(vb[i,j]) for j in range(ng)]
        leverage = ((sum(norms,arb(0)).sqrt()+positive_leak)**2/umin).upper()
        if not leverage < 1:
            raise ArithmeticError(f"Positive leverage at cell {i} is unresolved.")
        kii = diagonal_free + sum((sq(native.g.sensor_ball[j,i])
                                    for j in range(native.g.sensor_ball.nrows())),arb(0))/mu
        kplus = kii + sum((-d*a for d,a in zip(diagonal,norms,strict=True)),arb(0)) + correction
        weight_upper.append(high(kplus/(1-leverage)**2))
    details = dict(ritz_positive_leak_upper=high(positive_leak),
        ritz_residual_F_upper=high(f_norm(residual)), Ritz_gram_min_lower=low(gmin),
        negative_projected_gram_min_lower=low(umin), negative_K_norm_upper=high(max_negative),
        negative_action_approximation_error_upper=high(correction),
        green_norm_upper=high(max_green), max_scalar_weight_upper=max(weight_upper), **times(*start))
    return vb, av, positive_leak, umin, max_negative, max_green, np.array(weight_upper), details


def interval_endpoint_rows(native, pixels):
    g = native.g
    p,n = g.pad,g.side
    slide = [acb.dft([value.conjugate() for value in row],inverse=True) for row in g.spec_ball]
    kernel = [acb.dft([slide[i][j] for i in range(p)],inverse=True) for j in range(p)]
    cols = []
    cap = arb(native.w.cap)
    for i in range(g.n):
        iy,ix = divmod(i,n)
        first = [acb(1 if i == h else 0) for h in pixels]
        second = []
        for h, v in zip(pixels, first,strict=True):
            hy,hx = divmod(int(h),n)
            second.append(v-cap*kernel[(ix-hx)%p][(iy-hy)%p])
        cols.extend(first+second)
    return acb_mat(g.n,2*len(pixels),cols)


def endpoint_proof(native, state, saved, vb, leak, umin, max_negative):
    start = perf_counter(), process_time()
    x = ball(saved["endpoint_trial"])
    r = interval_endpoint_rows(native,saved["pixel_index"])
    kx = native.exact_action(x)
    error = r-kx
    remove = ball(mid(star(vb)*error))
    error_plus = (f_norm(error-vb*remove)+leak*f_norm(remove)).upper()
    negative_trial = ((f_norm(star(vb)*x)+leak*f_norm(x))/umin.sqrt()).upper()
    linear_shift = (2*f_norm(r)*negative_trial).upper()
    form = star(r)*x+star(x)*r-star(x)*kx
    midpoint = mid(form)
    midpoint = (midpoint+midpoint.conj().T)/2
    midpoint_shift = f_norm(form-ball(midpoint))
    e_upper = error_plus**2/arb(state["positive_gap_lower"])
    # Form after true projection differs by a linear term and a nonpositive
    # negative-space energy. The remainder on the true positive space is PSD.
    upper = high(midpoint_shift+linear_shift+e_upper)
    lower = high(midpoint_shift+linear_shift+max_negative*negative_trial**2)
    diff = f_norm(form-ball(saved["endpoint_gram"]))
    info = dict(rows=r.ncols(), upper_scalar_shift=upper, lower_scalar_shift=lower,
                projected_residual_upper=high(error_plus), trial_negative_content_upper=high(negative_trial),
                outward_form_to_cg_defect_upper=high(diff), **times(*start))
    return midpoint, info


def clipped_best(a,b,cap):
    nrm = np.sum(abs(b)**2,axis=1)
    return np.clip(np.real(np.sum(a.conj()*b,axis=1))/np.maximum(nrm, np.finfo(float).tiny),0.,cap)


def obstruction(native, state, saved, vb, av, leak, max_green, weights, source_rows):
    start = perf_counter(), process_time()
    eig, vec, gn = [saved[k] for k in ("negative_eig","negative_vec","negative_gn")]
    nu = len(eig)
    roots = 1/np.sqrt(-eig)
    n = vec*roots
    cap = native.w.cap
    # Saved candidates span either the deepest Ritz axes or the bottom homogeneous-row basis.
    candidates = [np.eye(nu),saved["homogeneous_subspace"]]
    limits = sorted({nu,2*nu,4*nu})
    remaining = set(limits)
    blocks, archive = [], {}
    gv = native.g.exact_green(vb)
    vbkv = star(vb)*av
    max_w = arb(float(weights.max())).sqrt()
    for d in range(nu,0,-1):
        if not remaining:
            break
        scores = []
        for index,basis in enumerate(candidates):
            sub = basis[:,:d]
            a,b = n@sub,gn@sub
            pix = clipped_best(a,b,cap)
            score = weights*np.sum(abs(a-pix[:,None]*b)**2,axis=1)
            scores.append((top_sum(score,max(remaining)),index,sub,pix))
        _, candidate, sub, pix = min(scores,key=lambda x:x[0])
        coordinate = ball(roots[:,None]*sub)
        vx,gx = vb*coordinate,gv*coordinate
        xx = [arb(float(x)) for x in pix]
        ax = vx-rowscale(gx,xx)
        aw = [arb(float(x)) for x in weights]
        power = np.array([high(aw[i]*sum((sq(ax[i,j]) for j in range(d)),arb(0)))
                         for i in range(native.g.n)])
        # In exact negative projection K-form only becomes smaller, because the
        # removed part is positive. Bound drift of the selected residual rows.
        negative_energy = -gershgorin(star(coordinate)*vbkv*coordinate, upper=True)
        proj_drift = (max_w*(1+arb(cap)*max_green)*leak*f_norm(coordinate)).upper()
        margins = {}
        for k in sorted(remaining):
            total = ((top_sum(power,k).sqrt()+proj_drift)**2).upper()
            margin = (negative_energy-total).lower()
            margins[str(k)] = dict(weighted_trace_upper=high(total), margin_lower=low(margin))
        active = [k for k in remaining if margins[str(k)]["margin_lower"]>0]
        if not active:
            continue
        remaining -= set(active)
        archive[f"contraction_pixels_x_{d}"] = pix
        archive[f"contraction_basis_{d}"] = roots[:,None]*sub
        archive[f"contraction_row_power_upper_{d}"] = power
        blocks.append(dict(contraction_dimension=d, max_compensated_dimension=nu-d,
            candidate=candidate, verified_cardinalities=sorted(active),
            negative_energy_lower=low(negative_energy), norm_projection_error_upper=high(proj_drift),
            margins=margins))
    result = dict(contractions=blocks, unresolved_cardinalities=sorted(remaining), **times(*start))

    start = perf_counter(), process_time()
    center, field, hs = [saved[k] for k in ("center","center_field","hs")]
    numeric = np.linalg.eigh(hs)
    mu = arb(native.w.scale)
    source_results = []
    for block in blocks:
        m = block["contraction_dimension"]
        for k in block["verified_cardinalities"]:
            for i,(j,e,h,proof) in enumerate(zip(center,field,hs,source_rows,strict=True)):
                q = numeric[1][i,:,:len(numeric[0][i])-m]
                ev = numeric[0][i,:q.shape[1]]
                jb,eb,qb = ball(j),ball(e),ball(q)
                qp,eq = j@q,e@q
                x = clipped_best(qp,eq,cap)
                exact = rowscale(eb*qb,[arb(float(z)) for z in x])-jb*qb
                row_power = np.array([high(mu*arb(float(w))*sum((sq(exact[t,l])
                                      for l in range(q.shape[1])),arb(0)))
                                      for t,w in enumerate(weights)])
                gram = star(qb)*qb
                qmin,qmax = gershgorin(gram,upper=False),gershgorin(gram,upper=True)
                if not qmin>0:
                    raise ArithmeticError("Source compression is ill-posed.")
                form = star(qb)*ball(h)*qb
                for index, v in enumerate(ev):
                    form[index,index] -= arb(float(v))
                herror = arb(proof["source_form_saved_to_trial_F_error_upper"]) + arb(
                         proof["absolute_gap_defect_correction_upper"])
                shift = (f_norm(form)+herror*qmax).upper()
                base = sum((max(arb(0),(arb(float(v))+shift).upper()) for v in ev),arb(0))
                center_error = arb(proof["center_error_frobenius_upper"])
                jerror = center_error+arb(proof["center_saved_to_trial_F_error_upper"])
                eerror = max_green*center_error+arb(proof["total_center_field_saved_to_trial_F_error_upper"])
                drift = (mu.sqrt()*max_w*qmax.sqrt()*(jerror+arb(cap)*eerror)).upper()
                power = (top_sum(row_power,k).sqrt()+drift)**2
                radius = arb(proof["radius"])
                upper = ((base+power)/qmin).sqrt()
                source_results.append(dict(id=proof["id"], row_budget=k, free_negative_directions=m,
                    compensation_dimension_upper=nu-m, source_base_squared_upper=high(base/qmin),
                    center_gain_squared_upper=high(power/qmin), center_drift_upper=high(drift),
                    distance_upper=high(upper), ratio_upper=high(upper/radius),
                    crossing_impossible_for_theor_2_4=bool(upper<radius)))
                archive[f"source_x_{m}_{i}"] = x
                archive[f"source_power_{m}_{i}"] = row_power
    result.update(signed_center_source_upper=source_results,
                  outward_source_cardinality_cost=times(*start))
    return result, archive


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("case", choices=CASES)
    p.add_argument("--skip-endpoint",action="store_true")
    args=p.parse_args()
    start=perf_counter(),process_time()
    geo,native,data,_=configure(args.case)
    base=HERE/args.case
    previous=load_json(base.with_suffix(".json"))
    state=previous["verified_state"]
    a=np.load(base.with_suffix(".npz"),allow_pickle=False)
    saved={k:a[k] for k in a.files}
    new=HERE/f"{args.case}_conditioned.npz"
    if new.exists():
        with np.load(new,allow_pickle=False) as normalized:
            saved.update({k:normalized[k] for k in normalized.files})
    vec,eig=saved["negative_vec"],saved["negative_eig"]
    print(args.case,"verified projector and pointwise weights",flush=True)
    vb,av,leak,umin,kmin,gmax,wcap,pc=safe_projector(native,state,vec,eig)
    archive={"frame_weight_upper":wcap}
    result=dict(case=args.case,verified_state=state,projector=pc)
    out=HERE/f"{args.case}_strict"
    if not args.skip_endpoint:
        print(args.case,"outward endpoint Gram",flush=True)
        gm,gp=endpoint_proof(native,state,saved,vb,leak,umin,kmin)
        archive["endpoint_gram_midpoint"]=gm
        result["endpoint_gram"]=gp
        save_json(out.with_suffix(".json"),result)
    rows=[dict(**r["outward_baseline"],id=r["id"],radius=r["padded_radius"]) for r in previous["rows"]]
    initial=HERE/f"{args.case}_strict_initial_nonorthogonal_full.json"
    if "center_saved_to_trial_F_error_upper" not in rows[0] and initial.exists():
        cached=load_json(initial)
        rows=cached["outward_centers"]
        result["outward_centers"]=rows
        result["outward_center_fit_reused_from"]=initial.name
        result["outward_center_fit_cost_in_prior_replay"]=cached["outward_center_fit_cost"]
    if "center_saved_to_trial_F_error_upper" not in rows[0]:
        print(args.case,"outward complete center fit",flush=True)
        t=perf_counter(),process_time()
        rr=outward_sources(native,state,data,saved["source_response"],saved["center"],
                           saved["center_field"],saved["hs"])
        rows=[dict(**r,id=m["id"],radius=m["padded_radius"]) for (m,_),r in zip(data,rr,strict=True)]
        result["outward_centers"]=rows
        result["outward_center_fit_cost"]=times(*t)
    print(args.case,"outward sparse exclusion-of-certificate",flush=True)
    comp,a_extra=obstruction(native,state,saved,vb,av,leak,gmax,wcap,rows)
    archive.update(a_extra)
    result.update(**comp,total=times(*start),operations=dict(geo.cost))
    save_json(out.with_suffix(".json"),result)
    np.savez_compressed(out.with_suffix(".npz"),**archive)
    peak=max((x["ratio_upper"] for x in comp["signed_center_source_upper"]),default=None)
    print(args.case,"done",result["total"],"largest theorem ratio upper",peak,flush=True)


if __name__=="__main__":
    ctx.prec,ctx.threads=112,2
    with threadpool_limits(1):
        main()
