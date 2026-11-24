"""Outward lower bounds on uncentered negative current in the fixed u=1 box."""
from time import perf_counter, process_time

from flint import acb_mat, arb, ctx
import numpy as np
from scipy.linalg import eigh
from threadpoolctl import threadpool_limits

from outward import f_norm
from replay import HERE, ball, configure, gershgorin, herm, high, load_json, low, mid, save_json, star, times

CASE = "u100_best"


def trace_real(a):
    return sum((a[i, i].real for i in range(a.nrows())), arb(0))


def pick_source_tests(hplus, lneg, radius, a_norm):
    """Choose replayable admissible projections from a small numerical grid."""
    covariance = herm(lneg.conj().T @ lneg)
    h_scale = max(float(np.linalg.norm(hplus, ord=2)), np.finfo(float).tiny)
    c_scale = max(float(np.linalg.norm(covariance, ord=2)), np.finfo(float).tiny)
    gammas = [0.0] + list(h_scale/c_scale*np.logspace(-5, 2, 29))
    tests = []
    for gamma in gammas:
        e, v = eigh(hplus-gamma*covariance, check_finite=False)
        q = v[:, e > 0.]
        if not q.shape[1]:
            continue
        h = float(np.trace(q.conj().T@hplus@q).real)
        b = float(np.linalg.norm(lneg@q))
        d = h-radius**2
        root = d/(b+np.sqrt(b*b+a_norm*d)) if d>0 else 0.
        tests.append((root, float(gamma), q))
    if not tests:
        raise ArithmeticError("No positive source subspace was located in Hplus.")
    best = max(tests,key=lambda x:x[0])
    return [tests[0]] if best[1]==0. else [tests[0],best]


def current_bound(native, baseline, center, h_saved, measurement, vb, vals,
                  leakage, projected_gram_min, negative_action_error, a_upper):
    start=perf_counter(),process_time()
    mu=arb(native.w.scale)
    radius=arb(baseline["radius"])
    j=ball(center)
    h=ball(h_saved)
    vbj=star(vb)*j
    weighted=acb_mat(vbj.nrows(),vbj.ncols(), [vals[i]*vbj[i,jj]
                      for i in range(vbj.nrows()) for jj in range(vbj.ncols())])
    hplus_trial=h-mu*(star(vbj)*weighted)
    delta_h=arb(baseline["source_form_saved_to_trial_F_error_upper"])+arb(
               baseline["absolute_gap_defect_correction_upper"])
    center_error=arb(baseline["center_error_frobenius_upper"])+arb(
                     baseline["center_saved_to_trial_F_error_upper"])
    jnorm=f_norm(j)
    # Physical K_- equals mu times the checked negative action of K/mu.
    # The same norm upper bound covers the true K_- in the center perturbation.
    hplus_delta=(delta_h+mu*negative_action_error*jnorm**2+
                 a_upper*(2*jnorm*center_error+center_error**2)).upper()
    hplus_num=herm(mid(hplus_trial))
    hplus_storage_delta=f_norm(hplus_trial-ball(hplus_num))
    hplus_matrix_error=(hplus_delta+hplus_storage_delta).upper()
    eigen, u=eigh(hplus_num,check_finite=False)
    ub=ball(u)
    trial_diag=star(ub)*ball(hplus_num)*ub
    for i,value in enumerate(eigen):
        trial_diag[i,i]-=arb(float(value))
    eig_drift=f_norm(trial_diag)
    gram=star(ub)*ub
    glo,ghi=gershgorin(gram,upper=False),gershgorin(gram,upper=True)
    if not glo>0:
        raise ArithmeticError("The represented Hplus eigenspace proposal is singular.")
    lower,upper=[],[]
    for value in eigen[::-1]:
        lo=(arb(float(value))-eig_drift).lower()
        hi=(arb(float(value))+eig_drift).upper()
        lower.append(low(lo/(ghi if lo>0 else glo)-hplus_matrix_error))
        upper.append(high(hi/(glo if hi>0 else ghi)+hplus_matrix_error))

    geo=native.g
    # L in the original operator units, formed from the frozen
    # receiving/incident matrices and the frozen binary64 observation.
    _,z=native.exact_scaled_coefficients()
    zb=acb_mat(j.nrows(),j.ncols(),[mu*z[i]*geo.incident_ball[i,k]
                              for i in range(j.nrows()) for k in range(j.ncols())])
    l=star(geo.sensor_ball)*ball(measurement)+zb
    vbl=star(vb)*l
    lnorm=f_norm(l)
    nominal_lneg=mid(vbl)
    proposals=pick_source_tests(hplus_num,nominal_lneg,float(radius),float(a_upper))
    checked=[]
    for j_num,gamma,q in proposals:
        qb=ball(q)
        qgram=star(qb)*qb
        qnorm_hi=arb(max(1.,high(gershgorin(qgram,upper=True))))
        # Z=Q Q*/qnorm_hi is exactly PSD and <= I even for a rounded Q.
        h_lower=(trace_real(star(qb)*hplus_trial*qb)/qnorm_hi-
                 hplus_delta*trace_real(qgram)/qnorm_hi).lower()
        # U=Pi_- V spans the exact negative space. U*U >= g I,
        # ||Pi_- LQ|| <= ||U*LQ||/sqrt(g), and U=V-Pi_+V.
        b_upper=((f_norm(vbl*qb)+leakage*qnorm_hi.sqrt()*lnorm)/
                   (projected_gram_min*qnorm_hi).sqrt()).upper()
        delta=(h_lower-radius**2).lower()
        directed=(delta/(b_upper+(b_upper**2+a_upper*delta).sqrt())).lower() if delta>0 else arb(0)
        checked.append(dict(negative_current_norm_lower=max(0.,low(directed)),
                      h_lower=low(h_lower),b_upper=high(b_upper),a_upper=high(a_upper),
                      h_minus_radius_squared_lower=low(delta),
                      relative_no_negative_distance_lower=low(h_lower.sqrt()/radius) if h_lower>0 else 0.,
                      represented_source_rank=q.shape[1],represented_source_gram_max_upper=high(qnorm_hi),
                      numerical_penalty=gamma,numerical_root_proposal=j_num))
    best=max(checked,key=lambda x:x["negative_current_norm_lower"])
    return dict(id=baseline["id"],radius=float(radius), best=best, checked_projectors=checked,
                hplus_loewner_error_upper=high(hplus_matrix_error),
                source_form_error_upper=high(delta_h),
                center_error_F_upper=high(center_error),
                hplus_source_eigenvalues_lower=lower,hplus_source_eigenvalues_upper=upper,
                actual_negative_current_formula="||Pi_-(K) J||_F",
                **times(*start))


def main():
    start=perf_counter(),process_time()
    geo,native,data,_=configure(CASE)
    record=load_json(HERE/f"{CASE}_strict.json")
    basis_proof=record["projector"]
    modes=np.load(HERE/f"{CASE}_conditioned.npz",allow_pickle=False)
    trial=np.load(HERE/f"{CASE}.npz",allow_pickle=False)
    v=modes["negative_vec"]
    vals=[arb(float(value)) for value in -modes["negative_eig"]]
    vb=ball(v)
    leakage=arb(basis_proof["ritz_positive_leak_upper"])
    gmin=arb(basis_proof["negative_projected_gram_min_lower"])
    error=arb(basis_proof["negative_action_approximation_error_upper"])
    vnorm=gershgorin(star(vb)*vb,upper=True)
    a=(arb(native.w.scale)*(max(vals)*vnorm+error)).upper()
    rows=record.get("outward_centers")
    if rows is None:
        replayed=load_json(HERE/f"{CASE}.json")["rows"]
        rows=[dict(**row["outward_baseline"],id=row["id"],radius=row["padded_radius"])
              for row in replayed]
    result=[]
    for row,center,h,(metadata,y) in zip(rows,trial["center"],trial["hs"],data,strict=True):
        if row["id"]!=metadata["id"]:
            raise ArithmeticError("The saved center and the frozen observation disagree.")
        result.append(current_bound(native,row,center,h,y,vb,vals,leakage,gmin,error,a))
        print(metadata["id"],result[-1]["best"],flush=True)
    save_json(HERE/"required_currents.json",dict(case=CASE,profile=native.w.as_dict(),
       frozen_cohort="six fixed native 32x32, 2% observations",frozen_observation_count=len(result),
       positive_necessary_current_bounds=sum(r["best"]["negative_current_norm_lower"]>0. for r in result),
       bound_scope="Every physical current at cap u=1 fitting the original padded radius; no claim that such current exists.",
       physical_negative_K_norm_upper=high(a),negative_basis_energy_error_scaled_upper=high(error),
       prerequisite_receipts=["u100_best_strict.json","u100_best.json","u100_best.npz","u100_best_conditioned.npz"],
       precision_bits=ctx.prec,rows=result,operations=dict(geo.cost),total=times(*start)))


if __name__=="__main__":
    ctx.prec,ctx.threads=112,2
    with threadpool_limits(1):
        main()
