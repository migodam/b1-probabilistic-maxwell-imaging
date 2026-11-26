"""Outward verified support-function separation of two genuine material boxes."""
import argparse
import json
from pathlib import Path
from time import perf_counter, process_time
from tempfile import TemporaryDirectory
import subprocess

import numpy as np
from scipy.linalg import cholesky, eigvalsh, lu_factor, lu_solve
from scipy.optimize import lsq_linear
from flint import acb, arb, ctx

from model import R
from native_verify import ExplicitDots, exact_bare_representative, gamma, tiny
from verify_witness import ball, frobenius, outward

HERE = Path(__file__).resolve().parent
WIDTHS = (".0001", ".001", ".005", ".01", ".02", ".05")
ETA = arb(".005")


def adiff(computed, first, second, dots):
    return dots.norm(computed) + gamma(6)*(dots.norm(first)+dots.norm(second)) + 8*computed.size*tiny


def spectral_from_gram(matrix, dots, tight=False):
    gram, error = dots.product(matrix, matrix.conj().T)
    n = len(gram)
    # Real and imaginary absolute row sums avoid assumptions about libm hypot.
    sums = np.sum(abs(gram.real), axis=1)+np.sum(abs(gram.imag), axis=1)
    row_upper = (arb(float(np.max(sums)))+8*n*tiny)/(1-gamma(4*n+8))
    if tight:
        gram_h=0.5*(gram+gram.conj().T)
        herm_error=gamma(6)*dots.norm(gram)+8*gram.size*tiny
        spectral_hint=float(eigvalsh(gram_h,subset_by_index=(n-1,n-1),check_finite=False)[0])
        shift=max(0.,spectral_hint)+max(1e-9,abs(spectral_hint)*1e-7)
        diag=shift*np.eye(n)
        L=cholesky(diag-gram_h,lower=True,check_finite=False)
        product,p_error=dots.product(L,L.conj().T)
        target=diag-gram_h
        remainder=adiff(target-product,target,product,dots)
        shadow=gamma(6)*(dots.norm(diag)+dots.norm(gram_h))+8*gram.size*tiny
        row_upper=min(row_upper,arb(shift)+herm_error+remainder+shadow+p_error)
    return min(dots.norm(matrix), (row_upper+error).sqrt().upper())


def row_format(matrix):
    return np.array([complex(v) for v in matrix.entries()]).reshape(matrix.nrows(), matrix.ncols())


def physical_arrays(data, dots, force_dense=True):
    """One enclosure covers released binary64 dressed and physical port formulas."""
    origin = perf_counter()
    n = len(data["E0"])
    inverse = ball(data["M_nominal"]).solve(ball(np.eye(len(data["M_nominal"]))))
    e, p = ball(data["E0_bare"]), ball(data["S_bare"])
    sqrt_r = arb(R).sqrt()
    # Use exact 2 sqrt(R) and sqrt(R) in both arrays.
    all_B = 2*acb(sqrt_r)*(e*inverse)
    all_S = -acb(0,sqrt_r)*(inverse*p)
    B, S = row_format(all_B), row_format(all_S)
    e_b = frobenius(all_B-ball(B)).upper()
    e_s = frobenius(all_S-ball(S)).upper()
    k, e_k = exact_bare_representative(data["spectrum"])
    y, x = np.indices(tuple(data["grid_shape"]))
    py, px = data["pad"]
    y, x = y.ravel(), x.ravel()
    lag = ((y[:,None]-y[None,:]) % py)*px + ((x[:,None]-x[None,:]) % px)
    G_bare = k[lag]
    w_exact = acb(0,1)*(e*inverse)
    w = row_format(w_exact)
    w_error = frobenius(w_exact-ball(w)).upper()
    series, series_error = dots.product(w, data["S_bare"])
    G_approx = G_bare+series
    G_error = (n*e_k+w_error*dots.norm(data["S_bare"])+series_error
               +gamma(5)*(dots.norm(G_bare)+dots.norm(series))+4*n*tiny)
    if force_dense and "G" in data:
        G = data["G"]
        G_error += adiff(G-G_approx, G, G_approx, dots)
    else:
        G = G_approx
        release_lift=data["low_rank_left"]
        G_error+=(adiff(release_lift-w,release_lift,w,dots)+w_error)*dots.norm(data["S_bare"])
    q = data["E0"].shape[1]
    trainB = data["E0"]
    trainS = data["S"]
    from_allB = B[:,data["source_index"]]/np.sqrt(q)
    from_allS = S[data["receive_index"]]
    train_errorB = (e_b/arb(q).sqrt()+adiff(trainB-from_allB, trainB, from_allB, dots)
                    +gamma(10)*dots.norm(from_allB))
    train_errorS = e_s+adiff(trainS-from_allS, trainS, from_allS, dots)
    matrix_inv = R*inverse
    matrix_inv_mid = row_format(matrix_inv)
    matrix_inv_error = frobenius(matrix_inv-ball(matrix_inv_mid)).upper()
    return dict(G=G, errorG=G_error, inverse_empty=matrix_inv_mid, error_inverse_empty=matrix_inv_error,
                train_index=data["receive_index"], source_index=data["source_index"], B=trainB, errorB=train_errorB, S=trainS, errorS=train_errorS,
                fullB=B, error_fullB=e_b, fullS=S, error_fullS=e_s,
                seconds_operator_enclosure=perf_counter()-origin)


def witness_matrices(op, m, dots):
    """Original full-box inverse, derivative, and resolvent data at one center."""
    origin = perf_counter(), process_time()
    G, B, S = [op[k] for k in ("G","B","S")]
    eG, eB, eS = [op[k] for k in ("errorG","errorB","errorS")]
    n = len(m)
    mnorm = arb(float(max(m)))
    A = np.eye(n) - m[:,None]*G
    eA = mnorm*eG+gamma(6)*(arb(n).sqrt()+mnorm*dots.norm(G))+8*n*tiny
    lu = lu_factor(np.asfortranarray(A), overwrite_a=True, check_finite=False)
    inverse = lu_solve(lu, np.eye(n,dtype=complex), check_finite=False)
    _, iu = inverse.shape
    product, product_error = dots.product(A, inverse)
    defect = adiff(np.eye(n)-product, np.eye(n), product, dots)+product_error+eA*dots.norm(inverse)
    if not defect < 1:
        raise ArithmeticError(f"state right inverse defect exceeds 1: {defect}")
    inv2 = spectral_from_gram(inverse,dots,tight=True)
    inv_exact2 = inv2/(1-defect)
    e_inverse = inv_exact2*defect
    field_resolvent, resolvent_round = dots.product(G,inverse)
    resolvent_err = dots.norm(G)*e_inverse+eG*inv_exact2+resolvent_round
    tnorm_fro = dots.norm(field_resolvent)+resolvent_err
    tnorm = min(tnorm_fro, spectral_from_gram(field_resolvent,dots,tight=True)+resolvent_err)

    def center(b, eb, s, es):
        rhs = m[:,None]*b
        current = lu_solve(lu, rhs, check_finite=False)
        predicted_rhs, r_current = dots.product(A, current)
        r = (mnorm*eb+eA*dots.norm(current)+r_current
             +adiff(predicted_rhs-rhs,predicted_rhs,rhs,dots)+gamma(4)*mnorm*dots.norm(b))
        e_current = inv_exact2*r
        gj, r_field = dots.product(G,current)
        field = b+gj
        ef = eb+eG*(dots.norm(current)+e_current)+dots.norm(G)*e_current+r_field
        ef += gamma(5)*(dots.norm(b)+dots.norm(gj))+4*field.size*tiny
        port, e_port_product = dots.product(s,current)
        e_port = es*(dots.norm(current)+e_current)+dots.norm(s)*e_current+e_port_product
        sensor, e_left = dots.product(s,inverse)
        e_sensor = es*inv_exact2+dots.norm(s)*e_inverse+e_left
        return dict(F=port, eF=e_port, E=field, eE=ef, L=sensor, eL=e_sensor)
    rect = center(B,eB,S,eS)
    full = center(op["fullB"],op["error_fullB"],op["fullS"],op["error_fullS"])
    full_left_bound = spectral_from_gram(full["L"],dots,tight=True)+full["eL"]
    jac = (rect["L"][:,None,:]*rect["E"].T[None,:,:]).reshape(rect["F"].size,n)
    raw = dict(defect=outward(defect)[1], inverse_norm_upper=outward(inv_exact2)[1],
               resolvent_norm_upper=outward(tnorm)[1], initial_Frobenius_resolvent_upper=outward(tnorm_fro)[1],
               port_error=outward(rect["eF"])[1], full_transfer_error=outward(full["eF"])[1],
               seconds_wall=perf_counter()-origin[0],seconds_CPU=process_time()-origin[1], LU=1,
               solved_inverse_columns=iu, rectangular_rhs=B.shape[1], full_terminal_rhs=op["fullB"].shape[1])
    return dict(m=m, rect=rect, full=full, full_left_bound=full_left_bound,
                jac=jac, T=field_resolvent, eT=resolvent_err, Tnorm=tnorm, raw=raw,
                left_gram=dots.product(field_resolvent,field_resolvent.conj().T),
                right_gram=dots.product(field_resolvent.conj().T,field_resolvent))


def endpoints(m,w):
    return np.maximum(-float(w),-m),np.minimum(float(w),4.5-m)


def exact_endpoints(m,w):
    return [(max(-w,-arb(float(x))), min(w,arb(4.5)-arb(float(x)))) for x in m]


def optimize_direction(delta,a,b,w):
    start=perf_counter()
    lo,hi=zip(*(endpoints(c["m"],w) for c in (a,b)),strict=True)
    linear=np.concatenate((a["jac"],-b["jac"]),axis=1)
    flat=delta.ravel()
    real=np.vstack((linear.real,linear.imag))
    d=np.r_[flat.real,flat.imag]
    fit=lsq_linear(real,-d,bounds=(np.concatenate(lo),np.concatenate(hi)),
                   method="trf",tol=1e-10,lsq_solver="lsmr",max_iter=100)
    u=(flat+linear@fit.x).reshape(delta.shape)
    if np.linalg.norm(u) <= 1e-16:
        u=delta
    return u/max(np.linalg.norm(u),1e-300),dict(status=fit.status,iterations=fit.nit,numeric_affine_gap=np.linalg.norm(u),
                                              seconds_wall=perf_counter()-start)


def exact_projection(a,delta):
    return sum((acb(float(x.real),-float(x.imag))*acb(float(y.real),float(y.imag))).real
               for x,y in zip(a.ravel(),delta.ravel(),strict=True))


def upper_real_cube(matrix,error,dots):
    n=len(matrix)
    emax=float(eigvalsh(matrix,subset_by_index=(n-1,n-1),check_finite=False)[0])
    shift=max(0.,emax)+max(1e-9,abs(emax)*1e-7)
    diag=shift*np.eye(n)
    target=diag-matrix
    factor=cholesky(target,lower=True,check_finite=False)
    reconstructed,roundoff=dots.product(factor,factor.T)
    alpha=arb(shift)+adiff(target-reconstructed,target,reconstructed,dots)
    alpha+=gamma(6)*(dots.norm(diag)+dots.norm(matrix))+roundoff+error
    absolute_sum=arb(float(np.sum(abs(matrix))))
    entrywise=(absolute_sum+8*matrix.size*tiny)/(1-gamma(3*matrix.size+8))+n*error
    return min(n*alpha,entrywise)


def quadratic_enclosure(a,v,eV,E,eE,sign,dots):
    T=a["T"]
    product,p_round=dots.product(v.conj(),E.T)
    ep=eV*(dots.norm(E)+eE)+dots.norm(v)*eE+p_round
    raw=sign*np.real(T*product)
    sym=0.5*(raw+raw.T)
    eq=(a["eT"]*dots.norm(product)+(dots.norm(T)+a["eT"])*ep
        +gamma(12)*dots.norm(T)*dots.norm(product)+8*sym.size*tiny)
    return upper_real_cube(sym,eq,dots)


def restricted_operator_enclosure(a, X, eX, left, dots):
    """Sup norm of T*diag(h)X or T diag(h)X, with |h_i| <= 1."""
    gram, e_gram=a["left_gram" if left else "right_gram"]
    X_cross, e_cross=dots.product(X.conj(),X.T)
    Q=np.real(gram*X_cross)
    sym=0.5*(Q+Q.T)
    error=(e_gram*dots.norm(X_cross)+(dots.norm(gram)+e_gram)*e_cross
           +gamma(12)*dots.norm(gram)*dots.norm(X_cross)+8*sym.size*tiny)
    last=(a["eT"]*dots.norm(X)+a["Tnorm"]*eX)
    projected=upper_real_cube(sym,error,dots).sqrt()+last
    unconditional=a["Tnorm"]*(dots.norm(X)+eX)
    return min(projected,unconditional)


def direction_bound(a,u,w,sign,dots):
    """Bound signed directional motion over every point in this material box."""
    v, vp=dots.product(a["rect"]["L"].conj().T,u)
    eV=a["rect"]["eL"]*dots.norm(u)+vp
    eE=a["rect"]["eE"]
    E=a["rect"]["E"]
    d=np.real(np.sum(v.conj()*E,axis=1))*sign
    g_round=gamma(8*v.shape[1]+8)*dots.norm(v)*dots.norm(E)
    g_error=eV*(dots.norm(E)+eE)+dots.norm(v)*eE+g_round+8*v.size*tiny
    linear = arb(0)
    # The sign is fixed for this hypothesis. Clip exactly against [0,4.5].
    for slope,(lo,hi) in zip(d,exact_endpoints(a["m"],w),strict=True):
        linear+=arb(float(slope))*(hi if slope>=0 else lo)
    linear+=w*g_error
    denominator=1-w*a["Tnorm"]
    if not denominator>0:
        return dict(linear=linear,remainder=None,denominator=denominator)
    # Remainder uses local directional row V, rather than full port residual.
    Vnorm=dots.norm(v)+eV
    Enorm=dots.norm(E)+eE
    basic = w*w*Vnorm*a["Tnorm"]*Enorm/denominator
    quadratic = w*w*quadratic_enclosure(a,v,eV,E,eE,sign,dots)
    cubic_basic=(w*a["Tnorm"])*basic
    left_norm=restricted_operator_enclosure(a,v,eV,True,dots)
    right_norm=restricted_operator_enclosure(a,E,eE,False,dots)
    cubic=min(cubic_basic,w**3*left_norm*right_norm/denominator)
    return dict(linear=linear,remainder=min(basic,quadratic+cubic),denominator=denominator,
                quadratic=quadratic,cubic=cubic,initial_remainder=basic,initial_cubic=cubic_basic)


def calibration_bounds(c,op,w,q,eta,dots):
    scale=arb(q).sqrt()
    denominator=1-w*c["Tnorm"]
    if not denominator>0:
        return None
    full=c["full"]
    a_inverse=op["inverse_empty"]
    e_a=op["error_inverse_empty"]
    b_inverse=a_inverse-0.5*full["F"]
    e_b=e_a+full["eF"]/2+gamma(5)*(dots.norm(a_inverse)+dots.norm(full["F"])/2)
    a=spectral_from_gram(a_inverse,dots,tight=True)+e_a
    b_center=spectral_from_gram(b_inverse,dots,tight=True)+e_b
    baseline=dots.norm(full["F"])+full["eF"]
    movement=w*c["full_left_bound"]*(dots.norm(full["E"])+full["eE"])/denominator
    # K is the nominal port-current transfer, even if a bare impedance
    # representation is singular. Adding D gives (I+KD)^-1 K.
    b=b_center+movement/2
    if not (eta*a<1 and eta*b<1):
        return None
    a_true,b_true=a/(1-eta*a),b/(1-eta*b)
    full_true_upper=(1+eta*a_true)*(1+eta*b_true)*(baseline+movement)
    prior_coefficient=eta*(a+b)+eta**2*a*b
    measured_equivalent=prior_coefficient*full_true_upper/scale
    uniform_direct=2*eta*(a*a_true+b*b_true)
    return dict(augmentation=min(uniform_direct,measured_equivalent),full_true=full_true_upper,
                movement=movement,a=a,b=b,a_true=a_true,b_true=b_true,uniform_direct=uniform_direct)


def upper_nuclear(matrix,dots):
    u,diag,vh=np.linalg.svd(matrix,full_matrices=False)
    weighted=u*diag
    approx,recon_error=dots.product(weighted,vh)
    total=sum(dots.norm(weighted[:,i])*dots.norm(vh[i]) for i in range(len(diag)))
    n=arb(len(diag)).sqrt()
    return total+n*(recon_error+adiff(matrix-approx,matrix,approx,dots))


def projected_calibration(c,op,u,cal,eta,dots):
    full=c["full"]["F"]
    ef=c["full"]["eF"]
    a0=op["inverse_empty"]
    ea=op["error_inverse_empty"]
    b0=a0-0.5*full
    eb=ea+ef/2+gamma(5)*(dots.norm(a0)+dots.norm(full)/2)
    # Use unscaled C, then divide by sqrt(number of independent shots).
    direction=np.zeros(full.shape,dtype=complex)
    direction[np.ix_(op["train_index"],op["source_index"])]=u
    c_adj=direction.conj().T
    C_norm=dots.norm(c_adj)
    fv,r1=dots.product(full,c_adj)
    q1,r2=dots.product(fv,b0)
    ba,r3=dots.product(a0,c_adj)
    q2,r4=dots.product(ba,full)
    h=q1+q2
    q1_err=(r2+(ef*C_norm+r1)*(dots.norm(b0)+eb)+dots.norm(fv)*eb)
    q2_err=(r4+(ea*C_norm+r3)*(dots.norm(full)+ef)+dots.norm(ba)*ef)
    h_err=q1_err+q2_err+gamma(5)*(dots.norm(q1)+dots.norm(q2))
    nuclear=upper_nuclear(h,dots)+arb(len(h)).sqrt()*h_err
    norm_center=dots.norm(full)+ef
    move=cal["movement"]
    a,b,a_true,b_true=[cal[k] for k in ("a","b","a_true","b_true")]
    remainder=eta*((move/2+eta*b*b_true+eta*a*a_true)*norm_center+(a_true+b_true)*move)
    remainder+=eta**2*a_true*b_true*(norm_center+move)
    return (eta*nuclear+C_norm*remainder)/arb(u.shape[1]).sqrt()


def optional_measurement(archive,eta):
    """An independent synthetic full-port measurement with a nonzero error."""
    actual=archive["F_true"]
    actual_delta_norm = np.linalg.norm(archive["normalized_additive_error"],ord=2)
    # The fixture perturbation is generated from a spectral-normalized network shape.
    if abs(actual_delta_norm-float(eta))>1e-12:
        raise ValueError("Use the paired normalized perturbation for this fixture.")
    rng=np.random.default_rng(96117731)
    n=rng.normal(size=actual.shape)+1j*rng.normal(size=actual.shape)
    perturb=(0.01*np.linalg.norm(actual)/np.linalg.norm(n))*n
    measured=actual+perturb
    observed=frobenius(ball(measured)).upper()
    actual_abs= frobenius(ball(measured)-ball(actual)).upper()
    declared=(arb(".02")*frobenius(ball(actual))).upper()
    if not actual_abs<declared:
        raise ArithmeticError("auxiliary reading exceeds its stated noise")
    bound=observed+declared
    return dict(provenance="actual pressure-test transfer with synthetic independent all-port readout noise",
                full_transfer_observed_norm_upper=outward(observed)[1],
                declared_full_readout_error=outward(declared)[1],
                certified_fixture_readout_error_upper=outward(actual_abs)[1],
                full_transfer_norm_bound=outward(bound)[1],
                difference_omitted_error_at_eta_005=outward((2*eta+eta**2)*declared/arb(8).sqrt())[1]),measured


def run(op,materials,r,dirs,dots,observed=None):
    rec=[]
    a,b=[witness_matrices(op,m,dots) for m in materials]
    delta=a["rect"]["F"]-b["rect"]["F"]
    eA=adiff(delta, a["rect"]["F"], -b["rect"]["F"], dots) - dots.norm(delta)
    latter=(a["rect"]["eF"]+b["rect"]["eF"]+eA)
    Q=a["rect"]["F"].shape[1]
    for width in WIDTHS:
        w=arb(width)
        for kind in dirs:
            if kind=="center_difference":
                u=delta/np.linalg.norm(delta)
                optimizer=None
            else:
                u,optimizer=optimize_direction(delta,a,b,w)
            umag=dots.norm(u)
            # F_a-F_b loses separation when a decreases or b increases.
            x=[direction_bound(c,u,w,sign,dots) for c,sign in [(a,-1),(b,1)]]
            projected=exact_projection(u,delta)-umag*latter
            row=dict(width=width, direction=kind, projected_center_gap=outward(projected)[0],
                     projected_center_gap_over_two_r=outward(projected/(2*r*umag))[0],
                     linear_support_upper=[outward(v["linear"])[1] for v in x],
                     remainder_support_upper=[outward(v["remainder"])[1] if v["remainder"] is not None else None for v in x],
                     resolvent_contraction_margin_lower=[outward(v["denominator"])[0] for v in x],
                     quadratic_support_upper=[outward(v["quadratic"])[1] if "quadratic" in v else None for v in x],
                     cubic_support_upper=[outward(v["cubic"])[1] if "cubic" in v else None for v in x],
                     initial_cubic_upper=[outward(v["initial_cubic"])[1] if "initial_cubic" in v else None for v in x],
                     initial_remainder_upper=[outward(v["initial_remainder"])[1] if "initial_remainder" in v else None for v in x],
                     optimizer=optimizer)
            if any(v["remainder"] is None for v in x):
                row["status"]="contraction_unresolved"
            else:
                class_gap=projected-sum(v["linear"]+v["remainder"] for v in x)
                nominal=class_gap-2*r*umag
                masks=4*ETA
                passive=nominal-2*masks*umag
                prior=[calibration_bounds(c,op,w,Q,ETA,dots) for c in (a,b)]
                if any(v is None for v in prior):
                    row.update(status="calibrated_resolvent_unresolved",nominal_margin_lower=outward(nominal)[0],
                               gap_over_two_r=outward(class_gap/(2*r*umag))[0])
                    rec.append(row)
                    continue
                calibrated=nominal-umag*sum(v["augmentation"] for v in prior)
                directional_penalties=[projected_calibration(c,op,u,v,ETA,dots)
                                       for c,v in zip((a,b),prior,strict=True)]
                direction_calibrated=nominal-sum(min(p,umag*v["augmentation"]) for p,v in zip(directional_penalties,prior,strict=True))
                original_gap=projected-sum(v["linear"]+v["initial_remainder"] for v in x)
                row.update(status="separated_nominal" if nominal>0 else "nominal_unresolved",
                           gap_over_two_r=outward(class_gap/(2*r*umag))[0],
                           nominal_margin_lower=outward(nominal)[0],
                           masked_absolute_error_margin_lower=outward(passive)[0],
                           separated_masked_absolute_error=bool(passive>0),
                           prior_model_calibrated_margin_lower=outward(calibrated)[0],
                           separated_prior_model_calibrated=bool(calibrated>0),
                           original_margin_lower=outward(original_gap-2*r*umag)[0],
                           original_prior_calibrated_margin_lower=outward(original_gap-2*r*umag-umag*sum(v["augmentation"] for v in prior))[0],
                           directional_calibration_margins_lower=outward(direction_calibrated)[0],
                           separated_directional_calibration=bool(direction_calibrated>0),
                           worst_projected_network_calibration_upper=[outward(v)[1] for v in directional_penalties],
                           class_calibration_augmentation_upper=[outward(v["augmentation"])[1] for v in prior],
                           class_full_transfer_upper=[outward(v["full_true"])[1] for v in prior],
                           shared_full_probe_observed_plus_error_cap=max(0.,outward(max(arb(0), (class_gap/(2*umag)-r)*arb(Q).sqrt()/(2*ETA+ETA**2)))[0]))
            rec.append(row)
            print(width,kind,row.get("gap_over_two_r"),row.get("nominal_margin_lower"),
                  row.get("prior_model_calibrated_margin_lower"),row.get("directional_calibration_margins_lower"),flush=True)
    first_second=[]
    if observed is not None:
        for c in (a,b):
            Y=c["rect"]["F"]
            Y0=np.asarray(observed)
            numeric=frobenius(ball(Y)-ball(Y0))
            lower=(numeric-c["rect"]["eF"])/r
            upper=(numeric+c["rect"]["eF"])/r
            first_second.append(dict(observation_residual_over_r_interval=[max(0.,outward(lower)[0]),outward(upper)[1]],
                                     endpoint_in_observation=bool(upper<1),endpoint_excluded_by_observation=bool(lower>1)))
    return rec,[a["raw"],b["raw"]],first_second


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pair",type=Path,default=HERE/"results/mechanism_witness.npz")
    p.add_argument("--keys",nargs=2,default=("with_region","without_region"))
    p.add_argument("--radius",type=float)
    p.add_argument("--observed-key",default="Y")
    p.add_argument("--name",default="mechanism")
    p.add_argument("--count",type=int,choices=(6,8,12,16),default=8)
    p.add_argument("--band",type=int,choices=(400,460),default=400)
    p.add_argument("--power",type=int,choices=(16,48,144),default=48)
    p.add_argument("--directions",nargs="+",choices=["center_difference","affine_closest_points"],
                   default=["center_difference","affine_closest_points"])
    args=p.parse_args()
    ctx.prec,ctx.threads=128,1
    started=perf_counter(),process_time()
    operator_id=("n32_f400_a08" if (args.band,args.count,args.power)==(400,8,48)
                 else f"n32_f{args.band}_a{args.count:02}_c{args.power:03}")
    with np.load(HERE/"public/operators"/f"{operator_id}.npz") as d:
        data=dict(d)
    with np.load(args.pair) as w:
        material=[w[k].ravel() for k in args.keys]
        r=arb(float(args.radius)) if args.radius is not None else arb(3)*args.count*arb(float(w["reference_sigma"]))/2
        observed=w[args.observed_key] if args.observed_key in w else None
    assert all(np.min(m)>=0 and np.max(m)<=4.5 and len(m)==1024 for m in material)
    with TemporaryDirectory() as tmp:
        library=Path(tmp)/"certify.dylib"
        subprocess.run(["clang","-O3","-fno-fast-math","-ffp-contract=off","-dynamiclib",
                        "-o",str(library),str(HERE/"certify_dot.c")],check=True)
        dots=ExplicitDots(library,4)
        op=physical_arrays(data,dots)
        rows,costs,original_witnesses=run(op,material,r,args.directions,dots,observed)
    pth=HERE/"results"/f"{args.name}_class_separation.json"
    arb6=json.loads((HERE/"results/mechanism_witness_arb.json").read_text())["conditions"]["400_6"]
    with np.load(HERE/"public/calibration_example_e005.npz") as ex:
        finite_probe,measured=optional_measurement(ex,ETA)
        if args.name=="mechanism":
            np.savez_compressed(HERE/"public"/"calibration_example_e005_finite_probe.npz",
                                F_measured=measured, radius=finite_probe["declared_full_readout_error"],
                                receive_index=ex["receive_index"],source_index=ex["source_index"],
                                reference_ohms=ex["reference_ohms"],
                                error_eta=.005)
    pth.write_text(json.dumps(dict(pair_file=str(args.pair.relative_to(HERE.parent)) if args.pair.is_relative_to(HERE.parent) else str(args.pair),
         center_keys=args.keys, mode=args.name, operator_id=operator_id, radius=outward(r), calibration_eta=outward(ETA),
         actual_measurement_endpoint_witnesses=original_witnesses,measurement_key=args.observed_key,
         six_ports_center_pair_q6=arb6 if args.name=="mechanism" else None,
         operator_enclosure=outward(op["errorG"])[1], operator_enclosure_wall=op["seconds_operator_enclosure"],
         finite_full_probe=finite_probe, center_cost=costs, rows=rows, seconds_wall=perf_counter()-started[0],
         seconds_CPU=process_time()-started[1]),indent=2))


if __name__=="__main__":
    main()
