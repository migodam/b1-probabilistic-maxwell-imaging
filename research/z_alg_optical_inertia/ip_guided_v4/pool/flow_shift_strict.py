"""Strong one-sided proof of native signed L4 proposals with numerical zero inertia."""

import json
import sys

import numpy as np
from flint import arb, ctx
from scipy.linalg import cho_factor, cho_solve, lu_factor, lu_solve
from threadpoolctl import threadpool_limits

from box_certificates import ball, gershgorin, lower, midpoint, norm, spectral_upper, star, upper
from flow_audit import OUTPUT
from flow_verify import SignedNative, row_scale
from pool_common import clock, elapsed, write_arrays, write_json


def verify_l4(verifier, g, s, b, y, lo, hi, tau):
    started = clock()
    saved = {}
    m, r = .5*(lo+hi), .5*(hi-lo)
    p = -m[:, None]*g
    p.flat[::len(m)+1] += 1.
    z = lu_factor(p, check_finite=False)
    state_and_transfer = lu_solve(z, np.column_stack((m[:, None]*b, np.diag(r))), check_finite=False)
    jc, numz = state_and_transfer[:, :y.shape[1]], state_and_transfer[:, y.shape[1]:]
    numeric_s_z = s@numz
    numeric_time = elapsed(started)
    begin = clock()
    ep = [(arb(float(a)), arb(float(b))) for a,b in zip(lo, hi, strict=True)]
    ma, ra = [(a+b)/2 for a,b in ep], [(b-a)/2 for a,b in ep]
    _, _, _, inv_defect, inv, targeted = verifier.inverse(ma)
    ze = ball(numz)
    gz = verifier.g*ze
    defect = ze-row_scale(gz, ma)
    for i, v in enumerate(ra):
        defect[i, i] -= v
    transport_defect = norm(defect).upper()
    # G and P^{-1} define a complex symmetric dressed Green operator.
    if not np.array_equal(g, g.T):
        raise ArithmeticError("Disk cell coefficients are not exactly reciprocal.")
    kappa_hi = (norm(gz).upper()+verifier.g_norm*inv*transport_defect).upper()
    if not kappa_hi < 1:
        raise ArithmeticError("Global positive quadratic gap unresolved.")
    taua = arb(tau)
    gap = (taua*(1-kappa_hi*kappa_hi)/(inv*inv)).lower()
    transport_time = elapsed(begin)

    begin = clock()
    rg = r[:, None]*g
    core = p.conj().T@p-rg.conj().T@rg
    k = s.conj().T@s+tau*.5*(core+core.conj().T)
    lin = s.conj().T@y + tau*(p.conj().T@(m[:, None]*b)+rg.conj().T@(r[:, None]*b))
    x = cho_solve(cho_factor(k, lower=True, check_finite=False), lin, check_finite=False)
    saved["quadratic_approximate_minimizer"] = x
    xa, ya = ball(x), ball(y)
    field = verifier.b+verifier.g*xa
    constitutive = xa-row_scale(field, ma)
    rx = verifier.s*xa-ya
    rw = row_scale(field,ra)
    half_d = star(verifier.s)*rx + taua*(constitutive-star(verifier.g)*
                                          (row_scale(constitutive,ma)+row_scale(rw,ra)))
    quadratic = norm(rx)**2+taua*(norm(constitutive)**2-norm(rw)**2)
    correction = (norm(half_d).upper()**2/gap).upper()
    result = (quadratic-correction).lower()
    strict = max(0., lower(result.sqrt())) if result>0 else 0.
    shift_time = elapsed(begin)

    begin = clock()
    center = ball(jc)
    ec = verifier.b+verifier.g*center
    center_defect = norm(center-row_scale(ec,ma)).upper()
    current_err = (inv*center_defect).upper()
    e_bound = (norm(ec).upper()+verifier.g_norm*current_err).upper()
    yshift = ya-verifier.s*center
    stop = numeric_s_z
    sza = verifier.s*ze
    receiver_uncertainty = (verifier.s_norm*inv*current_err).upper()
    coefficient_uncertainty = (verifier.s_norm*inv*transport_defect).upper()
    standard_optical = (spectral_upper(star(sza),stop.conj().T)+coefficient_uncertainty).upper()
    base_lower = (norm(yshift).lower()-receiver_uncertainty-standard_optical*e_bound/(1-kappa_hi)).lower()
    basis, _, _ = np.linalg.svd(stop, full_matrices=True)
    projected_lower = []
    for keep in range(1, s.shape[0]+1):
        w = ball(basis[:, -keep:].conj().T)
        wnorm = gershgorin(w*star(w),upper=True).sqrt().upper()
        gain = norm(w*sza).upper()+wnorm*coefficient_uncertainty
        error = wnorm*receiver_uncertainty
        edge = (norm(w*yshift).lower()-error-gain*e_bound/(1-kappa_hi))/wnorm
        projected_lower.append(max(0.,lower(edge)))
    # Audited upper on the ideal centered standard Neumann screen. The omitted
    # factor 1/(1-kappa) is >=1 when contraction holds and can only enlarge radius.
    test_vector = ball(np.linalg.svd(stop,full_matrices=False)[0][:, :1])
    gain_lower = (norm(star(test_vector)*sza).lower()/norm(test_vector).upper()-coefficient_uncertainty).lower()
    ec_lower = (norm(ec).lower()-verifier.g_norm*current_err).lower()
    ideal_standard_lower_bound_upper = (norm(yshift).upper()+receiver_uncertainty
                                        -max(arb(0),gain_lower)*max(arb(0),ec_lower)).upper()
    neumann_time = elapsed(begin)
    return dict(negative_inertia=0, positive_gap_lower=lower(gap),
                strict_residual_lower=strict, tau=tau, optical=0.,
                quadratic_completion_correction=upper(correction),
                physical_center_inverse_bound=upper(inv), physical_center_inverse_defect=upper(inv_defect),
                strict_kappa_upper=upper(kappa_hi), transport_equation_defect=upper(transport_defect),
                used_targeted_inverse_defect=targeted,
                neumann_standard_residual_lower=max(0.,lower(base_lower)),
                neumann_ideal_standard_lower_bound_upper=max(0.,upper(ideal_standard_lower_bound_upper)),
                neumann_directional_residual_lower=max(projected_lower),
                costs=dict(numeric_center=numeric_time, exact_center_transport=transport_time,
                           quadratic_completion=shift_time, neumann_curves=neumann_time),
                counts=dict(dense_center_lu=1, center_state_rhs=y.shape[1],
                            center_transfer_rhs=len(m), auxiliary_cholesky=1, auxiliary_rhs=y.shape[1]),
                total=elapsed(started)), saved


def main():
    ctx.prec, ctx.threads = 96, 1
    total = clock()
    proposed = json.loads((OUTPUT/"screen_proposals.json").read_text())
    if proposed["status"] != "complete numeric proposals; exact proof required for any exclusion":
        raise ValueError("Numerical proposals not sealed.")
    with np.load(OUTPUT/"native_geometry.npz") as geometry:
        g, s, b = geometry["green"], geometry["sensor"][geometry["used_rx"]], geometry["incident"]
    verifier = SignedNative(g, s, b)
    receipts, saved = {}, {}
    for name, groups in proposed["cases"].items():
        radius = groups["radius"]
        with np.load(OUTPUT/f"{name}.npz", allow_pickle=False) as test:
            y, lo, hi = test["train"], test["lows"], test["highs"]
        strict = {}
        for box, item in groups["boxes"].items():
            pick = max(item["weights"], key=lambda x: x["numeric_lower"])
            row = dict(status="RETAIN_UNRESOLVED")
            if (pick["numeric_lower"] > radius and pick["optical"]==0.
                    and pick.get("numeric_negative_inertia")==0):
                try:
                    classical, payload = verify_l4(verifier,g,s,b,y,lo[int(box)],hi[int(box)],pick["tau"])
                    row.update(classical)
                    saved.update({f"{name}__box{box}__{key}":value for key,value in payload.items()})
                    row["status"] = "STRICTLY_EXCLUDED" if row["strict_residual_lower"] > radius else "RETAIN_UNRESOLVED"
                except (ArithmeticError, ValueError, ZeroDivisionError) as exc:
                    row["reason"] = str(exc)
            strict[box] = row
            if "strict_residual_lower" in row:
                print(f"{name} box {box}: L4/r {row['strict_residual_lower']/radius:.6f}; "
                      f"Neumann classical ideal upper/r {row['neumann_ideal_standard_lower_bound_upper']/radius:.6f}",
                      flush=True)
            receipts[name] = dict(radius=radius, boxes=strict)
            write_json(OUTPUT/"native_l4.json", dict(cases=receipts, setup=verifier.setup, total=elapsed(total),
                       scope="Strict interval proof over the exact stored original 1024-cell signed finite model."))
    write_arrays(OUTPUT/"native_l4_witnesses.npz",**saved)


if __name__ == "__main__":
    with threadpool_limits(limits=1):
        main()
