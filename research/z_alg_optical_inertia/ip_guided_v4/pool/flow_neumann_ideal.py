"""Upper-bound conventional and direction-projected Neumann candidate scores."""

import json

import numpy as np
from flint import arb, ctx
from scipy.linalg import lu_factor, lu_solve
from threadpoolctl import threadpool_limits

from box_certificates import ball, gershgorin, lower, norm, star, upper
from flow_audit import OUTPUT
from flow_verify import SignedNative, column_scale, row_scale
from pool_common import clock, elapsed, write_json


def ideal_upper(verifier, g, s, b, y, lo, hi, witness):
    begin = clock()
    kappa = arb(witness["strict_kappa_lower"])
    if kappa >= 1:
        return dict(classical_inapplicable=True, strict_contraction_lower=lower(kappa),
                    reason="Verified lower bound on operator norm is at least one.", **elapsed(begin))
    center = .5*(lo+hi)
    halfwidth = .5*(hi-lo)
    matrix = -center[:, None]*g
    matrix.flat[::len(center)+1] += 1.
    center_factor = lu_factor(matrix, check_finite=False)
    rj = lu_solve(center_factor, center[:, None]*b, check_finite=False)
    rt = lu_solve(center_factor, s.conj().T, trans=2, check_finite=False)
    transfer = rt.conj().T*halfwidth[None, :]
    vm, _, _ = np.linalg.svd(transfer, full_matrices=True)
    endpoints = [(arb(float(l)), arb(float(h))) for l,h in zip(lo, hi, strict=True)]
    m, d = [(l+h)/2 for l,h in endpoints], [(h-l)/2 for l,h in endpoints]
    _, _, _, _, inv, _ = verifier.inverse(m)
    j, t, ye = ball(rj), ball(rt), ball(y)
    ec = verifier.b+verifier.g*j
    defect_j = (inv*norm(j-row_scale(ec,m)).upper()).upper()
    defect_t = (inv*norm(t-star(verifier.g)*row_scale(t,m)-star(verifier.s)).upper()
                *max(v.upper() for v in d)).upper()
    err_c = (verifier.s_norm*defect_j).upper()
    ec_lo = max(arb(0.), (norm(ec).lower()-verifier.g_norm*defect_j).lower())
    residual = ye-verifier.s*j
    optic = column_scale(star(t),d)
    v = ball(vm[:, :1])
    gain = (norm(star(v)*optic).lower()/norm(v).upper()-defect_t).lower()
    standard = norm(residual).upper()+err_c-max(arb(0),gain)*ec_lo/(1-kappa)
    projected = []
    for keep in range(1, s.shape[0]+1):
        w = ball(vm[:,-keep:].conj().T)
        metric = w*star(w)
        w_hi = gershgorin(metric, upper=True).sqrt().upper()
        w_lo = max((metric[i,i].real.lower() for i in range(keep))).sqrt().lower()
        woptic = w*optic
        row_norms = []
        for i in range(keep):
            row = sum((abs(woptic[i,j])**2 for j in range(woptic.ncols())), arb(0))
            row_norms.append(row.sqrt().lower())
        s_lower = (max(row_norms) - w_hi*defect_t).lower()
        numerator = (norm(w*residual).upper()+w_hi*err_c-max(arb(0),s_lower)*ec_lo/(1-kappa)).upper()
        edge = (numerator/(w_hi if numerator < 0 else w_lo)).upper()
        projected.append(upper(edge))
    return dict(classical_inapplicable=False,
                ideal_standard_screen_upper=upper(standard), ideal_projected_screens_upper=projected,
                ideal_projected_max_upper=max(projected), inverse_norm_upper=upper(inv),
                center_residual_additive_uncertainty=upper(err_c), center_transfer_uncertainty=upper(defect_t),
                cost_counts=dict(lu=1, forward_rhs=b.shape[1], receiver_adjoint_rhs=s.shape[0]), **elapsed(begin))


def main():
    ctx.prec, ctx.threads = 96, 1
    total = clock()
    neumann = json.loads((OUTPUT/"native_neumann.json").read_text())["cases"]
    with np.load(OUTPUT/"native_geometry.npz", allow_pickle=False) as geometry:
        g, s, b = geometry["green"], geometry["sensor"][geometry["used_rx"]], geometry["incident"]
    checked = SignedNative(g,s,b)
    rows = {}
    for name, receipt in neumann.items():
        with np.load(OUTPUT/f"{name}.npz", allow_pickle=False) as data:
            lo, hi, y = data["lows"], data["highs"], data["train"]
        proven = {}
        for label, kappa in receipt["branches"].items():
            proof = ideal_upper(checked,g,s,b,y,lo[int(label)],hi[int(label)],kappa)
            proven[label] = proof
        rows[name] = dict(radius=receipt["radius"], rows=proven)
        audit = [v for v in proven.values() if not v["classical_inapplicable"]]
        print(name, "max classic and directional ideal scores divided by r",
              max(v["ideal_standard_screen_upper"] for v in audit)/receipt["radius"],
              max(v["ideal_projected_max_upper"] for v in audit)/receipt["radius"],flush=True)
        write_json(OUTPUT/"native_neumann_ideal.json", dict(cases=rows, setup=checked.setup,total=elapsed(total),
            scope="Ideal standard and the same fixed 16 trailing SVD-projected Neumann sphere scores. "
                  "True contraction when available can only enlarge the sphere compared with the checked lower bound on kappa."))


if __name__ == "__main__":
    with threadpool_limits(limits=1):
        main()
