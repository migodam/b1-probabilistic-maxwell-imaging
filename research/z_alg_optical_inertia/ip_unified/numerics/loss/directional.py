"""Directional absorption multipliers using only receiver Gram blocks."""
import argparse
import json

from flint import arb, ctx
import numpy as np
from scipy.linalg import eigh, svd
from threadpoolctl import threadpool_limits

from finite_states import pure_shift_score
from pure_shift import center_ball
from run import (HERE, adj, cap_coefficients, clock, elapsed, finite_geometry, herm, native_geometry,
                 receiver_setup, write_json)
from verify import FiniteExact, ball, blocks, exact_number, eye, gram_interval, hi, lo, norm, star


def compare_bases(y, rs, cb, fc, common_tau):
    """Hold one orthonormal basis fixed and optimize one price per direction."""
    setup = clock()
    d = y-fc
    rho, u = eigh(herm(rs))
    if rho[0] <= 0:
        raise ArithmeticError("The comparison requires positive receiver Gram directions.")
    projected = adj(u)@d
    _, _, v_right = svd(d, full_matrices=False)
    scaled = adj(projected)@(projected/(rho[:,None]+common_tau))-cb
    dominance = adj(projected)@(projected/rho[:,None])-cb
    right = (
        ("canonical", np.eye(d.shape[1])),
        ("right_singular", adj(v_right)),
        ("best_common_shift", eigh(herm(scaled))[1]),
        ("zero_price_derivative", eigh(herm(dominance))[1]),
    )
    methods = {}
    base_cost = elapsed(setup)
    for method, v in right:
        start = clock()
        contracts = np.diag(adj(v)@cb@v).real
        if min(contracts) <= 0:
            raise ArithmeticError("The comparison requires positive source Gram directions.")
        masses = abs(projected@v)**2
        active = (masses/rho[:,None]).sum(axis=0) > contracts
        upper = np.sqrt((rho[:,None]*masses).sum(axis=0)/contracts)
        lower = np.zeros(len(upper))
        for _ in range(60):
            t = (upper+lower)/2
            deriv = (rho[:,None]*masses/(rho[:,None]+t[None,:])**2).sum(axis=0) - contracts
            lower = np.where(deriv > 0, t, lower)
            upper = np.where(deriv > 0, upper, t)
        prices = np.where(active, (upper+lower)/2, 0.)
        objective = prices*((masses/(rho[:,None]+prices[None,:])).sum(axis=0)-contracts)
        value = np.sqrt(max(0.,objective.sum()))
        methods[method] = (v,prices,dict(distance=float(value),prices=prices.tolist(),
                                         priced_directions=int(np.count_nonzero(prices)), **elapsed(start)))
    return methods,base_cost


def outward_basis(geo,y,v,prices,trial,error):
    """Certify the candidate with the exact unitary polar factor of a stored matrix."""
    start = clock()
    rs, cb, fc = blocks(geo,trial,error)
    D = ball(y)-center_ball(fc,error)
    vb = ball(v)
    e = norm(star(vb)*vb-eye(v.shape[1]))
    if not e < 1:
        raise ArithmeticError("Basis unitary defect was not controlled.")
    polar_error = e/(1+(1-e).sqrt())
    basis_correction = arb(0)
    d_sq_upper = (norm(D)**2).upper()
    cb_upper = norm(cb).upper()
    objective = arb(0)
    for i,t in enumerate(prices):
        if t == 0:
            continue
        column = ball(v[:,[i]])
        dv = D*column
        t = arb(float(t))
        resolvent = rs+eye(geo.p)*t
        quadratic = ((star(dv)*resolvent.solve(dv))[0,0]-(star(column)*cb*column)[0,0]).real
        objective += t*quadratic
        # t(R+tI)^−1 ≤ I, so this bound avoids small denominator loss.
        basis_correction += (d_sq_upper+t*cb_upper)*(1+(1+e).sqrt())*polar_error
    certified = max(arb(0),(objective-basis_correction).lower())
    return dict(distance_lower=max(0.,lo(certified.sqrt())),
                basis_defect_upper=hi(e), objective_lower=lo(certified),
                polar_correction_upper=hi(basis_correction), **elapsed(start))


def compare(study, geo, only):
    start=clock()
    recv = json.loads((HERE/f"{study}.json").read_text())
    out = dict(study=study, profiles=[], mp_precision_bits=ctx.prec)
    with np.load(HERE/f"{study}_inputs.npz",allow_pickle=False) as src:
        for geo_row in recv["receivers"]:
            t,cap = geo_row["tangent"],geo_row["cap"]
            if cap not in only:
                continue
            _,alpha=cap_coefficients(cap,t)
            (rs,cb,fc,_),setup,_=receiver_setup(geo,alpha)
            block=dict(tangent=t,cap=cap,setup=setup,rows=[])
            trial=None
            for row in geo_row["rows"]:
                ident,radius=row["id"],row["rho"]
                if "_low" in ident:
                    continue
                y=src[ident]
                # Small exact-operator Gram setup above is shared across profiles;
                # this reference numerical optimization is charged separately.
                v, u = rs.shape[0], cb.shape[0]
                gram=np.block([[4*cb,adj(-2j*fc)],[-2j*fc,rs]])
                common = pure_shift_score(geo,y,ball(gram))
                bases,cost = compare_bases(y,rs,cb,fc,common["tau"])
                methods={method:stats for method,(_,_,stats) in bases.items()}
                result=dict(id=ident,rho=radius,methods=methods,common_scalar=common,
                            receiver_reference=row["ellipsoid"],shared_per_object=cost)
                if study=="finite_0" and any(stats["distance"] > radius for stats in methods.values()):
                    if trial is None:
                        exact = FiniteExact(geo,.6)
                        trial,error,_,relation = gram_interval(geo,exact,exact_number(cap),exact_number(t))
                        block["exact_setup"]=relation
                    for method,(v,prices,stats) in bases.items():
                        stats["verified"]=outward_basis(geo,y,v,prices,trial,error)
                block["rows"].append(result)
                print(study,t,cap,ident,[(name,round(stats["distance"]/radius,6))
                                       for name,stats in methods.items()],flush=True)
            out["profiles"].append(block)
    write_json(HERE/f"{study}_directional.json",{**out,**elapsed(start)})


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("study", choices=["finite","native","all"], default="all", nargs="?")
    args = parser.parse_args()
    ctx.prec,ctx.threads=128,1
    with threadpool_limits(1):
        if args.study in ("finite","all"):
            compare("finite_0",finite_geometry(.6),{.75,1.,1.25})
        if args.study in ("native","all"):
            compare("native",native_geometry(),{.75,1.,1.25,1.5,2.,3.})
