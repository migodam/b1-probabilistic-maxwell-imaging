"""Outward receiver-placement perturbation bound for the saved g=1 witness."""
import json
from math import log2
from time import perf_counter, process_time

from flint import acb, acb_mat, arb

from continuous_band import Physical, dist, hankels, lower, norm, pt, rowsum, up
from data14 import ROOT, load_cases
from geometry_bounds import real_jacobian, spectral_bounds

OUTPUT = ROOT/"results"
WITNESS = OUTPUT/"continuous_band/fresh6-16_known_coherent.json"
# The final displacement must lie in this preliminary radius, in geometry units.
PILOT_RADIUS = arb("0.005")


def reconstruct_cover(saved):
    band = arb(str(saved["half_band"]))
    pieces = []
    for c,w,_,coherent_upper in saved["terminal_intervals"]:
        depth = round(log2(saved["half_band"]/w))
        ratio = 1 << depth
        index = round(c*ratio/saved["half_band"])
        assert index % 2 and abs(index)<ratio
        pieces.append((index, depth, pt(coherent_upper)))
    depth_max = max(d for _,d,_ in pieces)
    integer_pieces = sorted(((i-1)<<(depth_max-d),(i+1)<<(depth_max-d)) for i,d,_ in pieces)
    assert integer_pieces[0][0] == -(1<<depth_max)
    assert integer_pieces[-1][1] == (1<<depth_max)
    assert all(x[1]==y[0] for x,y in zip(integer_pieces,integer_pieces[1:]))
    return [(band*i/(1<<d),band/(1<<d),v) for i,d,v in pieces]


def diagonal(material):
    return acb_mat([[1/pt(m) if i==j else 0 for j in range(len(material))]
                    for i,m in enumerate(material)])


def receiver_bounds(phy, d, k0, center, width, rho):
    """S_i Lipschitz constants and bounds on any real unit-direction second derivative."""
    high = (1+center+width).upper()
    low = (1+center-width).lower()
    kr = phy.self_x*(1+center+arb(0,width.upper()))
    _, j1 = hankels(phy.self_x,center,width)
    prefactor = (arb.pi()/2*kr*abs(j1.real)).upper()
    lipschitz, hessian = [], []
    for distance in d:
        assert distance>rho and low>0
        scaled_width = (width+high*rho/distance).upper()
        hk0,hk1 = hankels(k0*distance,center,scaled_width)
        argument = k0*distance*(1+center+arb(0,scaled_width))
        assert argument>0
        fprime = (prefactor*(k0*high)*abs(hk1)).upper()
        fsecond = (prefactor*(k0*high)**2*abs(hk0-hk1/argument)).upper()
        lipschitz.append(fprime)
        hessian.append(max(fsecond,(fprime/(distance-rho)).upper()))
    combined = lambda x: sum(v**2 for v in x).upper().sqrt().upper()
    return combined(lipschitz),combined(hessian)


def old_analysis(phy, diagonals, case, k0, receiver, d):
    G,B,S,_ = phy.kernels(arb(0),arb(0))
    L,H = receiver_bounds(phy,d,k0,arb(0),arb(0),PILOT_RADIUS)
    eps = pt(case.additive_radius)
    measured = acb_mat([[acb(pt(z.real),pt(z.imag)) for z in case.old]])
    _,j1 = hankels(phy.self_x,arb(0),arb(0))
    prefactor = acb(0,arb.pi()/2)*phy.self_x*j1.real
    s_gradient = []
    for axis in range(2):
        entries = []
        for di,cell in zip(d,case.geometry["centers"],strict=True):
            _,h1 = hankels(k0*di,arb(0),arb(0))
            entries.append(-prefactor*k0*h1*(pt(receiver[axis])-pt(cell[axis]))/di)
        s_gradient.append(acb_mat([entries]))
    results = []
    for diag in diagonals:
        current = (diag-G[0]).inv()*B[0]
        current_norm = norm(current)
        residual = S[0]*current-measured
        res0 = norm(residual)
        derivative = acb_mat([[z[0,i] for z in [x*current for x in s_gradient]]
                              for i in range(current.ncols())])
        _,singular_max = spectral_bounds(real_jacobian(derivative),False)
        projection = residual.conjugate()*derivative
        slope = norm(acb_mat([[z.real for z in projection.entries()]]))
        results.append((res0,current_norm*L,slope,pt(singular_max),current_norm*H/2))
    return eps,results,L,H


def new_analysis(phy, diagonals, pieces, k0, d):
    new = []
    for center,width,nominal in pieces:
        G,B,_,P = phy.kernels(center,width)
        midpoint = phy.kernels(center,arb(0))[3]
        # The nominal certificate already used this center with a positive Gram.
        inv_power = rowsum(midpoint.inv())
        eta = inv_power*rowsum(P-midpoint)
        if not eta<1:
            raise ArithmeticError("Insufficient outward bound for inverse power Gram")
        Q = [(z-G[0]).inv() for z in diagonals]
        J = [r*B[0] for r in Q]
        delta = diagonals[1]-diagonals[0]
        current_diff = min(norm(J[0]-J[1]),norm(Q[0]*delta*J[1]),norm(Q[1]*delta*J[0]))
        L,_ = receiver_bounds(phy,d,k0,center,width,PILOT_RADIUS)
        coefficient = (L*current_diff*(3*inv_power/(1-eta)).sqrt()).upper()
        new.append((center,width,nominal,coefficient,current_diff,inv_power/(1-eta)))
    return new


def safe_limit(check):
    lo,hi = 0.,float(PILOT_RADIUS.mid())
    for _ in range(45):
        mid = (lo+hi)/2
        if check(pt(mid)): lo=mid
        else: hi=mid
    return lo


def old_second_order(shift,row):
    res,_,linear_projection,lin_norm,hess = row
    return (res**2+2*shift*linear_projection+(shift*lin_norm)**2).sqrt() + hess*shift**2


def main():
    wall,cpu = perf_counter(),process_time()
    saved = json.loads(WITNESS.read_text())
    assert saved["certified_continuum"] and saved["gains_zeta_ri"]==[[0.,0.],[0.,0.]]
    case = next(c for c in load_cases() if c.case_id==saved["case_id"])
    dgs = [diagonal(m) for m in saved["materials"]]
    receiver = case.geometry["receivers"][0]
    d = [dist(receiver,x) for x in case.geometry["centers"]]
    k0 = pt(case.geometry["k"])
    lam = 2*arb.pi()/k0
    phi = Physical(case)
    eps,old,L,H = old_analysis(phi,dgs,case,k0,receiver,d)
    checked = new_analysis(phi,dgs,reconstruct_cover(saved),k0,d)
    new_eps = pt(saved["threshold_2r"])
    first_order_max = min(lower((eps-r)/l) for r,l,*_ in old)
    refined = [safe_limit(lambda shift: old_second_order(shift,row)<eps) for row in old]
    new_margins = [lower((new_eps-r)/l) for _,_,r,l,*_ in checked]
    argmin = min(range(len(new_margins)),key=new_margins.__getitem__)
    new_max = new_margins[argmin]
    limit = min(*refined,new_max,float(PILOT_RADIUS.mid()))
    # Use a compact decimal radius strictly inside all proven limits.
    rounded = float(f"{.999*limit:.3g}")
    if rounded>=limit:
        rounded = float(f"{.997*limit:.3g}")
    rho = pt(rounded)
    assert rho>0 and rho<PILOT_RADIUS
    assert all(old_second_order(rho,row)<eps for row in old)
    assert all(r+rho*l<new_eps for _,_,r,l,*_ in checked)
    old_plain = [
        dict(residual_upper=up(row[0]),margin_lower=lower(eps-row[0]),
             norm_product_coefficient_upper=up(row[1]),
             residual_direction_derivative_norm_upper=up(row[2]),
             response_bilinear_linear_norm_upper=up(row[3]),
             response_quadratic_remainder_coefficient_upper=up(row[4]),
             refined_radius=refined[i],
             shifted_residual_upper=up(old_second_order(rho,row)))
        for i,row in enumerate(old)]
    new_at_rho = max((entry[2]+rho*entry[3]).upper() for entry in checked)
    result = dict(case=saved["case_id"], source=str(WITNESS.relative_to(ROOT)),
                  arbitrary_direction=True, nominal_receiver=receiver, wavelength=up(lam),
                  unit="archived geometry length unit", fixed_gain=1., pilot_radius=up(PILOT_RADIUS),
                  certified_receiver_radius=rounded, certified_radius_over_wavelength=lower(rho/lam),
                  old_norm_product_radius_lower=first_order_max, old_refined_radius_lower=min(refined),
                  new_only_radius_lower_without_subdividing=min(new_max,lower(PILOT_RADIUS)),
                  smallest_receiver_cell_distance=lower(min(d)), old=old_plain,
                  old_operator_lipschitz_upper=up(L),old_operator_hessian_upper=up(H),
                  old_noise_budget=case.additive_radius,
                  new=dict(half_band=saved["half_band"],interval_count=len(checked),
                      largest_shifted_coherent_gap_upper=up(new_at_rho),
                      largest_shifted_gap_over_threshold_upper=up(new_at_rho/new_eps),
                      threshold_lower=saved["threshold_2r"],
                      nominal_upper_ratio=saved["upper_ratio"],
                      interval_limiting_new_radius=dict(center=float(checked[argmin][0]),
                          half_width=up(checked[argmin][1]),gap_upper=up(checked[argmin][2]),
                          perturbation_per_unit_length_upper=up(checked[argmin][3]),
                          current_difference_frobenius_upper=up(checked[argmin][4]),
                          inverse_gram_norm_upper=up(checked[argmin][5])),
                      interval_operators=phi.kernel_calls, full_current_matrix_inverses=2+2*len(checked)),
                  runtime_seconds=dict(wall=perf_counter()-wall,cpu=process_time()-cpu))
    (OUTPUT/"receiver_shift_check.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print(json.dumps(result,indent=2))


if __name__=="__main__":
    main()
