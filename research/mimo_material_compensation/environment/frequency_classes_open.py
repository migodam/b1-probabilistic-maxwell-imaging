"""Transfer each certified stock-operator material box to physical open circuitry."""
from frequency_classes import (HERE, AuditedState, VerifiedWave, down, elapsed, explicit_difference,
                               fro, material, now, np, opnorm, rd, released_operator, up, write_json)
from flint import arb
import json


def empty_network_coefficients(op):
    """Enclose the bare 24-port law exact in the stored input coefficients."""
    z=op["Z_AA"]
    R=float(op["reference_ohms"])
    n=len(z)
    if n!=24 or R!=50.:
        raise ValueError("Expected the fixed 24-port nominal network")
    m=z+R*np.eye(n)
    material_arithmetic=up(rd(2)*(fro(z)+R*np.sqrt(n)))
    q=np.linalg.solve(m,np.eye(n))
    qnorm=opnorm(q)
    qr=m@q
    defect=up(fro(np.eye(n)-qr)+rd(n+3)*(np.sqrt(n)+fro(m)*fro(q))+
              material_arithmetic*qnorm)
    if defect>=1:
        raise ArithmeticError("Uncertified open empty-network inverse")
    qerr=up(qnorm*defect/(1.-defect))
    e,read=op["E0_bare"],op["S_bare"]
    s_hat=q@read
    s_hat_norm=fro(s_hat)
    e_hat=e@q
    en=opnorm(e)
    rn=opnorm(read)
    h=5.*np.eye(24,8)
    if not np.array_equal(h,op["H"]):
        raise ValueError("Drive should be exact 5 V and 1 W available")
    root=R**.5
    rooterr=up(float(abs(arb(R).sqrt()-arb(root)).upper()))
    x=np.arange(8,16)
    left=1j*e_hat
    right=read
    b_hat=e_hat@h
    c_hat=-1j*root*s_hat[x]
    l_err=up(en*qerr+rd(n)*fro(e)*fro(q))
    delta_left=up(fro(left-op["low_rank_left"])+l_err+rd(3)*(fro(left)+fro(op["low_rank_left"])))
    r_err=up(fro(right-op["low_rank_right"])+rd(3)*(fro(right)+fro(op["low_rank_right"]))) if not np.array_equal(right,op["low_rank_right"]) else 0.
    gerr=up(delta_left*rn+opnorm(op["low_rank_left"])*r_err)
    berr=up(fro(b_hat-op["E0"])+l_err*opnorm(h)+rd(n)*fro(e_hat)*fro(h)+
            rd(3)*(fro(b_hat)+fro(op["E0"])))
    read_rounding=up(qerr*fro(read)+rd(n)*fro(q)*fro(read))
    cerr=up(fro(c_hat-op["S"])+(root+rooterr)*read_rounding+rooterr*s_hat_norm+
            rd(4)*(fro(c_hat)+fro(op["S"])+root*s_hat_norm))
    # One 8 x 8 empty offset. H is exactly sqrt(8)-scaled in power waves.
    # Prefer -sqrt(R) (QH)_rx so the 5 V drive stays explicitly identical.
    port=q@h
    empty=-root*port[x]
    empty_err=up((root+rooterr)*(qerr*fro(h)+rd(n)*fro(q)*fro(h))+
                 rooterr*fro(port)+rd(4)*root*fro(port))
    return dict(g_delta_op_upper=gerr,incident_b_delta_Fro_upper=berr,
                receiver_s_delta_op_upper=cerr,empty_total_Fro_error_upper=empty_err,
                open_empty_inverse_difference_norm_upper=qerr,open_empty_inverse_defect_upper=defect,
                open_empty_matrix_rounding_error_upper=material_arithmetic,
                stored_nominal_vs_bare_empty_matrix_diagnostic=fro(op["M_nominal"]-m),
                low_rank_right_exactly_shared=(r_err==0),one_watt_drive_exact=True),empty


def uniform_wave_error(state,k,beta,operator_error):
    """One-sided response difference, uniformly over all material coordinates."""
    m=up(min(4.5,float(state.x.max())+beta))
    contraction=up(k*beta)
    if contraction>=1:
        raise ArithmeticError("The material neighborhood lacks a verified inverse")
    dressed=up(k/(1.-contraction))
    inverse=up(1.+m*dressed)
    total_field=up((fro(state.e)+state.e_error)/(1.-contraction))
    receive_lift=up((fro(state.u)+state.u_error)/(1.-contraction))
    j=up(m*total_field)
    g_delta=operator_error["g_delta_op_upper"]
    b_delta=operator_error["incident_b_delta_Fro_upper"]
    s_delta=operator_error["receiver_s_delta_op_upper"]
    q=up(inverse*m*g_delta)
    if q>=1:
        raise ArithmeticError("The open operator perturbation lacks a verified inverse")
    j_open=up((j+inverse*m*b_delta)/(1.-q))
    error=up(receive_lift*m*(b_delta+g_delta*j_open)+s_delta*j_open)
    return dict(response_error_fro_upper=error, material_resolvent_contraction_upper=contraction,
                dressed_k_upper=dressed, stored_inverse_upper=inverse,
                whole_class_stored_total_incident_field_Fro_upper=total_field,
                whole_class_stored_receiver_lift_Fro_upper=receive_lift,
                stored_current_Fro_upper=j,open_current_Fro_upper=j_open,
                open_operator_perturbation_contraction_upper=q)


def phase_rows(side):
    prior=json.loads((HERE/f"frequency_classes_{side}.json").read_text())
    stock=prior["frequency"]["460"]
    main=[(b["beta"],b["certified_lower_over_two_r"],"native standard directional")
          for b in stock["budget"] if b["beta"]>0 and b["status"]=="verified_excluded"]
    if side==32:
        more=json.loads((HERE/"frequency_classes_32_refine.json").read_text())
        main.extend((b["beta"],b["lower_over_two_r"],"verified signed quadratic")
                    for b in more["rows"] if b["status"]=="verified_excluded")
    return prior,stock,main


def main():
    start=now()
    out=[]
    arrays={}
    for side in (32,96):
        work=now()
        t,src,r,obs,ix,roi,a,b=material(side)
        prior,stock,targets=phase_rows(side)
        perband={}
        for band in (400,460):
            oid,op=released_operator(side,band)
            coefficient_cert,empty=empty_network_coefficients(op)
            with np.load(HERE/f"frequency_classes_{side}.npz",allow_pickle=False) as raw:
                stock_empty=raw[f"p_{band}_empty_total"]
            shift_diagnostic=fro(stock_empty-empty)
            enclosure=up(coefficient_cert["empty_total_Fro_error_upper"]+
                         rd(48)*(fro(stock_empty)+fro(empty)+1.))
            if shift_diagnostic>enclosure:
                raise ArithmeticError("Primed empty-field translation is inconsistent")
            arrays[f"empty_open_{side}_{band}"]=empty
            perband[str(band)]=dict(operator=oid,coefficient_certificate=coefficient_cert,
                 empty_numerical_shift_norm=shift_diagnostic,empty_shift_compatibility_guard=enclosure,
                 constant_cancels_between_two_open_total_fields=True)
            if band==400:
                cover=prior["frequency"]["400"]
                # Use the independent full 25-port open response uncertainty.
                # A common exact pair midpoint lies in both r-balls as soon
                # as their independently enclosed distance is less than 2r.
                open_upper=up(cover["independent_open_band_center_distance_over_two_r"]+
                              cover["independent_open_error_over_two_r"]+
                              rd(160)*cover["independent_open_band_center_distance_over_two_r"])
                if open_upper>=1:
                    raise ArithmeticError("Open 400 MHz midpoint was not certified")
                perband[str(band)].update(open_400_common_midpoint_residual_over_r_upper=open_upper,
                                         zero_offset_and_open_total_contract_share_pair_distance=True)
                continue
            wave=VerifiedWave(op)
            states=[AuditedState(wave,op["S"],op["E0"],x,require_k=False) for x in (a,b)]
            k=[stock[tag]["k_total_upper"] for tag in ("plus","minus")]
            budgets=[]
            for beta,lower,method in sorted(targets):
                terms=[uniform_wave_error(s,kappa,beta,coefficient_cert) for s,kappa in zip(states,k)]
                perturbation=up(sum(q["response_error_fro_upper"] for q in terms)/(2*r))
                actual=max(0.,down(lower-perturbation-rd(4)*(lower+perturbation)))
                budgets.append(dict(beta=beta,stock_lower_over_two_r=lower,uniform_transfer_penalty_over_two_r=perturbation,
                                    open_lower_over_two_r=actual,stock_certificate=method,endpoint_terms=terms,
                                    verified_open_excluded=actual>1))
                print(side,beta,actual,perturbation,flush=True)
            perband[str(band)].update(budgets=budgets,
                      largest_transferred_budget=max(row["beta"] for row in budgets if row["verified_open_excluded"]),
                      uniform_field_validation_material_LU=2,
                      uniform_field_validation_inverse_right_hand_sides=sum(len(s.active) for s in states))
        out.append(dict(side=side,reference=src,radius=r,radius_source_id=obs,core_cells=int(roi.sum()),
                        contract="All pixels in both reference neighborhoods. 25th physical port open and fixed original 24-port network. Fixed nominal bare coefficients; shared real, frequency-independent contrast.",
                        bands=perband,**elapsed(work)))
    with (HERE/"frequency_classes_open.npz").open("wb") as f:
        np.savez_compressed(f,**arrays)
    write_json(HERE/"frequency_classes_open.json",dict(records=out,**elapsed(start)))


if __name__=="__main__":
    main()
