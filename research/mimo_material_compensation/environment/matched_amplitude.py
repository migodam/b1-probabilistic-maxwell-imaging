"""Native 96² equal-core-amplitude control under two whole-image box priors."""
import json

from frequency_classes import (HERE, HARDWARE, LAYOUT, STUDY, AuditedState, VerifiedWave,
                               check_direction, choose_direction, down, elapsed, explicit_difference,
                               fro, now, np, rd, released_operator, up, write_json)
from frequency_classes_independent import evaluate
from frequency_classes_open import empty_network_coefficients, uniform_wave_error
from passive_replay import load_cover
from scipy.ndimage import label


SIDE=96
AMPLITUDE=2.
STEM="wave_96_2.0_1"
REGISTERED_BETAS=(.001,.003,.01)


def read_pair():
    with np.load(HERE/(STEM+".npz"),allow_pickle=False) as f:
        ix=f["ix"]
        base=f["group"]
        material=[]
        for key in ("plus","minus"):
            x=np.zeros(SIDE**2)
            x[ix]=f[key]
            material.append(x)
        roi=np.zeros(SIDE**2,bool)
        roi[ix[base==0]]=True
        r=float(f["r"])
    with np.load(STUDY/"ports/public/cases"/"n96_c01.npz",allow_pickle=False) as f:
        oid=str(f["observation_id"].item())
        if r!=float(f["radius"]):
            raise ArithmeticError("Frozen radius disagreement")
    if not (np.all(material[0][roi]==AMPLITUDE) and np.all(material[1][roi]==0.)
            and int(roi.sum())==164 and label(roi.reshape(SIDE,SIDE))[1]==1):
        raise ValueError("Unexpected registered core")
    if any(np.min(x)<0. or np.max(x)>4.5 for x in material):
        raise ValueError("Reference material outside physical box")
    return material,ix,roi,r,oid


def point_replay(states,band,r,raw,ix,materials):
    source=HERE/f"verified_responses_{SIDE}_{AMPLITUDE}_{STEM}_{band}.npz"
    with np.load(source,allow_pickle=False) as q:
        vs=(q["s1"],q["s0"])
        es=(float(q["e1"]),float(q["e0"]))
    opens=[]
    audits=[]
    for key,x,st,v,e in zip(("plus","minus"),materials,states,vs,es):
        p,j,s,audit=evaluate(SIDE,band,ix,x)
        # The bounded 25-port receipt encloses the exact open point response.
        value,err,_=load_cover(v,e,1.,0.,0.)
        exact_bound=up(err/np.sqrt(8)+rd(4)*fro(value))
        diff=up(fro(p-st.p))
        print("point",band,key,diff,"guard",up(exact_bound+st.response_error),flush=True)
        if diff>up(2.*exact_bound+st.response_error+rd(8)*(fro(p)+fro(st.p))):
            raise ArithmeticError("Bare/loaded point check did not agree")
        opens.append((value[8:16,:8]/np.sqrt(8),exact_bound))
        raw.update({f"stock_{band}_{key}":st.p,f"open_{band}_{key}":p,
                    f"currents_bare_{band}_{key}":j,f"scattering_{band}_{key}":s})
        audits.append(dict(material=key,difference_frozen_to_independent_Fro=diff,
                           independent_open_point_error_upper=exact_bound,**audit))
    delta=opens[0][0]-opens[1][0]
    pair_bound=up((fro(delta)+opens[0][1]+opens[1][1]+rd(4)*(fro(opens[0][0])+fro(opens[1][0])))/(2*r))
    if band==400:
        raw["midpoint_400_total"]=(opens[0][0]+opens[1][0])*.5
        raw["midpoint_400_subtracted"]=(raw["open_400_plus"]+raw["open_400_minus"])*.5
    return dict(independent_open_point_pair_upper_over_two_r=pair_bound,
                independent_open_point_pair_numeric_over_two_r=float(np.linalg.norm(delta)/(2*r)),
                bare_point_audits=audits,bounded_source=source.name)


def clip_budget(beta,states,jacs,d,ed,operr,r):
    j=[p[0] for p in jacs]
    x=[s.x for s in states]
    initial=d/np.linalg.norm(d)
    candidate,spent=choose_direction(d,j,*x,beta)
    possibilities=[check_direction(d,ed,l,states,jacs,beta) for l in (initial,candidate)]
    point,lam,projection,terms=max(possibilities,key=lambda q:q[0])
    stock=max(0.,down(point/(2*r)))
    if all(q["k_times_max_edit"]<1. for q in terms):
        errors=[uniform_wave_error(s,s.k_total,beta,operr) for s in states]
        penalty=up(sum(a["response_error_fro_upper"] for a in errors)/(2*r))
        physical=max(0.,down(stock-penalty-rd(5)*(stock+penalty)))
    else:
        errors=[]
        penalty=None
        physical=0.
    entry=dict(beta=beta,registered_beta=beta in REGISTERED_BETAS,stock_lower_over_two_r=stock,
               projection_lower=projection,optimized_direction_evaluations=spent,
               empirical_only_frozen_lower_over_two_r=max(0.,down(possibilities[0][0]/(2*r))),
               open_lower_over_two_r=physical,transfer_penalty_over_two_r=penalty,
               endpoint_stock_terms=terms,endpoint_open_uniform_terms=errors,
               core_mean_gap_lower=down(AMPLITUDE-2*beta),
               status="all_pairs_excluded" if physical>1. else "unresolved")
    print("beta",beta,"bounds",stock,physical,"transfer",penalty,flush=True)
    return entry,lam


def main():
    start=now()
    material,ix,roi,r,obs=read_pair()
    raw=dict(plus=material[0],minus=material[1],core=roi,active_union=ix,radius=r)
    cover_file=HERE/f"cover_{SIDE}_{AMPLITUDE}_{STEM}.json"
    cover=json.loads(cover_file.read_text())["outputs"]["400"]
    if not cover["arbitrary_wave_upper_over_2r"]<1:
        raise ArithmeticError("The independent 400 MHz all-load witness is missing")
    result=dict(side=SIDE,nominal_core_mean_gap=AMPLITUDE,core_cells=int(roi.sum()),
                core_area_m2=float(roi.sum())*(2/SIDE)**2,
                sample_file=STEM+".npz",radius_source_id=obs,radius_sqrtW=r,hardware=HARDWARE,
                contract="Both independently varied all-pixel boxes of common L-infinity radius, intersected with [0,4.5]. Single A layout. Real frequency-independent material contrast and fixed nominal instrument.",
                outputs={})
    for band in (400,460):
        work=now()
        tag,op=released_operator(SIDE,band)
        wave=VerifiedWave(op)
        states=[AuditedState(wave,op["S"],op["E0"],x,require_k=band==460) for x in material]
        opp,empty=empty_network_coefficients(op)
        raw[f"open_empty_offset_{band}"]=empty
        point=point_replay(states,band,r,raw,ix,material)
        audit=[]
        for state in states:
            info=dict(state.log,native_norm_audit=getattr(state,"norm_audit",None))
            if band==400:
                info["k_total_upper"]=None
            audit.append(info)
        report=dict(operator_id=tag,open_coefficient_transfer=opp,anchor_audits=audit,**point)
        if band==400:
            if point["independent_open_point_pair_upper_over_two_r"]>=1:
                raise ArithmeticError("400 MHz common midpoint could not be certified")
            report.update(open_common_midpoint_residual_over_r_upper=point["independent_open_point_pair_upper_over_two_r"],
                          independent_all_16_coherent_and_passive_load_upper_over_two_r=cover["arbitrary_wave_upper_over_2r"],
                          independent_all_load_receipt=cover_file.name)
        else:
            jacs=[s.jacobian(LAYOUT) for s in states]
            d,base_error=explicit_difference(states[0].p,states[1].p)
            ed=up(base_error+sum(s.response_error for s in states))
            rows=[]
            for beta in (0.,)+REGISTERED_BETAS:
                row,lam=clip_budget(beta,states,jacs,d,ed,opp,r)
                key=f"lambda_{len(rows):02d}"
                row["lambda_key"]=key
                raw[key]=lam
                rows.append(row)
            if rows[-1]["status"]=="unresolved":
                for beta in (.004,.006,.008,.009):
                    row,lam=clip_budget(beta,states,jacs,d,ed,opp,r)
                    key=f"lambda_{len(rows):02d}"
                    row["lambda_key"]=key
                    raw[key]=lam
                    rows.append(row)
            report.update(budgets=sorted(rows,key=lambda z:z["beta"]),
                          largest_verified_beta=max(z["beta"] for z in rows if z["status"]=="all_pairs_excluded"))
        report["cost"]=dict(**elapsed(work),material_LU_frozen=2,material_LU_bare=2,
                            frozen_inverse_rhs=sum(len(s.active) for s in states),
                            frozen_forward_rhs=16,bare_forward_rhs=50)
        result["outputs"][str(band)]=report
    result.update(**elapsed(start))
    write_json(HERE/"matched_amplitude.json",result)
    with (HERE/"matched_amplitude.npz").open("wb") as f:
        np.savez_compressed(f,**raw)


if __name__=="__main__":
    main()
