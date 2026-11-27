"""Frozen one-shot complex scalar designs and one-sided whole-class bounds."""
from datetime import datetime, timezone
from itertools import combinations,product
from pathlib import Path
import argparse
import json
import os
import sys

for variable in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","VECLIB_MAXIMUM_THREADS","MKL_NUM_THREADS"):
    os.environ[variable]="1"
sys.dont_write_bytecode=True

HERE=Path(__file__).resolve().parent
ENV=HERE.parent
sys.path.insert(0,str(ENV))
from frequency_classes import (LAYOUT, AuditedState, VerifiedWave, check_direction, down, elapsed,
                               explicit_difference, fro, material, now, np, rd, released_operator, up, write_json)
from frequency_classes_open import empty_network_coefficients, uniform_wave_error
from matched_amplitude import read_pair
from flint import arb


BAND=460
BETA=.01
PHASE=np.array([1.,1.j,-1.,-1.j])
N_KEEP=40


def stamp():
    return datetime.now(timezone.utc).isoformat(timespec="microseconds")


def lower_ball(x):
    return float(np.nextafter(float(x.lower()),-np.inf))


def upper_ball(x):
    return float(np.nextafter(float(x.upper()),np.inf))


def upper(x):
    return float(np.nextafter(x+rd(3)*abs(x),np.inf))


def real(m):
    return np.r_[m.real.ravel(),m.imag.ravel()]


def matrix(v):
    return v[:64].reshape(8,8)+1j*v[64:].reshape(8,8)


def descriptor(z):
    return [[float(y.real),float(y.imag)] for y in z]


def distinct_phase_choices():
    vectors=[]
    for support in combinations(range(8),4):
        for phases in product(range(4),repeat=3):
            w=np.zeros(8,complex)
            w[np.asarray(support)]=PHASE[[0,*phases]]*.5
            vectors.append(w)
    return np.asarray(vectors)


def top_pairs(weights,signal):
    """Rank all canonical aperture pairs by their public complex projection."""
    left=weights.conj()@signal
    best=[]
    chunk=120
    for low in range(0,len(weights),chunk):
        value=left[low:low+chunk]@weights.T
        max_bound=np.maximum(abs(value.real),abs(value.imag))
        top=np.argpartition(max_bound.ravel(),-N_KEEP)[-N_KEEP:]
        ij=np.unravel_index(top,max_bound.shape)
        for a,b in zip(*ij):
            p=value[a,b]
            phase=int(np.argmax((p*PHASE).real))
            best.append((max_bound[a,b],int(a+low),int(b),phase))
    return sorted(best,reverse=True, key=lambda x:x[0])[:N_KEEP]


def packet(side):
    if side==32:
        t,filename,r,obs,ix,roi,p,m=material(32)
        inputs=json.loads((ENV/"frequency_classes_32.json").read_text())['frequency']['460']
        record=next(y for y in inputs["budget"] if y["beta"]==BETA)
        data=ENV/"frequency_classes_32.npz"
        source="frequency_classes_32.json"
        xs=[p,m]
    else:
        xs,ix,roi,r,obs=read_pair()
        filename="wave_96_2.0_1.npz"
        inputs=json.loads((ENV/"matched_amplitude.json").read_text())["outputs"]["460"]
        record=next(y for y in inputs["budgets"] if y["beta"]==BETA)
        data=ENV/"matched_amplitude.npz"
        source="matched_amplitude.json"
    with np.load(data,allow_pickle=False) as raw:
        direction=raw[record["lambda_key"]]
    return xs,ix,roi,r,obs,filename,direction,source


def scalar_bound(u,v,states,jacs,dp,derr,errors,empty,empty_error,r):
    # Each real/imaginary entry is exactly 0, signed 1/4 or signed 1.
    if sum(abs(u)**2)!=1 or sum(abs(v)**2)!=1:
        raise ValueError("Unphysical launch or receive norm")
    z=real(u[:,None]*v[None,:].conj())
    val,w,projection,terms=check_direction(dp,derr,z,states,jacs,BETA)
    kk=np.flatnonzero(z)
    ratio=float(w[kk[0]]/z[kk[0]])
    if not 0.<ratio<=1 or not np.array_equal(w,z*ratio):
        raise ArithmeticError("The verified direction did not contract by a common factor")
    model_penalty=up(sum(errors))
    # Because inner q has norm at most one and ratio<=1, a positive lower
    # bound only becomes larger on undoing the internal contraction.
    positive=max(0.,down(val-model_penalty-rd(5)*(val+model_penalty)))
    physical0=max(0.,lower_ball(arb(positive)*arb(8).sqrt()))
    raw_factor=arb(8).sqrt()/arb(ratio)
    center_scores=[]
    for n,state in enumerate(states):
        point=real(state.p+empty)
        arithmetic=up(rd(140)*fro(w)*fro(point)+rd(4)*(fro(state.p)+fro(empty)))
        t=terms[n]
        if t["nonlinear"] is None:
            raise ArithmeticError("No material resolvent contraction")
        loss=up(t["safe_support"]+t["derivative_term_error"]+t["nonlinear"]+arithmetic+
                errors[n]+state.response_error+empty_error)
        value=float(w@point)
        if n==0:
            center_scores.append(lower_ball(raw_factor*arb(down(value-loss))))
        else:
            center_scores.append(upper_ball(raw_factor*arb(upper(value+loss))))
    high,low=center_scores
    sep_endpoint=max(0.,down(high-low-rd(4)*(abs(high)+abs(low))))
    lower=max(physical0,sep_endpoint)
    hi_r=down(high-r-rd(4)*(abs(high)+r))
    lo_r=upper(low+r+rd(4)*(abs(low)+r))
    threshold=(lo_r+hi_r)/2 if hi_r>lo_r else None
    if threshold is not None and not lo_r<threshold<hi_r:
        raise ArithmeticError("Invalid rounded threshold")
    sufficient=upper_ball((2*arb(r)/arb(lower))**2) if lower>0 else None
    return dict(
        u=descriptor(u),v=descriptor(v),nonzero_rx=int(np.count_nonzero(u)),nonzero_tx=int(np.count_nonzero(v)),
        native_unit_norm=True,internal_verifier_contraction=ratio,stock_directional_lower=val,
        open_model_transfer_subtraction_old_units=model_penalty,projected_center_lower_old_units=projection,
        endpoint_terms_stock_units=terms,
        physical_lower_by_scaled_difference=physical0,
        physical_lower_by_separate_thresholds=sep_endpoint,
        scalar_gap_lower=lower,scalar_gap_over_two_r=max(0.,down(lower/(2*r))),
        high_noiseless_score_lower=high,low_noiseless_score_upper=low,
        high_noisy_score_lower=hi_r,low_noisy_score_upper=lo_r,frozen_threshold=threshold,
        sufficient_normalized_energy_upper=sufficient,
        passes_fixed_unit_energy=lower>2*r and threshold is not None)


def build(side):
    clock=now()
    fname=HERE/f"proposal_{side}.json"
    if fname.exists() or (HERE/f"readings_{side}.npz").exists():
        raise FileExistsError("An existing design or acquisition is sealed")
    xs,ix,roi,r,obs,source,dual,old_json=packet(side)
    operator,op=released_operator(side,BAND)
    wave=VerifiedWave(op)
    state_start=now()
    states=[AuditedState(wave,op["S"],op["E0"],x) for x in xs]
    jac=[s.jacobian(LAYOUT) for s in states]
    state_cost=elapsed(state_start)
    opened,empty=empty_network_coefficients(op)
    errors=[uniform_wave_error(s,s.k_total,BETA,opened) for s in states]
    open_errors=[e["response_error_fro_upper"] for e in errors]
    diff,delta=explicit_difference(states[0].p,states[1].p)
    delta=up(delta+sum(s.response_error for s in states))
    signal=states[0].p-states[1].p
    weights=distinct_phase_choices()
    rank_timer=now()
    nominations={}
    for name,q in (("public center difference",signal),("prior certified direction",matrix(dual))):
        called=top_pairs(weights,q)
        for value,u,v,phase in called:
            nominations.setdefault((u,v,phase),[]).append(name)
    ranking_cost=elapsed(rank_timer)
    verify_timer=now()
    trials=[]
    starter=(states,jac,diff,delta,open_errors,empty,opened["empty_total_Fro_error_upper"],r)
    for key,reasons in nominations.items():
        u,v,phase=key
        candidate=scalar_bound(weights[u],weights[v]*PHASE[phase],*starter)
        trials.append(dict(origin=reasons,**candidate))
    picks=[]
    eye=np.eye(8,dtype=complex)
    for rr,tx in product(range(8),repeat=2):
        for phase in range(4):
            result=scalar_bound(eye[rr],eye[tx]*PHASE[phase],*starter)
            picks.append(dict(receive_index=rr+8,transmit_index=tx,quadrature=phase,**result))
    verify_cost=elapsed(verify_timer)
    def winning(q):
        return max(q,key=lambda z:z["scalar_gap_lower"])
    best4=winning(trials)
    best1=winning(picks)
    # Retain the best scalar of each hardware arrangement. Independent
    # acquisition is permitted for any arrangement whose threshold is frozen.
    post=dict(four_by_four=best4,single_pair=best1)
    print(side,"four",best4["scalar_gap_over_two_r"],best4["sufficient_normalized_energy_upper"],
          "single",best1["scalar_gap_over_two_r"],best1["sufficient_normalized_energy_upper"],flush=True)
    if not any(x["passes_fixed_unit_energy"] for x in post.values()):
        print(side,"No positive unit-power contract; search results remain unresolved",flush=True)
    receipt=dict(side=side,frequency_MHz=BAND,operator=operator,known_material_input=source,
                 prior_direction_receipt=old_json,calibration_tolerance_per_pixel=BETA,
                 core_area_m2=float(np.sum(roi))*(2/side)**2,core_cells=int(np.sum(roi)),nominal_core_mean_gap=2,
                 observation_radius_source=obs,absolute_scalar_noise_radius=r,
                 design_access_scope="Public references and previously saved class directions only. No realized specimen readings.",
                 search=dict(unit_vectors_per_side=len(weights),canonical_pair_projections_for_each_objective=len(weights)**2,
                             fully_certified_four_port_directions=len(trials),single_port_pairs=64,
                             fully_certified_single_port_directions=len(picks),rank_top_per_objective=N_KEEP,
                             ranking_cost=ranking_cost,verification_cost=verify_cost),
                 physical_contract=dict(unit_energy=1.,coherent_launches_per_action=1,complex_scalar_readouts_per_action=1,
                                        noise="Fixed absolute complex modulus ball",
                                        normalization="Unscaled full port S gives physical u* S v; legacy 8-column blocks were each normalized by sqrt(1/8).",
                                        open_extra_port=True,held_termination_ohms=50.),
                 states=[s.log for s in states],open_coefficient_errors=opened,uniform_endpoint_transfer=errors,
                 best=post,matrix_state_cost=state_cost,freeze_utc=stamp(),**elapsed(clock))
    write_json(HERE/f"trials_{side}.json",dict(side=side,four_by_four=trials,single_pairs=picks))
    with (HERE/f"reference_responses_{side}.npz").open("xb") as f:
        np.savez_compressed(f,p_plus=states[0].p,p_minus=states[1].p,empty=empty,core=roi,
                            reference_plus=xs[0],reference_minus=xs[1],dual=dual)
    write_json(fname,receipt)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("side",type=int,choices=(32,96),nargs="+")
    a=p.parse_args()
    for side in a.side:
        build(side)
