"""Verify Neon’s same-box virtual illumination separator from stored inputs."""
from itertools import product
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[1]
sys.path.insert(0,str(HERE.parent/"A"))
from common import digest, forward, herm, load, save, star, write_json
from bounds import (acb, acb_mat, arb, ball, bdiag, bt, ctx, frob, high, low,
                    midpoint, propose_frame, select, shift, validate_frame)
import numpy as np
from scipy.linalg import eigh, solve

ctx.prec = 192
ctx.threads = 1

WEIGHTS = (0.0001,0.001,0.0001)
LAMBDAS = ((0.,0.,-0.0002),(0.0005,-0.0005,0.0005))
OLD = (0,4,8)
VIRTUAL = 12
T_SQUARED = .0001
MATERIAL = (2.2,.2,.2)
L, U = (0.,0.,0.),(1.5,1.5,1.5)


def plain_matrix(g,s,b,w,lam):
    """Expanded local Hermitian blocks on [J,v,y], independent of Neon’s Q formula."""
    n = len(g)
    m = b.shape[1]
    j = np.hstack((np.eye(n),np.zeros((n,m+s.shape[0]))))
    e = np.hstack((g,b,np.zeros((n,s.shape[0]))))
    p = np.hstack((s,np.zeros((s.shape[0],m)),-np.eye(s.shape[0])))
    jb, eb, pb = map(ball,(j,e,p))
    wb = [arb(float(v)) for v in w]
    lb = [arb(float(v)) for v in lam]
    edges = tuple(zip(map(arb,L),map(arb,U)))
    a = bt(jb)*bdiag([acb(-wval*(lo+up)/2,-lv/2)
                     for wval,lv,(lo,up) in zip(wb,lb,edges,strict=True)],eb)
    out = bt(pb)*pb+bt(jb)*bdiag(wb,jb)+a+bt(a)
    out += bt(eb)*bdiag([wval*lo*up for wval,(lo,up) in zip(wb,edges,strict=True)],eb)
    return out


def numeric_matrix(g,s,b,w,lam):
    n = len(g)
    m = b.shape[1]
    j = np.hstack((np.eye(n),np.zeros((n,m+s.shape[0]))))
    e = np.hstack((g,b,np.zeros((n,s.shape[0]))))
    p = np.hstack((s,np.zeros((s.shape[0],m)),-np.eye(s.shape[0])))
    w,lam = np.array(w),np.array(lam)
    a = star(j)@((-.75*w-.5j*lam)[:,None]*e)
    return herm(star(p)@p+star(j)@(w[:,None]*j)+a+star(a))


def signed(matrix):
    d,vec=eigh(herm(midpoint(matrix)),check_finite=False)
    gap=float(min(abs(d))*.7)
    if gap<=0:
        raise ArithmeticError("Eigenvalue interval can include zero")
    ident=ball(np.eye(len(d)))
    counts={}
    for sign,label in ((1,"positive"),(-1,"negative")):
        ids=np.flatnonzero(sign*d>0)
        counts[label]=len(ids)
        if len(ids):
            frame=vec[:,ids]/np.sqrt(abs(d[ids]))
            strict,_=validate_frame(sign*matrix-arb(gap)*ident,frame)
            if not strict>0:
                raise ArithmeticError("Interval signed inertia check failed")
    return dict(positive=counts["positive"],negative=counts["negative"],zero=0,
                min_absolute_eigenvalue_lower=gap)


def schur(exact,n=3):
    above,below=range(n),range(n,exact.nrows())
    K=select(exact,above,above)
    ab=select(exact,above,below)
    return K,select(exact,below,below)-bt(ab)*K.solve(ab)


def candidate(g,s,b,z,w,lam):
    Z=ball(z)
    Qform=plain_matrix(g,s,b,w,lam)
    K,Q=schur(Qform)
    sK=signed(K)
    H=bt(Z)*Q*Z
    hsign=signed(H)
    ed = eigh(herm(midpoint(H)),eigvals_only=True)[::-1]
    rank = hsign["positive"]-sK["negative"]
    assert rank >= 0
    proposed = ed[sK["negative"]:hsign["positive"]]
    slack = T_SQUARED-float(sum(proposed))
    assert slack > 0
    tail_upper = arb(0)
    upper_records = []
    for i,eig in enumerate(proposed,sK["negative"]+1):
        endpoint = float(eig+.05*slack/max(1,rank))
        shifted = signed(H-arb(endpoint)*ball(np.eye(z.shape[1])))
        assert shifted["positive"] < i
        tail_upper += arb(endpoint)
        upper_records.append(dict(eigenvalue_rank_descending=i,verified_upper=endpoint,
                                  shifted_inertia = shifted))
    t=arb(T_SQUARED)
    assert tail_upper < t
    Po=H-t*ball(np.eye(z.shape[1]))
    condition=signed(Po)
    assert condition["positive"]==sK["negative"]
    QZ=Q*Z
    Qc=Q-QZ*Po.solve(bt(QZ))
    Ry=Qc[4,4].real
    assert Ry>0
    center=-Qc[4,3]/Ry
    beta=abs(Qc[4,3])**2/Ry-Qc[3,3].real
    assert beta>0
    radius=(beta/Ry).sqrt()
    out=dict(weights=list(w),optical=list(lam),K_inertia=sK,H_old_inertia=hsign,
             H_old_tail_positive_part_upper=high(tail_upper),tail_eigenvalue_upper_proofs=upper_records,
             radius_old_squared=T_SQUARED,
             P_old_inertia=condition,saturated=True,R_y_enclosure=[low(Ry),high(Ry)],
             center_real=[low(center.real),high(center.real)],
             center_imaginary=[low(center.imag),high(center.imag)],
             beta_enclosure=[low(beta),high(beta)],
             radius_enclosure=[low(radius),high(radius)])
    return out,center,radius


def ordinary_candidate(g,s,b,z,key):
    w,lam=key
    A=numeric_matrix(g,s,b,w,lam)
    K,ab=A[:3,:3],A[:3,3:]
    nu=int(np.count_nonzero(eigh(K,eigvals_only=True)<0.))
    H=herm(star(z)@(A[3:,3:]-star(ab)@solve(K,ab,assume_a="her"))@z)
    evals=eigh(H,eigvals_only=True)[::-1][nu:]
    return (float(sum(np.maximum(evals,0.))),w,lam,nu,evals)


def strict_ordinary(g,s,b,z,key):
    numeric,w,lam,nu,evals=key
    if numeric<=T_SQUARED:
        return None
    A=plain_matrix(g,s,b,w,lam)
    lift = ball(np.block([[np.eye(3),np.zeros((3,3))],[np.zeros((5,3)),z]]))
    same_data=bt(lift)*A*lift
    total=arb(0)
    terms=[]
    for j,v in enumerate(evals,1):
        if v<=0.:
            continue
        t=float(v*(1-1.e-5))
        Aq=shift(same_data,3,t)
        frame,edge=propose_frame(midpoint(Aq),3+j)
        if frame is None:
            return None
        gap,_=validate_frame(Aq,frame)
        if not gap>0:
            return None
        total+=arb(t)
        terms.append(dict(j=j,t=t,positive_count=3+j,pullback_gap_lower=low(gap)))
    if not total>arb(T_SQUARED):
        return None
    return dict(weights=list(w),optical=list(lam),negative_K_count=nu,
                verified_tail_lower=low(total),verified_ratio_lower=low((total/arb(T_SQUARED)).sqrt()),
                verified_thresholds=terms)


def main():
    inp=PAPER/"inputs/Z_ALG_Research_Package/checks/active_ambiguity_closed_loop_operators.npz"
    raw=load(inp)
    g,s,b0=4*raw["G"],raw["S"],raw["Bd"]
    b=b0[:,[*OLD,VIRTUAL]]
    cold=forward(g,s,b[:,:3],np.array(MATERIAL))[0]
    z=np.vstack((np.eye(4)[:,:3],cold))
    chi=list(map(arb,MATERIAL))
    M=ball(np.eye(3))-bdiag(chi,ball(g))
    exact_old=ball(s)*M.solve(bdiag(chi,ball(b[:,:3])))
    old_residual=frob(exact_old-ball(cold))
    assert old_residual**2<arb(T_SQUARED)
    out = [candidate(g,s,b,z,WEIGHTS,lam) for lam in LAMBDAS]
    original_centers=[val[1] for val in out]
    direction=complex(original_centers[1]-original_centers[0])
    direction/=abs(direction)
    v=acb(direction.real,direction.imag)
    projected=[]
    for record,center,rad in out:
        proj=(v.conjugate()*center).real
        a,bb=proj-abs(v)*rad,proj+abs(v)*rad
        record["projected_outer_interval"]=[low(a),high(bb)]
        projected.append((a,bb))
    gap=projected[1][0]-projected[0][1]
    assert gap>0

    multipliers=[]
    for mu,spatial,opt in product((1e-4,1e-3,1e-2,1e-1,1.),
            ((1,1,1),(10,1,1),(1,10,1),(1,1,10),(10,.1,.1),(.1,10,.1),(.1,.1,10)),
            ((0,0,0),(2,0,0),(-2,0,0),(0,2,0),(0,-2,0),(0,0,2),(0,0,-2),(5,-5,5),(-5,5,-5))):
        multipliers.append((tuple(mu*x for x in spatial),tuple(mu*x for x in opt)))
    highest = sorted((ordinary_candidate(g,s,b,z,key) for key in multipliers),
                     key=lambda item:item[0],reverse=True)
    unconstrained=next((result for trial in highest if (result:=strict_ordinary(g,s,b,z,trial)) is not None),None)
    uniform_trials=[r for r in highest if len(set(r[1]))==1 and not any(r[2])]
    uniform=next((result for trial in uniform_trials if (result:=strict_ordinary(g,s,b,z,trial)) is not None),None)
    report=dict(status="PASS",precision_bits=ctx.prec,
                model="Neon three-cell algebraic perturbation of archived scalar-TM data. G=4G_original; S and B retained. No calibrated continuum interpretation.",
                operator_path=str(inp.relative_to(PAPER)),operator_digest_blake2b=digest(inp),
                distinction="Existing-support intersection from old data only. No superiority claim over an ordinary multiplier search.",
                old_columns=OLD,virtual_column=VIRTUAL,material_box=[L,U],synthetic_generator=MATERIAL,
                old_measured_residual_upper=high(old_residual),old_error_radius_squared=T_SQUARED,
                future_response_used=False,projected_direction_binary64=[v.real.__float__(),v.imag.__float__()],
                projected_gap_lower=low(gap),weights_checked=[row[0] for row in out],
                ordinary_weight_bank_size=len(multipliers),
                ordinary_bank_numeric_exclusion_count=sum(row[0]>T_SQUARED for row in highest),
                ordinary_best_strict=unconstrained,ordinary_uniform_spatial_strict=uniform,
                uniform_search_numeric_squared_tails=[row[0] for row in uniform_trials])
    write_json(HERE/"result.json",report)
    save(HERE/"inputs_and_projection.npz",G=g,S=s,B4=b,Y_old=cold,
         old_columns=np.array(OLD),virtual_column=np.array(VIRTUAL),projected_direction=np.array(direction),
         ell=np.array(L),upp=np.array(U),weights=np.array(WEIGHTS),optical=np.array(LAMBDAS),
         r_old_squared=np.array(T_SQUARED))
    print(json.dumps(report,indent=2))


if __name__=="__main__":
    main()
