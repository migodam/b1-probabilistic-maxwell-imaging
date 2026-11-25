"""Fixed-budget, no-truth screen of five unsettled public Phase-A nodes."""
from itertools import combinations
import json
from time import perf_counter

from verify import HERE, signed, schur
from common import REGION, PUBLIC, RESULTS, digest, herm, load, star, write_json
from bounds import (arb, ball, bt, ctx, high, low, sector, sector_ball, select,
                    weight_vector)
import numpy as np
from scipy.linalg import eigh, solve
from run import PARAMS

ctx.prec = 192


def system(g,s,b,ell,upp,weights,optical,exact=False):
    n,m = b.shape
    p = len(s)
    current=np.hstack((np.eye(n),np.zeros((n,m+p))))
    field=np.hstack((g,b,np.zeros((n,p))))
    port=np.hstack((s,np.zeros((p,m)),-np.eye(p)))
    if exact:
        port,cur,fld=map(ball,(port,current,field))
        return bt(port)*port+sector_ball(cur,fld,ell,upp,weights,optical)
    return herm(star(port)@port+sector(current,field,ell,upp,weights,optical))


def numeric(g,s,b,ell,upp,w,opt,z,radius):
    n,m=b.shape
    radius2=radius**2
    A=system(g,s,b,ell,upp,w,opt)
    K,ab=A[:n,:n],A[:n,n:]
    Kd=eigh(K,eigvals_only=True)
    if min(abs(Kd))<1.e-11:
        return None,"numeric_K_singular"
    Q=herm(A[n:,n:]-star(ab)@solve(K,ab,assume_a="her"))
    H=herm(star(z)@Q@z)
    Hd=eigh(H,eigvals_only=True)[::-1]
    nu=int(np.count_nonzero(Kd<0.))
    if np.maximum(Hd[nu:],0.).sum()>radius2:
        return None,"ordinary_multiplier_already_excludes_float"
    Po=H-radius2*np.eye(z.shape[1])
    Pd=eigh(Po,eigvals_only=True)
    if min(abs(Pd))<1.e-11 or np.count_nonzero(Pd>0.)!=nu:
        return None,"not_saturated_float"
    Qz=Q@z
    Qc=herm(Q-Qz@solve(Po,star(Qz),assume_a="her"))
    Ry=Qc[m:,m:]
    if eigh(Ry,eigvals_only=True)[0]<=1.e-10:
        return None,"receiver_not_positive_float"
    cross=Qc[m:,8:m]
    centers=-solve(Ry,cross,assume_a="pos")
    beta=np.real(np.diag(-star(cross)@centers-Qc[8:m,8:m]))
    if np.min(beta)<0.:
        return None,"negative_ellipse_beta_float"
    return dict(w=w,opt=opt,centers=centers,beta=beta,Ry_inv=solve(Ry,np.eye(len(Ry)),
                                                                 assume_a="pos")),"eligible_float"


def exact_envelopes(g,s,b,ell,upp,weights,optical,z,radius,v):
    n,m=b.shape
    Z=ball(z)
    r2=arb(radius)**2
    K,Q=schur(system(g,s,b,ell,upp,weights,optical,exact=True),n)
    Ksign=signed(K)
    H=bt(Z)*Q*Z
    Po=H-r2*ball(np.eye(z.shape[1]))
    Psign=signed(Po)
    if Ksign["negative"]!=Psign["positive"]:
        raise ArithmeticError("no strict saturation")
    Qz=Q*Z
    qc=Q-Qz*Po.solve(bt(Qz))
    rec=range(m,m+len(s))
    ry=select(qc,rec,rec)
    rsign=signed(ry)
    if rsign["negative"]!=0:
        raise ArithmeticError("no strictly positive receiver block")
    cross=select(qc,rec,[v])
    c=-ry.solve(cross)
    beta=-(bt(cross)*c)[0,0].real-qc[v,v].real
    if not beta>0:
        raise ArithmeticError("negative or unresolved radius")
    return c,beta,ry,dict(K=Ksign,Po=Psign,Ry=rsign)


def data_for_public_nodes(ids):
    all_a=load(PUBLIC)
    directory=RESULTS/"known_zero_control"
    r=[]
    for k in range(1,25):
        stem=f"pair_{k:02d}"
        candidate=json.loads((directory/f"{stem}.json").read_text())
        rows=json.loads((RESULTS/"raw"/f"{stem}.json").read_text())["rows"]
        for q,s in zip(candidate["nodes"],rows,strict=True):
            if q["C_known_zero"]=="UNRESOLVED":
                r.append((s["pair_id"],s["target_id"],s["condition_id"],k-1,s["box_positions"]))
    seen=set()
    selected=[]
    for info in sorted(r):
        key=info[0]
        if key not in seen:
            selected.append(info)
            seen.add(key)
    return [(info,all_a) for info in selected[:ids]],len(r)


def run_one(info, data):
    pair_id,target_id,condition_id,pair,(target,bg,box_width,slot)=info
    ell_all,upp_all=data["boxes"][pair,bg,box_width,slot]
    nonzero=upp_all>0.
    g=data["G"][np.ix_(nonzero,nonzero)]
    s=data["S"][:,nonzero]
    angles=(2*np.arange(8)+1)*np.pi/8
    virtual=np.exp(2j*np.pi*data["xy"]@np.vstack((np.cos(angles),np.sin(angles))))
    b=np.hstack((data["B"],virtual))[nonzero]
    Z=np.vstack((np.eye(16)[:,:8],data["Y"][pair,target,bg]))
    radius=float(data["radii"][pair,target,bg])
    ell,upp=ell_all[nonzero],upp_all[nonzero]
    attempts={}
    eligible=[]
    for pars in PARAMS:
        w,opt=[f[nonzero] for f in weight_vector(len(nonzero),*pars)]
        match,status=numeric(g,s,b,ell,upp,w,opt,Z,radius)
        attempts[status]=attempts.get(status,0)+1
        if match is not None:
            eligible.append(match)
    winner=None
    for a,bb in combinations(eligible,2):
        for k in range(8):
            direction=bb["centers"][:,k]-a["centers"][:,k]
            norm=np.linalg.norm(direction)
            if norm==0.:
                continue
            direction/=norm
            r1=np.sqrt(a["beta"][k]*np.vdot(direction,a["Ry_inv"]@direction).real)
            r2=np.sqrt(bb["beta"][k]*np.vdot(direction,bb["Ry_inv"]@direction).real)
            gap=norm-r1-r2
            if winner is None or gap>winner[0]:
                winner=gap,direction,k,a,bb
    proof=None
    if winner and winner[0]>0.:
        gap,v,k,a,bb=winner
        try:
            e1=exact_envelopes(g,s,b,ell,upp,a["w"],a["opt"],Z,radius,k+8)
            e2=exact_envelopes(g,s,b,ell,upp,bb["w"],bb["opt"],Z,radius,k+8)
            V=ball(v[:,None])
            stats=[]
            for center,beta,ry,signs in (e1,e2):
                cent=(bt(V)*center)[0,0].real
                rad=(beta*(bt(V)*ry.solve(V))[0,0].real).sqrt()
                stats.append((cent-rad,cent+rad,signs))
            L=stats[1][0]-stats[0][1]
            proof=dict(certified=bool(L>0.),gap_lower=low(L),
                       projected_intervals=[[low(x),high(y)] for x,y,_ in stats],
                       inertia=[r[2] for r in stats])
        except (ArithmeticError,ZeroDivisionError,ValueError) as error:
            proof={"certified":False,"reason":str(error)}
    return dict(pair_id=pair_id,target_id=target_id,condition_id=condition_id,
                position=info[-1],remaining_state_size=int(sum(nonzero)),
                method_trials=attempts,best_gap_float=None if winner is None else float(winner[0]),
                virtual_angle_degrees=None if winner is None else float(np.rad2deg(angles[winner[2]])),
                independent_certificate=proof)


def main():
    t0=perf_counter()
    cohort,pool=data_for_public_nodes(5)
    rows=[]
    for info,a in cohort:
        row=run_one(info,a)
        rows.append(row)
        print(json.dumps(row),flush=True)
    write_json(HERE/"phase_a_followup.json",dict(
        scope="Exploratory mechanism check after primary Phase A freeze; selects five distinct pairs by first publicly ordered unresolved node",
        selection="Uses publicly saved C-known-vacuum UNRESOLVED states, sorted by pair_id, target_id, condition_id; no sealed labels or virtual responses read",
        actual_predecessor_unresolved=pool,source_caps="8 new virtual directions shifted by 22.5 degrees with disjoint stored columns from observed directions",
        cert_class="Kernel-tested Schur envelopes per ordinary C multiplier; only positive proposal triggers Arb replay",
        rows=rows,certified=sum(bool(row["independent_certificate"] and
                                     row["independent_certificate"]["certified"]) for row in rows),
        elapsed_seconds=perf_counter()-t0))


if __name__=="__main__":
    main()
