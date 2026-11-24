"""Finite native response to retaining a few weak positive modes explicitly."""
from time import perf_counter, process_time

import numpy as np
from flint import acb_mat, arb, ctx
from scipy.linalg import eigh
from threadpoolctl import threadpool_limits

from outward import clipped_best, f_norm, rowscale
from replay import HERE, ball, configure, gershgorin, herm, high, load_json, low, mid, save_json, star, times

CASE = "u100_best"
STEPS = (0, 2, 4, 8)
ETAS = (.05, .1, .2, .35, .5, .8, 1.2, 2.)
MATERIAL_SAMPLES = np.linspace(0., 1., 121)


def block(a, left, right):
    return acb_mat(a.nrows(), right-left, [a[i,j] for i in range(a.nrows()) for j in range(left,right)])


def compressed_upper(a):
    _, vec=eigh(herm(mid(a)),check_finite=False)
    v=ball(vec)
    g=star(v)*v
    q=gershgorin(star(v)*a*v,upper=True)
    return (q/gershgorin(g,upper=q<0)).upper()


def numerical_frame_candidates(positive, gp, q, cap):
    p,g=positive[:,q:],gp[:,q:]
    rb=p-cap*g
    midpoint=p-cap*g/2
    rng=np.random.default_rng(802+q)
    seeds=np.column_stack((np.eye(p.shape[1])[:,:7],rng.normal(size=(p.shape[1],3))))
    best=(0.,None,None)
    for coeff in seeds.T:
        for _ in range(8):
            x=np.where(abs(rb@coeff)>abs(p@coeff),cap,0.)
            a=p-x[:,None]*g
            score,rot=eigh(herm(a.conj().T@a),subset_by_index=(p.shape[1]-1,p.shape[1]-1),
                            check_finite=False)
            right=rot[:,-1]
            if score[0]>best[0]:
                best=(float(score[0]),right,x)
            if abs(np.vdot(coeff,right))>.999999:
                break
            coeff=right
    actual=dict(denom=best[0],current=p@best[1],materials=best[2])
    scores=[]
    for eta in ETAS:
        gram=(1+eta)*(midpoint.conj().T@midpoint)+(1+1/eta)*(cap/2)**2*(g.conj().T@g)
        _,b=eigh(herm(gram),subset_by_index=(p.shape[1]-1,p.shape[1]-1),check_finite=False)
        coeff=b[:,-1]
        denom=(np.linalg.norm(midpoint@coeff)+cap/2*np.linalg.norm(g@coeff))**2
        scores.append((denom,eta,coeff))
    value,eta,coef=max(scores,key=lambda x:x[0])
    young=dict(denom=float(value),current=p@coef,proposal_eta=eta)
    return dict(all_box=actual,young=young)


def strict_frame_upper(native, state, vneg, p, kp, gp, pl, kpl, materials=None):
    leak=arb(state["ritz_positive_leak_upper"])
    gmin=arb(state["negative_projected_gram_min_lower"])
    kneg=arb(state["negative_K_norm_upper"])
    gnorm=arb(state["green_norm_upper"])
    def negative_leak(t):
        return ((f_norm(star(vneg)*t)+leak*f_norm(t))/gmin.sqrt()).upper()
    p_neg=negative_leak(p)
    pl_neg=negative_leak(pl) if pl.ncols() else arb(0)
    p_drift=p_neg
    if pl.ncols():
        lower=gershgorin(star(pl)*kpl,upper=False)
        if not lower>0:
            raise ArithmeticError("The projected weak basis is not certified injective in the positive metric.")
        cross=f_norm(star(pl)*kp)+kneg*pl_neg*p_neg
        p_drift=(p_neg+f_norm(pl)*cross/lower).upper()
    energy=((star(p)*kp)[0,0].real+kneg*p_neg**2).upper()
    cap=arb(native.w.cap)
    if materials is None:
        # For every eta>0, Young envelope energy is at least
        # (||[I-uG/2]p_r|| + (u/2)||Gp_r||)^2 (AM-GM).
        left=max(arb(0),(f_norm_lower(p-cap/2*gp)-(1+cap/2*gnorm)*p_drift).lower())
        right=max(arb(0),(f_norm_lower(gp)-gnorm*p_drift).lower())
        denom=((left+cap/2*right)**2).lower()
    else:
        a=p-rowscale(gp,[arb(float(x)) for x in materials])
        denom=(max(arb(0),(f_norm_lower(a)-(1+cap*gnorm)*p_drift).lower())**2).lower()
    if not denom>0:
        raise ArithmeticError("Unable to resolve a positive response of the remaining-frame witness.")
    return dict(weight_upper=high(energy/denom),witness_energy_upper=high(energy),
                witness_frame_norm_squared_lower=low(denom),weak_frame_test_projection_drift_upper=high(p_drift),
                weak_columns_negative_leak_upper=high(pl_neg))


def f_norm_lower(v):
    from replay import frobenius_squared
    return frobenius_squared(v).sqrt().lower()


def augmented_schur(n, gn, p, gp, x, weight):
    if isinstance(x,np.ndarray):
        x=x[:,None]
    a=n-x*gn
    b=p-x*gp
    neg=weight*a.conj().T@a
    cross=weight*b.conj().T@a
    weak=np.eye(p.shape[1])+weight*b.conj().T@b
    l=-np.linalg.solve(weak,cross)
    return herm(neg+cross.conj().T@l)-np.eye(n.shape[1]),l


def graph_candidates(n,gn,p,gp,weight,u,allow_independent=True):
    rows=[]
    for x in u*MATERIAL_SAMPLES:
        mat,l=augmented_schur(n,gn,p,gp,float(x),weight)
        ev,v=eigh(mat,check_finite=False)
        rows.append((int((ev<0).sum()),ev,v,float(x),l))
    def order(c):
        d,ev,*_=c
        return -d,(ev[d-1] if d else np.inf)
    best=min(rows,key=order)
    if not allow_independent:
        return best,best
    d_top=best[0]
    for d in range(max(1,d_top),min(n.shape[1],max(d_top+2,9))+1):
        if d_top==n.shape[1]:
            break
        x=best[3]
        for _ in range(10):
            mat,l=augmented_schur(n,gn,p,gp,x,weight)
            ev,v=eigh(mat,check_finite=False)
            cur=(int((ev<0).sum()),ev,v,x,l)
            if order(cur)<order(best):
                best=cur
            v=v[:,:d]
            a=n@v+p@(l@v)
            b=gn@v+gp@(l@v)
            new=clipped_best(a,b,u)
            if isinstance(x,np.ndarray) and np.linalg.norm(x-new)<1.e-7:
                break
            x=new
    return min(rows,key=order),best


def outward_graph(native,state,vb,gnb,kvv,negative_eig,pl,gpl,kpl,
                  graphs,omega,pl_neg):
    d,ev,v,x,l=graphs
    leak=arb(state["ritz_positive_leak_upper"])
    gnorm=arb(state["green_norm_upper"])
    kneg=arb(state["negative_K_norm_upper"])
    v=v[:,:d]
    coord=ball(v/np.sqrt(-negative_eig)[:,None])
    weak_coord=ball(l@v)
    if not gershgorin(star(coord)*coord,upper=False)>0:
        raise ArithmeticError("A candidate negative graph has unresolved rank.")
    zv,zg=vb*coord,gnb*coord
    h=star(coord)*kvv*coord
    v_drift=leak*f_norm(coord)
    if weak_coord.nrows():
        zv+=pl*weak_coord
        zg+=gpl*weak_coord
        h+=star(weak_coord)*(star(pl)*kpl)*weak_coord
        q_drift=arb(pl_neg)*f_norm(weak_coord)
        v_drift+=q_drift
        kshift=kneg*q_drift**2
    else:
        kshift=arb(0)
    if isinstance(x,np.ndarray):
        residual=zv-rowscale(zg,[arb(float(a)) for a in x])
        magnitude=max(x)
    else:
        residual=zv-arb(x)*zg
        magnitude=x
    row_drift=((1+arb(float(magnitude))*gnorm)*v_drift).upper()
    fshift=arb(omega)*(2*f_norm(residual)*row_drift+row_drift**2)
    h+=arb(omega)*(star(residual)*residual)
    upper=(compressed_upper(h)+kshift+fshift).upper()
    return dict(negative_graph_dimension=d,augmented_extra_negative_count_lower=d if upper<0 else None,
                maximum_compensated_dimension=negative_eig.size-d if upper<0 else None,
                largest_graph_eigenvalue_upper=high(upper),negative_graph_operator_drift_upper=high(kshift+fshift),
                numerical_augmented_schur_eigenvalues=ev.tolist(),
                material="independent native cells" if isinstance(x,np.ndarray) else float(x))


def joint_rows(q,n,gn,p,gp,gram,pixels,cut):
    result=[]
    for checked in cut["sets"]:
        k=checked["rows"]
        x=checked["material"]
        i=pixels[:k]
        c=np.vstack(((1-x)*np.eye(k),x*np.eye(k)))
        dims=np.r_[np.arange(k), len(pixels)+np.arange(k)]
        gamma=herm(c.conj().T@gram[np.ix_(dims,dims)]@c)
        weak=(p-x*gp)[i]
        neg=(n-x*gn)[i]
        rem=herm(gamma-weak@weak.conj().T)
        rem_lam=eigh(rem,eigvals_only=True,check_finite=False)
        weight=1/rem_lam[-1]
        cross=weight*weak.conj().T@neg
        eff=weight*neg.conj().T@neg-cross.conj().T@np.linalg.solve(
            np.eye(q)+weight*weak.conj().T@weak,cross)
        result.append(dict(rows=k,material=x,remaining_frame_uniform_weight_upper_at_point_numeric=weight,
            compressed_negative_response_norm_numeric=float(eigh(herm(eff),eigvals_only=True)[-1]),
            reused_full_signed_row_gram_lower=checked["homogeneous_min_low"],
            strict_any_split_rank_obstruction=bool(checked["strictly_positive"])))
    return result


def source_controls(current_case,optimistic,trial,cap,comp,weight,mu):
    rows=[]
    for meta,old,j,e,h in zip(current_case["rows"],optimistic["rows"],trial["center"],
                                trial["center_field"],trial["hs"],strict=True):
        m=current_case["verified_state"]["exact_nu"]-comp
        e_h,q_h=eigh(herm(h),check_finite=False)
        q=q_h[:,:len(e_h)-m]
        jj,ee=j@q,e@q
        x=clipped_best(jj,ee,cap)
        power=mu*weight*np.linalg.norm(jj-x[:,None]*ee)**2
        rows.append(dict(id=meta["id"],
            optimistic_psd_center_ratio_upper=old["optimistic_center_free_ratio_upper_by_c"][comp],
            signed_center_compression_upper_numerical=float(np.sqrt(
                np.maximum(e_h[:len(e_h)-m],0.).sum()+power)/meta["padded_radius"])))
    return rows


def main():
    start=perf_counter(),process_time()
    geo,native,_,_=configure(CASE)
    prior=load_json(HERE/f"{CASE}_strict.json")
    cur=load_json(HERE/f"{CASE}.json")
    joint=load_json(HERE/f"{CASE}_joint_cut.json")
    opt=load_json(HERE/"OUTWARD_OPTIMISTIC_CONTROL.json")[CASE]
    np_input=lambda name:np.load(HERE/name,allow_pickle=False)
    trial=np_input(f"{CASE}.npz")
    raw=np_input(f"{CASE}_conditioned.npz")
    positive=np_input(f"{CASE}_all_pixel_probe.npz")
    gp=positive["positive_trial_green"]
    p=positive["positive_trial"]
    eig=raw["negative_eig"]
    v=raw["negative_vec"]
    n=v/np.sqrt(-eig)
    gn=raw["negative_gn"]
    cap=native.w.cap
    qmax=max(STEPS)
    print("Frame witness proposals",flush=True)
    tests={q:numerical_frame_candidates(p,gp,q,cap) for q in STEPS}
    tag=[(q,k) for q in STEPS for k in tests[q]]
    proposed=np.column_stack([v,p[:,:qmax]]+[tests[q][k]["current"] for q,k in tag])
    print("Outward native actions",flush=True)
    abc=ball(proposed)
    kabc=native.exact_action(abc)
    gabc=geo.exact_green(abc)
    vb,gnb=block(abc,0,len(eig)),block(gabc,0,len(eig))
    kvv=star(vb)*block(kabc,0,len(eig))
    gram=np_input(f"{CASE}_strict.npz")["endpoint_gram_midpoint"]
    rows=[]
    archive={}
    for q in STEPS:
        t0=perf_counter(),process_time()
        pl=block(abc,len(eig),len(eig)+q)
        kpl=block(kabc,len(eig),len(eig)+q)
        gpl=block(gabc,len(eig),len(eig)+q)
        r=dict(weak_positive_modes=q,
               fixed_row_groups=joint_rows(q,n,gn,p[:,:q],gp[:,:q],gram,raw["pixel_index"],joint))
        for case,c in tests[q].items():
            ii=len(eig)+qmax+tag.index((q,case))
            proof=strict_frame_upper(native,prior["projector"],vb,block(abc,ii,ii+1),
                block(kabc,ii,ii+1),block(gabc,ii,ii+1),pl,kpl,c.get("materials"))
            w=proof["weight_upper"]
            uniform,best=graph_candidates(n,gn,p[:,:q],gp[:,:q],w,cap)
            g=outward_graph(native,prior["projector"],vb,gnb,kvv,eig,pl,gpl,kpl,best,w,
                            proof["weak_columns_negative_leak_upper"])
            c_bound=g["maximum_compensated_dimension"]
            r[case]=dict(frame=proof,**g,uniform_only_negative_graph_dimension=uniform[0],
                         uniform_only_material=uniform[3],numerical_frame_capacity_preview=c["denom"])
            if c_bound is not None:
                r[case]["sources"]=source_controls(cur,opt,trial,cap,c_bound,w,native.w.scale)
            archive[f"{q}_{case}_negative_test_material"]=best[3]
            archive[f"{q}_{case}_negative_graph_basis"]=best[2][:,:best[0]]
            archive[f"{q}_{case}_weak_graph_basis"]=best[4]@best[2][:,:best[0]]
            if "materials" in c:
                archive[f"{q}_{case}_positive_test_material"]=c["materials"]
            archive[f"{q}_{case}_positive_witness"]=c["current"]
        rows.append(dict(**r,**times(*t0)))
        print(q, [(c,r[c]["frame"]["weight_upper"], r[c]["maximum_compensated_dimension"],
                   r[c]["largest_graph_eigenvalue_upper"]) for c in tests[q]],flush=True)
    save_json(HERE/"deflated.json",dict(case=CASE,profile=native.w.as_dict(),results=rows,
       hypotheses=dict(all_box="any spatially uniform nonnegative weight satisfying the exact remaining frame for all box points",
                       young="any spatially uniform nonnegative weight certified by the stated Young envelope, any eta>0, fixed weak basis",
                       weak_basis="exact positive projections of the represented weakest K_scaled Ritz directions, remaining space orthogonal in positive energy",
                       source_envelope="optimistic strict for a PSD center debit; signed center values are numerical upper proposals only and require pure-source containment"),
       source_center_optimization_performed=False,
       operations=dict(geo.cost),total=times(*start)))
    np.savez_compressed(HERE/"deflated.npz",**archive)


if __name__=="__main__":
    ctx.prec,ctx.threads=112,2
    with threadpool_limits(1):
        main()
