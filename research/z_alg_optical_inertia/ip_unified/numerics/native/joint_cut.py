"""Uniform-material signed Gram obstruction for the frozen nested row groups."""
import argparse
from time import perf_counter,process_time

from flint import acb_mat,arb,ctx
import numpy as np
from scipy.linalg import eigh
from threadpoolctl import threadpool_limits

from outward import f_norm,interval_endpoint_rows
from replay import CASES,HERE,ball,configure,gershgorin,herm,high,load_json,low,mid,save_json,star,times


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("case",choices=CASES)
    args=p.parse_args()
    start=perf_counter(),process_time()
    case=args.case
    strict=load_json(HERE/f"{case}_strict.json")
    checked=np.load(HERE/f"{case}_strict.npz",allow_pickle=False)
    source=HERE/f"{case}_conditioned.npz"
    if not source.exists():
        source=HERE/f"{case}.npz"
    raw=np.load(source,allow_pickle=False)
    geo,native,_,_=configure(case)
    v,eig,pixels=[raw[k] for k in ("negative_vec","negative_eig","pixel_index")]
    flat_gram=checked["endpoint_gram_midpoint"]
    vb=ball(v)
    ident=ball(np.eye(len(eig)))
    gv=star(vb)*vb
    gram_error=f_norm(gv-ident)
    gmax=gershgorin(gv,upper=True)
    projector=strict["projector"]
    leak=arb(projector["ritz_positive_leak_upper"])
    resid=arb(projector["ritz_residual_F_upper"])
    gneg=arb(strict["verified_state"]["negative_gap_lower"])
    dmin=arb(float(-max(eig)))
    neg_inv_error=((gram_error+leak**2+resid*gmax.sqrt()/dmin)/gneg+
                   2*gmax.sqrt()*leak/dmin).upper()
    rhs=interval_endpoint_rows(native,pixels)
    vv=star(rhs)*vb
    gg=ball(flat_gram)
    signed_data=[]
    scans=[]
    levels=sorted({8,2*len(eig),4*len(eig)})
    for k in levels:
        i=np.r_[np.arange(k),len(pixels)+np.arange(k)]
        proj=mid(vv)[i]
        gn=herm(flat_gram[np.ix_(i,i)]-(proj/(-eig))@proj.conj().T)
        grid=np.linspace(0.,1.,121)
        minimum=[]
        for fraction in grid:
            c=np.vstack(((1-fraction)*np.eye(k),fraction*np.eye(k)))
            minimum.append(eigh(herm(c.conj().T@gn@c),subset_by_index=(0,0),eigvals_only=True)[0])
        f=float(grid[int(np.argmax(minimum))])
        fraction=arb(f)
        rr=acb_mat(k,len(eig),[(1-fraction)*vv[a,l]+fraction*vv[a+len(pixels),l]
                   for a in range(k) for l in range(len(eig))])
        rr_inv=acb_mat(k,len(eig),[rr[a,l]/arb(float(-eig[l]))
                       for a in range(k) for l in range(len(eig))])
        gram_at=acb_mat(k,k,[
            (1-fraction)**2*gg[a,b]+fraction**2*gg[a+len(pixels),b+len(pixels)]
            +fraction*(1-fraction)*(gg[a,b+len(pixels)]+gg[a+len(pixels),b])
            for a in range(k) for b in range(k)])
        exact_trial=gram_at-rr_inv*star(rr)
        mid_trial=herm(mid(exact_trial))
        eigen,rotation=eigh(mid_trial)
        rot=ball(rotation)
        mlo=gershgorin(star(rot)*exact_trial*rot,upper=False)
        unit=gershgorin(star(rot)*rot,upper=mlo>0)
        if not unit>0:
            raise ArithmeticError("Signed Gram Ritz basis was singular.")
        g_error=arb(strict["endpoint_gram"]["lower_scalar_shift"])*((1-fraction)**2+fraction**2)
        action=1+arb(native.w.cap)*fraction*arb(projector["green_norm_upper"])
        exact_lower=(mlo/unit-g_error-action**2*neg_inv_error).lower()
        signed_data.append(dict(rows=k,fraction=f,material=float(f*native.w.cap),
            homogeneous_min_low=low(exact_lower),strictly_positive=bool(exact_lower>0),
            numerical_homogeneous_min=float(eigen[0]),
            scope="This fixed nested pixel set; all nonnegative frame weights (and all PSD block weights)."))
        scans.append(minimum)
        print(case,signed_data[-1],flush=True)
    out=HERE/f"{case}_joint_cut"
    save_json(out.with_suffix(".json"),dict(case=case,negative_inverse_perturbation_upper=high(neg_inv_error),
              sets=signed_data,operations=dict(geo.cost),**times(*start)))
    np.savez_compressed(out.with_suffix(".npz"),pixels=pixels,material_fraction=grid,
                        signed_gram_min_numeric=np.array(scans))


if __name__=="__main__":
    ctx.prec,ctx.threads=112,2
    with threadpool_limits(1):
        main()
