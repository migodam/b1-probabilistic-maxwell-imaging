"""Localized finite positive caps using the complete exterior/source Schur loss."""
import argparse
import json
import time
import numpy as np
from scipy.linalg import eigh
from threadpoolctl import threadpool_limits

from sector import HERE, PUBLIC, UPPER, Green, adj, blocks, herm, tail_from_blocks
from run import profile
from verify import (acb, acb_mat, adj as badj, arb, ball, ctx, diagonal, fhi, interval_spectrum, lower,
                    norm, row_scaled, source_tail, sub, verify_signs)


def numeric_upper_re(green):
    left=adj(green.left)@green.left
    right=green.right@adj(green.right)
    product=np.sqrt(max(np.abs(left).sum(axis=1))*max(np.abs(right).sum(axis=1)))
    return float(green.spectrum.real.max()+product)


def rigorous_upper_re(green):
    left,right=ball(green.left),ball(green.right)
    def dia_bound(g):
        return max(sum((abs(g[i,j]).upper() for j in range(g.ncols())),arb(0)).upper()
                   for i in range(g.nrows()))
    kl,kr=dia_bound(badj(left)*left),dia_bound(right*badj(right))
    return (max(z.real.upper() for z in interval_spectrum(green).entries())+(kl*kr).sqrt()).upper()


def prep(green,s,b,y,active,weight,optical,pilot):
    region=np.setdiff1d(np.arange(green.n),active)
    p=profile(np.full(len(active),UPPER),s,w=weight,optical=optical)
    e=s[:,active]
    r=s[:,region]
    k,h,c=blocks(green.matrix(active,active),e,b[active],y,p)
    z0=np.vstack((e.conj().T@r-p.mix[:,None]*green.matrix(active,region),-y.conj().T@r))
    z1=np.vstack((-green.matrix(region,active).conj().T,-b[region].conj().T))
    full=(float(max(np.abs(green.spectrum).ravel())+np.linalg.norm(green.left)*np.linalg.norm(green.right))
          if pilot else numeric_upper_re(green))
    terms=(z0@adj(z0),z0@adj(z1)+z1@adj(z0),z1@adj(z1))
    return np.block([[k,-h],[-adj(h),c]]),terms,full


def choose(green,s,b,y,r,active,weight,optical,pilot):
    base,terms,gnorm=prep(green,s,b,y,active,weight,optical,pilot)
    m=len(active)
    choices=[]
    numcap,numbeta=(32,35) if pilot else (24,25) if green.side==32 else (15,14)
    for cap in np.geomspace(1e-8,0.85/gnorm,numcap):
        best=None
        for beta in np.geomspace(2e-3,1e5,numbeta):
            mix=beta*cap/2
            delta=(terms[0]+mix*terms[1]+mix**2*terms[2])/(beta*(1-cap*gnorm))
            az=base-delta
            dist,nu,gap,_=tail_from_blocks(az[:m,:m],-az[:m,m:],az[m:,m:])
            row=dict(cap=float(cap),roi_weight=float(beta),distance=dist,nu=nu,gap=gap,margin=dist-r)
            if best is None or dist>best["distance"]:
                best=row
        choices.append(best)
    return choices


def exact_matrix(green,rows,cols,kernel):
    cells=list(rows)
    other=list(cols)
    g=acb_mat(len(cells),len(other),[
        kernel[(i%green.side-j%green.side)%green.pad][(i//green.side-j//green.side)%green.pad]
        for i in cells for j in other])
    return g + ball(green.left[rows])*ball(green.right[:,cols])


def exact_kernel(green):
    spec=interval_spectrum(green)
    kernel=[acb.dft([spec[i,j] for j in range(green.pad)],inverse=True) for i in range(green.pad)]
    return [acb.dft([kernel[j][i] for j in range(green.pad)],inverse=True) for i in range(green.pad)]


def one_sided(green,s,b,y,r,active,weight,optical,candidate,kernel,pilot):
    m=len(active)
    rr=np.setdiff1d(np.arange(green.n),active)
    se,be,yb=map(ball,(s[:,active],b[active],y))
    bs=[arb(float(w)) for w in weight]
    mix=[acb(w*arb(UPPER),arb(float(t)))/2 for w,t in zip(bs,optical,strict=True)]
    g=exact_matrix(green,active,active,kernel)
    c,w=arb(candidate["cap"]),arb(candidate["roi_weight"])
    me=w*c/2
    if pilot:
        full=max(abs(z).upper() for z in interval_spectrum(green).entries())
        full+=(norm(ball(green.left))*norm(ball(green.right))).upper()
    else:
        full=rigorous_upper_re(green)
    floor_rr=(w*(1-c*full)).lower()
    if not floor_rr>0:
        return dict(proved=False,reason="ROI pivot positivity unresolved")
    total=m+yb.ncols()
    down=acb_mat(total,total)
    for start in range(0,len(rr),192):
        local=rr[start:start+192]
        sr,br=ball(s[:,local]),ball(b[local])
        ge_r=exact_matrix(green,active,local,kernel)
        gr_e=exact_matrix(green,local,active,kernel)
        z0=badj(se)*sr-row_scaled(ge_r,mix)-me*badj(gr_e)
        z1=-badj(yb)*sr-me*badj(br)
        patch=acb_mat(total,len(local),z0.entries()+z1.entries())
        down+=patch*badj(patch)
    down/=floor_rr
    mg=row_scaled(g,mix)
    ii,ss=range(m),range(m,m+yb.ncols())
    k=badj(se)*se+diagonal(bs)-mg-badj(mg)-sub(down,ii,ii)
    l=badj(se)*yb+row_scaled(be,mix)+sub(down,ii,ss)
    h0=badj(yb)*yb-sub(down,ss,ss)
    signs=verify_signs(k)
    nu=signs["negative_at_least"]
    if signs["unresolved"] or nu>=yb.ncols():
        return dict(proved=False,reason="no surviving relaxed source directions",signs=signs)
    h=h0-badj(l)*k.solve(l)
    dist,tails=source_tail(h,nu)
    threshold=arb(r)+arb("1e-8")*norm(yb)
    return dict(proved=bool(dist>threshold),roi_peak_lower=candidate["cap"],
                lower_distance=lower(dist),noise_radius_upper=fhi(threshold),
                lower_margin=lower(dist-threshold),nminus_exterior_envelope=nu,
                roi_pivot_floor_lower=lower(floor_rr),source_eigenvalue_lowers=tails,
                outside_material_interval=[0.,UPPER],green_hermitian_upper=fhi(full),precision_bits=ctx.prec)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case",required=True)
    parser.add_argument("--suffix",default="main")
    parser.add_argument("--region",required=True,action="append")
    parser.add_argument("--pilot-norm",action="store_true")
    args=parser.parse_args()
    ctx.prec,ctx.threads=112,2
    clock=time.perf_counter()
    case=args.case
    sid=int(case.split("_")[0][1:])
    data=np.load(PUBLIC/"cases"/f"{case}.npz")
    op=PUBLIC/"operators"/f'{str(data["operator_id"])}.npz'
    p=np.load(op)
    g=Green.read(op,sid)
    s,b,y,r=p["S"],p["E0"],data["Y"],float(data["radius"])
    arrays=np.load(HERE/"results"/f"{case}_{args.suffix}.npz")
    out=dict(case=case,mask_version=args.suffix,regions={})
    kernel=exact_kernel(g)
    for site in args.region:
        if site+"_active" not in arrays:
            continue
        site_clock=time.perf_counter()
        active,weight,optical=(arrays[site+"_"+k] for k in ("active","weight","optical"))
        proposals=choose(g,s,b,y,r,active,weight,optical,args.pilot_norm)
        tested=[]
        for cand in sorted((a for a in proposals if a["distance"]>r), key=lambda a:-a["cap"]):
            audit=one_sided(g,s,b,y,r,active,weight,optical,cand,kernel,args.pilot_norm)
            tested.append(dict(candidate=cand,verification=audit))
            if audit["proved"]:
                break
        out["regions"][site]=dict(proposals=proposals,external_cells=len(active),
                                  strict_tests=tested,wall_seconds=time.perf_counter()-site_clock)
        print(case,site,tested[-1] if tested else None,flush=True)
    out["wall_seconds"]=time.perf_counter()-clock
    tag="" if args.pilot_norm else "_v2"
    (HERE/"results"/f"{case}_{args.suffix}_refined_cap{tag}.json").write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    with threadpool_limits(1):
        main()
