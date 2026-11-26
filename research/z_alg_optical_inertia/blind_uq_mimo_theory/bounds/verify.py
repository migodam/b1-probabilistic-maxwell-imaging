"""One-sided Arb replay of native loaded Green operators and residual thresholds."""
import argparse
import json
import sys
import time

from flint import acb, acb_mat, arb, ctx
import numpy as np
from scipy.linalg import eigh
from threadpoolctl import threadpool_limits

from sector import HERE, PUBLIC, UPPER, Green, herm

historical = HERE.parents[1]/"ip_closed_loop"/"reduction"
sys.path.insert(0, str(historical if historical.exists() else HERE.parent/"vendor"/"inertia"))
from interval_backend import adj, ball, diagonal, floor, lower, row_scaled, sub, verify_signs


def norm(a):
    # Abs upper is nonnegative even when a tiny interval encloses exact zero.
    squares = sum((abs(z).upper()**2 for z in a.entries()), arb(0))
    return squares.sqrt().upper() if squares > 0 else arb(0) if squares.is_zero() else arb.nan()


def fhi(v):
    return float(np.nextafter(float(v.upper()), np.inf))


def hstack(arrays):
    nrow = arrays[0].nrows()
    entries = [a[i,j] for i in range(nrow) for a in arrays for j in range(a.ncols())]
    return acb_mat(nrow, sum(a.ncols() for a in arrays), entries)


def fft_apply(green, values, spectrum):
    """Independent 112-bit replay without any dense native state operator."""
    n, pad = green.side, green.pad
    unset = acb(0)
    cols = []
    for c in range(values.ncols()):
        rows = [acb.dft([values[n*i+j,c] for j in range(n)]+[unset]*(pad-n)) for i in range(n)]
        filtered = []
        for j in range(pad):
            transformed = acb.dft([rows[i][j] for i in range(n)]+[unset]*(pad-n))
            filtered.append(acb.dft([transformed[i]*spectrum[i,j] for i in range(pad)],inverse=True)[:n])
        output = [acb.dft([filtered[j][i] for j in range(pad)],inverse=True)[:n] for i in range(n)]
        cols.extend(z for r in output for z in r)
    return acb_mat(values.ncols(),green.n,cols).transpose()


def interval_spectrum(green):
    return ball(green.spectrum)


def verify_dark(green, s, current, active, upper, weight, optical, spectrum=None):
    clock = time.perf_counter()
    active = np.asarray(active)
    s_act = ball(s[:,active])
    x = ball(current[active])
    correction = adj(s_act)*((s_act*adj(s_act)).solve(s_act*x))
    dark = x-correction
    apart = np.full(green.n,-1,int)
    apart[active] = np.arange(len(active))
    x = acb_mat(green.n,dark.ncols(),[dark[apart[i],j] if apart[i] >= 0 else acb(0)
                 for i in range(green.n) for j in range(dark.ncols())])
    g = fft_apply(green, x, interval_spectrum(green) if spectrum is None else spectrum)
    g += ball(green.left)*(ball(green.right)*x)
    m = [acb(arb(float(w))*arb(float(u)), arb(float(l)))/2
         for w,u,l in zip(weight,upper,optical,strict=True)]
    w = [arb(float(v)) for v in weight]
    mixed = adj(x)*row_scaled(g,m)
    matrix = mixed+adj(mixed)-adj(x)*row_scaled(x,w)
    counts = verify_signs(matrix)
    nneg = counts["positive_at_least"]
    return dict(n_minus_at_least=nneg, inertia_tail_surviving_columns_at_most=max(0,len(s)-nneg),
                negative_margin_lower=counts["floors"][0], receiver_null_by_exact_right_inverse=True,
                receiver_gram_solved_in_interval=True, correction_norm_upper=fhi(norm(correction)),
                orthogonality_bound_not_required=True, precision_bits=ctx.prec,
                represented_operator="all native FFT frequencies plus exact binary64 low-rank port loading",
                wall_seconds=time.perf_counter()-clock)


def source_tail(h, nu):
    size = h.nrows()
    mid = np.array([complex(z) for z in h.entries()]).reshape(size,size)
    _, v = eigh(herm(mid), check_finite=False)
    limits = []
    for j in range(nu,size):
        basis = ball(v[:,-(j+1):])
        proj = floor(adj(basis)*h*basis)
        g = adj(basis)*basis
        upper = max((g[i,i].real + sum((abs(g[i,k]) for k in range(j+1) if k != i),arb(0))).upper()
                    for i in range(j+1))
        if not proj > 0:
            break
        limits.append((proj/upper).lower())
    total = sum(limits,arb(0))
    return total.sqrt().lower() if total > 0 else arb(0), [lower(v) for v in limits]


def exact_support(green,active):
    kern = acb_mat(green.pad,green.pad,[z for row in
                   [acb.dft([interval_spectrum_value for interval_spectrum_value in r],inverse=True)
                    for r in green.spectrum.tolist()] for z in row])
    ker_rows = [acb.dft([kern[j,i] for j in range(green.pad)],inverse=True) for i in range(green.pad)]
    cells=list(active)
    g=acb_mat(len(cells),len(cells),[
        ker_rows[(i%green.side-j%green.side)%green.pad][(i//green.side-j//green.side)%green.pad]
        for i in cells for j in cells])
    return g + ball(green.left[active])*ball(green.right[:,active])


def verify_support(green,s,b,y,radius,active,weight,optical):
    clock=time.perf_counter()
    g=exact_support(green,active)
    w=[arb(float(x)) for x in weight]
    m=[acb(v*arb(UPPER),arb(float(t)))/2 for v,t in zip(w,optical,strict=True)]
    sb,bb,yb=map(ball,(s[:,active],b[active],y))
    mg=row_scaled(g,m)
    k=adj(sb)*sb+diagonal(w)-mg-adj(mg)
    rhs=adj(sb)*yb+row_scaled(bb,m)
    signs=verify_signs(k)
    nu=signs["negative_at_least"]
    if signs["unresolved"] or nu>=y.shape[1]:
        return dict(floor_signs=signs, proved=False, wall_seconds=time.perf_counter()-clock)
    sol=k.solve(rhs)
    h=adj(yb)*yb-adj(rhs)*sol
    distance, tails=source_tail(h,nu)
    threshold=arb(float(radius))+arb("1e-8")*norm(yb)
    out=dict(floor_signs=signs, nu=nu, distance_lower=lower(distance),
             padded_noise_radius_upper=fhi(threshold), safe_margin_lower=lower(distance-threshold),
             proved=bool(distance>threshold), eigenvalue_lowers=tails, precision_bits=ctx.prec)
    if out["proved"]:
        out["strict_positive_cap"] = sector_continuation(
             green,k,sol,s,b,y,active,weight,optical,signs,tails,threshold)
        out["resolvent_positive_cap_control"] = positive_continuation(green,g,s,b,active,distance-threshold)
    out["wall_seconds"]=time.perf_counter()-clock
    return out


def sector_continuation(green,k,sol,s,b,y,active,weight,optical,signs,tails,threshold):
    kpos=signs["positive_at_least"]
    center=np.array([complex(z) for z in k.entries()]).reshape(k.nrows(),k.ncols())
    _,v=eigh(herm(center),driver="evd",check_finite=False)
    basis=ball(v[:,-kpos:])
    # The numerical vectors are exact binary64 constants. Their Gram upper
    # converts positivity of the compressed form to a Euclidean gap.
    gram=adj(basis)*basis-diagonal([arb(1)]*kpos)
    kval=arb(signs["floors"][0])/(1+norm(gram).upper()) if kpos else arb(1)
    source=[arb(float(a)) for a in tails]
    fraction=((1+threshold**2/sum(source,arb(0)))/2).upper()
    if not fraction < 1:
        return dict(proved=False,reason="source thresholds unresolved")
    thresholds=[fraction*t for t in source]
    hgap=min((t-z).lower() for t,z in zip(source,thresholds,strict=True))
    eps=min(kval.lower(),hgap)/(1+norm(sol).upper())**2
    full=max(abs(v).upper() for v in interval_spectrum(green).entries())
    full+=(norm(ball(green.left))*norm(ball(green.right))).upper()
    local=np.setdiff1d(np.arange(green.n),active)
    nr=norm(ball(s[:,local])).upper()
    ne=norm(ball(s[:,active])).upper()
    yr=norm(ball(y)).upper()
    br=norm(ball(b[local])).upper()
    largest_mix=max(abs(acb(arb(float(w))*arb(UPPER),arb(float(t)))/2).upper()
                    for w,t in zip(weight,optical,strict=True))
    a=(nr*(ne+yr)+largest_mix*full).upper()
    bnorm=((full+br)/2).upper()
    beta=(32*a**2/eps).upper()
    cap=min((1/(2*full)).lower(),(a/(beta*bnorm)).lower())
    loss=(2*(a+beta*cap*bnorm)**2/beta).upper()
    if not (eps>loss and cap>0 and sum(thresholds,arb(0))>threshold**2):
        return dict(proved=False,reason="uniform cap extension unresolved",eps=lower(eps),
                    loss=fhi(loss),cap=lower(cap),hgap=lower(hgap),
                    threshold_squared_upper=fhi(threshold**2),
                    source_threshold_total_lower=lower(sum(thresholds,arb(0))),
                    a=fhi(a),b=fhi(bnorm),beta=fhi(beta))
    return dict(proved=True,roi_peak_lower=lower(cap),exterior_upper=UPPER,
                roi_sector_weight=fhi(beta),border_positive_gap_lower=lower(eps),
                border_loss_upper=fhi(loss),state_positive_gap_lower=lower(kval),
                source_thresholds_lower=[lower(t) for t in thresholds],
                strict_noise_ball_margin_lower=lower(sum(thresholds,arb(0)).sqrt()-threshold),
                global_green_norm_upper=fhi(full), frozen_same_noise_radius=True)


def positive_continuation(green,g,s,b,active,margin):
    """Uniform exterior resolvent, full-space coupling majorants, and positive ROI cap."""
    n=len(active)
    m0=green.matrix(active,active)
    base=np.eye(n)/UPPER-herm(m0)
    wave=(m0-m0.conj().T)/(2j)
    scored=[]
    for t in (0.,0.125,0.25,0.5,1.,2.,4.,8.,16.,32.,64.):
        gap=eigh(base+t*wave,subset_by_index=(0,0),check_finite=False,eigvals_only=True)[0]
        scored.append((float(gap/np.hypot(1.,t)),float(gap),t))
    _,proposal,ratio=max(scored)
    if proposal<=0:
        return dict(proved=False,reason="no rotated accretivity gap on exterior with tested ratios",
                    best_numeric_gap=proposal, ratio=ratio)
    gamma=arb(proposal)/2
    mult=acb(1,ratio)
    up=arb(UPPER)
    rotated=mult*g
    coercivity=diagonal([1/up-gamma]*n)-(rotated+adj(rotated))/2
    sign=verify_signs(coercivity)
    if sign["positive_at_least"]!=n:
        return dict(proved=False,reason="exterior positive gap not verified", verified=sign)
    inverse=(abs(mult)/gamma).upper()
    spectrum=interval_spectrum(green)
    full=max(abs(v).upper() for v in spectrum.entries())
    full+=(norm(ball(green.left))*norm(ball(green.right))).upper()
    local=np.setdiff1d(np.arange(green.n),active)
    b0=norm(ball(b[local])).upper()+full*inverse*norm(ball(b[active])).upper()
    feedback=full+full**2*inverse
    response=norm(ball(s[:,local])).upper()+norm(ball(s[:,active])).upper()*inverse*full
    # Half of the threshold root keeps the final margin strictly positive.
    cap=(margin/(2*(b0*response+margin*feedback))).lower()
    loss=(b0*response*cap/(1-cap*feedback)).upper()
    if not (cap>0 and margin-loss>0):
        return dict(proved=False, reason="unresolved finite continuation margin")
    return dict(proved=True, roi_peak_lower=lower(cap), exterior_upper=UPPER,
                rotated_gap_lower=lower(gamma), ratio=ratio, resolvent_norm_upper=fhi(inverse),
                full_green_operator_norm_upper=fhi(full), trace_current_affine=fhi(b0),
                feedback_upper=fhi(feedback), receiver_current_upper=fhi(response),
                worst_added_residual_upper=fhi(loss), safety_margin_lower=lower(margin-loss),
                frozen_same_noise_radius=True)


class ExactUniform:
    def __init__(self,green,s,b,y):
        self.green=green
        self.y=ball(y)
        self.b=ball(b)
        self.s=ball(s)
        self.lys=hstack([ball(green.left),ball(green.right.conj().T),adj(self.s),self.b])
        self.other_g=adj(self.lys)*self.lys
        self.roots=[acb(arb(2)*k/green.pad).exp_pi_i() for k in range(green.pad)]
        n,p=green.side,green.pad
        self.dirichlet=[sum((self.roots[k*j%p] for j in range(n)),acb(0)) for k in range(p)]
        self.spectrum=interval_spectrum(green)
        self.cached={}

    def gram(self, q):
        ids=tuple(q)
        if ids in self.cached:
            return self.cached[ids]
        g=self.green
        n,p=g.side,g.pad
        needed=sorted(set(y for y,_ in ids)|set(x for _,x in ids))
        index={a:i for i,a in enumerate(needed)}
        f=acb_mat(len(needed),n,[self.roots[-k*j%p] for k in needed for j in range(n)])
        columns=[]
        for j in range(self.lys.ncols()):
            a=acb_mat(n,n,[self.lys[i,j] for i in range(g.n)])
            ft=f*a*f.transpose()/p
            columns.extend(ft[index[y],index[x]] for y,x in ids)
        a=acb_mat(self.lys.ncols(),len(ids),columns).transpose()
        dirichlet=self.dirichlet
        top=acb_mat(len(ids),len(ids),[dirichlet[(yv-yu)%p]*dirichlet[(xv-xu)%p]/p**2
                                     for yu,xu in ids for yv,xv in ids])
        right=adj(a)
        m=len(ids)+self.other_g.nrows()
        self.cached[ids]=acb_mat(m,m,[z for i in range(m) for z in
                             ([top[i,j] for j in range(len(ids))]+[a[i,j] for j in range(a.ncols())]
                              if i<len(ids) else [right[i-len(ids),j] for j in range(len(ids))]
                              +[self.other_g[i-len(ids),j] for j in range(self.other_g.ncols())])])
        return self.cached[ids]

    def proof(self,upper,ratio,weight,radius):
        clock=time.perf_counter()
        m,w,t=arb(upper),arb(weight),arb(ratio)
        alpha=w/2
        keep,co=[],[]
        omitted=arb(1000)
        g=self.green
        for i in range(g.pad):
            for j in range(g.pad):
                z=self.spectrum[i,j]
                d=w*(1-m*z.real+t*z.imag)-alpha
                if d<0:
                    keep.append((i,j))
                    co.append(d)
                elif d>=0:
                    omitted=min(omitted,d.lower())
                else:
                    raise ArithmeticError("Unresolved Fourier envelope sign")
        mix=acb(w*m,w*t)/2
        ff=len(keep)
        rank=g.left.shape[1]
        rec,ns=self.y.nrows(),self.y.ncols()
        size=ff+2*rank+rec
        gram=self.gram(keep)
        short=range(size)
        M=sub(gram,short,short)/alpha
        for j,c in enumerate(co):
            M[j,j]+=1/c
        for j in range(rank):
            M[ff+j,ff+rank+j]-=1/mix.conjugate()
            M[ff+rank+j,ff+j]-=1/mix
        for j in range(size-rec,size):
            M[j,j]+=1
        signs=verify_signs(M)
        if signs["unresolved"]:
            return dict(proved=False,floor_signs=signs,wall_seconds=time.perf_counter()-clock)
        nu=signs["positive_at_least"]-rank-rec
        if not 0 <= nu < ns:
            return dict(proved=False,surrogate_nu=nu,floor_signs=signs,wall_seconds=time.perf_counter()-clock)
        factors=acb_mat(rec+ns,ns,[self.y[i,j] if i<rec else mix if j==i-rec else acb(0)
                                  for i in range(rec+ns) for j in range(ns)])
        last=range(size-rec,size+ns)
        a=sub(gram,last,last)
        F=sub(gram,short,last)*factors
        schur=adj(self.y)*self.y-adj(factors)*a*factors/alpha+adj(F)*M.solve(F)/alpha**2
        dist,tails=source_tail(schur,nu)
        threshold=arb(radius)+arb("1e-8")*norm(self.y)
        return dict(distance_lower=lower(dist),padded_noise_radius_upper=fhi(threshold),
                    safe_margin_lower=lower(dist-threshold),proved=bool(dist>threshold),
                    upper=upper,ratio=ratio,weight=weight,surrogate_nu=nu,
                    retained_fourier_modes=len(keep),reduced_dimension=size,
                    omitted_symbol_slack_lower=lower(omitted),floor_signs=signs,
                    uncertainty_scope="exact binary64 stored native Fourier symbol, port loading, and data",
                    eigenvalue_lowers=tails,precision_bits=ctx.prec,wall_seconds=time.perf_counter()-clock)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case",action="append",required=True)
    parser.add_argument("--suffix",default="main")
    parser.add_argument("--dark",action="store_true")
    parser.add_argument("--support",action="store_true")
    parser.add_argument("--support-key",action="append")
    parser.add_argument("--global-cap",action="store_true")
    parser.add_argument("--dark-key",action="append")
    parser.add_argument("--ratio-family",action="store_true")
    args=parser.parse_args()
    ctx.prec,ctx.threads=112,2
    for case in args.case:
        sid=int(case.split("_")[0][1:])
        numeric_path=HERE/"results"/f"{case}_{args.suffix}"
        data=json.loads(numeric_path.with_suffix(".json").read_text())
        obs=np.load(PUBLIC/"cases"/f"{case}.npz")
        ids=data["operator"]
        operator=PUBLIC/"operators"/f"{ids}.npz"
        g=Green.read(operator,sid)
        p=np.load(operator)
        s,b,y,r=p["S"],p["E0"],obs["Y"],float(obs["radius"])
        sections=data["sections"]
        out=dict(case=case,operator=ids,source_radius=r,section={})
        if args.global_cap:
            exact=ExactUniform(g,s,b,y)
            caps=sections.get("global_cap",{})
            for cap,candidates in sorted(caps.items(),key=lambda a: -float(a[0])):
                cand=candidates[0]
                if cand["distance"] <= r:
                    continue
                proof=exact.proof(float(cap),cand["optical"],cand["weight"],r)
                out["section"]["global_cap"]=proof
                if proof["proved"]:
                    break
        store=np.load(numeric_path.with_suffix(".npz"))
        if args.dark:
            keys=["global_prior"]
            keys.extend(k for k in sections.get("roi",{}) if k.startswith(("roi","background")))
            for key in keys:
                if args.dark_key and key not in args.dark_key:
                    continue
                active=np.arange(g.n) if key=="global_prior" else np.flatnonzero(~store[key+"_mask"].ravel())
                base=sections["global_at_prior_cap_dark"] if key=="global_prior" else sections["roi"][key]["obstruction"]
                names={"real_uniform", base["profile_with_smallest_dark_count"]["method"]}
                if key=="roi_01" and key+"_optimized_pixel_currents" in store:
                    names.add("optimized_pixel")
                out["section"][key] = {}
                for name in names:
                    pref=key+"_"+name
                    if pref+"_currents" not in store:
                        continue
                    out["section"][key][name]=verify_dark(g,s,store[pref+"_currents"],active,
                                 store[pref+"_upper"],store[pref+"_weight"],store[pref+"_optical"])
                if args.ratio_family and key in ("global_prior","roi_01","roi_union"):
                    family_pref=key+"_optical_uniform_256"
                    if family_pref+"_currents" in store:
                        frame=store[family_pref+"_currents"]
                    else:
                        from sector import dark_trial, dark_compression
                        from run import profile
                        basis,scattered=dark_trial(g,s,active,rank=104)
                        orth=profile(np.full(g.n,UPPER),s,optical=256)
                        comp=dark_compression(basis,scattered,s,orth)
                        _,vec=eigh(comp,subset_by_index=(0,min(len(comp),len(s))-1),check_finite=False)
                        frame=basis@vec
                    ends=[]
                    for t in (0.,256.):
                        ends.append(verify_dark(g,s,frame,active,np.full(g.n,UPPER),
                                               np.ones(g.n),np.full(g.n,t)))
                    out["section"][key]["continuous_exterior_ratio_family"] = dict(
                        ratio_interval=[0.,256.],same_exact_null_frame=True,endpoint_proofs=ends,
                        all_exterior_positive_sector_scales=True,all_psd_receiver_weights=True,
                        arbitrary_inside_roi_sector_and_optical_weights=True,
                        s_negatives_verified_for_every_ratio=all(a["n_minus_at_least"]>=len(s) for a in ends))
        if args.support:
            for key in sections.get("roi",{}):
                if args.support_key and key not in args.support_key:
                    continue
                if key+"_active" not in store:
                    continue
                out["section"][key+"_zero"]=verify_support(g,s,b,y,r,store[key+"_active"],
                                                            store[key+"_weight"],store[key+"_optical"])
        output=HERE/"results"/f"{case}_{args.suffix}_verified.json"
        output.write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
        print(case, {k:v.get("proved","dark") for k,v in out["section"].items()},flush=True)


if __name__=="__main__":
    with threadpool_limits(1):
        main()
