"""Padded FFT fixed-trace sector certificates for a moving local class.

The padded unknown current is an outer relaxation. Every physical current
extends by zero and satisfies the padded background sector.
"""
import argparse
from collections import defaultdict
from fractions import Fraction
import json
import os

for variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[variable] = "1"

from flint import acb, arb, ctx
import numpy as np
from scipy.linalg import eigh
from threadpoolctl import threadpool_limits

from certified_map import ARITH, fnorm, local_models
from common import FREEZE, HERE, Geometry, clock, elapsed, geom_path, save_json
from operator_bounds import ball, inverse_dft_table


def hermitian(x):
    return (x + x.conj().swapaxes(-1, -2)) / 2


def lag_tables(g, b):
    cache = HERE / "results" / f"sector_arb_lags_{g.side}_{b:g}.npz"
    if cache.exists():
        with np.load(cache, allow_pickle=False) as a:
            return a["blocks"], float(a["entry_error"])
    started = clock()
    ctx.prec = 128
    atom = [[ball(a) for a in row] for row in g.symbol]
    denominator = [[1 / (1 - arb(b)*z.real) for z in row] for row in atom]
    lookup = [
        [[acb(d) for d in row] for row in denominator],
        [[d*z.conjugate() for d, z in zip(dr, ar)] for dr, ar in zip(denominator, atom)],
        [[d*(z.real**2+z.imag**2) for d,z in zip(dr, ar)] for dr, ar in zip(denominator, atom)],
    ]
    kernel = []
    err = 0.
    for table in lookup:
        values = []
        for row in inverse_dft_table(table):
            vals = [complex(z) for z in row]
            radius = max(float(abs(z-ball(v)).upper()) for z, v in zip(row, vals))
            err = max(err, radius)
            values.append(vals)
        kernel.append(values)
    blocks = np.asarray(kernel)
    err = float(np.nextafter(err, np.inf))
    pending=cache.with_name(cache.stem+f".{os.getpid()}.npz")
    np.savez_compressed(pending,blocks=blocks,entry_error=err)
    pending.replace(cache)
    print(f"Arb lag tables {g.side} {b:g}", err, elapsed(started), flush=True)
    return blocks, err


class SectorGeometry:
    """Receiver, incident and translated trace covariances."""
    def __init__(self, g, b):
        self.g = g
        self.pad = g.symbol.shape[0]
        self.b = b
        self.blocks, self.lag_error = lag_tables(g, b)
        self.started = clock()
        # Validate a direct DFT root table with Arb. Matrix arithmetic
        # then uses standard complex dot-product error bounds.
        ctx.prec = 128
        roots = [acb(0, -2*arb.pi()*i/self.pad).exp() for i in range(self.pad)]
        self.root_error = max(float(abs(v-ball(complex(v))).upper()) for v in roots)
        roots = np.asarray([complex(v) for v in roots])
        idx = np.arange(self.pad)
        self.dft = roots[(idx[:, None]*idx[None, :]) % self.pad]
        self.full_g = max(abs(z) for z in g.symbol.flat) + ARITH
        self.kappa = 1 - b*self.full_g
        assert self.kappa > 0
        self.gamma = 1/(1 - b*g.symbol.real)
        self.normS = fnorm(g.S)
        self.normB = fnorm(g.B)

        self.base_B = self.apply(b/2*g.B, self.gamma)
        self.base_GB = self.apply(b/2*g.B, self.gamma*g.symbol)
        ks = self.apply(g.S.conj().T, self.gamma)
        kgs = self.apply(g.S.conj().T, self.gamma*g.symbol)
        self.cross0 = ks.conj().T
        self.cross1 = kgs.conj().T
        self.p0 = hermitian(g.S @ ks)
        self.y0 = g.S @ self.base_B
        self.t0 = hermitian(b/2*g.B.conj().T @ self.base_B)
        # Spectral/frobenius errors of covariances, accounting for
        # both direct DFT factors and every product with exact S/B.
        self.rel_fft_bound = (64.0*self.pad*2.0**-53 +
                              32.0*np.sqrt(self.pad)*self.root_error) / self.kappa**2
        gscale = 1+self.full_g
        err = 2*self.rel_fft_bound + ARITH
        self.cross_error = err*max(1., self.normS*gscale)
        self.v_error = err*max(1., b/2*self.normB*gscale)
        self.p_error = err*max(1., self.normS**2/self.kappa)
        self.y_error = err*max(1., b/2*self.normS*self.normB/self.kappa)
        self.t_error = err*max(1., (b/2*self.normB)**2/self.kappa)

    def apply(self, x, symbol):
        side, pad = self.g.side, self.pad
        z = np.zeros((x.shape[1], pad, pad), dtype=np.complex128)
        z[:, :side, :side] = x.reshape(side, side, -1).transpose(2, 0, 1)
        U = self.dft
        transformed = U @ z @ U.T
        cropped = (U.conj() @ (transformed * symbol) @ U.conj().T)[:, :side, :side]
        return cropped.transpose(1,2,0).reshape(x.shape) / pad**2

    def covariance(self, a):
        G0 = self.g.restrict(a, self.blocks[0])
        G1 = self.g.restrict(a, self.blocks[1])
        G2 = self.g.restrict(a, self.blocks[2])
        return hermitian(np.block([[G0, G1], [G1.conj().T, G2]]))

    def selected(self, a):
        cross = np.concatenate((self.cross0[:, a], self.cross1[:, a]), axis=2).transpose(1, 0, 2)
        V = np.concatenate((self.base_B[a], self.base_GB[a] + self.g.B[a]), axis=1)
        return cross, V


def small_update(ell, high, b, weight):
    d = len(ell)
    x = weight-1.
    y = -weight*(ell + high)/2 + b/2
    z = weight*ell*high
    det = x*z-y*y
    if np.any(np.abs(det) <= ARITH*(1+np.abs(x*z)+y*y)):
        return None
    n_pos = int(np.count_nonzero(det < 0) +
                2*np.count_nonzero((det > 0)&(x+z > 0)))
    diagonal = np.arange(d)
    D = np.zeros((2*d, 2*d))
    D[diagonal, diagonal] = z/det
    D[diagonal+d, diagonal+d] = x/det
    D[diagonal+d, diagonal] = D[diagonal, diagonal+d] = -y/det
    # All endpoints and weights enter as exact binary64 scalars.
    error_sq = 0.
    for i in range(d):
        w = Fraction(float(weight))
        l, u, b_ = [Fraction(float(q)) for q in (ell[i], high[i], b)]
        a, c, f = w-1, -(l+u)*w/2+b_/2, w*l*u
        delta = a*f-c*c
        aerr = float(abs(Fraction(D[i,i])-f/delta))
        berr = float(abs(Fraction(D[i+d,i+d])-a/delta))
        cerr = float(abs(Fraction(D[i,i+d])+c/delta))
        error_sq += aerr*aerr + berr*berr + 2*cerr*cerr
    return D, n_pos, np.nextafter(np.sqrt(error_sq)*(1+matrix_epsilon(d)),np.inf)


def matrix_epsilon(n):
    return (64*n+128)*2.0**-53


def verified_inverse(A, error, separate=False):
    """Use checked approximate eigenvectors and inverse residual.
    
    Return null if the spectral sign or inverse cannot be verified.
    """
    if separate:
        v, Q = eigh(A, check_finite=False, driver="evd")
    else:
        v, Q = np.linalg.eigh(A)
    qt = Q.conj().swapaxes(-1, -2)
    identity = np.eye(A.shape[-1])
    f = fnorm(A, axis=(-1, -2))
    # Each approximate Q is a basis if its checked Gram is nonsingular.
    eta = matrix_epsilon(A.shape[-1])
    gram = fnorm(qt @ Q - identity, axis=(-1, -2))
    gram += eta*(1.+A.shape[-1])
    diag_error = fnorm(A @ Q - Q*v[..., None, :], axis=(-1, -2))
    test = (np.sqrt(1+gram)*diag_error + gram*np.max(np.abs(v),axis=-1)
            + (1+gram)*error + eta*(1+f)*A.shape[-1])
    distance = np.min(np.abs(v), axis=-1)
    valid = (gram < 1.) & (distance > test)
    v_safe = np.where(np.abs(v)>1e-100,v,1e-100)
    X = hermitian((Q/v_safe[...,None,:]) @ qt)
    xnorm = ((1.+gram)/np.min(np.abs(v_safe),axis=-1)) + ARITH
    rho = fnorm(identity-X@A,axis=(-1,-2)) + xnorm*error + (
        eta*max(1., np.sqrt(A.shape[-1]))*(1.+xnorm*f)
    )
    valid &= rho < 1.
    err = np.where(valid, xnorm*rho/np.maximum(1.-rho, 1e-100), np.inf)
    return valid, X, xnorm, err, np.count_nonzero(v < 0, axis=-1)


def local_certificate(sector, a, profile, lower, upper, weight):
    update = small_update(lower*profile, upper*profile, sector.b, weight)
    if update is None:
        return None
    Dinv, n_pos, d_error = update
    G = sector.covariance(a)
    d = len(Dinv)
    error = (d*sector.lag_error + d_error +
             16*2.0**-53*max(1.,fnorm(G)+fnorm(Dinv)))
    A = G + Dinv
    ok, X, xnorm, enorm, n_neg = verified_inverse(A, error, separate=True)
    if not ok:
        return None
    nu = d - n_neg - n_pos
    if nu < 0 or nu >= sector.g.S.shape[0] + sector.g.B.shape[1]:
        return None
    return A, error, X, xnorm+enorm, int(nu)


def response_reduction(sector, cells, certificate):
    A, ae, X, xnorm, nu = certificate
    C, V = sector.selected(cells)
    ct = C.conj().transpose(0,2,1)
    def batch_product(matrix,right):
        batch,d,cols=right.shape
        merged=right.transpose(1,0,2).reshape(d,batch*cols)
        return (matrix@merged).reshape(d,batch,cols).transpose(1,0,2)

    combined=np.concatenate((ct,V),axis=2)
    rays=batch_product(X,combined)
    propagated=batch_product(A,rays)
    outputs=C.shape[1]
    XC,XV=rays[:,:,:outputs],rays[:,:,outputs:]
    AC,AV=propagated[:,:,:outputs],propagated[:,:,outputs:]
    RC = ct - AC
    RV = V - AV
    n = lambda x: fnorm(x,axis=(-1,-2))
    cn,vn,ac,av = n(C),n(V),n(XC),n(XV)
    dc,dv = sector.cross_error,sector.v_error
    eta = matrix_epsilon(A.shape[-1])
    na = fnorm(A)
    erc = n(RC)+dc+ae*ac+eta*(1.+cn+na*ac)
    erv = n(RV)+dv+ae*av+eta*(1.+vn+na*av)

    def block(left,xl,rl,el,dl,nl,right,xr,rr,er,dr,nr):
        xlstar = xl.conj().transpose(0,2,1)
        leftstar = left.conj().transpose(0,2,1)
        # Exact at the chosen trial vectors up to the quadratic residual.
        # This avoids multiplying a large inverse norm by irrelevant
        # directions in the full trace-coordinate inverse.
        value = leftstar@xr + xlstar@rr
        arithmetic = eta*(1.+n(left)*nr+nl*(n(right)+na*nr))
        err = dl*nr+dr*nl+ae*nl*nr+xnorm*el*er+arithmetic
        return value,err

    pp,pe=block(ct,XC,RC,erc,dc,ac,ct,XC,RC,erc,dc,ac)
    yy,ye=block(ct,XC,RC,erc,dc,ac,V,XV,RV,erv,dv,av)
    tt,te=block(V,XV,RV,erv,dv,av,V,XV,RV,erv,dv,av)
    P=hermitian(sector.p0-pp)
    yc=sector.y0-yy
    theta=hermitian(sector.t0-tt)
    return P,yc,theta,sector.p_error+pe,sector.y_error+ye,sector.t_error+te,nu


def margin_distances(prepared, Y, r, scale_list, require_psd=False):
    P, yc, theta, pe, ye, te, nu = prepared
    e = Y-yc
    ef = fnorm(e,axis=(-1,-2))
    nt = len(P)
    s = Y.shape[1]
    margins = np.full(nt,-np.inf)
    ranks = np.full(nt,-1,dtype=np.int16)
    scales = np.zeros(nt)
    for scale in scale_list:
        port = np.eye(Y.shape[0]) + P/scale
        error_port = pe/scale + ARITH*(1.+fnorm(P,axis=(-1,-2))/scale)
        oks, port_inv, pinv_norm, pinv_err, neg = verified_inverse(port, error_port)
        n = nu-neg
        allow_rank=(n==0) if require_psd else ((n>=0)&(n<s))
        selected=np.flatnonzero(oks&allow_rank)
        if not len(selected):
            continue
        v = e[selected]
        vt = v.conj().transpose(0,2,1)
        H = hermitian(vt @ port_inv[selected] @ v - scale*theta[selected])
        dH = (pinv_err[selected]*(ef[selected]+ye[selected])**2 +
              pinv_norm[selected]*(2*ef[selected]+ye[selected])*ye[selected] +
              scale*te[selected] + ARITH*(1.+fnorm(H,axis=(-1,-2))+
                               pinv_norm[selected]*ef[selected]**2 +
                               scale*fnorm(theta[selected],axis=(-1,-2))))
        # Small source eigensystem gets an independent congruence check.
        sv, sq = np.linalg.eigh(H)
        st = sq.conj().transpose(0,2,1)
        orthogonal = fnorm(st@sq - np.eye(s),axis=(-1,-2)) + ARITH*s
        relation = fnorm(H@sq-sq*sv[:,None,:],axis=(-1,-2))
        correction = (np.sqrt(1.+orthogonal)*relation + orthogonal*(
                 np.max(np.abs(sv),axis=1) + fnorm(H,axis=(-1,-2))) +
                ARITH*s*(1.+fnorm(H,axis=(-1,-2))))
        valid = orthogonal < .5
        sv_lo = sv - (dH + correction)[:,None]
        sv_lo = np.where(np.arange(s)[None,:] < (s-n[selected])[:,None],sv_lo,0.)
        gap = np.sqrt(np.maximum(0.,sv_lo).sum(axis=1)) - r - ARITH*(1+r)
        change = (gap > margins[selected]) & valid
        margins[selected[change]] = gap[change]
        ranks[selected[change]] = n[selected[change]]
        scales[selected[change]] = scale
    return margins, ranks, scales


def run(side, shape, label, lo, hi, only_interior=False, broad=False,
        library=False, fast=False, require_psd=False, whole_panel=False, reuse_interior=False, blind=False):
    start = clock()
    g = Geometry(geom_path(side))
    b = FREEZE["model"]["background_upper"]
    setup = clock()
    sector = SectorGeometry(g, b)
    setup_time = elapsed(setup)
    if blind:
        public=HERE/"public"/"gallery"
        receipt=HERE/"results"/"gallery_sealed"
        receipt.mkdir(exist_ok=True)
    else:
        public=HERE/"public"
        receipt=HERE/"results"
    acq = json.loads((public/"acquisition.json").read_text())
    observations = {}
    for key, meta in acq["observations"].items():
        if meta["side"] != side:
            continue
        if blind:
            pass
        elif library:
            # The four preselected IDs give the panel membership only.
            # Their morphology never enters the online scoring.
            if "__weak__" not in key or (not whole_panel and "01__" not in key):
                continue
        elif meta["shape"] != shape:
            continue
        with np.load(public / f"{key}.npz", allow_pickle=False) as f:
            observations[key] = (f["Y"], float(f["r"]))
    output = {oid: dict(margin_sector=np.full(g.n,-np.inf),
                         negative_index=np.full(g.n,-1,dtype=np.int16),
                         weight=np.zeros(g.n), sector_scale=np.zeros(g.n))
              for oid in observations}
    if reuse_interior:
        for oid in output:
            stored=np.load(HERE/"results"/f"{oid}__b{b:g}__{label}_interior_sector.npz",
                           allow_pickle=False)
            for item in output[oid]:
                output[oid][item]=stored[item].copy()
    costs = defaultdict(float)
    nclip = 0
    time_local = clock()
    if fast:
        weights=(0.125,0.5,1.,2.,8.) if side==32 else (0.125,0.03125,0.5)
        scales=(.5,1.,.25,2.,.125,4.,.0625,.03125,.015625)
    else:
        weights = (0.125, 0.5, 1., 2., 8.)
        scales = tuple(2.**np.arange(-12,5,1))
    clips = sorted(local_models(g,shape),key=lambda a:len(a[0]),reverse=True)
    for ids,cells,profile in clips:
        if only_interior and nclip >= 1:
            break
        nclip += 1
        if reuse_interior and nclip==1:
            continue
        for weight in weights:
            work = (np.any([v["margin_sector"][ids]<=0 for v in output.values()],axis=0)
                    if fast else np.ones(len(ids),dtype=bool))
            if not np.any(work):
                break
            if broad:
                lower = np.where(profile == 1., lo/profile,0.)
                upper = hi/profile
            else:
                lower,upper = lo,hi
            one = clock()
            local = local_certificate(sector,cells[0],profile,lower,upper,weight)
            for k,v in elapsed(one).items():
                costs["state_"+k] += v
            if local is None:
                continue
            selected_ids,selected_cells=ids[work],cells[work]
            for offset in range(0,len(selected_ids),128):
                one = clock()
                sl = slice(offset,offset+128)
                row = selected_ids[sl]
                reduced = response_reduction(sector,selected_cells[sl],local)
                for oid,(Y,r) in observations.items():
                    margin,ranks,scale = margin_distances(reduced,Y,r,scales,require_psd)
                    entry = output[oid]
                    chosen = margin > entry["margin_sector"][row]
                    loc = row[chosen]
                    entry["margin_sector"][loc] = margin[chosen]
                    entry["negative_index"][loc] = ranks[chosen]
                    entry["weight"][loc] = weight
                    entry["sector_scale"][loc] = scale[chosen]
                for k,v in elapsed(one).items():
                    costs["online_"+k] += v
            if nclip == 1:
                print("sector weight",side,shape,label,weight,elapsed(start),
                      [int(np.count_nonzero(v["margin_sector"]>0)) for v in output.values()],
                      flush=True)
        if nclip % 32 == 1:
            print("sector",side,shape,label,nclip,len(clips),elapsed(start),
                  [int(np.count_nonzero(v["margin_sector"]>0)) for v in output.values()],flush=True)
    costs["evaluation"] = elapsed(time_local)
    summary = {}
    for oid,arrays in output.items():
        prefix = f"{oid}__b{b:g}__{label}_sector.npz"
        if library or blind:
            prefix=f"{shape}__{prefix}"
        np.savez_compressed(receipt/prefix,**arrays)
        summary[oid] = {
            "excluded": int(np.count_nonzero(arrays["margin_sector"]>0)),
            "negative_state_excluded": int(np.count_nonzero((arrays["margin_sector"]>0)&
                                                             (arrays["negative_index"]>0))),
            "verified_margin_available": int(np.count_nonzero(np.isfinite(arrays["margin_sector"]))),
        }
    costs["setup"] = setup_time
    save_json(receipt/f"sector_{side}_{shape}_{label}.json",
              dict(centers=g.n, clip_classes=nclip, total_clips=len(clips), interval=label,
                   weights=weights, scales=scales, footprint_broad=broad, interior_only=only_interior,
                   costs=dict(costs), summary=summary, elapsed=elapsed(start), fast=fast,
                   require_psd=require_psd, whole_panel=whole_panel, reused_interior_clip=reuse_interior,
                   opaque_public_ids=blind))
    print("sector result",side,shape,label,summary,elapsed(start),flush=True)


if __name__=="__main__":
    a=argparse.ArgumentParser()
    a.add_argument("side",type=int)
    a.add_argument("shape")
    a.add_argument("label")
    a.add_argument("lower",type=float)
    a.add_argument("upper",type=float)
    a.add_argument("--interior",action="store_true")
    a.add_argument("--broad",action="store_true")
    a.add_argument("--library",action="store_true")
    a.add_argument("--fast",action="store_true")
    a.add_argument("--psd",action="store_true")
    a.add_argument("--whole-panel",action="store_true")
    a.add_argument("--reuse-interior",action="store_true")
    a.add_argument("--blind",action="store_true")
    args=a.parse_args()
    with threadpool_limits(limits=1,user_api="blas"):
        run(args.side,args.shape,args.label,args.lower,args.upper,args.interior,args.broad,args.library,args.fast,args.psd,args.whole_panel,args.reuse_interior,args.blind)
