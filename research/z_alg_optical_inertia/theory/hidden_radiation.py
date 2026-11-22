"""Upper bound on the ability of a global optical multiplier to remove invisible currents.

The certificate reads only the original discrete Green, sources, receiver and noisy fields.
Candidate projections and eigenpairs use SciPy; all inequalities marked verified use Arb.
"""
from pathlib import Path
import json
import sys
import time
import numpy as np
from scipy.linalg import eigh, svd
from scipy.sparse.linalg import LinearOperator, eigsh
from threadpoolctl import threadpool_limits
from flint import arb, acb, ctx

STUDY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(STUDY / 'core'))
from certificate import ball, star, frobenius_squared, safe_float, gershgorin
sys.path.insert(0, str(STUDY.parent / 'z_alg_secant_tasks/highdim'))
from replay import load_operators, SOURCE

CASES = {'r10_twod_0001': 5., 'austria_eps4_400MHz': 1.,
         's2d_007__pw360_10x20': 1., 's2d_007__pw180_10x20': 1.,
         's2d_010__pw180_10x20': 1.}
MU_MIN = .0001

def herm(a):
    return (a + a.conj().T) / 2

def upper(x):
    return safe_float(x, upper=True)

def lower(x):
    return safe_float(x, upper=False)

def neg_patch(green):
    side = int(np.sqrt(len(green)))
    for width in (7, 10, 16, side):
        ids = np.arange(len(green)).reshape(side, side)[:width, :width].ravel()
        sub = green[np.ix_(ids, ids)]
        if np.any(sub != sub.T):
            raise ValueError('Patch witness expects symmetric stored Green.')
        arr = sub.imag.copy()
        _, vectors = eigh(arr, driver='evd')
        R = ball(arr)
        for vec in vectors[:, :min(55, len(ids))].T:
            code = np.round((2.**32) * vec).astype('int64')
            v = ball(code[:, None])
            q = (star(v) * (R * v))[0, 0].real / frobenius_squared(v)
            if q < 0:
                return dict(status='strict negative Rayleigh quotient on stored Green', patch_width=width,
                            indices=ids.tolist(), integer_vector=code.tolist(),
                            normalized_upper=upper(q))
    return dict(status='no strictly negative patch witness found')

def receiver_bound(g, S, abg, abs_, u):
    """For SV=0, V*(R-S*X-X*S)V=V*RV for every X."""
    R, H = herm(g / 1j), herm(g)
    left, sig, vh = svd(S, full_matrices=False)
    v = vh.conj().T
    def proj(z):
        return z-v@(v.conj().T@z)
    n = len(g)
    v0 = np.random.default_rng(1950).normal(size=n)
    A = LinearOperator((n,n), lambda z:proj(H@proj(z)), dtype=complex)
    _, vec = eigsh(A, k=1, which='LA', tol=1e-12, v0=proj(v0), ncv=28)
    proposal = ball(proj(vec))
    st = star(abs_)
    ctx.prec=256
    exact = proposal-st*(abs_*st).solve(abs_*proposal)
    a = frobenius_squared(exact)
    response=(star(exact)*(abg*exact))[0,0]
    b,c = response.real, response.imag
    ctx.prec=96
    witnesses = dict(projected_vector_norm2_lower=lower(a), reactive_quadratic_form_lower=lower(b),
                     radiation_to_reactive_lower=lower(c/b), radiation_to_reactive_upper=upper(c/b),
                     cap_intercept_upper=upper(a/b))
    if c > 0 and arb(u)*b > a:
        witnesses['minimum_ratio_lambda_over_mu_to_remove_this_vector']=lower((arb(u)*b-a)/c)
    elif c <= 0 and arb(u)*b > a:
        witnesses['invisible_vector_stays_negative_for_all_positive_lambda']=True

    # A truncated fit need not invert S or solve any constraint exactly.
    # Its residual still proves an upper bound on the exact receiver nullspace.
    fits = []
    for cut in (1e-6,1e-10):
        mask=sig>sig[0]*cut
        V=v[:,mask]
        vr=V.conj().T@R
        X=(left[:,mask]/sig[mask])@(vr-0.5*(vr@V)@V.conj().T)
        computed=np.linalg.norm(R-S.conj().T@X-X.conj().T@S)
        fits.append((computed,X,int(mask.sum()),cut))
    _, X, rank, cut = min(fits,key=lambda f:f[0])
    sx=st*ball(X)
    Rball=(abg-star(abg))*acb(0,-.5)
    bound=frobenius_squared(Rball-sx-star(sx)).sqrt()
    return dict(kernel_radiation_operator_norm_upper=upper(bound), rank_used_in_fit=rank,
                fit_relative_singular_cutoff=cut, hidden_threshold=witnesses)

def cutoff(g,S,B,Y,u,abg,abs_):
    """Construct X=i alpha B, verify T=-Qim>0 and E*E/mu0+Qre-t*T<0."""
    # Scale the columns to make the proof insensitive to unequal incident norms.
    _,binv=np.linalg.qr(B,mode='reduced')
    W=np.linalg.inv(binv)
    b=B@W
    R=herm(g/1j)
    rr=eigh(herm(b.conj().T@R@b),eigvals_only=True)[-1]
    alpha=0.5/max(float(rr),1e-8)
    X=1j*alpha*b
    E=b+g@X
    gap=S@X-Y@W
    t= -herm((X.conj().T@E)/1j)
    d=herm(gap.conj().T@gap)/MU_MIN+herm(X.conj().T@X)-u*herm(X.conj().T@E)
    if eigh(t,eigvals_only=True)[0] <=0: return {'status':'non-positive trial'}
    est=max(0.,float(eigh(d,t,eigvals_only=True)[-1]))
    test=np.ceil(1.01*est+1.0)
    # W and X can be arbitrary. Treat these represented entries as exact.
    Wb=ball(W); Xb=ball(X)
    b=ball(B)*Wb; yy=ball(Y)*Wb
    e=b+abg*Xb
    xe=star(Xb)*e
    T=(xe-star(xe))*acb(0,.5)
    gap=abs_*Xb-yy
    D=(star(gap)*gap)/arb(MU_MIN)+star(Xb)*Xb-arb(u)*(xe+star(xe))/2
    # Verify via a real lower precision rotation proposal and strict Gershgorin.
    def rotate_positive(candidate,actual):
        qq=eigh(herm(candidate))[1]; Q=ball(qq)
        smallest=gershgorin(star(Q)*actual*Q,upper=False)
        gram=gershgorin(star(Q)*Q,upper=False)
        return lower(smallest),lower(gram)
    td,tg=rotate_positive(t,T)
    dd,dg=rotate_positive(test*t-d,arb(float(test))*T-D)
    wgram=lower(gershgorin(star(Wb)*Wb,upper=False))
    # A large Gershgorin disk in W*W can be avoided using its eigenbasis.
    if wgram<=0: wgram,_=rotate_positive(W.conj().T@W,star(Wb)*Wb)
    return dict(status='Arb verified' if min(td,tg,dd,dg,wgram)>0 else 'unresolved',
                mu_min=MU_MIN, cap=u, t_at_which_bound_vanishes=float(test),
                trial_alpha=alpha, min_rotated_T_gershgorin=td,
                min_rotated_tT_minus_D_gershgorin=dd, certified_congruence_lower=wgram)

def run():
    ctx.prec,ctx.threads=96,1
    results={}; cache={}
    for case,u in CASES.items():
        before=time.perf_counter()
        a=load_operators(case);g,s,b=(a[k] for k in ('G','S','E0'))
        with np.load(SOURCE/(case+'.npz')) as f:y=f['Y_archive']
        G,S=ball(g),ball(s)
        family='24-ring' if len(g)==576 else '32-400MHz'
        if family not in cache:cache[family]=neg_patch(g)
        result=receiver_bound(g,s,G,S,u)
        result['cutoff']=cutoff(g,s,b,y,u,G,S)
        result['green_family']=family
        result['wall_seconds']=time.perf_counter()-before
        results[case]=result
        print(case, result['kernel_radiation_operator_norm_upper'], result['hidden_threshold'],
              result['cutoff'],flush=True)
        out=STUDY/'results/theory/hidden_radiation.json'
        out.write_text(json.dumps(dict(cases=results,stored_power_witness=cache),indent=2)+'\n')
if __name__=='__main__':
    with threadpool_limits(1): run()
