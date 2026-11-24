"""Verify a hard-support branch by exact active-current elimination on the native model."""

import numpy as np
from flint import acb, acb_mat, arb
from scipy.fft import ifft2
from scipy.linalg import cholesky, eigh

from branches import expand
from native_box import (ball, clock, elapsed, frobenius_squared, gershgorin, herm, hi, identity,
                        lo, mid, star)


def exact_kernel(geometry):
    if not hasattr(geometry,'active_kernel'):
        start=clock()
        g=geometry.g
        f=[acb.dft([acb(float(v.real),float(v.imag)) for v in row],inverse=True) for row in g]
        columns=[acb.dft([row[j] for row in f],inverse=True) for j in range(len(f))]
        geometry.active_kernel=[columns[j][i] for i,j in np.ndindex(g.shape)]
        geometry.active_kernel_cost=elapsed(start)
    return geometry.active_kernel


def numeric_box(S,B,G,u,y,t,alpha):
    H=(G+G.conj().T)/2
    R=(G-G.conj().T)/(2j)
    A=S.conj().T@S+alpha*(np.diag(1/u)-H+t*R)
    L=S.conj().T@y+(.5*alpha*(1+1j*t))*B
    ev,v=eigh(herm(A),check_finite=False,driver='evd')
    nu=int(np.count_nonzero(ev<0))
    lens=y.shape[1]
    if nu>=lens or np.min(np.abs(ev))<1e-12:
        return 0.,nu,None
    x=v@((v.conj().T@L)/ev[:,None])
    F=herm(y.conj().T@y-L.conj().T@x)
    feasibility=float(np.sqrt(np.maximum(0,eigh(F,eigvals_only=True)[::-1][nu:]).sum()))
    return feasibility,nu,(A,ev,v,x)


def interval_box(G,S,B,u,y,t,alpha,proposal):
    A,ev,v,x=proposal
    n=len(u)
    at,aa=arb(t),arb(alpha)
    Z=acb(aa,aa*at)/2
    k=star(S)*S-Z*G-star(Z*G)
    for i,val in enumerate(u):
        k[i,i]+=aa/arb(float(val))
    L=star(S)*y+Z*B
    nu=int(np.count_nonzero(ev<0))
    sigma=float(np.min(np.abs(ev))*.45)
    stabilizer=v[:,:nu]*np.sqrt(-2*ev[:nu])
    F=ball(cholesky(herm(A+stabilizer@stabilizer.conj().T)-sigma*np.eye(n),lower=True,check_finite=False))
    w=ball(stabilizer)
    lifted=k+w*star(w)-F*star(F)-arb(sigma)*identity(n)
    gap=(arb(sigma)-frobenius_squared(lifted).sqrt().upper()).lower()
    if not gap>0:
        raise ArithmeticError('Active shifted Cholesky residual did not prove a positive gap.')
    if nu:
        neg=ball(v[:,:nu])
        if not gershgorin(star(neg)*k*neg,upper=True)<0:
            raise ArithmeticError('Active negative subspace did not verify.')
    x=ball(x)
    kx=k*x
    correction=(frobenius_squared(kx-L)/gap).upper()
    h=star(y)*y-star(L)*x-star(x)*L+star(x)*kx-correction*identity(y.ncols())
    _,vec=eigh(herm(mid(h)),check_finite=False)
    terms=[]
    for j in range(nu+1,y.ncols()+1):
        vs=ball(vec[:,-j:])
        low=gershgorin(star(vs)*h*vs,upper=False)
        gram=star(vs)*vs
        if not low>0:
            break
        if not gershgorin(gram,upper=False)>0:
            raise ArithmeticError('Unresolved source congruence.')
        terms.append((low/gershgorin(gram,upper=True)).lower())
    total=sum(terms,arb(0))
    return (max(0.,lo(total.sqrt())) if total>0 else 0.),dict(
        exact_active_inertia=nu,positive_gap_lower=lo(gap),elimination_correction_upper=hi(correction),
        source_tail_eigenvalue_lowers=list(map(lo,terms)))


def score_hard(geometry,y,radius,branch,*,t_grid=(0.,.75,2.,4.,8.),
               alpha_grid=(2.**-17,2.**-14,2.**-11,2.**-8,2.**-6),max_active=1152,max_proofs=2):
    start=clock()
    material=expand(branch,geometry.side)
    active=np.flatnonzero(material)
    n=len(active)
    if n==0:
        return geometry.score(y,radius,branch)
    result=dict(branch=branch.name,decision='RETAIN',residual_lower=0.,radius=radius,
                rigor='no_pruning_certificate',active_native_cells=n,native_side=geometry.side,
                method='exact active support',screen=[],proofs=[])
    if np.all(material>0) or n>max_active:
        return {**result,'status':'active_support_budget',**elapsed(start)}
    jj,ii=np.divmod(active,geometry.side)
    p=geometry.p
    index=((jj[:,None]-jj[None,:])%p)*p+(ii[:,None]-ii[None,:])%p
    green=ifft2(geometry.g).ravel()[index]
    ports=geometry.s[:,active]
    incident=geometry.b[active]
    vals=material[active]
    proposed=[]
    for t in t_grid:
        for alpha in alpha_grid:
            value,nu,ev=numeric_box(ports,incident,green,vals,y,t,alpha)
            result['screen'].append(dict(t=t,alpha=alpha,nu=nu,numeric_residual=value))
            if value>radius:
                proposed.append((value,t,alpha))
    verify_start=clock()
    prepared=hasattr(geometry,'active_kernel')
    for _,t,alpha in sorted(proposed,reverse=True)[:max_proofs]:
        started=clock()
        try:
            kernel=exact_kernel(geometry)
            g=acb_mat(n,n,[kernel[int(k)] for k in index.ravel()])
            _,_,proposal=numeric_box(ports,incident,green,vals,y,t,alpha)
            lower,diagnostic=interval_box(g,ball(ports),ball(incident),vals,ball(y),t,alpha,proposal)
            result['proofs'].append(dict(t=t,alpha=alpha,strict_residual_lower=lower,**diagnostic,**elapsed(started)))
            result['residual_lower']=max(result['residual_lower'],lower)
            if lower>radius:
                result['decision']='PRUNE'
                result['rigor']='strict_native_exact_support'
                break
        except (ArithmeticError,ZeroDivisionError,ValueError) as error:
            result['proofs'].append(dict(t=t,alpha=alpha,error=str(error),**elapsed(started)))
    if hasattr(geometry,'active_kernel'):
        result['kernel_reused']=prepared
        result['kernel_cost']=geometry.active_kernel_cost
    return {**result,'interval_time':elapsed(verify_start),**elapsed(start)}
