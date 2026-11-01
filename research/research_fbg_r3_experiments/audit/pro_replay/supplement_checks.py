"""Algebra checks for A17 dressing equivalence and invariant soft current release.
Synthetic complex symmetric VIE; not a Green-grid imaging benchmark.
"""
import json
from pathlib import Path
import numpy as np

rng=np.random.default_rng(2026103107)
def cr(shape): return rng.normal(size=shape)+1j*rng.normal(size=shape)
def realify(a): return np.block([[a.real,-a.imag],[a.imag,a.real]])
def isqrt(a):
    d,v=np.linalg.eigh(a)
    return (v/d**.5)@v.T
def rel(a,b):return float(np.linalg.norm(a-b)/max(1.,np.linalg.norm(a),np.linalg.norm(b)))

n,p,r,k=7,2,3,2
G0=cr((n,n));G=.04*(G0+G0.T)
m=np.linspace(.15,.6,n);L=np.eye(n)-m[:,None]*G
E=cr((n,p));j=np.linalg.solve(L,m[:,None]*E);Q=E+G@j
S=cr((r,n));rho=1.3
qall=np.linalg.qr(cr((n,k+1)))[0];U=qall[:,:k];C=qall[:,k:]
R_U=U@np.linalg.solve(U.conj().T@L@U,U.conj().T)
Cbar=(np.eye(n)-R_U@L)@C
Cscaled=3.7*Cbar
def block_current(B):
    a=np.vstack((S@B,np.sqrt(rho)*L@B))
    return realify(np.kron(np.eye(p),a))
Av=block_current(U);Ac=block_current(C);Ab=block_current(Cbar);At=block_current(Cscaled)
P=Av@np.linalg.pinv(Av)
Z=(np.eye(P.shape[0])-P)@Ac
Zbar=(np.eye(P.shape[0])-P)@Ab
Zscaled=(np.eye(P.shape[0])-P)@At
bc=[];Bc=[]
for s in range(p):
    bc.extend(np.r_[.2*cr((r,)),np.zeros(n,dtype=complex)])
    Bc.append(np.vstack((np.zeros((r,n)), -np.sqrt(rho)*np.diag(Q[:,s]))))
bc=np.asarray(bc);Bc=np.vstack(Bc)
b=np.r_[bc.real,bc.imag];B=np.vstack((Bc.real,Bc.imag))
H=B.T@(np.eye(P.shape[0])-P)@B+.08*np.eye(n)
g=B.T@(np.eye(P.shape[0])-P)@b
po=-np.linalg.solve(H,g)

def release(Z,gamma,metric='residual'):
    gram=Z.T@Z
    ridge=gamma*gram if metric=='residual' else gamma*np.eye(gram.shape[0])
    Cmat=gram+ridge
    D=isqrt(Cmat)@Z.T@B;v=isqrt(Cmat)@Z.T@b
    Hn=H-D.T@D;gn=g-D.T@v
    pn=-np.linalg.solve(Hn,gn)
    eta=float(np.linalg.norm(D@isqrt(H),2)**2)
    return Hn,gn,pn,eta,D,v

old=release(Z,0);dressed=release(Zbar,0);scaled=release(Zscaled,0)
gamma=.8
soft=release(Z,gamma);soft_scaled=release(Zscaled,gamma)
naive=release(Z,gamma,'coeff');naive_scaled=release(Zscaled,gamma,'coeff')
eta0=old[3];cap=.05;chosen_gamma=max(0.,eta0/cap-1.)
chosen=release(Z,chosen_gamma)

# Normal signatures must also be invariant under retained-space dressing and
# invertible coefficient changes. Use a canonical current-space quotient and
# whole-block whitening, not independent column normalization.
Ur=realify(np.kron(np.eye(p),U))
Cr=realify(np.kron(np.eye(p),C))
Lr=realify(np.kron(np.eye(p),L))
Mcur=Lr.T@Lr
Pcur=Ur@np.linalg.solve(Ur.T@Mcur@Ur,Ur.T@Mcur)
Bq_complex=np.vstack([np.diag(Q[:,s]) for s in range(p)])
Bq=np.vstack((Bq_complex.real,Bq_complex.imag))
Nop=(np.eye(2*n*p)-Bq@np.linalg.pinv(Bq))@Lr
Tm=rng.normal(size=(Cr.shape[1],Cr.shape[1]))
Tm=Tm+3*np.eye(Tm.shape[0])
Km=rng.normal(size=(Ur.shape[1],Cr.shape[1]))
Cmixed=Cr@Tm+Ur@Km
def normal_score(current_columns,quotient):
    Cq=(np.eye(Mcur.shape[0])-Pcur)@current_columns if quotient else current_columns
    W=Cq@isqrt(Cq.T@Mcur@Cq)
    return float(np.linalg.norm(Nop@W,'fro')**2/W.shape[1])
nu=normal_score(Cr,True);nu_mixed=normal_score(Cmixed,True)
nu_raw=normal_score(Cr,False);nu_raw_mixed=normal_score(Cmixed,False)
result={
 'label':'SYNTHETIC_COMPLEX_VIE_ALGEBRA_NOT_IMAGING',
 'seed':2026103107,'n':n,'sources':p,'receivers':r,
 'dressing_projected_Z_error':rel(Z,Zbar),
 'hard_dressing_H_error':rel(old[0],dressed[0]),
 'hard_dressing_step_error':rel(old[2],dressed[2]),
 'hard_scaled_H_error':rel(old[0],scaled[0]),
 'invariant_soft_scaled_H_error':rel(soft[0],soft_scaled[0]),
 'invariant_soft_scaled_step_error':rel(soft[2],soft_scaled[2]),
 'naive_coefficient_ridge_scaled_H_difference':rel(naive[0],naive_scaled[0]),
 'naive_coefficient_ridge_scaled_step_difference':rel(naive[2],naive_scaled[2]),
 'soft_downdate_scaling_error':rel(H-soft[0],(H-old[0])/(1+gamma)),
 'soft_eta_scaling_error':float(abs(soft[3]-eta0/(1+gamma))),
 'eta_hard':eta0,'eta_cap':cap,'chosen_gamma':chosen_gamma,'eta_after_soft_release':chosen[3],
 'soft_step_identity_error':rel(soft[2]-po,np.linalg.solve(soft[0],soft[4].T@(soft[5]+soft[4]@po)))
 ,'quotient_normal_score':nu
 ,'quotient_normal_score_reparameterization_error':abs(nu-nu_mixed)
 ,'raw_normal_score_retained_dressing_difference':abs(nu_raw-nu_raw_mixed)
 ,'separable_group_normal_trace_prediction':1.-1./(2*p)
 ,'separable_group_normal_trace_error':abs(nu-(1.-1./(2*p)))
}
Path(__file__).with_name('supplement_checks.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
