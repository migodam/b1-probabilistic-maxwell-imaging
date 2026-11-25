"""Algebraic checks for fixed certificate and constitutive-trace reductions.

These checks validate identities in binary64; they are not interval certificates
and are not an electromagnetic performance benchmark.
"""
from pathlib import Path
import json
import numpy as np
from scipy.linalg import null_space, orth

rng = np.random.default_rng(261126)
def crand(*shape):
    return (rng.standard_normal(shape) + 1j*rng.standard_normal(shape))/np.sqrt(2)
def herm(a): return (a+a.conj().T)/2
def inertia(a, tol=1e-8):
    ev=np.linalg.eigvalsh(herm(a))
    return [int(np.sum(ev < -tol)), int(np.sum(abs(ev)<=tol)), int(np.sum(ev>tol))]

n, k, d, s = 34, 10, 8, 5
q, _ = np.linalg.qr(crand(n,n))
lam=np.r_[np.linspace(0.8,4,23), -np.linspace(0.5,2,7), np.zeros(4)]
K=herm((q*lam)@q.conj().T)
T=crand(k,n)/np.sqrt(n)
E=crand(k,d)/np.sqrt(d)
L=crand(n,d)/np.sqrt(n)
C=herm(crand(d,d))
U=orth(np.c_[T.conj().T, q[:,23:]])
W=null_space(U.conj().T)
r=U.shape[1]
Cpp=herm(W.conj().T@K@W)
Qb=np.block([[K,-L],[-L.conj().T,C]])
Q_ret=np.block([[U.conj().T@K@U,-U.conj().T@L],[-L.conj().T@U,C]])
Cross=np.c_[W.conj().T@K@U, -W.conj().T@L]
Q_red=herm(Q_ret-Cross.conj().T@np.linalg.solve(Cpp,Cross))
F=np.c_[T@U,E]
Z=crand(d,s)
Pfull=np.block([[np.eye(n),np.zeros((n,s))],[np.zeros((d,n)),Z]])
Pred=np.block([[np.eye(r),np.zeros((r,s))],[np.zeros((d,r)),Z]])
cases=[]
for it in range(15):
    D=herm(crand(k,k))*(0 if it==0 else .3)
    Q=Qb+np.c_[T,E].conj().T@D@np.c_[T,E]
    QR=Q_red+F.conj().T@D@F
    for threshold in [0.,.3,1.1]:
        A=Pfull.conj().T@Q@Pfull
        Ar=Pred.conj().T@QR@Pred
        A[n:,n:]-=threshold*np.eye(s)
        Ar[r:,r:]-=threshold*np.eye(s)
        ia, ir=inertia(A),inertia(Ar)
        ok=ia==[ir[0],ir[1],ir[2]+n-r]
        assert ok,(it,threshold,ia,ir)
        cases.append(ok)

# Exact trace/energy completion, including a singular Gram from repeated rows.
pdim, rows, sources=19,9,4
F0=crand(rows,pdim)
F0=np.r_[F0,F0[:2]]
v=crand(pdim,sources)
b=F0@v
Sigma=herm(F0@F0.conj().T)
vmin=np.linalg.pinv(F0)@b
min_gram=herm(b.conj().T@np.linalg.pinv(Sigma)@b)
trace_error=float(np.linalg.norm(F0@vmin-b))
gram_error=float(np.linalg.norm(vmin.conj().T@vmin-min_gram))
slack_min=float(np.linalg.eigvalsh(herm(v.conj().T@v-min_gram)).min())
assert trace_error<1e-10 and gram_error<1e-10 and slack_min>-1e-10

# Sequential compensation innovation, with pseudoinverse/range handling.
j=5
S11,S12=Sigma[:j,:j],Sigma[:j,j:]
S21,S22=Sigma[j:,:j],Sigma[j:,j:]
b1,b2=b[:j],b[j:]
S11i=np.linalg.pinv(S11)
conditional=herm(S22-S21@S11i@S12)
innov=b2-S21@S11i@b1
sequential=herm(b1.conj().T@S11i@b1+innov.conj().T@np.linalg.pinv(conditional)@innov)
sequential_error=float(np.linalg.norm(sequential-min_gram))
assert sequential_error<1e-10

# One-sided truncation bound after whitening the omitted positive block.
Rdim, Edim, Ydim=8,14,11
Ee=crand(Rdim,Edim)*.015
Ff=crand(Rdim,Ydim)*.2
D=herm(crand(Rdim,Rdim))*.5
eps2=float(np.linalg.norm(Ee,2)**2)
dmax=float(np.linalg.norm(D,2))
alpha=1-dmax*eps2
assert alpha>0
actual=herm(Ff.conj().T@D@Ff-Ff.conj().T@D@Ee@np.linalg.solve(np.eye(Edim)+Ee.conj().T@D@Ee,Ee.conj().T@D@Ff))
lower=herm(Ff.conj().T@D@Ff-(eps2/alpha)*(Ff.conj().T@D@D@Ff))
lower_gap=float(np.linalg.eigvalsh(herm(actual-lower)).min())
assert lower_gap>-1e-10

result={
 "seed":261126,"kind":"binary64 algebraic checks; not rigorous interval certificates",
 "full_state_dimension":n,"retained_state_dimension":r,
 "positive_eliminated_dimension":n-r,
 "positive_complement_min_eigenvalue":float(np.linalg.eigvalsh(Cpp).min()),
 "max_update_leakage":float(np.linalg.norm(T@W)),
 "exact_inertia_cases_passed":len(cases),
 "trace_completion_residual":trace_error,
 "minimum_gram_identity_error":gram_error,
 "unused_positive_energy_min_eigenvalue":slack_min,
 "nested_trace_energy_identity_error":sequential_error,
 "truncation_complement_gap":alpha,
 "truncation_lower_bound_psd_min_eigenvalue":lower_gap,
}
out=Path(__file__).with_name("fixed_space_checks.json")
out.write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
