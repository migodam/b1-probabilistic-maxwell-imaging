"""Small finite-dimensional checks for FBG-R3-T; not a spatial imaging benchmark."""
from pathlib import Path
import json
import numpy as np
from numpy.linalg import solve, norm, eigh, eigvalsh
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT=Path(__file__).parent/'root_checks'
OUT.mkdir(exist_ok=True)
rng=np.random.default_rng(31102026)
n,p,r=6,3,4
z=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
G=.028*(z+z.T)+.08j*np.eye(n)
E=rng.normal(size=(n,p))+1j*rng.normal(size=(n,p))
S=rng.normal(size=(r,n))+1j*rng.normal(size=(r,n))
mref=np.full(n,.4)
mtarget=rng.uniform(.8,2.0,n)
K=rng.normal(size=(n,n))*.15

def forward(m,Es=E):
    L=np.eye(n)-m[:,None]*G
    j=solve(L,m[:,None]*Es)
    Q=Es+G@j
    return j,Q,L
Y=S@forward(mtarget)[0]

def objective(m,mode='plain',Es=E,Ys=Y,eps=0.,C=None):
    j,Q,L=forward(m,Es)
    phi=np.zeros(n) if mode=='plain' else K@(m-mref)
    V=np.exp(phi)[:,None]*Q
    H=V.conj().T@V
    if eps: H=H+eps*(np.eye(p) if C is None else C)
    R=S@j-Ys
    return .5*np.real(np.trace(R@solve(H,R.conj().T)))

def gradient(m,mode):
    j,Q,L=forward(m)
    phi=np.zeros(n) if mode=='plain' else K@(m-mref)
    scale=np.exp(phi)
    V=scale[:,None]*Q
    H=V.conj().T@V
    R=S@j-Y
    W=solve(H,R.conj().T).conj().T
    D=W.conj().T@W
    ans=np.zeros(n)
    for a in range(n):
        h=np.eye(n)[a]
        T=solve(L,h[:,None]*Q)
        dV=scale[:,None]*(G@T)
        if mode!='plain': dV+=scale[:,None]*K[:,a,None]*Q
        ans[a]=np.real(np.vdot(W,S@T)-np.vdot(V@D,dV))
    return ans

g0=gradient(mref,'plain'); g1=gradient(mref,'exp')
j,Q,L=forward(mref); R=S@j-Y
QR=R@solve(Q.conj().T@Q,Q.conj().T)
d=np.sum(np.abs(QR)**2,axis=0)
identity=norm((g1-g0)+K.T@d)
step=2e-6
fd=np.array([(objective(mref+step*np.eye(n)[a],'exp')-objective(mref-step*np.eye(n)[a],'exp'))/(2*step) for a in range(n)])
fd_error=norm(fd-g1)/max(1.,norm(g1))

C=np.array([[1.,.5+.2j,0.],[0.,3.,.2j],[.1,0.,.25]],complex)
invariant=abs(objective(mref,'exp')-objective(mref,'exp',E@C,Y@C))
ridge=.3
ridge_identity_discrepancy=abs(objective(mref,'exp',eps=ridge)-objective(mref,'exp',E@C,Y@C,eps=ridge))
ridge_covariant_error=abs(objective(mref,'exp',eps=ridge)-objective(mref,'exp',E@C,Y@C,eps=ridge,C=C.conj().T@C))

alpha=.7
mprobe=rng.uniform(.1,3.0,n)
_,Qprobe,_=forward(mprobe)
phi=alpha*np.tanh(K@(mprobe-mref)/alpha)
H0=Qprobe.conj().T@Qprobe
e,U=eigh(H0); Hminushalf=(U*(1/np.sqrt(e)))@U.conj().T
H1=(np.exp(phi)[:,None]*Qprobe).conj().T@(np.exp(phi)[:,None]*Qprobe)
relative_eigen=eigvalsh(Hminushalf@H1@Hminushalf)

beta=.2;b=.1;kappa=beta/(1+beta*b);mstar=1.;g=.5+.1j
xx=np.linspace(.02,20.,3000)
weights={
    'CIE':1+kappa*(xx-b),
    'Exponential':np.exp(kappa*(xx-b)),
    'Bounded (alpha=1)':np.exp(np.tanh(kappa*(xx-b)))
}
fig,ax=plt.subplots(figsize=(8.4,4.5))
for label,w in weights.items():
    obj=(xx-mstar)**2/((b-mstar)**2*w*w)
    ax.semilogy(xx,np.maximum(obj,1e-7),label=label,lw=2)
barrier=mstar+1/kappa
ax.axvline(mstar,color='black',ls='--',lw=1)
ax.axvline(barrier,color='#b45309',ls=':',lw=1.5)
ax.set(xlabel='Real material m',ylabel='Objective / common reference value',title='Same initial value and gradient do not determine finite search')
ax.text(barrier+.3,.045,'Exp: descent moves toward upper box',fontsize=9,color='#b45309')
ax.grid(alpha=.2);ax.legend(frameon=False);fig.tight_layout()
fig.savefig(OUT/'matched_metric_counterexample.png',dpi=200)
plt.close(fig)

def deriv_exp(x):
    return 2*(x-mstar)*np.exp(-2*kappa*(x-b))*(1-kappa*(x-mstar))
def deriv_cie(x):
    w=1+kappa*(x-b)
    return 2*(x-mstar)*(1+kappa*(mstar-b))/w**3
def deriv_bounded(x):
    phi=np.tanh(kappa*(x-b)); dp=kappa/np.cosh(kappa*(x-b))**2
    return 2*(x-mstar)*np.exp(-2*phi)*(1-(x-mstar)*dp)

summary={
 'scope':'Finite complex algebra only. No new 24x24 Maxwell images.',
 'seed':31102026,
 'gradient_identity_absolute_error':float(identity),
 'analytic_gradient_finite_difference_relative_error':float(fd_error),
 'nonunitary_source_mixing_absolute_error':float(invariant),
 'fixed_identity_ridge_breaks_covariance_difference':float(ridge_identity_discrepancy),
 'covariant_ridge_absolute_error':float(ridge_covariant_error),
 'bounded_gram_relative_eigenvalues':relative_eigen.tolist(),
 'bounded_gram_guaranteed_interval':[float(np.exp(-2*alpha)),float(np.exp(2*alpha))],
 'one_cell_model':{'G_real':g.real,'G_imag':g.imag,'E':1.,'S':1.,'true_material':mstar,'reference':b,'beta':beta,'kappa':kappa,'box':[.02,20.], 'exp_false_upper_basin_threshold':barrier,'derivative_at_m18':{'CIE':deriv_cie(18.),'exp':deriv_exp(18.),'bounded':deriv_bounded(18.)}},
 'scalar_K_field_with_G_0p1_plus_0p2i_at_m1':float(np.real((.1+.2j)/(1-(.1+.2j)))),
 'scalar_projector_derivative':0.
}
assert identity<1e-10 and fd_error<1e-7
assert invariant<1e-10 and ridge_covariant_error<1e-10
assert relative_eigen.min()>np.exp(-2*alpha) and relative_eigen.max()<np.exp(2*alpha)
assert deriv_exp(18.)<0<deriv_cie(18.) and deriv_bounded(18.)>0
(OUT/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps(summary,indent=2))
