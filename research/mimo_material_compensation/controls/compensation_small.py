"""Reproducible passive reciprocal six-cell VIE network mechanism check.
This is a dissipative algebraic VIE model, not a discretized free-space Maxwell array.
The .json records floating-point checks; theorem inequalities are analytic,
but these floating-point norm evaluations are not interval certificates.
"""
import numpy as np
from scipy.optimize import root
import json
from pathlib import Path
np.set_printoptions(precision=12,suppress=True)
n=6
pos=np.array([[0.,0.],[.5,0.],[1.,0.],[0.,.6],[.5,.6],[1.,.6]])
dist=np.linalg.norm(pos[:,None,:]-pos[None,:,:],axis=-1)
# Rational symmetric ring coupling: this is a finite passive network model.
ring=np.minimum(np.abs(np.arange(n)[:,None]-np.arange(n)[None,:]),n-np.abs(np.arange(n)[:,None]-np.arange(n)[None,:]))
H=np.choose(ring,[.08,.06,.025,.01])
B=.19*np.array([[1.,.2,.8],[.8,1.,.1],[.2,.7,1.],[-.5,.9,.7],[.6,-.4,1.],[-.7,-.5,.3]])
gamma=.24
Gamma=B@B.T+gamma*np.eye(n)
G=H+1j*Gamma
x0=np.array([1.7,1.5,2.0,1.3,1.8,1.6])
Btx=B[:,:2]
S=2j*B[:,2:3].T
I=np.eye(n)

def pack(z):
    z=np.ravel(z)
    return np.r_[z.real,z.imag]

def forward(x):
    return S@np.linalg.solve(I-np.diag(x)@G,np.diag(x)@Btx)

def jac(x):
    R=np.linalg.inv(I-np.diag(x)@G)
    U=S@R
    E=np.linalg.solve(I-G@np.diag(x),Btx)
    return np.column_stack([pack(U[:,i:i+1]@E[i:i+1,:]) for i in range(n)])
J=jac(x0)
A=J[:,1:]
C=np.zeros((n,4))
C[1:]=np.linalg.pinv(A)
h=np.zeros(n);h[0]=1.
v=h-C@J@h
R=np.linalg.inv(I-np.diag(x0)@G)
U=S@R
K=G@R
E=np.linalg.solve(I-G@np.diag(x0),Btx)
k=float(np.linalg.norm(K,2))
D0=float(np.linalg.norm(U,2)*np.linalg.norm(E,'fro'))
c=float(np.max(np.linalg.norm(C,axis=1)))
vnorm=float(np.max(np.abs(v)))
eta=float(np.linalg.norm(np.eye(4)-J@C,2))
zeta=float(np.linalg.norm(J@v))
margins=np.minimum(x0,4.5-x0)
rows=np.linalg.norm(C,axis=1)

# A nonzero explicit certified-under-exact-norm-bounds segment.
def tube(a,b):
    d=a*vnorm+b*c
    if k*d>=1: return None
    q=k*d
    kap=eta+c*D0*((1-q)**(-2)-1)
    ga=a*zeta+D0*k*(a*vnorm)**2/(1-k*a*vnorm)
    legal=bool(np.all(a*np.abs(v)+b*rows<=margins))
    return dict(a=a,b=b,d=d,q=q,kappa=kap,forcing_bound=ga,
                contraction_budget=ga+kap*b,legal_box=legal,
                success=bool(kap<1 and ga+kap*b<b and legal))
best=None
for a in np.geomspace(1e-8,.5,4000):
    # Radius need not be tied to its natural O(a^2) scale if optimization improves it.
    scale=D0*k*(a*vnorm)**2/max(1-k*a*vnorm,1e-15)
    for mul in (1.1,1.25,1.5,2.,3.,4.,8.):
        z=tube(float(a),float(mul*scale))
        if z and z['success'] and (best is None or a>best['a']): best=z
if best is None: raise RuntimeError('no nonzero tube found')
a=best['a'];b=best['b']
Y=forward(x0)
path=[]
for t in np.linspace(-a,a,41):
    f=lambda z:pack(forward(x0+t*v+C@z)-Y)
    sol=root(f,np.zeros(4),tol=1e-11)
    xx=x0+t*v+C@sol.x
    path.append(dict(t=float(t),x=xx.tolist(),z_norm=float(np.linalg.norm(sol.x)),
                     residual=float(np.linalg.norm(f(sol.x))),
                     task=float(xx[0]),min_x=float(xx.min()),max_x=float(xx.max())))
# Longer non-certified predictor-corrector continuation: exact data constraint,
# one ROI coordinate; reset least-norm compensation at each step.
long=[]
for sgn in (-1,1):
    xx=x0.copy()
    for it in range(300):
        dt=sgn*.005
        JJ=jac(xx);AA=JJ[:,1:]
        CC=np.zeros((n,4));CC[1:]=np.linalg.pinv(AA)
        vv=h-CC@JJ@h
        xp=xx+dt*vv
        fun=lambda z:pack(forward(xp+CC@z)-Y)
        sol=root(fun,np.zeros(4),tol=1e-10)
        xn=xp+CC@sol.x
        res=float(np.linalg.norm(fun(sol.x)))
        if not (res<1e-10 and xn.min()>.02 and xn.max()<4.48): break
        xx=xn
        long.append(dict(t=float(xx[0]-x0[0]),x=xx.tolist(),residual=res,
                         exterior_min_singular=float(np.linalg.svd(jac(xx)[:,1:],compute_uv=False)[-1])))
Tfull=np.eye(3)+2j*B.T@np.linalg.solve(np.diag(1/x0)-G,B)
report=dict(model='six-cell dissipative reciprocal algebraic VIE; not free-space Maxwell',
            positive_material_box=[0.,4.5],x0=x0.tolist(),
            H=H.tolist(),B=B.tolist(),gamma=gamma,
            min_eigenvalue_ImG=float(np.linalg.eigvalsh(G.imag)[0]),
            full_scattering_max_singular=float(np.linalg.svd(Tfull,compute_uv=False)[0]),
            external_jacobian_singular_values=np.linalg.svd(A,compute_uv=False).tolist(),
            k=k,D0=D0,c=c,v=v.tolist(),eta=eta,zeta=zeta,
            tube=best,floating_tube_candidate_width=2*a,
            max_path_residual=max(p['residual'] for p in path),
            max_path_z_norm=max(p['z_norm'] for p in path),
            path=path,long_continuation=long,
            long_task_range=[min(p['t'] for p in long),max(p['t'] for p in long)],
            max_long_residual=max(p['residual'] for p in long),
            rigor_status='analytic passive/stability identities; floating-point numerical inequalities, no outward-rounded interval arithmetic')
out=Path(__file__).with_suffix('.json');out.write_text(json.dumps(report,indent=2))
print(json.dumps({k:report[k] for k in ['model','min_eigenvalue_ImG','full_scattering_max_singular','external_jacobian_singular_values','k','D0','c','v','eta','zeta','tube','floating_tube_candidate_width','max_path_residual','max_path_z_norm','long_task_range','max_long_residual','rigor_status']},indent=2))
