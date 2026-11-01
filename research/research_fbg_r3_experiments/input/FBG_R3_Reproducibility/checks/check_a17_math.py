import json
from pathlib import Path
import numpy as np
from numpy.linalg import norm, solve, eigvalsh
from numpy.polynomial import Polynomial as Poly

rng=np.random.default_rng(20261031)
def cr(shape): return rng.normal(size=shape)+1j*rng.normal(size=shape)
def rel(a,b): return float(norm(a-b)/max(1.,norm(a),norm(b)))
def invsqrt(a):
    w,v=np.linalg.eigh(a)
    return (v/w**.5)@v.conj().T
def realify(a): return np.block([[a.real,-a.imag],[a.imag,a.real]])
def grad(fun,m,eps=1e-6):
    out=[]
    for i in range(len(m)):
        h=np.zeros_like(m);h[i]=eps
        out.append((fun(m+h)-fun(m-h))/(2*eps))
    return np.array(out)

n,p,r=6,3,4
G0=cr((n,n));G=.035*(G0+G0.T)
E=cr((n,p));S=cr((r,n))
m0=np.linspace(.2,.7,n);K=.13*rng.normal(size=(n,n));c=1.7
def state(m):
    L=np.eye(n)-m[:,None]*G
    j=solve(L,m[:,None]*E);Q=E+G@j
    return L,j,Q
L,j,Q=state(m0);Y=S@state(m0+.12*rng.normal(size=n))[1]+.1*cr((r,p));R=S@j-Y
def fm(m):
    _,jj,qq=state(m);rr=S@jj-Y
    vv=np.exp(K@(m-m0))[:,None]*qq
    return float(c/2*np.trace(rr@solve(vv.conj().T@vv,rr.conj().T)).real)
def f0(m):
    _,jj,qq=state(m);rr=S@jj-Y
    return float(c/2*np.trace(rr@solve(qq.conj().T@qq,rr.conj().T)).real)
H=Q.conj().T@Q;Qdag=solve(H,Q.conj().T)
d=np.sum(abs(R@Qdag)**2,axis=0)
gdiff=grad(fm,m0)-grad(f0,m0)
results={'gradient_identity_relative_error':rel(gdiff,-c*K.T@d)}

def kfield(E,Q):
    Z=Q@invsqrt(Q.conj().T@Q);out=np.zeros((n,n));wh=invsqrt(Q.conj().T@Q)
    for a in range(n):
        rhs=np.zeros_like(Q);rhs[a]=Q[a]
        T=solve(L,rhs)
        out[:,a]=np.real(np.sum(Z.conj()*(G@T@wh),axis=1))/(np.sum(abs(Z)**2,axis=1)+1e-5)
    return out
C=cr((p,p))+3*np.eye(p)
results['source_mix_Kfield_error']=rel(kfield(E,Q),kfield(E@C,Q@C))
results['source_mix_objective_error']=float(abs(c/2*np.trace((R@C)@solve((Q@C).conj().T@(Q@C),(R@C).conj().T)).real-f0(m0)))

B=.3*np.eye(n)+.03*rng.normal(size=(n,n));V=Q+B@j
A=V.conj().T@V;W=R@np.linalg.inv(A);D=W.conj().T@W
h=rng.normal(size=n);T=solve(L,h[:,None]*Q)
z=S.conj().T@W-(G+B).conj().T@V@D
pred=c*np.vdot(z,T).real
def fixed(m):
    _,jj,qq=state(m);rr=S@jj-Y;vv=qq+B@jj
    return float(c/2*np.trace(rr@solve(vv.conj().T@vv,rr.conj().T)).real)
eps=1e-6;fd=(fixed(m0+eps*h)-fixed(m0-eps*h))/(2*eps)
results['fixed_B_direction_error']=float(abs(pred-fd))

# A complex VIE, one source, real shared material and realified current coefficients.
rho=2.3;reg=.2*np.eye(n);rr=R[:,0]
Q1=Q[:,0];inj=np.diag(Q1)
Uall=np.linalg.qr(cr((n,3)))[0];U=Uall[:,:2];Ca=Uall[:,2:]
def stacked(U):return np.vstack((realify(S@U),rho**.5*realify(L@U)))
Av=stacked(U);Ac=stacked(Ca)
b0=np.r_[rr.real,rr.imag,np.zeros(2*n)]
BB=np.vstack((np.zeros((2*r,n)),-rho**.5*np.vstack((inj.real,inj.imag))))
Pv=Av@np.linalg.pinv(Av);Z=(np.eye(len(b0))-Pv)@Ac
Umat=invsqrt(Z.T@Z)@Z.T@BB
v=invsqrt(Z.T@Z)@Z.T@b0
Ho=BB.T@(np.eye(len(b0))-Pv)@BB+reg;go=BB.T@(np.eye(len(b0))-Pv)@b0
An=np.column_stack((Av,Ac));Pn=An@np.linalg.pinv(An)
Hn=BB.T@(np.eye(len(b0))-Pn)@BB+reg;gn=BB.T@(np.eye(len(b0))-Pn)@b0
po=-solve(Ho,go);pn=-solve(Hn,gn)
dp=solve(Hn,Umat.T@(v+Umat@po))
eta=float(norm(Umat@invsqrt(Ho),2)**2)
actual=float(np.sqrt((pn-po)@Ho@(pn-po)))
bound=float(eta**.5/(1-eta)*norm(v+Umat@po))
results.update(profile_H_error=rel(Hn,Ho-Umat.T@Umat),profile_g_error=rel(gn,go-Umat.T@v),profile_step_error=rel(pn-po,dp),profile_eta=eta,profile_step_Hnorm=actual,profile_step_bound=bound)

# False stiffness counterexample: no current vs one free current, scalar G=0.
rr0=1.;rho0=1.;epsS=1e-4;lam=1e-8
results['false_stiffness_example']={'H_old':rho0+lam,'H_new':rho0*epsS**2/(rho0+epsS**2)+lam,'eta':1-(rho0*epsS**2/(rho0+epsS**2)+lam)/(rho0+lam)}

# Exact two-cell-support finite rational response, checked on dense complex VIE.
inds=np.array([1,4]);UU=np.eye(n)[:,inds];dh=np.array([.4,-.2]);DD=np.diag(dh)
P=solve(L,UU);HH=UU.T@G@P;O=S@P;BB0=UU.T@Q;AA=DD@HH
detpoly=np.array([1.,-np.trace(AA),np.linalg.det(AA)],dtype=complex)
qcoeff=[R,detpoly[1]*R+O@DD@BB0,detpoly[2]*R+O@(AA-np.trace(AA)*np.eye(2))@DD@BB0]
NN=np.zeros(5,dtype=complex)
for ii in range(3):
    for jj in range(3):NN[ii+jj]+=np.vdot(qcoeff[ii],qcoeff[jj])
DP=np.convolve(detpoly.conj(),detpoly).real
NP=NN.real
poly=Poly(NP).deriv()*Poly(DP)-Poly(NP)*Poly(DP).deriv()
roots=poly.roots();ts=[float(t.real) for t in roots if abs(t.imag)<1e-7 and -.4<t.real<.4]+[-.4,.4]
def ray(t):return .5*norm(R+t*O@solve(np.eye(2)-t*AA,DD@BB0))**2
errs=[]
for t in [-.4,-.17,.09,.4]:
    md=m0+UU@(t*dh)
    full=S@state(md)[1]-Y
    approx=R+t*O@solve(np.eye(2)-t*AA,DD@BB0)
    errs.append(rel(full,approx))
results['finite_response_max_error']=max(errs)
results['ray_stationarity_polynomial_degree']=poly.trim(1e-11).degree()
results['ray_global_interval_minimum']={'t':min(ts,key=ray),'loss':min(map(ray,ts)),'candidate_count':len(ts)}

# CIE vs CSI readout: off-manifold bias enters at order two, tangent bias at four.
beta=.7;hh=rng.normal(size=n);TT=solve(L,hh[:,None]*Q)
VV=cr((n,p));orders={}
def project_real(q,v):
    cc=np.sum(abs(q)**2,axis=1)
    aa=np.sum(q.conj()*v,axis=1).real/cc
    return v-aa[:,None]*q
for name,V0,order in [('arbitrary_current',VV,2),('physical_tangent',TT,4)]:
    normal=project_real(Q,L@V0)
    coeff=beta*np.sum(abs(normal)**2,axis=1)/(np.sum(abs(Q)**2,axis=1)*(1+beta*m0))
    if order==4:
        coeff=beta*hh**2*np.sum(abs(project_real(Q,G@V0))**2,axis=1)/(np.sum(abs(Q)**2,axis=1)*(1+beta*m0))
    rows=[]
    for t in [.04,.02,.01,.005]:
        qt=Q+t*G@V0;jt=j+t*V0
        aa=np.sum(abs(jt)**2,axis=1);bb=np.sum(qt.conj()*jt,axis=1).real;cc=np.sum(abs(qt)**2,axis=1)
        mcsi=bb/cc;mbeta=(bb+beta*aa)/(cc+beta*bb)
        ee=jt-mcsi[:,None]*qt
        exact=beta*np.sum(abs(ee)**2,axis=1)/(cc+beta*bb)
        rows.append({'t':t,'bias_norm':float(norm(exact)),'scaled_coefficient_relative_error':float(norm(exact/t**order-coeff)/norm(coeff)),'finite_identity_error':float(norm(mbeta-mcsi-exact))})
    orders[name]={'expected_order':order,'observed_last_halving_order':float(np.log2(rows[-2]['bias_norm']/rows[-1]['bias_norm'])),'rows':rows}
results['CIE_bias_orders']=orders

out=Path(__file__).with_name('a17_math_checks.json');out.write_text(json.dumps(results,indent=2))
print(json.dumps(results,indent=2))
