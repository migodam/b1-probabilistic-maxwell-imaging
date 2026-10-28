"""Material-normalized geometry and discrete Maxwell jets, all real coordinates.

Taylor coefficients have NO factorials: x(h)=x0+h*v+h²*x2+h³*x3.
The decoder is the imported original P1 decoder; formulas below differentiate it.
These are local, piecewise-analytic jets, not a global C3 seam assertion.
"""
from __future__ import annotations
import numpy as np
from scipy.linalg import eigh, null_space
from scipy.optimize import brentq
from scipy.special import ndtr


def real_material(a):
    a = np.asarray(a)
    return np.concatenate((a.real, a.imag), axis=0) / 24.


def decode(prior, x):
    q = np.asarray(x, float).copy()
    q[64:] = ndtr(q[64:])
    return prior.decode(q)


def material_derivative(prior, x):
    q = np.asarray(x, float).copy()
    q[64:] = ndtr(q[64:])
    chain = np.r_[np.ones(64), np.exp(-.5*x[64:]**2)/np.sqrt(2*np.pi)]
    return prior.derivative(q)*chain[None, :]


def decoder_coefficients(prior, x, v, x2=None, x3=None):
    x2 = np.zeros(66) if x2 is None else np.asarray(x2)
    x3 = np.zeros(66) if x3 is None else np.asarray(x3)
    material, info = decode(prior, x)
    s = info['field']/prior.temperature
    inside = (abs(s)<1) & prior.mask
    hp = (15./16)*(1-s*s)**2*inside
    hpp = (15./4)*(s**3-s)*inside
    hppp = (15./4)*(3*s*s-1)*inside
    s1, s2, s3 = [prior.A@t[:64]/prior.temperature for t in (v,x2,x3)]
    h0 = info['fraction']
    h1 = hp*s1
    h2 = hp*s2 + .5*hpp*s1*s1
    h3 = hp*s3 + hpp*s1*s2 + hppp*s1**3/6
    widths = np.array([np.diff(prior.real_interval)[0],np.diff(prior.imag_interval)[0]])
    xi = x[64:]
    pp = widths*np.exp(-.5*xi*xi)/np.sqrt(2*np.pi)
    d1 = pp*v[64:]
    d2 = pp*(x2[64:]-.5*xi*v[64:]**2)
    d3 = pp*(x3[64:]-xi*v[64:]*x2[64:]+(xi*xi-1)*v[64:]**3/6)
    c = lambda t: complex(t[0],t[1])
    d0 = info['inclusion_increment']
    return (material, d0*h1+c(d1)*h0,
            d0*h2+c(d1)*h1+c(d2)*h0,
            d0*h3+c(d1)*h2+c(d2)*h1+c(d3)*h0)


def geometry(prior, x, J, seed):
    D = material_derivative(prior,x)
    M = real_material(D)
    U,s,Vh = np.linalg.svd(M,full_matrices=False)
    good = s > max(s[0]*1e-8,1e-10)
    if not good.any():
        return dict(rank=0,singular_material=s,material_cutoff=float(max(s[0]*1e-8,1e-10)))
    T = Vh[good].T/s[good][None,:]
    Q = M@T
    K = J@T
    A=K.T@K+T.T@T
    ev,ew=eigh(A)
    w=ew[:,0]
    # Fix direction gauge by largest absolute native component, independent of GT.
    v=T@w
    if v[np.argmax(abs(v))]<0: w=-w;v=-v
    _,ks,kVh=np.linalg.svd(K,full_matrices=False)
    strong=[]
    for candidate in kVh:
        c=candidate-w*(w@candidate)
        for z in strong:c-=z*(z@c)
        if np.linalg.norm(c)>1e-9:strong.append(c/np.linalg.norm(c))
        if len(strong)>=min(8,len(s[good])-1):break
    W=np.column_stack(strong) if strong else np.zeros((len(w),0))
    latent_U=T@W
    # A normal prior displacement, projected off the SAME material strong frame.
    rng=np.random.default_rng(seed)
    # Project a NATIVE prior displacement using the true material constraints.
    # Projecting in Q-coordinates then lifting through T spuriously amplifies
    # nearly null latent coordinates and makes an artificially weak random arm.
    N=null_space((Q@W).T@M,rcond=1e-10)
    rv=N@(N.T@rng.normal(size=66))
    rv/=max(np.linalg.norm(M@rv),1e-300)
    rw=Q.T@(M@rv)
    if rv[np.argmax(abs(rv))]<0:rv=-rv;rw=-rw
    groups=[slice(0,256),slice(256,512),slice(None)]
    spectra=[];frames=[]
    for sl in groups:
        _,sv,vh=np.linalg.svd(K[sl],full_matrices=False)
        spectra.append(sv)
        frames.append(Q@vh[-min(4,len(vh)):].T)
    angles=np.linalg.svd(frames[0].T@frames[1],compute_uv=False)
    _,wa=eigh(K[:256].T@K[:256]+T.T@T)
    va=T@wa[:,0]
    if va[np.argmax(abs(va))]<0:va=-va
    return dict(D=D,M=M,T=T,Q=Q,K=K,singular_material=s,
                rank=int(good.sum()),material_cutoff=float(max(s[0]*1e-8,1e-10)),v_group_A=va,
                material_frame=U[:,good],v=v,random_v=rv,strong_U=latent_U,
                weak_material=real_material(D@v),random_material=real_material(D@rv),
                generalized_joint_singular=ks,prior_aware_eigenvalues=ev,
                source_spectra=np.array(spectra),source_weak_principal_cosines=angles)


def response_action(target, base, dmaterial):
    """Cached receiver transfer action; no new solve, count arithmetic explicitly."""
    if base.T is None:target.jacobian(base)
    target.book.counts['cached_material_receiver_actions']+=1
    return target.white(base.T@(dmaterial[:,None]*base.state.field))


def common_compensation(A,b):
    # Real least squares ONLY. All illumination blocks share the same columns.
    z=-np.linalg.lstsq(A,b,rcond=1e-10)[0]
    return z,b+A@z


def jets(target, base, v, strong_U):
    x=base.x
    m0,m1,m2raw,m3raw=decoder_coefficients(target.prior,x,v)
    G=target.backend().G
    J=target.jacobian(base)
    A=J@strong_U
    j1=target.solve_anchor(m1[:,None]*base.state.field,base,role='jet1')
    f1=target.white(target.backend().S@j1)
    propagation2=target.white(base.T@(m1[:,None]*(G@j1)))
    decoder2=response_action(target,base,m2raw)
    b2=propagation2+decoder2
    z2,c2=common_compensation(A,b2)
    x2=strong_U@z2
    _,_,m2,_=decoder_coefficients(target.prior,x,v,x2)
    j2=target.solve_anchor(m2[:,None]*base.state.field+m1[:,None]*(G@j1),base,role='jet2')
    # Third-order fixed term contains the already shared second-order material.
    _,_,_,m3fixed=decoder_coefficients(target.prior,x,v,x2)
    propagation3=target.white(base.T@(m1[:,None]*(G@j2)+m2[:,None]*(G@j1)))
    decoder3=response_action(target,base,m3fixed)
    b3=propagation3+decoder3
    z3,c3=common_compensation(A,b3)
    x3=strong_U@z3
    _,_,_,m3=decoder_coefficients(target.prior,x,v,x2,x3)
    j3=target.solve_anchor(m3[:,None]*base.state.field+m1[:,None]*(G@j2)+m2[:,None]*(G@j1),base,role='jet3')
    checks=dict(chain_rule_relative=float(np.linalg.norm(f1-J@v)/max(np.linalg.norm(f1),1e-300)),
                recurrence2_relative=float(np.linalg.norm(target.white(target.backend().S@j2)-c2)/max(np.linalg.norm(c2),1e-300)),
                recurrence3_relative=float(np.linalg.norm(target.white(target.backend().S@j3)-c3)/max(np.linalg.norm(c3),1e-300)))
    group={}
    zz=[];cc=[];HH=[]
    for name,sl in [('A',slice(0,256)),('B',slice(256,512)),('JOINT',slice(None))]:
        z,c=common_compensation(A[sl],b2[sl])
        group[name]=dict(z2=z,c2=c,A=A[sl],b2=b2[sl],rank=int(np.linalg.matrix_rank(A[sl])))
        if name!='JOINT':zz.append(z);cc.append(c);HH.append(A[sl].T@A[sl])
    zbar=group['JOINT']['z2']
    within=sum(float(c@c) for c in cc)
    conflict=sum(float((z-zbar)@H@(z-zbar)) for z,H in zip(zz,HH))
    group['identity']=dict(joint_squared=float(c2@c2),within=within,conflict=conflict,
                           absolute_defect=abs(float(c2@c2)-within-conflict),
                           shared_lower_order=True,real_coordinate_compensation=True)
    zd2,_=common_compensation(A,decoder2)
    xd2=strong_U@zd2
    md3=decoder_coefficients(target.prior,x,v,xd2)[3]
    zd3,_=common_compensation(A,response_action(target,base,md3))
    xd3=strong_U@zd3
    M=real_material(material_derivative(target.prior,x))
    uu,ss,vh=np.linalg.svd(M,full_matrices=False)
    keep=ss>max(ss[0]*1e-8,1e-10)
    r2=real_material(m2raw)
    material_coordinate2=uu[:,keep]@(uu[:,keep].T@r2)
    material_normal2=r2-material_coordinate2
    # Unit material speed; Maxwell first order is generally NONZERO.
    return dict(v=v,strong_U=strong_U,x2=x2,x3=x3,m1=m1,m2raw=m2raw,m2=m2,
                m3raw=m3raw,m3fixed=m3fixed,m3=m3,j1=j1,j2=j2,j3=j3,
                c1=f1,c2=c2,c3=c3,b2=b2,b3=b3,decoder2=decoder2,decoder3=decoder3,
                decoder_only_x2=xd2,decoder_only_x3=xd3,
                material_coordinate2=material_coordinate2,material_normal2=material_normal2,
                propagation2=propagation2,propagation3=propagation3,A=A,
                checks=checks,groups=group)


def polynomial(x,v,h,x2=None,x3=None):
    return x+h*v+(0 if x2 is None else h*h*x2)+(0 if x3 is None else h*h*h*x3)


def match_distance(prior,x0,proposal,distance,max_latent=4.):
    """Find the FIRST positive bracketed decoder distance, never use GT/Maxwell."""
    m0=decode(prior,x0)[0]
    def fun(h):return np.linalg.norm(decode(prior,proposal(h))[0]-m0)/24-distance
    lo=0.; hi=max(distance/8,1e-8)
    for _ in range(100):
        x=proposal(hi)
        if not np.isfinite(x).all():
            return None,dict(status='NONFINITE_DECODER_PROPOSAL',last_h=hi)
        if np.linalg.norm(x-x0)>max_latent:
            # An overshooting search bracket is not itself an infeasibility proof.
            hi=brentq(lambda h:np.linalg.norm(proposal(h)-x0)-max_latent,lo,hi)
            if fun(hi)<0:
                return None,dict(status='INFEASIBLE_LATENT_RADIUS',last_h=hi,
                                 maximum_bracket_RMS=float(fun(hi)+distance))
        if fun(hi)>=0:
            h=brentq(fun,lo,hi,xtol=1e-14,rtol=1e-12)
            xx=proposal(h)
            return xx,dict(status='MATCHED',h=h,actual_RMS=float(fun(h)+distance),latent_step=float(np.linalg.norm(xx-x0)))
        lo=hi;hi*=1.35
    return None,dict(status='DISTANCE_BRACKET_NOT_FOUND')


def retract_distance(prior,x0,x,distance):
    delta=x-x0
    return match_distance(prior,x0,lambda h:x0+h*delta,distance)


def tasks(prior,material):
    mask=prior.mask
    right=mask & (prior.points[:,0]>=0)
    mean=material[mask].mean()
    fraction=np.mean(material[mask].real>=.75)
    return np.array([mean.real,mean.imag,np.mean(abs(material[right]-prior.host)**2),fraction,float(fraction>=.25)])


def feasible_progress(prior,x0,x,unit,distance):
    a=real_material(decode(prior,x)[0]-decode(prior,x0)[0])
    rms=np.linalg.norm(a)
    return abs(rms-distance)<=.02*distance and float(unit@a)>=.5*distance


def gn_step(prior,x0,x,v,J,residual,trust=.5,material_limit=.01):
    """All effective material coordinates, signed-progress tangent constraint.

    The distance is retracted with the REAL decoder after this subproblem.
    Identity prior precision is original Gaussianized prior, never an extra prior.
    """
    M=real_material(material_derivative(prior,x))
    M0=real_material(material_derivative(prior,x0))
    unit=M0@v;unit/=max(np.linalg.norm(unit),1e-300)
    row=unit@M
    _,sv,vh=np.linalg.svd(M,full_matrices=False)
    keep=sv>max(sv[0]*1e-8,1e-10)
    effective=vh[keep].T
    Z=effective@null_space((row@effective)[None,:],rcond=1e-10)
    augmented=np.vstack((J@Z,Z))
    r=np.r_[residual,x]
    u=-np.linalg.lstsq(augmented,r,rcond=1e-10)[0]
    step=Z@u
    factor=min(1.,trust/max(np.linalg.norm(step),1e-300),material_limit/max(np.linalg.norm(M@step),1e-300))
    step*=factor
    predicted=.5*float(r@r-np.linalg.norm(r+np.r_[J@step,step])**2)
    return step,predicted
