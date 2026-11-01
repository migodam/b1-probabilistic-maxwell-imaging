"""FBG-R3-T: algebra checks only; not a Maxwell imaging benchmark."""
import json
from pathlib import Path
import numpy as np

def check_block():
    rng = np.random.default_rng(402)
    n,k,p,r = 9,3,4,5
    H = rng.normal(size=(n,n))*.025
    B = rng.normal(size=(n,n))*.02
    G = (H+H.T)/2+1j*(B@B.T+.01*np.eye(n))
    E = rng.normal(size=(n,p))+1j*rng.normal(size=(n,p))
    S = rng.normal(size=(r,n))+1j*rng.normal(size=(r,n))
    m = rng.uniform(.3,1,n)
    U = np.eye(n)[:,[0,3,7]]
    L = np.eye(n)-m[:,None]*G
    j = np.linalg.solve(L,m[:,None]*E)
    Q = E+G@j
    V = np.linalg.solve(L,U)
    C,q,O = U.T@G@V,U.T@Q,S@V
    d = np.array([.2,.4,-.1])
    A = np.eye(k)-d[:,None]*C
    eta = np.linalg.solve(A,d[:,None]*q)
    eta_tilde = eta+.015*(rng.normal(size=(k,p))+1j*rng.normal(size=(k,p)))
    b = q+C@eta_tilde
    dfit = np.real(np.sum(b.conj()*eta_tilde,axis=1))/np.sum(abs(b)**2,axis=1)
    defect = eta_tilde-dfit[:,None]*b
    Ad = np.eye(k)-dfit[:,None]*C
    mf = m+U@dfit
    jf = np.linalg.solve(np.eye(n)-mf[:,None]*G,mf[:,None]*E)
    corrected = S@j+O@eta_tilde-O@np.linalg.solve(Ad,defect)
    err = np.linalg.norm(S@jf-corrected)/np.linalg.norm(S@jf)
    # Derivative of material-profile defect at a physically realizable current.
    b0 = q+C@eta
    h = rng.normal(size=(k,p))+1j*rng.normal(size=(k,p))
    v = A@h
    coeff = np.real(np.sum(b0.conj()*v,axis=1))/np.sum(abs(b0)**2,axis=1)
    predicted = v-coeff[:,None]*b0
    def profile(x):
        bx=q+C@x
        dx=np.real(np.sum(bx.conj()*x,axis=1))/np.sum(abs(bx)**2,axis=1)
        return x-dx[:,None]*bx
    eps=1e-6
    fd=(profile(eta+eps*h)-profile(eta-eps*h))/(2*eps)
    derivative=np.linalg.norm(fd-predicted)/np.linalg.norm(predicted)
    return {'physical_output_identity_relative_error':float(err),'profile_normal_derivative_relative_error':float(derivative)}

def two_roots():
    # Rational decimal coefficients, symmetric G with positive Hermitian imaginary part.
    G=np.array([[-.06+.02j,.065+.007j],[.065+.007j,-.15+.02j]])
    E=np.array([1.,.7]); S=np.array([1.,1.3]); m=np.array([2.,7.])
    def state(x,e=E):
        L=np.eye(2)-x[:,None]*G
        j=np.linalg.solve(L,x*e)
        return S@j,L,e+G@j
    y,_,_=state(m)
    a,b,d=G[0,0],G[0,1],G[1,1]
    A=S[0]*E[0]+y*a
    B=S[1]*E[1]+y*d
    C=b*(S[0]*E[1]+S[1]*E[0])-d*S[0]*E[0]-a*S[1]*E[1]-y*(a*d-b*b)
    poly=np.array([np.imag(-A*np.conj(C)),np.imag(y*np.conj(C)-A*np.conj(B)),np.imag(y*np.conj(B))])
    roots=[]
    for x in np.roots(poly):
        z=(y-A*x)/(B+C*x)
        v=np.array([x.real,z.real]);f,L,Q=state(v)
        J=(S@np.linalg.inv(L))*Q
        JR=np.vstack([J.real,J.imag])
        fresh,_,_=state(v,np.array([.2,1.]))
        roots.append({'m':v.tolist(),'residual':float(abs(f-y)), 'imaginary_material':float(abs(z.imag)), 'jacobian_singular_values':np.linalg.svd(JR,compute_uv=False).tolist(),'L_minimum_singular_value':float(np.linalg.svd(L,compute_uv=False)[-1]),'second_source_response':[float(fresh.real),float(fresh.imag)]})
    m0=np.array(roots[0]['m']);m1=np.array(roots[1]['m'])
    j0=np.linalg.solve(np.eye(2)-m0[:,None]*G,m0*E)
    j1=np.linalg.solve(np.eye(2)-m1[:,None]*G,m1*E)
    q0=E+G@j0;q1=E+G@j1;t=.5
    jt=(1-t)*j0+t*j1;qt=E+G@jt;mt=(1-t)*m0+t*m1
    defect=jt-mt*qt;chord=t*(1-t)*(m1-m0)*(q1-q0)
    profile_m=np.real(qt.conj()*jt)/abs(qt)**2
    profile_defect=jt-profile_m*qt
    angle_formula=t*t*(1-t)**2*(m1-m0)**2*(abs(q0)**2*abs(q1)**2-np.real(q0.conj()*q1)**2)/abs(qt)**2
    bridge={'linear_material_defect_norm':float(np.linalg.norm(defect)),'chord_identity_absolute_error':float(np.linalg.norm(defect-chord)),'profile_material':profile_m.tolist(),'profile_defect_norm':float(np.linalg.norm(profile_defect)),'profile_angle_identity_relative_error':float(np.linalg.norm(abs(profile_defect)**2-angle_formula)/np.linalg.norm(abs(profile_defect)**2))}
    return {'G':[[[float(z.real),float(z.imag)] for z in row] for row in G], 'E':E.tolist(),'S':S.tolist(),'y':[float(y.real),float(y.imag)],'quadratic_coefficients':poly.tolist(),'roots':roots,'bridge':bridge,'claim_scope':'Exact rational algebraic VIE demonstration; not claimed to come from the user Hankel disk-pulse geometry.'}

result={'scope':'Algebra checks, no image reconstruction claim','block':check_block(),'isolated_roots':two_roots()}
path=Path(__file__).with_name('branches_checks.json')
path.write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
