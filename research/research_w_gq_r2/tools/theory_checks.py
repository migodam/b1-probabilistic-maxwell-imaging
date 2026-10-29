"""Only new R2 formula wiring; independent dense algebra, zero PDE."""
from pathlib import Path
import json,numpy as np
ROOT=Path(__file__).resolve().parents[1]
def main():
    rng=np.random.default_rng(28102941)
    def rand(*shape):return rng.normal(size=shape)+1j*rng.normal(size=shape)
    n=7;m=5;s=2;I=np.eye(n);C0=.02*rand(n,n);H=.03*rand(n,n);a=rand(n,s);B=rand(m,n);d=rand(m,s);y=rand(m,s);V=.1*rand(n,n);V0=.07*rand(n,n);E=.1*rand(n,n)
    def state(H,V):
        T=V@np.linalg.inv(I-(C0+H)@V);q=T@a;f=d+B@q;r=y-f;p=T.conj().T@B.conj().T@r
        return T,q,f,r,p
    T,q,f,r,p=state(H,V);T0,q0,f0,r0,p0=state(H,V0)
    def contrast(H):
        *_,r,p=state(H,V);*_,r0,p0=state(H,V0)
        return -.5*np.sum(abs(r)**2)+.5*np.sum(abs(r0)**2)
    xi=p@q.conj().T-p0@q0.conj().T;analytic=np.vdot(xi,E).real;h=1e-5;fd=(contrast(H+h*E)-contrast(H-h*E))/(2*h)
    freeT=V@np.linalg.inv(I-C0@V);split=freeT@np.linalg.inv(I-H@freeT)
    # Scattered field: b(W,X) is affine but not constant; mixed derivative retains both loads.
    A0=3*np.eye(n)+.1*rand(n,n);Mw=.05*rand(n,n);Mx=.04*rand(n,n);b0=rand(n,s);bw=rand(n,s);bx=rand(n,s)
    w=.3;x=.4;L=A0-w*Mw-x*Mx;u=np.linalg.solve(L,b0+w*bw+x*bx);uw=np.linalg.solve(L,bw+Mw@u);ux=np.linalg.solve(L,bx+Mx@u);uwx=np.linalg.solve(L,Mw@ux+Mx@uw)
    def sol(w,x):return np.linalg.solve(A0-w*Mw-x*Mx,b0+w*bw+x*bx)
    hh=1e-3;mix=(sol(w+hh,x+hh)-sol(w+hh,x-hh)-sol(w-hh,x+hh)+sol(w-hh,x-hh))/(4*hh**2)
    # Exact two-background contrast identity includes residual projection.
    f10=rand(m,s);f20=rand(m,s);D1=rand(m,s);D2=rand(m,s)
    c1=np.vdot(y-f10,D1).real-.5*np.sum(abs(D1)**2);c2=np.vdot(y-f20,D2).real-.5*np.sum(abs(D2)**2)
    identity=np.vdot(y-.5*(f10+f20+D1+D2),D1-D2).real-np.vdot(f10-f20,.5*(D1+D2)).real
    rows=dict(T_self_split_relative=float(np.linalg.norm(T-split)/np.linalg.norm(T)),contrast_cotangent_relative=float(abs(analytic-fd)/max(abs(analytic),1e-12)),material_load_mixed_relative=float(np.linalg.norm(uwx-mix)/np.linalg.norm(uwx)),two_background_identity_error=float(abs(c1-c2-identity)),PDE_calls=0,scope='dense algebra formula wiring only')
    rows['status']='PASS' if rows['T_self_split_relative']<1e-12 and rows['contrast_cotangent_relative']<1e-7 and rows['material_load_mixed_relative']<1e-5 and rows['two_background_identity_error']<1e-10 else 'FAIL'
    (ROOT/'results/THEORY_CHECKS.json').write_text(json.dumps(rows,indent=2));print(json.dumps(rows))
    if rows['status']!='PASS':raise RuntimeError(rows)
if __name__=='__main__':main()
