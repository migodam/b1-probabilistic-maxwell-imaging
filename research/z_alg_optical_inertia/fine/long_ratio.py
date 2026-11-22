"""Find and prove a native fine-grid dark current over a range of optical weights."""
import json, time
import numpy as np
from scipy.linalg import svd
from scipy.sparse.linalg import eigsh,LinearOperator
from flint import acb,arb,ctx
from threadpoolctl import threadpool_limits
from hidden_current import DATA,PUBLIC,angular_energy,ball,star,frobenius_squared,safe_float

def run():
    t=10000.;cap=.75
    start=time.perf_counter()
    with np.load(PUBLIC/'fine_geometry.npz') as f:s=f['PW360_10x20__S'];g=f['spectrum']
    p=len(g);side=96;n=side**2
    coord=(np.arange(side)+.5)*2/side-1
    yy,xx=np.meshgrid(coord,coord,indexing='ij')
    ang=2*np.pi*np.arange(100)/100
    field=np.exp(1j*(2*np.pi*4e8/3e8)*(xx.ravel()[:,None]*np.cos(ang)+yy.ravel()[:,None]*np.sin(ang)))
    V=svd(s,full_matrices=False)[2].conj().T
    ext=field-V@(V.conj().T@field);W=svd(ext,full_matrices=False)[0]
    trials=[]
    for num in (6,10,14,18,22,26):
        basis=np.column_stack((V,W[:,:num]))
        def project(z):return z-basis@(basis.conj().T@z)
        def react(z):
            f=np.fft.fft2(project(z).reshape(side,side),s=(p,p))
            return project(np.fft.ifft2(f*g.real)[:side,:side].ravel())
        _,v=eigsh(LinearOperator((n,n),react,dtype=complex),which='LA',k=1,ncv=22,tol=2e-10,
                  v0=project(np.random.default_rng(7393).normal(size=n)))
        v=project(v).ravel();z=np.abs(np.fft.fft2(v.reshape(side,side),s=(p,p))/p)**2
        score=(1-cap*(z*g.real).sum()+t*(z*g.imag).sum())
        trials.append((score,num,v));print(num,score,flush=True)
    score,num,v=min(trials,key=lambda x:x[0])
    ctx.prec=160;ctx.threads=1
    ss=ball(s);Vb=ball(v[:,None]);w=Vb-star(ss)*(ss*star(ss)).solve(ss*Vb)
    z=angular_energy(w,side,p);b=arb(0);c=arb(0)
    for x,row in enumerate(z):
        for y,value in enumerate(row):
            b+=arb(float(g[y,x].real))*value;c+=arb(float(g[y,x].imag))*value
    a=frobenius_squared(w);q=safe_float(a-arb(cap)*b+arb(t)*c,upper=True)
    receipt=dict(ratio=t,cap=cap,masked_plane_wave_modes=num,verified_form_upper=q,
                 radiation_lower=safe_float(c,upper=False),reactive_lower=safe_float(b,upper=False),
                 all_smaller_ratios_have_a_negative_hidden_direction=bool(q<0 and c>0),
                 numeric_proposals=[[float(item[0]),item[1]] for item in trials],wall_seconds=time.perf_counter()-start)
    (DATA/'fine_long_ratio.json').write_text(json.dumps(receipt,indent=2)+'\n')
if __name__=='__main__':
    with threadpool_limits(1):run()
