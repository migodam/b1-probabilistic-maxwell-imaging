"""Test negative directions of the actual (not Fourier-truncated) 96² Green."""
import json
import time
import sys
import numpy as np
from scipy.linalg import svd
from scipy.sparse.linalg import LinearOperator,eigsh
from flint import arb,acb,ctx
from threadpoolctl import threadpool_limits

from fourier_bound import PUBLIC,DATA,STUDY
sys.path.insert(0,str(STUDY/'core'))
from certificate import ball,star,frobenius_squared,safe_float

def angular_energy(v,side,p):
    z=acb(0)
    rows=[acb.dft([v[i*side+j,0] for j in range(side)]+[z]*(p-side)) for i in range(side)]
    out=[]
    for x in range(p):
        col=acb.dft([rows[y][x] for y in range(side)]+[z]*(p-side))
        out.append([(w.real*w.real+w.imag*w.imag)/(p*p) for w in col])
    return out

def run():
    start=time.perf_counter();ctx.prec,ctx.threads=144,1
    with np.load(PUBLIC/'fine_geometry.npz') as data:
        spec=data['spectrum']; p=len(spec);side=96;n=side**2
        out={}
        for name in ('PW360_10x20','PW180_10x20'):
            S=data[name+'__S']; _,_,V=svd(S,full_matrices=False);V=V.conj().T
            def project(x):return x-V@(V.conj().T@x)
            def reactive(z):
                ext=np.fft.fft2(project(z).reshape(side,side),s=(p,p))
                return project(np.fft.ifft2(ext*spec.real)[:side,:side].ravel())
            _,vec=eigsh(LinearOperator((n,n),reactive,dtype=complex),k=1,which='LA',tol=2e-11,
                        ncv=30,v0=project(np.random.default_rng(3126).normal(size=n)))
            v=ball(project(vec));s=ball(S);v=v-star(s)*(s*star(s)).solve(s*v)
            power=angular_energy(v,side,p)
            a=frobenius_squared(v);b=arb(0);c=arb(0)
            for x,line in enumerate(power):
                for y,w in enumerate(line):
                    b+=arb(float(spec[y,x].real))*w
                    c+=arb(float(spec[y,x].imag))*w
            slopes=(a/b,c/b);rows={}
            for cap in (.5,.75):
                for t in (0.,.3,.4,.55,.75,1.):
                    upper=safe_float(a-arb(cap)*b+arb(t)*c,upper=True)
                    rows[f'{cap}_{t}']=dict(cap=cap,t=t,quadratic_upper=upper,
                                             exact_K_has_negative_direction_at_all_mu=upper<0)
            out[name]=dict(cap_intercept_upper=safe_float(slopes[0],upper=True),
                           optical_slope_upper=safe_float(slopes[1],upper=True),
                           cases=rows,wall_seconds=time.perf_counter()-start)
            print(name,out[name]['cap_intercept_upper'],out[name]['optical_slope_upper'],flush=True)
    (DATA/'fine_hidden.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':
    with threadpool_limits(1):run()
