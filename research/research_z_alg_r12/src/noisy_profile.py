"""Local minimum counts with finite, fold-scale noise in an overdetermined full VIE."""
from functools import lru_cache
import json
from time import perf_counter

import numpy as np
from scipy.optimize import brentq

from freq_model import Fit, ri
from freq_fold import exact_profile
from mirror_spectrum import CENTER, STRONG, V, Mirror, ROOT
from spectral_moments import PATTERNS

def run():
    t0=perf_counter()
    model=Mirror([1,2]); _,jac,h=model.response(CENTER,hessian=True)
    N=np.linalg.svd(ri(jac)@STRONG,full_matrices=True)[0][:,1:]
    c=N.T@ri(np.einsum('kij,i,j->k',h,V,V))/2
    d=ri(model.response(CENTER,derivative=True)[3]@V)
    dN=N.T@d
    A=c@c
    theta_star=np.sqrt(A)/np.linalg.norm(d)  # New energy fraction 1/2.
    records=[]
    for a in [.2,.6]:
        basis=(CENTER,V,STRONG,a,True)
        for rep in range(24):
            side=-1 if rep%2==0 else 1
            truth=CENTER+side*a*V
            rng=np.random.default_rng(238471+rep)
            gaussian=rng.standard_normal((3,6))+1j*rng.standard_normal((3,6))
            for pat in ('opposite_edges','two_positive'):
                for bandwidth in (.85,1.15):
                    slots,pow=PATTERNS[pat]
                    powers=np.asarray(pow); w=powers/sum(powers)
                    theta=theta_star*bandwidth*np.array(slots)
                    delta=a*theta
                    eps=.45*np.sqrt(A)*gaussian/np.sqrt(powers[:,None])
                    e=np.stack([ri(q) for q in eps]); em=w@e; mu=w@theta
                    var=w@((theta-mu)**2)
                    K=2*side*c+mu*dN
                    B=c@K
                    C=K@K+var*(d@d)-2*c@(N.T@em)
                    T=K@(N.T@em)+np.einsum('j,j,jk,k',w,theta-mu,e-em,d)
                    def predicted(co,t):
                        roots=np.roots([2*A,3*B,co,-t])
                        return sorted([float(z.real+side) for z in roots if abs(z.imag)<1e-8 and -4<z.real+side<4])
                    predicted_roots=predicted(C,T)
                    central_prediction=predicted(K@K+var*(d@d),0.)
                    observations=[model.response(truth,dt)[0]+a*a*n for dt,n in zip(delta,eps)]
                    f=Fit(model,delta,observations,np.sqrt(powers))
                    @lru_cache(maxsize=100)
                    def prof(x):
                        _,g,hh,sg,sh=exact_profile(f,basis,x)
                        return g/(A*a**4),hh/(A*a**4),sg,sh
                    grid=np.linspace(-4.,4.,17)
                    hh=[prof(x)[1] for x in grid]
                    turning=[brentq(lambda x:prof(x)[1],l,r,xtol=1e-9) for l,r,h0,h1 in zip(grid,grid[1:],hh,hh[1:]) if h0*h1<0]
                    limits=[-4.,*turning,4.]
                    roots=[brentq(lambda x:prof(x)[0],l,r,xtol=1e-9) for l,r in zip(limits,limits[1:]) if prof(l)[0]*prof(r)[0]<0]
                    records.append(dict(a=a,rep=rep,side=side,pattern=pat,band_over_symmetric_threshold=bandwidth,
                        band=float(a*max(abs(theta))),coefficients_A_B_C_T=[A,B,C,T],predicted=predicted_roots,
                        noise_free_prediction=central_prediction,measured=roots,
                        second_derivatives=[prof(x)[1] for x in roots],minimum_strong_hessian=min(prof(x)[3] for x in roots),
                        noise_seed=238471+rep))
    result=dict(per_real_component_sigma_over_a2_c=.45,window=[-4,4],halfband_scale=theta_star,
                wall_s=perf_counter()-t0,factorizations=model.factorizations,runs=records)
    (ROOT/'results/noisy_profile.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':
    run()
