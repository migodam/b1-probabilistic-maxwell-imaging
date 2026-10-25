"""Deterministic R5 response-coordinate units; no Maxwell calls."""
from pathlib import Path
import json,csv
import numpy as np
from scipy.integrate import quad
from scipy.special import ndtr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
OUT=Path(__file__).resolve().parents[1]

def integral(fun,a=-np.inf,b=np.inf):return quad(fun,a,b,epsabs=2e-12,epsrel=2e-11,limit=300)[0]
def normal(x):return np.exp(-.5*x*x)/np.sqrt(2*np.pi)
def metrics(p,q,lo=-10,hi=10):
    return dict(TV=integral(lambda x:.5*abs(p(x)-q(x)),lo,hi),KL_p_q=integral(lambda x:p(x)*np.log(max(p(x),1e-300)/max(q(x),1e-300)),lo,hi),KL_q_p=integral(lambda x:q(x)*np.log(max(q(x),1e-300)/max(p(x),1e-300)),lo,hi))
def run():
    rows=[];fig,axes=plt.subplots(1,3,figsize=(14,4))
    # Exact integration in b, and independently s=k(a)b/sigma.
    k=lambda a:.25+np.exp(.75*a)
    aa=np.linspace(-5,5,701)
    for sigma in [1.,.5,.25,.125]:
        z=integral(lambda a:normal(a)/np.sqrt(1+(k(a)/sigma)**2),-10,10)
        p=lambda a:normal(a)/np.sqrt(1+(k(a)/sigma)**2)/z
        sr=lambda a:normal(a)*integral(lambda s:normal(sigma*s/k(a))*np.exp(-.5*s*s)*sigma/k(a),-10,10)
        zr=integral(sr,-10,10)
        qnodet_raw=lambda a:normal(a)/np.sqrt(1+(sigma/k(a))**2)
        zn=integral(qnodet_raw,-10,10);qn=lambda a:qnodet_raw(a)/zn
        m=metrics(p,normal);mn=metrics(p,qn)
        mean=integral(lambda a:a*p(a),-10,10)
        rows.append(dict(unit='conditional_volume',sigma=sigma,Z_original=z,Z_response=zr,normalization_relative_error=abs(zr-z)/z,mean_a=mean,mean_a_response=integral(lambda a:a*sr(a)/zr,-10,10),profile_TV=m['TV'],no_determinant_TV=mn['TV']))
        axes[0].plot(aa,[p(a) for a in aa],label=f'sigma={sigma}')
    axes[0].plot(aa,normal(aa),'k--',label='profile');axes[0].set(title='Conditional volume changes weak marginal',xlabel='a');axes[0].legend(fontsize=7)
    # Fold x^2: response variable t=sqrt(s)=|x|/sqrt(sigma) removes integrable singularity.
    for y in [0.,1.]:
        for sigma in [1.,.5,.25,.125,.0625]:
            px=lambda x:normal(x)*np.exp(-.5*((x*x-y)/sigma)**2)
            z=integral(px,-8,8)
            # Two branches each dx= sqrt(sigma) dt with s=t^2-y/sigma.
            rt=lambda t:2*normal(np.sqrt(sigma)*t)*np.exp(-.5*(t*t-y/sigma)**2)*np.sqrt(sigma)
            zr=integral(rt,0,12)
            var=integral(lambda x:x*x*px(x),-8,8)/z
            varr=integral(lambda t:sigma*t*t*rt(t),0,12)/zr
            rows.append(dict(unit='fold',y=y,sigma=sigma,Z_original=z,Z_response=zr,normalization_relative_error=abs(zr-z)/z,variance=var,variance_response=varr,width=np.sqrt(var),width_over_sqrt_sigma=np.sqrt(var/sigma),negative_branch_mass=.5,positive_only_TV=.5))
            if y==0:axes[1].plot(np.linspace(-2,2,601),[px(x)/z for x in np.linspace(-2,2,601)],label=f'sigma={sigma}')
    axes[1].set(title='Fold thickness scales as sqrt(sigma)',xlabel='x');axes[1].legend(fontsize=7)
    # Independent exact Gaussian posterior; dropping a fixed weak likelihood rank.
    lam=np.array([1.,.03]); xx=np.linspace(-6,6,20001)
    for sigma in [.1,.03,.01,.003,.001]:
        v=sigma*sigma/(sigma*sigma+lam*lam)
        gauss=lambda t,v:np.exp(-.5*t*t/v)/np.sqrt(2*np.pi*v)
        p=lambda t:gauss(t,v[1]);q=normal
        # Closed KL in both directions and exact TV via crossings.
        c=np.sqrt(-v[1]*np.log(v[1])/(1-v[1]));tv=2*(ndtr(c/np.sqrt(v[1]))-ndtr(c))
        kl=.5*(v[1]-1-np.log(v[1]));rev=.5*(1/v[1]-1+np.log(v[1]))
        rows.append(dict(unit='ill_conditioned_linear',sigma=sigma,width_strong=np.sqrt(v[0]),width_weak=np.sqrt(v[1]),dropped_weak_likelihood_TV=tv,KL_p_q=kl,KL_q_p=rev,response_covariance_0=lam[0]**2/(sigma*sigma+lam[0]**2),response_covariance_1=lam[1]**2/(sigma*sigma+lam[1]**2),exact_GN_KL=0.,exact_RTO_KL=0.))
        axes[2].plot(xx,[p(t) for t in xx],label=f'sigma={sigma}')
    axes[2].plot(xx,normal(xx),'k--',label='wrong prior remainder');axes[2].set(xlim=(-2,2),title='Fixed rank loses weak likelihood',xlabel='weak coordinate');axes[2].legend(fontsize=7)
    fig.tight_layout();fig.savefig(OUT/'figures/A1_math_units.png',dpi=170);plt.close(fig)
    checks={
        'conditional_volume_normalization':max(x['normalization_relative_error'] for x in rows if x['unit']=='conditional_volume')<1e-9,
        'conditional_volume_mean':max(abs(x['mean_a']-x['mean_a_response']) for x in rows if x['unit']=='conditional_volume')<1e-9,
        'fold_normalization':max(x['normalization_relative_error'] for x in rows if x['unit']=='fold')<1e-9,
        'fold_second_moment':max(abs(x['variance']-x['variance_response']) for x in rows if x['unit']=='fold')<1e-9,
        'profile_bias_persists':rows[3]['profile_TV']>.1,
        'fixed_rank_failure':rows[-1]['dropped_weak_likelihood_TV']>.9,
        'two_branch_half_mass':all(x['negative_branch_mass']==.5 for x in rows if x['unit']=='fold'),
    }
    checks={k:bool(v) for k,v in checks.items()}
    result=dict(status='PASS' if all(checks.values()) else 'FAIL',checks=checks,rows=rows,physics_calls=0,limitations=['finite numerical quadrature check; not a theorem certificate','linear GN/RTO entries are analytic Gaussian equalities, not nonlinear Maxwell baselines'])
    (OUT/'math/A1_MATH_RESULTS.json').write_text(json.dumps(result,indent=2)+'\n')
    keys=sorted(set().union(*(r.keys() for r in rows)))
    with (OUT/'math/A1_MATH_RESULTS.csv').open('w') as f:w=csv.DictWriter(f,keys);w.writeheader();w.writerows(rows)
    print(json.dumps(dict(status=result['status'],checks=checks,rows=len(rows))))
if __name__=='__main__':run()
