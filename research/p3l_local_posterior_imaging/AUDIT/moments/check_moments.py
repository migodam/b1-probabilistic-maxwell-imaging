"""Independent zero-physics audit; run from the repository root."""
import json
import os
from pathlib import Path
import sys
import time
import warnings

os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
os.environ.setdefault('MKL_NUM_THREADS','1')
import numpy as np
from scipy.integrate import quad
from scipy.special import ndtr
import mpmath as mp

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT/'src'))
from p3l.moments import (gaussian_saturation_moments, gaussian_material_moments,
                         linear_gaussian_field_marginals)
from p1_material.prior import saturation

OUT = Path(__file__).resolve().parent
SEED = 20261022


def independent_quad(mu,sd,power):
    """Direct original P1 decoder; no moment recurrence or production integrand."""
    if sd == 0:
        return float(saturation(mu/.5)**power)
    a,b = (-.5-mu)/sd,(.5-mu)/sd
    lower,upper = max(a,-38),min(b,38)
    q = 0.
    if lower < upper:
        q = quad(lambda z: float(saturation((mu+sd*z)/.5))**power*
                 np.exp(-z*z/2)/np.sqrt(2*np.pi), lower,upper,
                 points=[0.] if lower < 0 < upper else None,
                 epsabs=1e-310,epsrel=2e-11,limit=300)[0]
    return q+ndtr(-b)


def tiny_interior_polynomial(mu,sd):
    """400-digit full Gaussian polynomial integration (remote seams negligible)."""
    with mp.workdps(400):
        mu,sd=mp.mpf(str(mu)),mp.mpf(str(sd))
        c=[mp.mpf('.5'),mp.mpf('1.875'),0,-5,0,6]
        d=[sum(c[j]*c[k-j] for j in range(6) if 0 <= k-j < 6)
           for k in range(11)]
        def moment(n):
            return sum(mp.binomial(n,j)*mu**(n-j)*sd**j*
                       (mp.factorial(j)/(2**(j//2)*mp.factorial(j//2)))
                       for j in range(0,n+1,2))
        m=sum(c[k]*moment(k) for k in range(6))
        s=sum(d[k]*moment(k) for k in range(11))
        return float(m),float(s-m*m)


def high_precision_tail(mu,sd):
    with mp.workdps(90):
        # Preserve exact binary input values when the mean is near a seam.
        m,s=mp.mpf(mu),mp.mpf(sd)
        a,b=(-mp.mpf('.5')-m)/s,(mp.mpf('.5')-m)/s
        lo,hi=max(a,-38),min(b,38)
        shift=max(lo,0)
        def moment(power):
            def integrand(z):
                t=m+s*z+mp.mpf('.5')
                p=t**3*(10+t*(-15+6*t))
                return p**power*mp.exp(-(z*z-shift*shift)/2)/mp.sqrt(2*mp.pi)
            breaks=[lo]+([0] if lo < 0 < hi else [])+[hi]
            q=mp.quad(integrand,breaks)*mp.exp(-shift*shift/2)
            return q+mp.erfc(b/mp.sqrt(2))/2
        return [float(moment(p)) for p in (1,2)]


def main():
    started=time.perf_counter()
    cases=[(0,0),(-.5,0),(.5,0),(0,1e-150),(.123,1e-12),(-.5,1e-8),
           (.5,1e-8),(-.50000001,1e-8),(-.49999999,1e-8),
           (-.25,.05),(-.4,.25),(-.15,.1),(0,.25),(0,1),(0,1000),
           (-1,.1),(-2,.1),(-3,.1),(-4,.1),(-4.2,.1),(4,.1),
           (1e5,1e5),(-1e5,1e5)]
    rng=np.random.default_rng(SEED)
    cases += list(zip(rng.uniform(-2,2,160),10**rng.uniform(-3,1,160)))
    checks=[]
    quad_warnings=[]
    for mu,sd in cases:
        got=gaussian_saturation_moments(mu,sd)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            expected=[independent_quad(mu,sd,p) for p in (1,2)]
        if caught:
            quad_warnings.append(dict(mean=float(mu),sd=float(sd),
                                      warnings=[str(x.message) for x in caught]))
        errors=[abs(float(got.mean)-expected[0]),abs(float(got.second_moment)-expected[1])]
        assert max(errors)<3e-12,(mu,sd,errors)
        if mu < -.5 and expected[0] > 1e-308:
            assert errors[0]/expected[0]<2e-8,(mu,sd,errors,expected)
        assert np.isfinite(got.variance) and 0 <= got.variance <= .250000000001
        checks.append(dict(mean=float(mu),sd=float(sd),fraction_mean=float(got.mean),
                           fraction_second=float(got.second_moment),variance=float(got.variance),
                           quad_mean=expected[0],quad_second=expected[1],max_abs_error=max(errors)))
    tiny=[]
    for mu,sd in [(0,1e-150),(.123,1e-12),(-.27,1e-7)]:
        expected,variance=tiny_interior_polynomial(mu,sd)
        got=gaussian_saturation_moments(mu,sd)
        relative=abs(float(got.variance)/variance-1)
        assert relative<3e-12,(mu,sd,relative)
        tiny.append(dict(mean=mu,sd=sd,expected_variance=variance,
                         actual_variance=float(got.variance),relative_error=relative))
    high_precision=[]
    for mu,sd in [(-.5,1e-8),(-.50000001,1e-8),(-.49999999,1e-8),(-4.2,.1)]:
        expected=high_precision_tail(mu,sd)
        got=gaussian_saturation_moments(mu,sd)
        actual=[float(got.mean),float(got.second_moment)]
        rel=[abs(a/e-1) if e > 1e-308 else None for a,e in zip(actual,expected)]
        assert max(x for x in rel if x is not None)<3e-11,(mu,sd,rel)
        high_precision.append(dict(mean=mu,sd=sd,expected=expected,actual=actual,relative_errors=rel))
    # Correlated, rank-deficient 8-dimensional Gaussian mapped to 12 pixels.
    d,n=8,12
    root=rng.normal(size=(d,5))*.1
    cov=root@root.T
    matrix=rng.normal(size=(n,d))*.6
    latent_mean=rng.normal(size=d)*.2
    offset=rng.normal(size=n)*.15
    fm,fs=linear_gaussian_field_marginals(matrix,latent_mean,cov,offset=offset)
    delta=1.3+.105j
    mask=np.ones(n);mask[-2:]=0
    material=gaussian_material_moments(fm,fs,delta,mask=mask)
    count=600000
    fields=offset+(rng.normal(size=(count,5))@root.T+latent_mean)@matrix.T
    fractions=saturation(fields/.5)*mask
    values=.1+.04j+delta*fractions
    mc_mean=values.mean(axis=0)
    mc_var=np.mean(abs(values-mc_mean)**2,axis=0)
    mc_second=np.mean(abs(values)**2,axis=0)
    se=np.sqrt(material.variance/count)
    mean_z=abs(mc_mean-material.mean)/np.maximum(se,1e-30)
    squared=abs(values-material.mean)**2
    var_se=squared.std(axis=0)/np.sqrt(count)
    var_z=abs(mc_var-material.variance)/np.maximum(var_se,1e-30)
    assert np.max(mean_z[:10])<6 and np.max(var_z[:10])<6
    assert np.max(abs(material.second_abs_moment-(abs(material.mean)**2+material.variance)))==0
    assert np.all(material.variance[-2:]==0)
    tail_mc=[]
    for mu,sd in [(-.5,1e-8),(-.50000001,1e-8),(-.8,.15)]:
        phi=saturation((mu+sd*rng.normal(size=count))/.5)
        target=gaussian_saturation_moments(mu,sd)
        zscores=[]
        for power,expected in [(1,float(target.mean)),(2,float(target.second_moment))]:
            samples=phi**power
            zscore=abs(samples.mean()-expected)/(samples.std()/np.sqrt(count))
            assert zscore<6,(mu,sd,power,zscore)
            zscores.append(float(zscore))
        tail_mc.append(dict(mean=mu,sd=sd,draws=count,moment_standard_errors=zscores))
    # Direct independent quadrature for complex modulus variance.
    complex_checks=[]
    for mu,sd in [(-.2,.15),(.4,.3),(0,1.)]:
        result=gaussian_material_moments(mu,sd,delta)
        cm=complex(result.mean)
        qv=quad(lambda z: abs(.1+.04j+delta*float(saturation((mu+sd*z)/.5))-cm)**2*
                np.exp(-z*z/2)/np.sqrt(2*np.pi),-12,12,
                points=[(-.5-mu)/sd,(.5-mu)/sd],epsabs=1e-13,epsrel=2e-12)[0]
        assert abs(qv-float(result.variance))<3e-12
        complex_checks.append(dict(mean=mu,sd=sd,variance=float(result.variance),quad_variance=qv))
    # A concrete posterior-like dependence example: bounded property U=Phi(Z).
    dep=quad(lambda z: (.8+ndtr(z))*float(saturation(z/.5))*
             np.exp(-z*z/2)/np.sqrt(2*np.pi),-12,12,
             points=[-.5,.5],epsabs=1e-13)[0]
    independence_wrong=1.3*.5
    assert abs(dep-independence_wrong)>.1
    invalid_checks=0
    for args in [(0,-1),(np.nan,1),(0,np.inf)]:
        try: gaussian_saturation_moments(*args)
        except ValueError: invalid_checks+=1
        else: raise AssertionError('invalid Gaussian accepted')
    for bad in [np.diag([1.,-1.]),np.array([[1.,.2],[0.,1.]])]:
        try: linear_gaussian_field_marginals(np.eye(2),np.zeros(2),bad)
        except ValueError: invalid_checks+=1
        else: raise AssertionError('invalid covariance accepted')
    stress=[]
    for mu,sd,h in [(1e308,1e308,1e-308),(-1e308,1e308,1e-308),
                    (0,1e308,.5),(0,1e307,.5),(0,5e-324,.5),(0,0,1e-308),
                    (1e308,0,1e-308)]:
        result=gaussian_saturation_moments(mu,sd,halfwidth=h)
        assert np.all(np.isfinite([result.mean,result.second_moment,result.variance]))
        stress.append(dict(mean=mu,sd=sd,halfwidth=h,result=[float(result.mean),float(result.second_moment),float(result.variance)]))
    broad=gaussian_saturation_moments(np.array([[-.2],[.1]]),np.array([0,.1,.7]),
                                     halfwidth=np.array([.3,.5,1.1]))
    assert broad.mean.shape==(2,3)
    for i,mu in enumerate([-.2,.1]):
        for j,(sd,h) in enumerate(zip([0,.1,.7],[.3,.5,1.1])):
            reference=gaussian_saturation_moments(mu*.5/h,sd*.5/h)
            assert abs(broad.mean[i,j]-reference.mean)<1e-15
    summary=dict(status='PASS_GAUSSIAN_REFERENCE_MOMENTS_ONLY',seed=SEED,
                 maxwell_solves=0,maxwell_rhs=0,training_runs=0,quadrature_cases=len(checks),
                 maximum_quadrature_absolute_error=max(x['max_abs_error'] for x in checks),
                 independent_quadrature_warnings=quad_warnings,
                 tiny_variance_checks=tiny,high_precision_tail_checks=high_precision,tail_monte_carlo=tail_mc,
                 correlated_gaussian_mc=dict(draws=count,latent_dimension=d,
                 covariance_rank=5,pixels=n,max_mean_standard_errors=float(np.max(mean_z[:10])),
                 max_variance_standard_errors=float(np.max(var_z[:10])),
                 max_abs_second_moment_difference=float(np.max(abs(mc_second-material.second_abs_moment)))),
                 complex_variance_checks=complex_checks,invalid_input_checks=invalid_checks,
                 broadcast_and_nondefault_halfwidth='PASS',
                 dependent_property_counterexample=dict(correct_mean_real=dep,independence_wrong_mean_real=independence_wrong),
                 floating_point_stress=stress,elapsed_seconds=time.perf_counter()-started,
                 limitations=['Supplied Gaussian reference only; no posterior-approximation validation.',
                              'No full spatial covariance or local-region mass computation.',
                              'Quadrature is floating point, not an interval certificate.',
                              'Subnormal-underflow moments necessarily round to zero.'])
    (OUT/'quadrature_cases.json').write_text(json.dumps(checks,indent=2)+'\n')
    (OUT/'results.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    main()
