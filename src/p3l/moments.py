"""P1 decoder moments under a supplied Gaussian conditional reference.

The field is Gaussian; inclusion properties are fixed *conditional on the
outer state*. This does not replace a Maxwell posterior, account for a local
hard truncation, or assume posterior independence of properties and field.
Only pixel marginal moments are returned, never a spatial covariance.

See Research P3-L report section 4.5 and probability appendix (10)--(13).
No physics module is imported and no Maxwell solve/RHS is performed.
"""
from dataclasses import dataclass
from math import comb, exp, pi, sqrt

import numpy as np
from scipy.integrate import quad
from scipy.special import ndtr


@dataclass(frozen=True)
class SaturationMoments:
    mean: np.ndarray
    second_moment: np.ndarray
    variance: np.ndarray


@dataclass(frozen=True)
class MaterialMoments:
    mean: np.ndarray
    variance: np.ndarray
    second_abs_moment: np.ndarray
    fraction: SaturationMoments


def _poly(t):
    return t*t*t*(10+t*(-15+6*t))


def _recurrence(mu, sd):
    """Long-double truncated normal recurrence in a well-conditioned region."""
    # Here the field is divided by halfwidth; transition is [-1,1].
    a, b = (-1-mu)/sd, (1-mu)/sd
    ld = np.longdouble
    t = np.zeros(11, dtype=ld)
    t[0] = ndtr(b)-ndtr(a)
    pa, pb = ld(exp(-a*a/2)/sqrt(2*pi)), ld(exp(-b*b/2)/sqrt(2*pi))
    t[1] = pa-pb
    for k in range(2, 11):
        t[k] = ld(a)**(k-1)*pa-ld(b)**(k-1)*pb+(k-1)*t[k-2]
    raw = np.array([sum(ld(comb(n,k))*ld(mu)**(n-k)*ld(sd)**k*t[k]
                        for k in range(n+1)) for n in range(11)], dtype=ld)
    c = np.array([.5, 15/16, 0, -5/8, 0, 3/16], dtype=ld)
    m = ld(ndtr(-b))+c@raw[:6]
    s = ld(ndtr(-b))+np.convolve(c,c)@raw
    v = s-m*m
    return float(m), float(s), float(v)


def _quadrature(mu, sd):
    """Positive, scaled Gaussian integration; mu <= 0, halfwidth = 1.

    Centered polynomial differences retain tiny variances even if mean and
    second moment round to their deterministic limits. Tail scaling preserves
    representable rare-event moments instead of subtracting nearly equal CDFs.
    The ±38 Gaussian integration cut has absolute tail mass < 6e-316.
    """
    with np.errstate(over='ignore',divide='ignore'):
        a, b = (-1-mu)/sd, (1-mu)/sd
    if a >= 38:
        return 0., 0., 0.
    tail = float(ndtr(-b))
    norm = 1/sqrt(2*pi)
    if mu <= -1:
        # z=a+r/(a+1), t=sd*r/[2(a+1)] avoids narrow rare-event peaks.
        rate = a+1
        scale = min(1., (sd/rate)/2)
        limit = min((b-a)*rate, (38-a)*rate, 100.)
        if limit <= 0:
            return tail, tail, tail*(1-tail)
        prefactor = exp(-a*a/2)*norm/rate
        def integrand(r, power):
            t = (sd/rate)*r/2
            # p(t)/scale**3, evaluated without forming scale**3 first.
            pscaled = (t/scale)**3*(10+t*(-15+6*t))
            return pscaled**power*exp(-a*r/rate-.5*(r/rate)**2)
        j1 = quad(integrand, 0, limit, args=(1,), epsabs=0., epsrel=3e-13,
                  limit=150)[0]
        j2 = quad(integrand, 0, limit, args=(2,), epsabs=0., epsrel=3e-13,
                  limit=150)[0]
        scale3 = scale**3
        m = tail+(prefactor*j1)*scale3
        s = tail+((prefactor*j2)*scale3)*scale3
        return m, s, max(0., s-m*m)
    # Interior baseline; compute p(t0+d)-p(t0) algebraically, avoiding
    # subtraction of two similar fractions for nearly deterministic fields.
    t0 = (mu+1)/2
    base = _poly(t0)
    coeff = [30*t0*t0*(1-t0)**2,
             30*t0*(1-3*t0+2*t0*t0),
             10-60*t0+60*t0*t0, -15+30*t0, 6.]
    # At a seam p' vanishes quadratically; scaling by SD alone would leave
    # an integrand so tiny that an absolute quadrature tolerance loses its
    # relative accuracy. Scale by the local polynomial response instead.
    scale = (1. if sd >= 2 else
             min(1., max(abs(c)*(sd/2)**(k+1) for k,c in enumerate(coeff))))
    if scale == 0:
        return base, base*base, 0.
    def centered(z):
        d = sd*z/2
        return (d/scale)*(coeff[0]+d*(coeff[1]+d*(coeff[2]+d*(coeff[3]+d*coeff[4]))))
    lo, hi = max(-38., a), min(38., b)
    # Splitting at zero keeps adaptive quadrature effective for long intervals.
    cuts = sorted(set([lo, hi]+([0.] if lo < 0 < hi else [])))
    def integral(power):
        return sum(quad(lambda z: centered(z)**power*exp(-z*z/2)*norm,
                        left, right, epsabs=2e-14, epsrel=3e-13, limit=150)[0]
                   for left, right in zip(cuts[:-1],cuts[1:]))
    pleft = float(ndtr(a))
    dleft = -base/scale if pleft else 0.
    dright = (1-base)/scale if tail else 0.
    ed = integral(1)+pleft*dleft+tail*dright
    ed2 = integral(2)+pleft*dleft*dleft+tail*dright*dright
    m = base+scale*ed
    s = base*base+2*base*scale*ed+(scale*ed2)*scale
    v = (scale*max(0., ed2-ed*ed))*scale
    return m, s, v


def gaussian_saturation_moments(mean, sd, *, halfwidth=.5):
    """E[phi(F)], E[phi(F)^2], Var(phi(F)) for F ~ N(mean, sd²).

    ``phi(F)`` is exactly ``p1_material.prior.saturation(F/halfwidth)``.
    Inputs broadcast; SD must be nonnegative and all inputs finite. SD=0 is
    evaluated deterministically. An analytic truncated-moment recurrence is
    used where conditioned well; scaled quadrature handles difficult tails,
    broad Gaussians and near-zero SD. Numerical integration is not an interval
    certificate. Values below float64's subnormal range necessarily underflow.
    """
    mean, sd, width = np.broadcast_arrays(np.asarray(mean,float),
                                          np.asarray(sd,float),
                                          np.asarray(halfwidth,float))
    if (not np.all(np.isfinite(mean)) or not np.all(np.isfinite(sd)) or
            not np.all(np.isfinite(width)) or np.any(sd < 0) or np.any(width <= 0)):
        raise ValueError('finite mean, nonnegative finite SD, positive finite halfwidth required')
    m, s, v = (np.empty(mean.shape,float) for _ in range(3))
    for index in np.ndindex(mean.shape):
        value, spread, h = mean[index], sd[index], width[index]
        if spread == 0:
            p = 0. if value <= -h else 1. if value >= h else float(_poly((value/h+1)/2))
            m[index], s[index], v[index] = p, p*p, 0.
            continue
        with np.errstate(over='ignore', invalid='ignore', divide='ignore'):
            mu, sigma = -abs(value/h), spread/h
        if not np.isfinite(sigma) or not np.isfinite(mu):
            # Normalize in extended precision if a finite input ratio overflowed.
            with np.errstate(over='ignore',invalid='ignore',divide='ignore'):
                mu_ld = -abs(np.longdouble(value)/np.longdouble(h))
                sd_ld = np.longdouble(spread)/np.longdouble(h)
            if sd_ld > np.finfo(float).max or abs(mu_ld) > np.finfo(float).max:
                # The transition width is negligible on this scale. Its
                # probability is below representable relative precision.
                with np.errstate(over='ignore',divide='ignore'):
                    ratio = float(np.longdouble(value)/np.longdouble(spread))
                p = float(ndtr(ratio))
                m[index],s[index],v[index] = p,p,p*(1-p)
                continue
            mu,sigma = float(mu_ld),float(sd_ld)
        if sigma < 1e-162:
            # The Lipschitz bound Var(phi) <= (15*sigma/16)^2 is
            # below float64's smallest subnormal. Avoid underflow in a
            # normalized quadrature integrand for an unrepresentable variance.
            p = float(_poly(np.clip((mu+1)/2,0,1)))
            mm,ss,vv = p,p*p,0.
        elif -.8 <= mu <= 0 and .1 <= sigma <= .5:
            mm,ss,vv = _recurrence(mu,sigma)
        else:
            mm,ss,vv = _quadrature(mu,sigma)
        # Reflection preserves tiny variance near the upper plateau; never
        # recover variance by subtracting E[phi]^2 from E[phi²].
        if value > 0:
            mm,ss = 1-mm,1-2*mm+ss
        m[index],s[index],v[index] = np.clip(mm,0,1),np.clip(ss,0,1),max(0.,vv)
    return SaturationMoments(m,s,v)


def gaussian_material_moments(mean, sd, inclusion_increment, *, host=.1+.04j,
                              mask=1., halfwidth=.5):
    """Complex pixel moments conditional on fixed properties and a Gaussian field.

    ``inclusion_increment`` may broadcast over outer states; every entry is
    fixed for its supplied field distribution. To integrate uncertain properties
    dependent on the field, call at each conditioned outer state and combine
    with the law of total variance. No product-of-expectations assumption is
    made. ``mask`` is deterministic in [0,1] (P1 uses a Boolean ROI).

    ``variance`` = E[|chi-E chi|²], not E[(chi-E chi)²]. ``second_abs_moment``
    includes the host/inclusion cross term. Spatial covariance is not returned.
    """
    fraction = gaussian_saturation_moments(mean,sd,halfwidth=halfwidth)
    delta, background, roi, fm, fs, fv = np.broadcast_arrays(
        np.asarray(inclusion_increment,complex),np.asarray(host,complex),
        np.asarray(mask,float),fraction.mean,fraction.second_moment,fraction.variance)
    if (not np.all(np.isfinite(delta)) or not np.all(np.isfinite(background)) or
            not np.all(np.isfinite(roi)) or np.any(roi < 0) or np.any(roi > 1)):
        raise ValueError('finite fixed material properties and mask in [0,1] required')
    fm,fs,fv = roi*fm,roi*roi*fs,roi*roi*fv
    chi = background+delta*fm
    var = np.abs(delta)**2*fv
    second = np.abs(chi)**2+var
    return MaterialMoments(chi,var,second,SaturationMoments(fm,fs,fv))


def linear_gaussian_field_marginals(matrix, mean, covariance, *, offset=0.):
    """Marginal means/SDs of offset + matrix @ Z, including correlated Z.

    Singular PSD covariance is supported. Eigenvalues negative only at floating
    roundoff are clipped; a materially indefinite or asymmetric matrix raises.
    The returned scalar marginals do not assert that distinct pixels are
    independent. No field or material spatial covariance is claimed.
    """
    matrix, mean, covariance = (np.asarray(x,float) for x in (matrix,mean,covariance))
    if (matrix.ndim != 2 or matrix.shape[1] == 0 or mean.shape != (matrix.shape[1],) or
            covariance.shape != (matrix.shape[1],matrix.shape[1]) or
            not all(np.all(np.isfinite(x)) for x in (matrix,mean,covariance))):
        raise ValueError('finite matrix (pixels,d), mean (d,), covariance (d,d) required')
    tolerance = 64*np.finfo(float).eps*max(np.finfo(float).tiny,
                                         np.linalg.norm(covariance,ord=2))
    if not np.allclose(covariance,covariance.T,atol=tolerance,rtol=0.):
        raise ValueError('covariance must be symmetric')
    eigenvalues,eigenvectors = np.linalg.eigh((covariance+covariance.T)/2)
    if eigenvalues.min(initial=0.) < -tolerance:
        raise ValueError('covariance must be positive semidefinite')
    root = eigenvectors*np.sqrt(np.maximum(eigenvalues,0.))
    field_mean = np.asarray(offset,float)+matrix@mean
    if field_mean.shape != (matrix.shape[0],) or not np.all(np.isfinite(field_mean)):
        raise ValueError('offset must broadcast to the pixel axis with finite result')
    field_sd = np.linalg.norm(matrix@root,axis=1)
    return field_mean,field_sd
