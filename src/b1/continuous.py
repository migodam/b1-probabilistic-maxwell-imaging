"""Continuous posterior summaries from already-paid polynomial responses.

This module has no physical solver or budget-ledger interface.  Its numerical
target is explicit: a frozen response polynomial, the original continuous
material prior, and the observation's fixed noise and gain model.
"""
from __future__ import annotations

from functools import lru_cache
from pathlib import Path
import heapq
import time
import numpy as np
from numpy.polynomial.legendre import leggauss, legvander2d
from scipy.optimize import brentq
from scipy.special import logsumexp
from scipy.stats import chi2
from .inference import log_likelihood


@lru_cache(None)
def unit_rule(order):
    x, w = leggauss(order)
    return (x + 1.) / 2., w / 2.


def coherence(mean, covariance, quantiles, widths):
    mean = np.asarray(mean); variance = np.diag(np.asarray(covariance))
    q = np.asarray(quantiles); widths = np.asarray(widths)
    lower = .05 * (np.maximum(mean-q[:, 0], 0.)**2 +
                   np.maximum(q[:, 2]-mean, 0.)**2)
    tol = 1e-12 * widths**2
    return dict(necessary_variance_lower_bound=lower.tolist(),
                reported_variance=variance.tolist(), floating_tolerance=tol.tolist(),
                coordinate_pass=(variance + tol >= lower).tolist(),
                passed=bool(np.all(variance + tol >= lower)),
                sufficient_for_accuracy=False)


class FrozenDensity:
    """A solver-free density over the original two-dimensional prior box."""
    def __init__(self, coefficients, bounds, observation, stats=None, label='poly7'):
        self.coeff = np.asarray(coefficients, float)
        self.bounds = np.asarray(bounds, float)
        self.lo = self.bounds[:, 0]; self.width = np.diff(self.bounds, axis=1)[:, 0]
        self.order = int(round(np.sqrt(self.coeff.shape[1])))
        self.y = np.asarray(observation['data'], float)
        self.sigma = float(observation['noise_real_sd'])
        self.gain_sd = float(observation['gain_sd'])
        self.prior = np.asarray(observation['topology_prior'], float)
        self.stats = stats if stats is not None else dict()
        self.label = label; self.shifts = []
        if self.coeff.shape != (3, self.order**2, len(self.y)):
            raise ValueError('frozen response dimensions')
        if self.sigma <= 0 or self.gain_sd < 0 or np.any(self.prior <= 0):
            raise ValueError('invalid declared likelihood/prior metadata')
        if abs(self.prior.sum()-1.) > 1e-12:
            raise ValueError('topology prior must already be normalized')
        z, _ = unit_rule(25); a, b = np.meshgrid(z, z, indexing='ij')
        points = np.column_stack((a.ravel(), b.ravel()))
        for h in range(3):
            ll, _ = self.loglik(h, points, 'density_setup')
            self.shifts.append(float(ll.max()))

    def response(self, h, z, kind):
        z = np.atleast_2d(z)
        if np.any(z < -1e-14) or np.any(z > 1+1e-14):
            raise ValueError('outside frozen prior support')
        v = legvander2d(2*z[:, 0]-1, 2*z[:, 1]-1,
                       [self.order-1, self.order-1]).reshape(-1, self.order**2)
        response = v @ self.coeff[int(h)]
        self.stats['polynomial_response_evaluations'] = self.stats.get('polynomial_response_evaluations', 0) + len(z)
        counts = self.stats.setdefault('by_kind', {})
        counts[kind] = counts.get(kind, 0) + len(z)
        return response

    def loglik(self, h, z, kind):
        response=self.response(h,z,kind)
        ll, mah = log_likelihood(response, self.y, self.sigma, self.gain_sd)
        if not np.all(np.isfinite(ll)):
            raise FloatingPointError('nonfinite declared density')
        return ll, mah

    def quadrature_bank(self,order):
        x,w=unit_rule(order);a,b=np.meshgrid(x,x,indexing='ij')
        z=np.column_stack((a.ravel(),b.ravel()))
        response=np.concatenate([self.response(h,z,'CLI_saved_bank') for h in range(3)])
        return dict(H=np.repeat(np.arange(3),len(z)),theta=np.tile(self.lo+z*self.width,(3,1)),
                    weights=np.tile(np.outer(w,w).ravel(),3),F=response,bounds=self.bounds)

    def rectangle(self, h, box, order, kind='tensor_integration'):
        box = np.asarray(box, float); lengths = box[:, 1]-box[:, 0]
        if np.any(lengths <= 0):
            return np.zeros(8)
        x, w = unit_rule(order)
        a = box[0, 0] + lengths[0]*x; b = box[1, 0] + lengths[1]*x
        aa, bb = np.meshgrid(a, b, indexing='ij')
        z = np.column_stack((aa.ravel(), bb.ravel()))
        ll, mah = self.loglik(h, z, kind)
        weight = np.outer(w, w).ravel()*np.prod(lengths)*np.exp(ll-self.shifts[h])
        features = np.column_stack((np.ones(len(z)), z[:, 0], z[:, 1],
                                    z[:, 0]**2, z[:, 1]**2, z[:, 0]*z[:, 1],
                                    mah, chi2.sf(mah, len(self.y))))
        return weight @ features

    def adaptive(self, h, box, relative_tolerance=2e-8, max_cells=1024):
        """Independent partitioned 5/9-point quadrature, not a certificate."""
        box = np.asarray(box, float)
        if np.any(box[:, 1] <= box[:, 0]):
            return np.zeros(8), np.zeros(8), dict(converged=True, cells=0)
        heap=[]; counter=0; value=np.zeros(8); error=np.zeros(8)
        def cell(rect):
            nonlocal counter
            a=self.rectangle(h, rect, 5, 'adaptive_low')
            b=self.rectangle(h, rect, 9, 'adaptive_high')
            e=np.abs(b-a);counter+=1
            return (-float(e[:6].max()), counter, rect, b, e)
        grids=[np.linspace(*b, 5) for b in box]
        for i in range(4):
            for j in range(4):
                c=cell(np.array([[grids[0][i],grids[0][i+1]],
                                 [grids[1][j],grids[1][j+1]]]))
                heapq.heappush(heap,c);value+=c[3];error+=c[4]
        def converged():
            return error[:6].max() <= 1e-13 + relative_tolerance*abs(value[0])
        while not converged() and len(heap)<max_cells:
            _,_,rect,b,e=heapq.heappop(heap);value-=b;error-=e
            axis=int(np.argmax(rect[:,1]-rect[:,0]));mid=float(rect[axis].mean())
            left=rect.copy();right=rect.copy();left[axis,1]=mid;right[axis,0]=mid
            for part in [left,right]:
                c=cell(part);heapq.heappush(heap,c);value+=c[3];error+=c[4]
        error=np.maximum(error,0.)
        self.stats['adaptive_integrals']=self.stats.get('adaptive_integrals',0)+1
        self.stats['adaptive_unresolved']=self.stats.get('adaptive_unresolved',0)+int(not converged())
        return value,error,dict(converged=bool(converged()),cells=len(heap),
                               error_estimate_is_certified=False)

    def moments(self, h, integral):
        norm=float(integral[0])
        if norm <= 0 or not np.isfinite(norm):
            raise FloatingPointError('density integral not positive finite')
        m=integral[1:3]/norm
        second=np.array([[integral[3],integral[5]], [integral[5],integral[4]]])/norm
        cov=(second-np.outer(m,m))*np.outer(self.width,self.width)
        if np.min(np.diag(cov)) < -1e-12*np.max(self.width**2):
            raise FloatingPointError('negative material variance')
        return dict(H=int(h),mean=(self.lo+self.width*m).tolist(),
                    covariance=cov.tolist(),sd=np.sqrt(np.maximum(np.diag(cov),0)).tolist(),
                    log_evidence=float(self.shifts[h]+np.log(norm)),
                    discrepancy_tail=float(integral[7]/norm),
                    expected_mahalanobis=float(integral[6]/norm))

    def summary(self, order=None, grid=None, adaptive=False):
        start=time.perf_counter();fullbox=np.array([[0.,1.],[0.,1.]])
        grid=np.linspace(0.,1.,129 if not adaptive else 33) if grid is None else np.asarray(grid)
        cond=[];curves=[];integrals=[];errors=[];statuses=[];cdf_cache={}
        def integrate(h,box):
            if adaptive:
                v,e,s=self.adaptive(h,box);statuses.append(s);return v,e
            return self.rectangle(h,box,int(order)),np.zeros(8)
        for h in range(3):
            v,e=integrate(h,fullbox);integrals.append(v);errors.append(e);cond.append(self.moments(h,v))
        logz=np.array([r['log_evidence'] for r in cond]);total=logsumexp(np.log(self.prior)+logz)
        probabilities=np.exp(np.log(self.prior)+logz-total)
        def cdf(h,axis,cut):
            if cut <= 0:return 0.
            if cut >= 1:return 1.
            key=(h,axis,float(cut))
            if key not in cdf_cache:
                box=fullbox.copy();box[axis,1]=cut;part,err=integrate(h,box)
                cdf_cache[key]=float(part[0]/integrals[h][0])
            return cdf_cache[key]
        def quantile(fun,p):return float(brentq(lambda t:fun(t)-p,0.,1.,xtol=2e-11))
        for h in range(3):
            curve=np.array([[cdf(h,j,float(t)) for t in grid] for j in range(2)])
            curves.append(curve)
            q=np.array([[quantile(lambda t:cdf(h,j,t),p) for p in [.05,.5,.95]] for j in range(2)])
            q=self.lo[:,None]+self.width[:,None]*q
            cond[h]['quantiles_05_50_95']=q.tolist()
            cond[h]['coherence']=coherence(cond[h]['mean'],cond[h]['covariance'],q,self.width)
            cond[h]['summary_measure']='same continuous polynomial-likelihood density'
            cond[h]['cdf_monotone']=bool(np.min(np.diff(curve,axis=1))>=-1e-10)
            cond[h]['cdf_range_valid']=bool(curve.min()>=-1e-10 and curve.max()<=1+1e-10)
        curves=np.asarray(curves);mcdf=np.sum(probabilities[:,None,None]*curves,axis=0)
        means=np.asarray([r['mean'] for r in cond]);mu=probabilities@means
        covariance=sum(probabilities[h]*(np.asarray(cond[h]['covariance'])+
                       np.outer(means[h]-mu,means[h]-mu)) for h in range(3))
        mq=np.array([[quantile(lambda t:sum(probabilities[h]*cdf(h,j,t) for h in range(3)),p)
                      for p in [.05,.5,.95]] for j in range(2)])
        mq=self.lo[:,None]+self.width[:,None]*mq
        marginal=dict(mean=mu.tolist(),covariance=covariance.tolist(),
                      sd=np.sqrt(np.maximum(np.diag(covariance),0)).tolist(),
                      quantiles_05_50_95=mq.tolist(),
                      coherence=coherence(mu,covariance,mq,self.width),
                      cdf_monotone=bool(np.min(np.diff(mcdf,axis=1))>=-1e-10),
                      summary_measure='same continuous polynomial-likelihood density')
        result=dict(topology_probabilities=probabilities.tolist(),log_evidence=logz.tolist(),
                    log_total_evidence=float(total),log_bayes_factors_vs_solid=(logz-logz[0]).tolist(),
                    conditional_material=cond,marginal_material=marginal,
                    posterior_predictive_tail_probability=float(probabilities@np.array([r['discrepancy_tail'] for r in cond])),
                    posterior_expected_whitened_residual_per_component=float(probabilities@np.array([r['expected_mahalanobis'] for r in cond])/len(self.y)),
                    target=f'B1-E0-continuous-{self.label}', integration_order=order,
                    method='partitioned_adaptive_5_9' if adaptive else 'tensor_Gauss_direct_CDF',
                    certified_bound=False,wall_seconds=time.perf_counter()-start,
                    adaptive_converged=all(s['converged'] for s in statuses),
                    estimated_integral_errors=[e.tolist() for e in errors])
        return result,dict(normalized_grid=grid,conditional_cdf=curves,marginal_cdf=mcdf,
                           physical_grid=self.lo[:,None]+self.width[:,None]*grid)


def compare_summaries(candidate,reference,ca,ra,widths):
    widths=np.asarray(widths); rows=[]
    a=list(candidate['conditional_material'])+[candidate['marginal_material']]
    b=list(reference['conditional_material'])+[reference['marginal_material']]
    for i,(x,y) in enumerate(zip(a,b)):
        mean=np.abs(np.asarray(x['mean'])-y['mean'])
        sd=np.abs(np.asarray(x['sd'])-y['sd'])
        q=np.abs(np.asarray(x['quantiles_05_50_95'])-y['quantiles_05_50_95'])
        xc=ca['conditional_cdf'][i] if i<3 else ca['marginal_cdf']
        yc=ra['conditional_cdf'][i] if i<3 else ra['marginal_cdf']
        # Reference grids can differ; compare only on their common tested grid.
        xc=np.array([np.interp(ra['normalized_grid'],ca['normalized_grid'],r) for r in xc])
        cdf=np.max(np.abs(xc-yc),axis=1)
        sd_tol=np.maximum(.01*np.asarray(y['sd']),1e-6*widths)
        passed=bool(np.all(cdf<=.01) and np.all(mean<=.001*widths) and
                    np.all(q<=.001*widths[:,None]) and np.all(sd<=sd_tol))
        rows.append(dict(event=i if i<3 else 'material_marginal',cdf_tested_grid_max=cdf.tolist(),
                         mean_difference=mean.tolist(),sd_difference=sd.tolist(),
                         quantile_difference=q.tolist(),passed=passed))
    return dict(rows=rows,passed=all(r['passed'] for r in rows),
                max_event_mass_difference=float(np.max(np.abs(np.asarray(candidate['topology_probabilities'])-reference['topology_probabilities']))),
                max_log_BF_difference=float(np.max(np.abs(np.asarray(candidate['log_bayes_factors_vs_solid'])-reference['log_bayes_factors_vs_solid']))),
                max_log_evidence_difference=float(np.max(np.abs(np.asarray(candidate['log_evidence'])-reference['log_evidence']))),
                numerical_comparison_not_uniform_certificate=True)
