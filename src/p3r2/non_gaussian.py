"""Frozen nonlinear mixture proposals and 2-D Full-likelihood quadrature.

Pure array code: never imports a solver or campaign Book. Exact MH corrections
establish invariant targets only; finite-time mixing remains an empirical issue.
"""
from dataclasses import dataclass
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import solve_triangular
from scipy.special import gammaln, logsumexp, ndtr


def _logtarget(value):
    value = float(value)
    if np.isnan(value) or np.isposinf(value):
        raise ValueError('target must return finite value or -inf')
    return value


@dataclass
class StudentMixture:
    """Normalized multivariate Student-t mixture with common degrees of freedom."""
    centers: np.ndarray
    cholesky: np.ndarray
    weights: np.ndarray
    df: float = 7.

    def __post_init__(self):
        self.centers = np.asarray(self.centers, float).copy()
        self.weights = np.asarray(self.weights, float).copy()
        k, d = self.centers.shape
        self.cholesky = np.broadcast_to(np.asarray(self.cholesky, float), (k,d,d)).copy()
        if self.weights.shape != (k,) or not np.isfinite(self.centers).all() or not np.isfinite(self.cholesky).all():
            raise ValueError('finite matching mixture arrays required')
        if self.df <= 0 or not np.isfinite(self.df) or np.any(self.weights <= 0):
            raise ValueError('positive mixture weights and degrees of freedom required')
        if any(not np.allclose(L, np.tril(L)) or np.any(np.diag(L) <= 0) for L in self.cholesky):
            raise ValueError('positive lower triangular Cholesky factors required')
        self.weights /= self.weights.sum()
        self.logdet = np.log(np.diagonal(self.cholesky, axis1=1, axis2=2)).sum(axis=1)
        self.constant = gammaln((self.df+d)/2)-gammaln(self.df/2)-d/2*np.log(self.df*np.pi)

    def sample(self, rng):
        k = int(rng.choice(len(self.weights), p=self.weights))
        t = rng.normal(size=self.centers.shape[1])/np.sqrt(rng.chisquare(self.df)/self.df)
        return self.centers[k]+self.cholesky[k]@t, k

    def logpdf(self, x):
        x = np.asarray(x, float)
        terms = []
        for center, L, ld, w in zip(self.centers,self.cholesky,self.logdet,self.weights):
            t = solve_triangular(L,x-center,lower=True)
            terms.append(np.log(w)+self.constant-ld-(self.df+len(x))/2*np.log1p(t@t/self.df))
        return float(logsumexp(terms))


@dataclass
class NonlinearChartMixture:
    """Triangular Student-t charts on Gaussian64 x probit2.

    In chart k, x_eta=c_eta+s*t_eta and
    x_z=c_z+L*t_z+B*(Phi(x_eta)-Phi(c_eta)). The nonlinear shear
    has unit determinant; total log determinant is logdet L + sum(log s).
    Each chart has full support and a directly evaluable normalized density.
    """
    centers: np.ndarray
    field_cholesky: np.ndarray
    property_scales: np.ndarray
    compensation: np.ndarray
    weights: np.ndarray
    df: float = 7.

    def __post_init__(self):
        self.centers=np.asarray(self.centers,float).copy()
        if self.centers.ndim != 2 or self.centers.shape[1] != 66:
            raise ValueError('66-dimensional chart centers required')
        k=len(self.centers)
        self.field_cholesky=np.broadcast_to(self.field_cholesky,(k,64,64)).copy()
        self.property_scales=np.broadcast_to(self.property_scales,(k,2)).copy()
        self.compensation=np.broadcast_to(self.compensation,(k,64,2)).copy()
        self.weights=np.asarray(self.weights,float).copy()
        if self.weights.shape != (k,) or np.any(self.weights<=0) or np.any(self.property_scales<=0):
            raise ValueError('positive chart weights/scales required')
        if not all(np.isfinite(a).all() for a in (self.centers,self.field_cholesky,self.property_scales,self.compensation,self.weights)):
            raise ValueError('finite chart arrays required')
        for L in self.field_cholesky:
            if not np.allclose(L,np.tril(L)) or np.any(np.diag(L)<=0):
                raise ValueError('positive lower triangular field factors required')
        if not np.isfinite(self.df) or self.df <= 0: raise ValueError('positive df required')
        self.weights/=self.weights.sum()
        self.logdet=np.log(np.diagonal(self.field_cholesky,axis1=1,axis2=2)).sum(axis=1)+np.log(self.property_scales).sum(axis=1)
        self.constant=gammaln((self.df+66)/2)-gammaln(self.df/2)-33*np.log(self.df*np.pi)

    def sample(self,rng):
        k=int(rng.choice(len(self.weights),p=self.weights));c=self.centers[k]
        t=rng.normal(size=66)/np.sqrt(rng.chisquare(self.df)/self.df)
        eta=c[64:]+self.property_scales[k]*t[64:]
        z=c[:64]+self.field_cholesky[k]@t[:64]+self.compensation[k]@(ndtr(eta)-ndtr(c[64:]))
        return np.r_[z,eta],k

    def logpdf(self,x):
        terms=[];x=np.asarray(x,float)
        for k,c in enumerate(self.centers):
            residual=x[:64]-c[:64]-self.compensation[k]@(ndtr(x[64:])-ndtr(c[64:]))
            t=np.r_[solve_triangular(self.field_cholesky[k],residual,lower=True),
                    (x[64:]-c[64:])/self.property_scales[k]]
            terms.append(np.log(self.weights[k])+self.constant-self.logdet[k]-(self.df+66)/2*np.log1p(t@t/self.df))
        return float(logsumexp(terms))


def independent_mh(logtarget, proposal, starts, *, rng, draws=1024, callback=None):
    """Frozen independent MH; retain every repeat and proposal correction.

    State IDs, ancestor iteration, proposal component, exact log targets and
    both normalized proposal densities are persisted through callback each step.
    """
    starts=np.asarray(starts,float)
    if starts.ndim != 2 or draws < 1: raise ValueError('chains x dimension starts and positive draws required')
    states=np.empty((len(starts),draws,starts.shape[1]));logs=np.empty((len(starts),draws))
    accept=np.zeros((len(starts),draws),bool);ancestors=np.empty((len(starts),draws),int)
    metadata=[]
    for chain,x0 in enumerate(starts):
        x=x0.copy();lp=_logtarget(logtarget(x));lq=proposal.logpdf(x);ancestor=-1
        if not np.isfinite(lp): raise ValueError('finite initial target required')
        if callback: callback(dict(event='CHAIN_INITIALIZED',chain=chain,x=x,logtarget=lp,logproposal=lq))
        for t in range(draws):
            xp,component=proposal.sample(rng);lpp=_logtarget(logtarget(xp));lqp=proposal.logpdf(xp)
            ratio=lpp-lp+lq-lqp
            yes=bool(np.log(max(rng.random(),np.finfo(float).tiny))<min(0.,ratio))
            row=dict(event='TRANSITION',chain=chain,draw=t,component=component,accepted=yes,
                     old_logtarget=lp,proposed_logtarget=lpp,old_logproposal=lq,
                     proposed_logproposal=lqp,log_mh_ratio=ratio,proposed_x=xp.copy())
            if yes: x,lp,lq=xp,lpp,lqp;ancestor=t
            states[chain,t]=x;logs[chain,t]=lp;accept[chain,t]=yes;ancestors[chain,t]=ancestor
            row.update(state_x=x.copy(),state_logtarget=lp,ancestor_draw=ancestor)
            metadata.append(row)
            if callback: callback(row)
    return dict(states=states,logtarget=logs,accepted=accept,ancestor_draw=ancestors,
                metadata=metadata,status='EXACT_INVARIANT_MH_FINITE_RUN_UNQUALIFIED',
                posterior_qualification=False)


def tensor_rule_2d(order, bounds=((0.,1.),(0.,1.))):
    nodes,weights=leggauss(order);axes=[];ws=[]
    for low,high in bounds:
        if not low<high: raise ValueError('nonempty bounds required')
        axes.append((high+low)/2+(high-low)/2*nodes);ws.append((high-low)/2*weights)
    a,b=np.meshgrid(*axes,indexing='ij');wa,wb=np.meshgrid(*ws,indexing='ij')
    return np.column_stack((a.ravel(),b.ravel())),(wa*wb).ravel()


def quadrature_moments(logtarget, order, *, bounds=((0.,1.),(0.,1.)), observable=None):
    """Original scalar target evaluated at every node; no surrogate likelihood."""
    nodes,qw=tensor_rule_2d(order,bounds);lp=[];obs=[]
    for x in nodes:
        result=logtarget(x)
        if isinstance(result,tuple): value,material=result
        else: value=result;material=x if observable is None else observable(x)
        lp.append(_logtarget(value));obs.append(np.asarray(material))
    lp=np.asarray(lp);obs=np.asarray(obs);lz=float(logsumexp(lp+np.log(qw)))
    if not np.isfinite(lz): raise ValueError('quadrature has no positive target mass')
    w=np.exp(lp+np.log(qw)-lz);mean=np.einsum('i,i...->...',w,obs)
    var=np.einsum('i,i...->...',w,np.abs(obs-mean)**2)
    return dict(nodes=nodes,quadrature_weights=qw,logtarget=lp,weights=w,observables=obs,
                mean=mean,variance=var,coordinate_mean=w@nodes,
                coordinate_covariance=(nodes-w@nodes).T@(w[:,None]*(nodes-w@nodes)),
                log_integral=lz,order=order,status='NUMERICAL_QUADRATURE_NOT_CERTIFIED')


def mixture_from_quadrature(reference, *, components=9, inflation=1.6,df=7.):
    """Freeze a non-Gaussian Student mixture from paid quadrature only.

    Deterministic weighted-quantile centers follow distinct nodes. Component
    spread is the quadrature posterior covariance; exact MH fixes this fit.
    """
    w=reference['weights'];nodes=reference['nodes'];order=np.lexsort((nodes[:,1],nodes[:,0]))
    indices=np.searchsorted(np.cumsum(w[order]),(np.arange(components)+.5)/components)
    centers=nodes[order[np.minimum(indices,len(order)-1)]]
    C=reference['coordinate_covariance'];ev=np.linalg.eigvalsh(C)
    C=C+max(1e-12,1e-9*float(ev.max()))*np.eye(2)
    return StudentMixture(centers,inflation*np.linalg.cholesky(C),np.ones(components),df=df)
