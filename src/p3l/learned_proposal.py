"""Normalized full-dimensional proposals from supplied arrays and a fixed density.

No files, physics, truth, trainer or resource ledger are accessed. Coordinates
are x=(z,eta), u=Phi(eta): all returned densities are with respect to dz d_eta.
The caller owns likelihood payment and exact MH correction. The Gaussian
strong conditional is a proposal reference, never a posterior certificate.
"""
from dataclasses import dataclass

import numpy as np
from scipy.linalg import cho_solve, solve_triangular
from scipy.special import ndtr, logsumexp

from .moments import linear_gaussian_field_marginals, gaussian_material_moments


def _finite(value, name):
    a = np.array(value, dtype=float, copy=True)
    if not np.isfinite(a).all():
        raise ValueError(f'{name} must be finite')
    return a


@dataclass(frozen=True)
class GaussianWeakDensity:
    """Pure NumPy normalized control implementing the single-context API.

    Coordinates are (a,eta). The RNG is caller-supplied independently from the
    proposal RNG, just as torch networks own their torch random stream.
    """
    mean: np.ndarray
    covariance: np.ndarray
    rng: object = None

    def __post_init__(self):
        mean = _finite(self.mean, 'mean')
        cov = _finite(self.covariance, 'covariance')
        if mean.ndim != 1 or cov.shape != (mean.size, mean.size):
            raise ValueError('Gaussian control dimensions')
        if not np.allclose(cov, cov.T, atol=1e-12, rtol=1e-12):
            raise ValueError('Gaussian control covariance must be symmetric')
        chol = np.linalg.cholesky(cov)
        object.__setattr__(self, 'mean', mean)
        object.__setattr__(self, 'covariance', cov)
        object.__setattr__(self, '_chol', chol)
        object.__setattr__(self, 'rng', np.random.default_rng() if self.rng is None else self.rng)

    def sample(self, count, scalar=None, image=None):
        return (self.mean + self.rng.standard_normal((count, self.mean.size)) @ self._chol.T)[None]

    def log_prob(self, values, scalar=None, image=None):
        values = np.asarray(values, float)
        flat = (values-self.mean).reshape(-1, self.mean.size)
        v = solve_triangular(self._chol, flat.T, lower=True).T
        result = -.5*(np.sum(v*v, axis=1) + self.mean.size*np.log(2*np.pi)) - np.log(np.diag(self._chol)).sum()
        return result.reshape(values.shape[:-1])


class LearnedFullProposal:
    """Fixed weak/probit density times a correlated strong Gaussian reference.

    ``network`` implements sample(count,scalar,image)->(1,count,k+p) and
    log_prob(values,scalar,image)->(1,count). Existing torch densities are
    adapted with no_grad; NumPy controls use the same shape convention.
    Supplied context must contain exactly one parent. Freeze/evaluate torch
    networks before passing them; a training-mode network is rejected.
    """
    def __init__(self, network, scalar, image, Vw, Vs, q0, J, residual, *, defensive_weight=.1):
        self.Vw, self.Vs = _finite(Vw, 'Vw'), _finite(Vs, 'Vs')
        self.q0, self.J, self.residual = (_finite(q0, 'q0'), _finite(J, 'J'), _finite(residual, 'residual'))
        if self.Vw.ndim != 2 or self.Vs.ndim != 2:
            raise ValueError('frames must be matrices')
        self.d, self.k = self.Vw.shape
        self.s = self.Vs.shape[1]
        if self.Vs.shape[0] != self.d or self.k+self.s != self.d or min(self.k,self.s) < 1:
            raise ValueError('complete nonempty weak/strong orthogonal frame required')
        frame = np.column_stack((self.Vw, self.Vs))
        if not np.allclose(frame.T@frame, np.eye(self.d), atol=1e-10, rtol=0):
            raise ValueError('frame must be complete orthogonal; chart absolute Jacobian one')
        if self.q0.ndim != 1 or self.q0.size <= self.d:
            raise ValueError('q0=(z0,u0) required')
        self.p = self.q0.size-self.d
        self.dimension = self.d+self.p
        if np.any((self.q0[self.d:] < 0) | (self.q0[self.d:] > 1)):
            raise ValueError('q0 property coordinates must lie in [0,1]')
        if self.J.ndim != 2 or self.J.shape[1] != self.dimension or self.residual.shape != (self.J.shape[0],):
            raise ValueError('J=(observations,d+p), residual=(observations,) required')
        if not np.isfinite(defensive_weight) or not 0 < defensive_weight < 1:
            raise ValueError('defensive weight strictly between zero and one required')
        self.defensive_weight = float(defensive_weight)
        self.Jw = self.J[:,:self.d]@self.Vw
        self.Jb = self.J[:,:self.d]@self.Vs
        self.Ju = self.J[:,self.d:]
        self.H = np.eye(self.s)+self.Jb.T@self.Jb
        self.chol_H = np.linalg.cholesky(self.H)
        self.covariance_b = cho_solve((self.chol_H,True),np.eye(self.s))
        self.network = network
        self._torch = None
        if hasattr(network, 'parameters'):
            import torch
            if network.training:
                raise ValueError('network must be fixed in evaluation mode')
            self._torch = torch
            parameter = next(network.parameters(), None)
            if parameter is None:
                parameter = next(network.buffers(), None)
            dtype = torch.float32 if parameter is None else parameter.dtype
            device = 'cpu' if parameter is None else parameter.device
            self.scalar = torch.as_tensor(scalar,dtype=dtype,device=device).detach().clone()
            self.image = torch.as_tensor(image,dtype=dtype,device=device).detach().clone()
            if self.scalar.ndim == 1:self.scalar=self.scalar[None]
            if self.image.ndim == 3:self.image=self.image[None]
            if not torch.isfinite(self.scalar).all() or not torch.isfinite(self.image).all():
                raise ValueError('context must be finite')
        else:
            self.scalar,self.image=_finite(scalar,'scalar'),_finite(image,'image')
            if self.scalar.ndim == 1:self.scalar=self.scalar[None]
            if self.image.ndim == 3:self.image=self.image[None]
        if self.scalar.ndim != 2 or self.image.ndim != 4 or self.scalar.shape[0] != 1 or self.image.shape[0] != 1:
            raise ValueError('single-parent scalar/image context required')

    def _network_call(self, method, value):
        if self._torch is None:
            return np.asarray(getattr(self.network,method)(value,self.scalar,self.image),float)
        torch=self._torch
        if self.network.training:
            raise ValueError('network changed to training mode')
        with torch.no_grad():
            if method == 'log_prob':
                value=torch.as_tensor(value,dtype=self.scalar.dtype,device=self.scalar.device)
            return getattr(self.network,method)(value,self.scalar,self.image).detach().cpu().numpy().astype(float)

    def to_chart(self, x):
        x=_finite(x,'x')
        if x.shape[-1:] != (self.dimension,):raise ValueError('x dimension')
        delta=x[...,:self.d]-self.q0[:self.d]
        return delta@self.Vw, delta@self.Vs, x[...,self.d:]

    def from_chart(self, a, b, eta):
        a,b,eta=(_finite(v,n) for v,n in ((a,'a'),(b,'b'),(eta,'eta')))
        if a.shape[-1:] != (self.k,) or b.shape[-1:] != (self.s,) or eta.shape[-1:] != (self.p,):
            raise ValueError('chart dimension')
        z=self.q0[:self.d]+a@self.Vw.T+b@self.Vs.T
        shape=np.broadcast_shapes(z.shape[:-1],eta.shape[:-1])
        return np.concatenate((np.broadcast_to(z,shape+(self.d,)),np.broadcast_to(eta,shape+(self.p,))),axis=-1)

    def strong_mean(self, a, u):
        a,u=_finite(a,'a'),_finite(u,'u')
        if a.shape[-1:] != (self.k,) or u.shape[-1:] != (self.p,) or np.any((u<0)|(u>1)):
            raise ValueError('weak/property coordinate dimensions and property box required')
        v=self.residual+a@self.Jw.T+(u-self.q0[self.d:])@self.Ju.T
        rhs=self.Vs.T@self.q0[:self.d]+v@self.Jb
        return -cho_solve((self.chol_H,True),rhs.reshape(-1,self.s).T).T.reshape(rhs.shape)

    def learned_log_prob(self, x):
        x=_finite(x,'x');a,b,eta=self.to_chart(x)
        values=np.concatenate((a,eta),axis=-1)
        shape=values.shape[:-1]
        lp=self._network_call('log_prob',values.reshape(1,-1,self.k+self.p))
        if lp.shape != (1,values.reshape(-1,self.k+self.p).shape[0]):
            raise ValueError('network log_prob output shape')
        if np.isnan(lp).any() or np.isposinf(lp).any():
            raise ValueError('network returned invalid log density')
        delta=b-self.strong_mean(a,ndtr(eta))
        quadratic=np.sum((delta@self.chol_H)**2,axis=-1)
        strong=-.5*(quadratic+self.s*np.log(2*np.pi))+np.log(np.diag(self.chol_H)).sum()
        return lp.reshape(shape)+strong

    def log_prob(self, x):
        x=_finite(x,'x')
        learned=self.learned_log_prob(x)
        prior=-.5*(np.sum(x*x,axis=-1)+self.dimension*np.log(2*np.pi))
        return logsumexp(np.stack((np.log1p(-self.defensive_weight)+learned,np.log(self.defensive_weight)+prior)),axis=0)

    def sample(self, count, *, rng=None):
        if not isinstance(count,(int,np.integer)) or count < 1:raise ValueError('positive integer sample count required')
        rng=np.random.default_rng() if rng is None else rng
        prior_mask=rng.random(count)<self.defensive_weight
        result=np.empty((count,self.dimension))
        n=int((~prior_mask).sum())
        if n:
            values=self._network_call('sample',n)
            if values.shape != (1,n,self.k+self.p) or not np.isfinite(values).all():
                raise ValueError('network sample output shape/finite contract')
            a,eta=values[0,:,:self.k],values[0,:,self.k:]
            noise=solve_triangular(self.chol_H.T,rng.standard_normal((n,self.s)).T,lower=False).T
            b=self.strong_mean(a,ndtr(eta))+noise
            result[~prior_mask]=self.from_chart(a,b,eta)
        result[prior_mask]=rng.standard_normal((int(prior_mask.sum()),self.dimension))
        return result

    def conditional_field_moments(self, a, u, matrix, *, offset=0.):
        """Scalar pixel marginals under b|a,u, retaining all b correlations."""
        a,u=_finite(a,'a'),_finite(u,'u')
        if a.shape != (self.k,) or u.shape != (self.p,):
            raise ValueError('one conditioned outer state required')
        matrix=_finite(matrix,'field matrix')
        if matrix.ndim != 2 or matrix.shape[1] != self.d:raise ValueError('field matrix=(pixels,d) required')
        base=self.q0[:self.d]+self.Vw@a
        return linear_gaussian_field_marginals(matrix@self.Vs,self.strong_mean(a,u),self.covariance_b,offset=np.asarray(offset)+matrix@base)

    def conditional_material_moments(self, a, u, matrix, inclusion_increment, *, offset=0., host=.1+.04j, mask=1., halfwidth=.5):
        """Reference moments with caller-supplied fixed property increment at u.

        Returns pixel marginals only; outer-state integration must use total
        variance and preserve a/u dependence. No posterior-moment claim.
        """
        mean,sd=self.conditional_field_moments(a,u,matrix,offset=offset)
        return gaussian_material_moments(mean,sd,inclusion_increment,host=host,mask=mask,halfwidth=halfwidth)
