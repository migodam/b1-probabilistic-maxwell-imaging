"""Normalized linear/PCA conditional Gaussian; no physics or data IO."""
import numpy as np
from scipy.linalg import solve_triangular

class SimpleConditionalGaussian:
    def __init__(self,arrays,*,rng=None):
        self.arrays={k:np.asarray(v,float).copy() for k,v in arrays.items()}
        self.rng=np.random.default_rng() if rng is None else rng
        a=self.arrays
        self.chol=np.linalg.cholesky(a['covariance'])
        if a['label_mean'].shape!=(10,) or np.any(a['label_scale']<=0):
            raise ValueError('raw weak8/probit2 label density required')

    def mean_standard(self,scalar):
        a=self.arrays
        c=(np.asarray(scalar,float)-a['scalar_mean'])/a['scalar_scale']
        return (c@a['projection'])@a['weights']

    def log_prob(self,values,scalar,image=None):
        a=self.arrays
        d=(np.asarray(values,float)-a['label_mean'])/a['label_scale']-self.mean_standard(scalar)[...,None,:]
        v=solve_triangular(self.chol,d.reshape(-1,10).T,lower=True).T
        lp=-.5*(np.sum(v*v,axis=-1)+10*np.log(2*np.pi))-np.log(np.diag(self.chol)).sum()-np.log(a['label_scale']).sum()
        return lp.reshape(d.shape[:-1])

    def sample(self,count,scalar,image=None):
        a=self.arrays
        mu=self.mean_standard(scalar)
        noise=self.rng.normal(size=(len(mu),count,10))@self.chol.T
        return a['label_mean']+a['label_scale']*(mu[:,None,:]+noise)

    @classmethod
    def load(cls,path,*,rng=None):
        with np.load(path,allow_pickle=False) as z:
            return cls({k:z[k].copy() for k in ['scalar_mean','scalar_scale','label_mean','label_scale','projection','weights','covariance']},rng=rng)
