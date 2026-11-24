"""Reusable original Paper 1 uniform .5 certificate with exact archived Fourier profile."""

import sys

import numpy as np

from native_box import OPTICS, clock, elapsed

sys.path.insert(0,str(OPTICS/'fine'))
from verify_fourier import ExactGeometry, one_y, profile_proof
from fourier_bound import grid


FIXED = {
    'PW180_10x20': (12,.5,.4,.001,.8125),
    'PW360_10x20': (12,.5,0.,.0025,.78125),
}


class UniformLow:
    """Cold verification is separate; repeated score calls only form source matrices."""

    def __init__(self,geometry,aperture):
        self.geometry=geometry
        self.aperture=aperture
        self.ready=False
        self.trial=FIXED[aperture]
        self.prepare_cost=None

    def prepare(self):
        if self.ready:
            return self.prepare_cost
        start=clock()
        g=self.geometry
        retained,cap,t,mu,c=self.trial
        if g.side!=96 or g.p!=192:
            raise ValueError('The frozen uniform profiles are native 96² synthetic profiles.')
        with np.load(OPTICS/'synthetic/public/fine_geometry.npz',allow_pickle=False) as f:
            same=np.array_equal(g.g,f['spectrum']) and np.array_equal(g.s,f[self.aperture+'__S']) and np.array_equal(g.b,f[self.aperture+'__B'])
        if not same:
            raise ValueError('Uniform profile requires its archived native transmit/receive geometry.')
        q=grid(retained,g.p)
        exact=ExactGeometry(g.s,g.b,q,g.p,retained)
        self.response,self.nu,self.proof=profile_proof(exact,g.g,q,g.nr,g.ns,cap,t,mu,c)
        self.ready=True
        self.prepare_cost=elapsed(start)
        return self.prepare_cost

    def score(self,y,radius):
        start=clock()
        warm=self.ready
        self.prepare()
        _,cap,t,mu,c=self.trial
        bound,terms=one_y(self.response,self.nu,y,cap,t,mu)
        return dict(branch='uniform_low_050',decision='PRUNE' if bound>radius else 'RETAIN',
                    radius=radius,residual_lower=bound,rigor='strict_native_archived_uniform_profile',
                    fixed_profile=self.trial,warm=warm,prepare=self.prepare_cost,proof=self.proof,
                    source_tail_eigenvalue_lowers=terms,**elapsed(start))
