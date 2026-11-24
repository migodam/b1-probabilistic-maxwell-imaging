"""Shared translated-tile energy bound followed by a two-sided port restriction."""

import numpy as np
from flint import acb, acb_mat, arb
from scipy.fft import ifft2
from scipy.linalg import cholesky, eigh

from branches import expand
from native_box import ball, clock, elapsed, frobenius_squared, gershgorin, herm, hi, identity, lo, star
from restricted import exact_kernel

T=32.
GAP_PROPOSAL=.009


class TileProjection:
    """The same nonnegative 24² tile and unitless upper4 are fixed throughout."""

    def __init__(self,geometry):
        self.geometry=geometry
        self.gap=None
        self.setup_cost=None

    def prepare(self):
        if self.gap is not None:
            return
        start=clock()
        g=self.geometry
        side=g.side//4
        y,x=np.indices((side,side)).reshape(2,-1)
        p=g.p
        i=((y[:,None]-y[None,:])%p)*p+(x[:,None]-x[None,:])%p
        n=i.shape[0]
        kernel=exact_kernel(g)
        G=acb_mat(n,n,[kernel[int(k)] for k in i.ravel()])
        numeric=ifft2(g.g).ravel()[i]
        z=acb(1,arb(T))/2
        zz=(1+1j*T)/2
        A=identity(n)/4-z*G-star(z*G)
        A0=.25*np.eye(n)-zz*numeric-(zz*numeric).conj().T
        F=ball(cholesky(herm(A0)-GAP_PROPOSAL*np.eye(n),lower=True,check_finite=False))
        defect=frobenius_squared(A-arb(GAP_PROPOSAL)*identity(n)-F*star(F)).sqrt().upper()
        lower=(arb(GAP_PROPOSAL)-defect).lower()
        if not lower>0:
            raise ArithmeticError('Stationary tile spectral gap did not verify.')
        self.gap=lower
        self.factor=(arb(1)+arb(T)**2).sqrt()/lower
        self.setup_cost=elapsed(start)

    def score(self,y,radius,branch):
        start=clock()
        g=self.geometry
        m= np.asarray(branch.caps)
        if m.shape!=(4,4) or np.count_nonzero(m==4.)!=1 or np.count_nonzero(m)!=1:
            raise ValueError('Translated tile energy test requires one active tile of cap 4.')
        cold=self.gap is None
        self.prepare()
        ix=np.flatnonzero(expand(branch,g.side))
        S,B=g.s[:,ix],g.b[ix]
        _,left=eigh(herm(S@S.conj().T),check_finite=False)
        _,right=eigh(herm(B.conj().T@B),check_finite=False)
        ss=np.sqrt(np.cumsum(np.sum(abs(left.conj().T@S)**2,axis=1)))
        bb=np.sqrt(np.cumsum(np.sum(abs(B@right)**2,axis=0)))
        projected=left.conj().T@y@right
        yy=np.sqrt(np.cumsum(np.cumsum(abs(projected)**2,axis=0),axis=1))
        screen=yy-float(self.factor.mid())*ss[:,None]*bb[None,:]
        qr,qs=np.unravel_index(int(screen.argmax()),screen.shape)
        v,w=map(ball,(left[:,:qr+1],right[:,:qs+1]))
        norm=lambda x:frobenius_squared(x).sqrt()
        v_norm=gershgorin(star(v)*v,upper=True).sqrt()
        w_norm=gershgorin(star(w)*w,upper=True).sqrt()
        left_s=norm(star(v)*ball(S)).upper()
        right_b=norm(ball(B)*w).upper()
        projected=norm(star(v)*ball(y)*w)
        lower=(projected.lower()-left_s*self.factor*right_b)/(v_norm*w_norm)
        value=max(0.,lo(lower))
        return dict(branch=branch.name,decision='PRUNE' if value>radius else 'RETAIN',radius=radius,
                    residual_lower=value,rigor='strict_native_translated_tile_port_bound',
                    left_modes=int(qr)+1,right_modes=int(qs)+1,cold=cold,prepare=self.setup_cost,
                    exact_tile_state_gap_lower=lo(self.gap),current_norm_factor_upper=hi(self.factor),
                    projected_incident_norm_upper=hi(right_b),projected_sensor_norm_upper=hi(left_s),
                    **elapsed(start))
