"""Fast directed full-wave checks using a verified floating inverse residual."""
import numpy as np
from flint import arb
from scipy.linalg import lu_factor, lu_solve

from beam_common import clocks, elapsed
from verify import ball, diagonal_times, tr


def roundoff(n):
    # The generous constants cover complex products expanded to real arithmetic
    # and any order of pairwise or sequential addition. All data are finite
    # binary64 normal-scale, and the absolute allowance below includes underflow.
    unit=arb(2)**-53
    u=arb(int(n))*unit
    return u/(1-u)


def norm_f64(matrix):
    if not np.all(np.isfinite(matrix)):
        raise ValueError('Residual certification requires finite binary64 coefficients.')
    count=matrix.size
    initial=arb(float(np.linalg.norm(matrix)))
    return (initial/(1-roundoff(16*count+64))+arb(2)**-490).upper()


def exact_norm(a):
    return sum((abs(z).upper()**2 for z in a.entries()),arb(0)).sqrt().upper()


class ResidualBoundBox:
    def __init__(self, exact, numeric, center, width, green_error):
        start=clocks()
        self.g,self.s,self.b=exact
        self.n=self.g.nrows()
        self.p=self.s.nrows()
        self.m=self.b.ncols()
        self.center=[arb(float(x)) for x in center]
        self.width=[arb(float(x)) for x in width]
        self.center_np=center
        self.width_np=width
        n=len(center)
        G,S,B=numeric
        Gf=norm_f64(G)
        exact_Gf=(Gf+green_error).upper()
        max_chi=arb(float(max(abs(center))))
        max_width=arb(float(max(width)))
        A=np.eye(n)-center[:,None]*G
        I=np.eye(n)
        factor=lu_factor(A,check_finite=False)
        U=lu_solve(factor,I,check_finite=False)
        alpha=roundoff(16*n+64)
        Az=norm_f64(A)
        Uz=norm_f64(U)
        Iz=arb(n).sqrt()
        eA=green_error*max_chi+roundoff(64)*(Iz+max_chi*Gf)+arb(2)**-490
        rau=(norm_f64(I-A@U)+alpha*(Iz+Az*Uz)+eA*Uz).upper()
        if not rau<1:
            raise ArithmeticError('Approximate native inverse residual is not contractive.')
        self.verified_center_inverse_residual=rau
        self.Ubound=(Uz/(1-rau)).upper()
        u_error=(self.Ubound*rau).upper()
        T=G@U
        t_error=(alpha*Gf*Uz+green_error*Uz+exact_Gf*u_error).upper()
        width_product_margin=(t_error*max_width+roundoff(64)*norm_f64(T)*max_width
                              +arb(2)**-490).upper()
        self.eta=(norm_f64(T*width[None,:])+width_product_margin).upper()
        self.eta_reverse=(norm_f64(width[:,None]*T)+width_product_margin).upper()
        if not self.eta<1 or not self.eta_reverse<1:
            raise ArithmeticError('Entire-box native inverse residual did not contract.')
        self.Sbound=exact_norm(self.s)
        d_error=(self.Sbound*(u_error+alpha*Uz)).upper()
        D=S@U
        self.D=ball(D)
        self.d_error=d_error
        self.current=ball(lu_solve(factor,center[:,None]*B,check_finite=False))
        self.E=self.b+self.g*self.current
        defect=diagonal_times(self.center,self.E)-self.current
        j_err=(self.Ubound*exact_norm(defect)).upper()
        self.e_error=(j_err*exact_Gf).upper()
        self.y_error=(j_err*self.Sbound).upper()
        self.Y=self.s*self.current
        self.time=elapsed(start)

    def support(self,a,w):
        wnorm=exact_norm(w)
        anorm=exact_norm(a)
        response=tr(w)*self.Y*a
        direction_weights=tr(diagonal_times(self.width,tr(tr(w)*self.D)))
        max_width=max(x.upper() for x in self.width)
        d=exact_norm(direction_weights)+wnorm*self.d_error*max_width
        reach=d*(exact_norm(self.E*a)+anorm*self.e_error)/(1-self.eta)
        reach+=wnorm*anorm*self.y_error
        center=response[0,0].real
        return (center-reach).lower(),(center+reach).upper()
