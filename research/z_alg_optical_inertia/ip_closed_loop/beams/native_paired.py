"""One-sided native Schur verification by thin original-operator actions."""
import numpy as np
from flint import acb, acb_mat, arb
from scipy.linalg import lu_factor, lu_solve

from beam_common import clocks, elapsed, herm
from interval_certify import (ball, column, diagonal_times, endpoints, eye, joined, lo,
                              midpoint, norm_upper, tr)
from native_residual import ResidualBoundBox
from quadratic import ResidualForm


def green_kernel_roundoff(grid):
    nominal=np.fft.ifft2(grid.geom.spectrum).ravel()
    return max(abs(a-acb(float(b.real),float(b.imag))).upper()
               for a,b in zip(grid.kernel,nominal,strict=True))


def pad_entries(a,radii):
    radii=list(radii)
    return acb_mat(a.nrows(),a.ncols(),
                   [z+acb(arb(0,rad),arb(0,rad))
                    for i,z in enumerate(a.entries()) for rad in (radii[i%a.ncols()],)])


def pad_all(a,radius):
    return pad_entries(a,[radius]*a.ncols())


class FastOpticalBox(ResidualBoundBox):
    """Lift only the thin residuals needed by the small output quadratic."""
    def __init__(self, exact, numeric, center, width, green_error):
        super().__init__(exact,numeric,center,width,green_error)
        t=clocks()
        G,S,B=numeric
        self.g_adj=tr(self.g)
        self.s_adj=tr(self.s)
        self.lu=lu_factor(np.eye(self.n)-center[:,None]*G,check_finite=False)
        T=lu_solve(self.lu,G.T,trans=1,check_finite=False).T
        D=lu_solve(self.lu,S.T,trans=1,check_finite=False).T
        j=lu_solve(self.lu,center[:,None]*B,check_finite=False)
        y=S@j
        e=B+G@j
        self.proposals=ResidualForm(width,e,y,T,D)
        self.width_sq=[z**2 for z in self.width]
        self.exact_current=self.inverse_action(diagonal_times(self.center,self.b))
        self.center_e=self.b+self.g*self.exact_current
        self.center_y=self.s*self.exact_current
        self.center_e_aug=joined(self.center_e,acb_mat(self.n,self.p))
        self.center_field_read=joined(self.center_y,-eye(self.p))
        self.time_extra=elapsed(t)

    def inverse_action(self,source,*,adjoint=False):
        """Validate each thin action against A or A*, with a common inverse bound."""
        mode=2 if adjoint else 0
        estimate=ball(lu_solve(self.lu,midpoint(source),trans=mode,check_finite=False))
        def defect(approx):
            if adjoint:
                computed=approx-self.g_adj*diagonal_times(self.center,approx)
            else:
                computed=approx-diagonal_times(self.center,self.g*approx)
            return source-computed

        # One defect-correction step in the high-precision working coordinates.
        # The binary64 correction is added as exact stored coefficients. It
        # removes avoidable first-order forward uncertainty, distinct from an
        # unverified value returned directly by a floating solver.
        correction=ball(lu_solve(self.lu,midpoint(defect(estimate)),trans=mode,check_finite=False))
        estimate=estimate+correction
        residual=defect(estimate)
        caps=((self.Ubound*norm_upper(column(residual,j))).upper() for j in range(source.ncols()))
        return pad_entries(estimate,caps)

    def base_action(self,x,mu):
        forward=self.inverse_action(x)
        receiver=self.s*forward
        field=self.g*forward
        adjoint_arg=self.s_adj*receiver-acb(mu)*(self.g_adj*diagonal_times(self.width_sq,field))
        return acb(mu)*x+self.inverse_action(adjoint_arg,adjoint=True)

    def exact_optical_columns(self,where):
        ei=acb_mat(self.n,1,[int(i==where) for i in range(self.n)])
        gi=acb_mat(self.n,1,[self.g[where,i].conjugate() for i in range(self.n)])
        ti=self.inverse_action(gi,adjoint=True)
        return joined(ei,ti)

    def stationary_proposals(self,z,cell,U,mult):
        mode=self.proposals
        first=np.pad(mode.base_n,((0,0),(0,self.p)))
        K0=mode.base_a*z+mode.d.conj().T@mode.d
        L=first*z+mode.readout_right
        if U is None:
            K=K0
        else:
            K=K0+U@mult@U.conj().T
            first=np.r_[mode.field[cell],np.zeros(self.p)]
            L=L+mult[0,1]*U[:,0,None]*first[None,:]
        Z=lu_solve(lu_factor(herm(K),check_finite=False),L,check_finite=False)
        Z0=None if U is None else lu_solve(lu_factor(herm(K0),check_finite=False), U,check_finite=False)
        return ball(Z),None if U is None else ball(Z0)

    def optical_form(self,form):
        t=clocks()
        mu=arb(float(form['mu']))
        state_floor=(mu*(1-self.eta_reverse**2)).lower()
        if not state_floor>0:
            return None,dict(state_negative=None)
        rows=self.center_field_read
        fields=self.center_e_aug
        wfields=diagonal_times(self.width,fields)
        const=tr(rows)*rows-acb(mu)*(tr(wfields)*wfields)
        right=self.inverse_action(self.s_adj*rows-acb(mu)*(self.g_adj*diagonal_times(self.width_sq,fields)),adjoint=True)
        loc=int(form['cell'])
        use_optical=loc>=0 and form['ratio']!=0
        U=None
        H=None
        if use_optical:
            U=self.exact_optical_columns(loc)
            delta=-acb(0,1)*acb(mu*arb(float(form['ratio']))/2)
            H=acb_mat(2,2,[0,delta,delta.conjugate(),0])
            field_row=acb_mat(1,fields.ncols(),[fields[loc,i] for i in range(fields.ncols())])
            right=right+delta*(column(U,0)*field_row)
        z,z0=self.stationary_proposals(float(mu),loc,None if U is None else midpoint(U),
                            None if H is None else midpoint(H))
        upper_inverse=1/state_floor
        correction=None
        inertia=0
        data={}
        if use_optical:
            applied=self.base_action(z0,mu)
            residual=U-applied
            rho=norm_upper(residual)
            u_gain=norm_upper(z0)+rho/state_floor
            # Stationary evaluation cancels first-order approximate-solve error.
            term=tr(U)*z0+tr(z0)*U-tr(z0)*applied
            correction=H.inv()+pad_all(term,(rho**2/state_floor).upper())
            det=(correction[0,0]*correction[1,1]-correction[0,1]*correction[1,0]).real
            if det>0 and correction[0,0].real>0:
                inertia=1
            elif det<0:
                inertia=0
            else:
                return None,dict(state_negative=None,woodbury_determinant=endpoints(det))
            upper_inverse+=u_gain**2*norm_upper(correction.inv())
            data['woodbury_determinant']=endpoints(det)
            data['woodbury_action_defect_upper']=endpoints(rho)[1]
        applied=self.base_action(z,mu)
        if use_optical:
            applied=applied+U*(H*(tr(U)*z))
        residual=right-applied
        rho=norm_upper(residual)
        anchor=const-tr(right)*z-tr(z)*right+tr(z)*applied
        # R* K^-1 R can have either sign. Its operator norm bounds every entry.
        drift=(upper_inverse*rho**2).upper()
        Q=pad_all(anchor,drift)
        return Q,dict(state_negative=inertia,normal_base_lower=lo(state_floor),
                      schur_action_defect_upper=endpoints(rho)[1],
                      schur_second_order_error_upper=endpoints(drift)[1],
                      native_inverse_residual_upper=endpoints(self.verified_center_inverse_residual)[1],
                      **data,**elapsed(t))
