"""Exact native-grid actions, numerical stationary solves and Arb projections."""

from collections import Counter
from dataclasses import dataclass
from time import perf_counter, process_time

import numpy as np
from scipy.fft import fft2, ifft2
from scipy.linalg import eigh, cho_factor, cho_solve
from scipy.sparse.linalg import LinearOperator, eigsh, gmres
from flint import acb, acb_mat, arb

from curve_inputs import STUDY
from search import h, times
import sys

sys.path.insert(0, str(STUDY))
from core.certificate import ball, frobenius_squared, gershgorin, safe_float, star


def high(x):
    return safe_float(x, upper=True)


def low(x):
    return safe_float(x, upper=False)


def mid(a):
    return np.array([complex(x) for x in a.entries()]).reshape(a.nrows(), a.ncols())


def mult_rows(m, scalars):
    return acb_mat(m.nrows(),m.ncols(), [scalars[i]*m[i,j]
                   for i in range(m.nrows()) for j in range(m.ncols())])


@dataclass(frozen=True, order=True)
class Weights:
    cap: float
    scale: float
    real: tuple
    optical: tuple

    @classmethod
    def global_(cls, u, mu, t):
        return cls(u, mu, (mu,)*4, (mu*t,)*4)

    @property
    def is_global(self):
        return len(set(self.real)) == len(set(self.optical)) == 1

    def as_dict(self):
        return dict(cap=self.cap, scale=self.scale, weights=list(self.real), optical=list(self.optical))


class Geometry:
    def __init__(self, spec, sensor, incident, *, block_side=2):
        self.spec, self.sensor, self.incident = spec, sensor, incident
        self.side = int(np.sqrt(incident.shape[0]))
        self.n = self.side**2
        self.s = incident.shape[1]
        self.pad = len(spec)
        self.sensor_star = sensor.conj().T
        self.sensor_ball = ball(sensor)
        self.incident_ball = ball(incident)
        self.spec_ball = [[acb(float(z.real),float(z.imag)) for z in row] for row in spec]
        if block_side <= 0 or self.side % block_side:
            raise ValueError('The block grid must divide each native image side.')
        self.block_side = block_side
        self.block_count = block_side**2
        ix = np.arange(self.side)
        self.blocks = ((block_side*ix[:, None]//self.side)*block_side
                       + block_side*ix[None,:]//self.side).ravel()
        self.cost = Counter()

    def numerical_green(self, x, *, adjoint=False):
        n,p = self.side, self.pad
        y = fft2(x.reshape(n,n,-1),s=(p,p),axes=(0,1))
        spectrum = self.spec.conj() if adjoint else self.spec
        z = ifft2(spectrum[:,:,None]*y, axes=(0,1))[:n,:n]
        self.cost['native_binary64_green_columns'] += x.shape[1]
        return z.reshape(x.shape)

    def exact_green(self, x, *, adjoint=False):
        n,p = self.side,self.pad
        zero = acb(0)
        spectrum = self.spec_ball
        cols = []
        for c in range(x.ncols()):
            rows = [acb.dft([x[n*i+j,c] for j in range(n)]+[zero]*(p-n)) for i in range(n)]
            transformed = []
            for j in range(p):
                d = acb.dft([rows[i][j] for i in range(n)]+[zero]*(p-n))
                filtered = [d[i]*(spectrum[i][j].conjugate() if adjoint else spectrum[i][j]) for i in range(p)]
                transformed.append(acb.dft(filtered, inverse=True)[:n])
            native_rows = [acb.dft([transformed[j][i] for j in range(p)],inverse=True)[:n] for i in range(n)]
            cols.extend(k for r in native_rows for k in r)
        self.cost['native_arb_green_columns'] += x.ncols()
        return acb_mat(x.ncols(), self.n, cols).transpose()


class Native:
    def __init__(self, geometry, weights):
        self.g, self.w = geometry, weights
        if len(weights.real) != geometry.block_count or len(weights.optical) != geometry.block_count:
            raise ValueError('Weights must match the row-major block grid.')
        self.wr = np.asarray(weights.real)[geometry.blocks]/weights.scale
        optical = np.asarray(weights.optical)[geometry.blocks]/weights.scale
        self.z = (weights.cap*self.wr + 1j*optical)/2
        self.operator = LinearOperator((geometry.n, geometry.n),
                                       matvec=self.numeric_action, matmat=self.numeric_action,
                                       dtype=np.complex128)
        eye = np.eye(len(geometry.sensor))
        self.pre_factor = cho_factor(h(geometry.sensor @ geometry.sensor_star)+weights.scale*eye)

    def precondition(self, x):
        g = self.g
        return x-g.sensor_star@cho_solve(self.pre_factor,g.sensor@x,check_finite=False)

    def numeric_action(self, x):
        g, w = self.g, self.w
        shape = x.shape
        x = x.reshape(g.n,-1)
        az = (self.wr[:,None]*x + g.sensor_star@(g.sensor@x)/w.scale
              - self.z[:,None]*g.numerical_green(x)
              - g.numerical_green(self.z.conj()[:,None]*x,adjoint=True))
        return az.reshape(shape)

    def exact_scaled_coefficients(self):
        cap,scale = arb(self.w.cap),arb(self.w.scale)
        w = [arb(x)/scale for x in self.w.real]
        z = [acb(cap*r,arb(l)/scale)/2 for r,l in zip(w,self.w.optical,strict=True)]
        return [w[i] for i in self.g.blocks],[z[i] for i in self.g.blocks]

    def exact_action(self, x):
        g = self.g
        w,z = self.exact_scaled_coefficients()
        return (mult_rows(x,w)+star(g.sensor_ball)*(g.sensor_ball*x)/arb(self.w.scale)
                -mult_rows(g.exact_green(x),z)
                -g.exact_green(mult_rows(x,[v.conjugate() for v in z]),adjoint=True))

    def modes(self, nneg_expected):
        g = self.g
        start,cpu=perf_counter(),process_time()
        rng = np.random.default_rng(54273)
        trial = rng.standard_normal(g.n)+1j*rng.standard_normal(g.n)
        num = max(5,min(14,nneg_expected+2))
        eigs,v = eigsh(self.operator,k=num,which='SA',ncv=max(32,2*num+4),
                       maxiter=650,tol=1e-9,v0=trial)
        order=np.argsort(eigs)
        self.eigs,self.eig_basis=eigs[order],v[:,order]
        self.nu=int(np.count_nonzero(self.eigs<0.))
        return dict(numerical_smallest_scaled_K=self.eigs.tolist(),
                    numerical_eigen_defect=float(np.linalg.norm(self.numeric_action(self.eig_basis)
                                                               -self.eig_basis*self.eigs)),
                    **times(start,cpu))

    def ensure_negative(self, n, *, return_gap=False):
        if n == 0:
            return None
        basis=ball(self.eig_basis[:,:n])
        negative=gershgorin(star(basis)*self.exact_action(basis),upper=True)
        if not negative < 0:
            raise ArithmeticError('Native outward negative subspace condition did not hold.')
        result = high(negative)
        if return_gap:
            gram_upper = gershgorin(star(basis)*basis,upper=True)
            return result,low(-negative/gram_upper)
        return result

    def numerical_distance(self, y, nu):
        g = self.g
        start,cpu=perf_counter(),process_time()
        rhs = g.sensor_star@y/self.w.scale+self.z[:,None]*g.incident
        pre=LinearOperator((g.n,g.n),matvec=self.precondition,dtype=np.complex128)
        trial=np.empty_like(rhs)
        iterations=[]
        for j, l in enumerate(rhs.T):
            trace=[]
            trial[:,j],info=gmres(self.operator,l,M=pre,restart=65,maxiter=10,
                                  rtol=2.e-12,atol=0.,callback=trace.append,callback_type='pr_norm')
            iterations.append(len(trace))
            if info:
                raise ArithmeticError(f'Native GMRES failed at source {j}: {info}.')
        vals,v=eigh(h(y.conj().T@y-self.w.scale*rhs.conj().T@trial), check_finite=False)
        vals,v=vals[::-1],v[:,::-1]
        d=float(np.sqrt(np.maximum(vals[nu:],0.).sum()))
        return trial, v, vals, dict(distance=d,
                                    source_form_numeric_eigenvalues=vals.tolist(),
                                    GMRES_iterations=iterations,**times(start,cpu))

    def ensure_source(self, y, nu, gap, trial, modes, eigs, *, negative_gap=None):
        g = self.g
        total_modes = nu+np.count_nonzero(eigs[nu:]>0)
        num=g.s if negative_gap is not None else total_modes
        if num<=nu:
            return dict(distance_lower=0., positive_source_eigenvalues_lower=[], correction_upper=0.)
        vec=ball(modes[:,:num])
        x=ball(trial@modes[:,:num])
        yv=ball(y)*vec
        bv=g.incident_ball*vec
        scale=arb(self.w.scale)
        _,z=self.exact_scaled_coefficients()
        rhs=star(g.sensor_ball)*yv/scale+mult_rows(bv,z)
        ax=self.exact_action(x)
        defect=scale*frobenius_squared(ax-rhs)
        correction=(defect/arb(gap)).upper()
        form=star(yv)*yv+scale*(star(x)*ax-star(rhs)*x-star(x)*rhs)
        gram=star(vec)*vec
        lower=[]
        for j in range(nu+1,total_modes+1):
            projected=acb_mat(j,j,[form[i,k] for i in range(j) for k in range(j)])
            gg=acb_mat(j,j,[gram[i,k] for i in range(j) for k in range(j)])
            num0=gershgorin(projected,upper=False)-correction
            # Use a positive numerator before dividing by the Gram upper bound.
            if num0 > 0:
                lower.append((num0/gershgorin(gg,upper=True)).lower())
        total=sum(lower,arb(0))
        result=dict(distance_lower=max(0.,low(total.sqrt())),
                    positive_source_eigenvalues_lower=[low(x) for x in lower],
                    correction_upper=high(correction))
        if negative_gap is not None:
            absolute_gap = min(gap,negative_gap) if nu else gap
            if not absolute_gap>0.:
                raise ArithmeticError('An upper source certificate needs the absolute spectral gap.')
            for i,val in enumerate(eigs):
                form[i,i]-=arb(float(val))
            perturb=(frobenius_squared(form).sqrt()+defect/arb(absolute_gap)).upper()
            gram_lower=gershgorin(gram,upper=False)
            if not gram_lower>0:
                raise ArithmeticError('Full source trial matrix is not certified invertible.')
            terms=[max(arb(0),(arb(float(v))+perturb).upper())/gram_lower for v in eigs[nu:]]
            result['distance_upper']=high(sum(terms,arb(0)).sqrt())
            result['two_sided_absolute_gap_lower']=absolute_gap
            result['upper_form_perturbation']=high(perturb)
        return result
