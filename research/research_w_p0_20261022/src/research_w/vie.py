"""Rectangular scalar-TM VIE with analytic outgoing Sommerfeld boundary.

Convention exp(-i omega t), g=i H0^(1)/4, chi=epsilon_r-1,
 j=diag(chi)(incident+G j). G includes k² source-cell integration.
No PML or artificial truncation boundary is introduced: material is exactly
reference medium outside D=[-2,2]x[0,3]. This does not certify continuum error.

Cell order is (y,x), flattened C order. Data arrays are (source,receiver),
flattened C order. Jacobian rows follow that source-major order; columns are
all Re(epsilon) cells followed by all Im(epsilon) cells. Source strength is an
explicit complex multiplier of the free-space line-source Green function.
The incident field is evaluated at collocation centers; receiver kernels are
source-cell integrals. The self-cell integral is over the actual square.

Adapted from Gaussian/code/a2/physics.py; direct matrix assembly uses the
*integrated* offset kernel, avoiding that source's dense_direct GL mismatch.
All actual LU/RHS actions reserve the supplied Book first. Timings are
informational events; the caller's job wrapper owns CPU charging.
"""
from dataclasses import dataclass
import time
import warnings

import numpy as np
from scipy.integrate import quad
from scipy.linalg import LinAlgWarning, lu_factor, lu_solve
from scipy.special import hankel1


def square_self_green(k, h):
    """Integral of g over a centered square (analytic radial integration)."""
    def radial(radius):
        return radius*hankel1(1, k*radius)/k + 2j/(np.pi*k*k)
    def value(theta):
        radius=h/(2*max(abs(np.cos(theta)), abs(np.sin(theta))))
        return .25j*radial(radius)
    result=0j
    for lo in (-np.pi/4, np.pi/4, 3*np.pi/4, 5*np.pi/4):
        result += quad(lambda t:value(t).real,lo,lo+np.pi/2,epsabs=2e-12,epsrel=2e-12)[0]
        result += 1j*quad(lambda t:value(t).imag,lo,lo+np.pi/2,epsabs=2e-12,epsrel=2e-12)[0]
    return result


def point_green(dst, src, k):
    distance=np.linalg.norm(np.asarray(dst)[:,None,:]-np.asarray(src)[None,:,:],axis=-1)
    if np.any(distance == 0):
        raise ValueError('Coincident point source/observation requires separate regularization')
    return .25j*hankel1(0,k*distance)


def source_integral(dst, centers, k, h, quadrature):
    nodes,weights=np.polynomial.legendre.leggauss(quadrature)
    out=np.zeros((len(dst),len(centers)),complex)
    for a,wa in zip(nodes,weights):
        for b,wb in zip(nodes,weights):
            out += wa*wb*(h/2)**2*point_green(dst,centers+np.array([a,b])*h/2,k)
    return out


@dataclass
class State:
    epsilon: np.ndarray
    active: np.ndarray
    incident: np.ndarray
    field: np.ndarray
    current: np.ndarray
    data: np.ndarray
    scattered: np.ndarray
    factor: object
    residual: float
    jacobian: object
    timings: dict
    source_major: bool = True

    @property
    def E(self):
        return self.field


class VIE:
    def __init__(self,nx,freq,tx,rx,book,phase='validation',quadrature=4,
                 source_normalization=1.,residual_tolerance=1e-10,source_cell_average=False):
        tic=time.perf_counter()
        if int(nx)!=nx or nx<4 or nx%4:
            raise ValueError('nx must be a positive multiple of four, ny=3nx/4')
        if not np.isfinite(freq) or freq<=0:
            raise ValueError('freq=f/f0 must be positive')
        if int(quadrature)!=quadrature or quadrature<2:
            raise ValueError('quadrature order must be integer >=2')
        if book is None or not all(hasattr(book,key) for key in ('reserve','event')):
            raise TypeError('A physical-budget Book is required')
        self.nx=int(nx);self.ny=3*self.nx//4;self.N=self.nx*self.ny
        self.freq=float(freq);self.k=2*np.pi*self.freq;self.h=4/self.nx
        self.book=book;self.phase=phase;self.quadrature=int(quadrature)
        self.source_normalization=complex(source_normalization)
        if not np.isfinite(self.source_normalization):raise ValueError('finite source normalization required')
        self.residual_tolerance=float(residual_tolerance)
        self.tx=self._antennas(tx,4,'tx');self.rx=self._antennas(rx,12,'rx')
        x=-2+(np.arange(self.nx)+.5)*self.h
        y=(np.arange(self.ny)+.5)*self.h
        xx,yy=np.meshgrid(x,y,indexing='xy');self.points=np.c_[xx.ravel(),yy.ravel()]
        dx=np.arange(1-self.nx,self.nx)*self.h
        dy=np.arange(1-self.ny,self.ny)*self.h
        ox,oy=np.meshgrid(dx,dy,indexing='xy');offsets=np.c_[ox.ravel(),oy.ravel()]
        # Replace the central offset by a harmless nonzero one during ordinary
        # quadrature; the singular cell is then replaced by its actual integral.
        middle=(self.ny-1)*(2*self.nx-1)+self.nx-1
        safe=offsets.copy();safe[middle]=[10*self.h,10*self.h]
        values=source_integral(safe,np.zeros((1,2)),self.k,self.h,self.quadrature)[:,0]
        values[middle]=square_self_green(self.k,self.h)
        self.kernel=(self.k**2*values).reshape(2*self.ny-1,2*self.nx-1)
        # Keep only the Toeplitz kernel. A dense all-domain Green matrix would
        # exceed 2 GB already at nx=128; exact active blocks are assembled below.
        self.S=self.k**2*source_integral(self.rx,self.points,self.k,self.h,self.quadrature)
        self.source_cell_average=bool(source_cell_average)
        self.incident=self.source_normalization*(
            source_integral(self.tx,self.points,self.k,self.h,self.quadrature).T/self.h**2
            if self.source_cell_average else point_green(self.points,self.tx,self.k))
        self.direct=self.source_normalization*point_green(self.rx,self.tx,self.k)
        self.setup_seconds=time.perf_counter()-tic
        self.book.event(event='backend_setup',phase=self.phase,nx=self.nx,ny=self.ny,
                        freq=self.freq,wall_seconds=self.setup_seconds,
                        boundary='analytic outgoing Sommerfeld; no truncation/PML',
                        charged_cpu_seconds=0,quadrature=self.quadrature,
                        source_cell_average=self.source_cell_average)

    def green_block(self,rows,cols):
        ry,rx=np.divmod(np.asarray(rows),self.nx)
        cy,cx=np.divmod(np.asarray(cols),self.nx)
        return self.kernel[ry[:,None]-cy[None,:]+self.ny-1,
                           rx[:,None]-cx[None,:]+self.nx-1]

    @property
    def G(self):
        if self.N>2048:raise MemoryError('Use green_block; full Green storage disabled')
        return self.green_block(np.arange(self.N),np.arange(self.N))

    @staticmethod
    def _antennas(value,count,name):
        a=np.asarray(value,float)
        if a.shape!=(count,2) or not np.isfinite(a).all():
            raise ValueError(f'{name} requires finite ({count},2) coordinates')
        if np.any((a[:,0]>=-2)&(a[:,0]<=2)&(a[:,1]>=0)&(a[:,1]<=3)):
            raise ValueError('Antennas must lie outside closed material domain')
        return a.copy()

    def _material(self,epsilon):
        a=np.asarray(epsilon,complex)
        if a.shape not in ((self.N,),(self.ny,self.nx)):
            raise ValueError('epsilon shape must be N or (ny,nx)')
        a=a.ravel().copy()
        if not np.isfinite(a).all() or np.any(a.real<1) or np.any(a.real>9):
            raise ValueError('Requires finite 1<=Re epsilon<=9')
        if np.any(a.imag<0) or np.any(a.imag*self.freq>.15):
            raise ValueError('Requires passive conductivity 0<=f Im epsilon<=0.15')
        # Shared conductivity across frequencies is a caller-level contract.
        return a

    def _solve(self,factor,rhs,*,trans=0,label='forward'):
        self.book.reserve(self.phase,label,rhs=rhs.shape[1],freq=self.freq)
        tic=time.perf_counter()
        result=lu_solve(factor,rhs,trans=trans,check_finite=False)
        elapsed=time.perf_counter()-tic
        if not np.isfinite(result).all():raise RuntimeError('Nonfinite VIE linear solve')
        self.book.event(event='backend_solve',phase=self.phase,label=label,
                        rhs=rhs.shape[1],wall_seconds=elapsed,charged_cpu_seconds=0)
        return result,elapsed

    def solve(self,epsilon,jacobian=False):
        started=time.perf_counter();eps=self._material(epsilon);chi=eps-1
        # Exact zero only: no approximation or hidden support threshold.
        active=np.flatnonzero(chi!=0);n=len(active)
        current=np.zeros_like(self.incident);factor=None
        factor_seconds=solve_seconds=adjoint_seconds=0.;residual=0.
        if n:
            if n>4500:raise MemoryError('Active dense-LU bound 4500 cells; no approximate fallback')
            # Assemble in small rows directly into Fortran storage; the previous
            # expression held several full complex matrices simultaneously.
            matrix=np.empty((n,n),dtype=np.complex128,order='F')
            for first in range(0,n,128):
                block=np.arange(first,min(first+128,n))
                matrix[block]=-chi[active[block],None]*self.green_block(active[block],active)
                matrix[block,block]+=1
            rhs=chi[active,None]*self.incident[active]
            self.book.reserve(self.phase,'material_frequency_LU',factorizations=1,
                              freq=self.freq,active_cells=n,total_cells=self.N)
            tic=time.perf_counter()
            with warnings.catch_warnings():
                warnings.simplefilter('error',LinAlgWarning)
                factor=lu_factor(matrix,check_finite=False,overwrite_a=True)
            factor_seconds=time.perf_counter()-tic
            self.book.event(event='backend_factor',phase=self.phase,active_cells=n,
                            wall_seconds=factor_seconds,charged_cpu_seconds=0)
            current[active],solve_seconds=self._solve(factor,rhs,label='forward_sources')
            residual_sq=0.
            for first in range(0,n,128):
                block=np.arange(first,min(first+128,n))
                defect=current[active[block]]-chi[active[block],None]*(self.green_block(active[block],active)@current[active])-rhs[block]
                residual_sq+=float(np.linalg.norm(defect)**2)
            residual=float(np.sqrt(residual_sq)/max(np.linalg.norm(rhs),np.finfo(float).tiny))
            if not np.isfinite(residual) or residual>self.residual_tolerance:
                raise RuntimeError(f'VIE backward residual {residual} exceeds tolerance')
        field=self.incident.copy()
        for first in range(0,self.N,256):
            rows=np.arange(first,min(first+256,self.N))
            field[rows]+=self.green_block(rows,active)@current[active]
        scattered=self.S[:,active]@current[active]
        J=None
        if jacobian:
            dressed=self.S.copy()
            if n:
                # Bilinear transpose (trans=1), *not* conjugate adjoint:
                # rows of S_A L_A^{-1}. Exactly 12 complex RHS.
                receive,adjoint_seconds=self._solve(factor,self.S[:,active].T,
                                                   trans=1,label='receiver_transpose')
                for first in range(0,self.N,256):
                    cols=np.arange(first,min(first+256,self.N))
                    dressed[:,cols]+=(receive.T*chi[active])@self.green_block(active,cols)
            derivative=np.einsum('rn,ns->srn',dressed,field).reshape(-1,self.N)
            # Complex data derivative wrt two real material coordinate blocks.
            J=np.concatenate([derivative,1j*derivative],axis=1)
        timings=dict(setup_seconds=self.setup_seconds,factor_seconds=factor_seconds,
                     forward_solve_seconds=solve_seconds,receiver_solve_seconds=adjoint_seconds,
                     state_seconds=time.perf_counter()-started)
        self.book.event(event='backend_state',phase=self.phase,active_cells=n,
                        total_cells=self.N,freq=self.freq,residual=residual,
                        jacobian=bool(jacobian),timings=timings,charged_cpu_seconds=0)
        return State(eps,active,self.incident.copy(),field,current,
                     (self.direct+scattered).T,scattered.T,factor,residual,J,timings)
