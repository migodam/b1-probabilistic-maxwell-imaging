"""Counted, fixed-grid, nonlinear scalar-TM disk-pulse VIE.

Material is real contrast: j=diag(m)(E0+G j), y=S j.  One retained LU
serves all incident sources, a receiver-transpose analytic Full-J, or a
source-adjoint small-Gram derivative. No finite differences or full SVD.
"""
from collections import Counter
from contextlib import contextmanager
import time
import numpy as np
from scipy.linalg import lu_factor, lu_solve, cho_factor, cho_solve, qr, solve_triangular
from scipy.special import hankel1, jv


def build_operators(centers, sources, receivers, radius, k):
    centers, sources, receivers = [np.asarray(x, float) for x in (centers,sources,receivers)]
    radius, k = float(radius), float(k)
    d = np.linalg.norm(centers[:,None]-centers[None,:],axis=2)
    if len(centers)>1 and np.min(d+np.eye(len(centers))*1e10)<=2*radius:
        raise ValueError('Overlapping pulse disks')
    def external(distance):
        return .5j*np.pi*k*radius*jv(1,k*radius)*hankel1(0,k*distance)
    np.fill_diagonal(d,1.)
    G = external(d)
    np.fill_diagonal(G,.5j*np.pi*k*radius*hankel1(1,k*radius)-1)
    E0 = .25j*hankel1(0,k*np.linalg.norm(centers[:,None]-sources[None,:],axis=2))
    S = external(np.linalg.norm(receivers[:,None]-centers[None,:],axis=2))
    return dict(G=G,E0=E0,S=S,H=S,U0=E0)


class Model:
    def __init__(self, public=None, **geometry):
        public = dict(public or {})
        public.update(geometry)
        self.counts,self.timings,self.cpu_timings = Counter(),Counter(),Counter()
        self.events=[]
        self.guard=None
        self.state=self.factor=None
        if 'G' not in public:
            with self.operation('GREEN_SETUP',Green_setups=1):
                public.update(build_operators(**{key:public[key] for key in ('centers','sources','receivers','radius','k')}))
        self.G=np.asarray(public['G'],complex)
        self.E0=np.asarray(public.get('E0',public.get('U0')),complex)
        self.S=np.asarray(public.get('S',public.get('H')),complex)
        self.n=self.G.shape[0]
        self.ns,self.nr=self.E0.shape[1],self.S.shape[0]
        if self.G.shape!=(self.n,self.n) or self.E0.shape[0]!=self.n or self.S.shape[1]!=self.n:
            raise ValueError('Invalid G/E0/S dimensions')
        if not all(np.all(np.isfinite(a)) for a in (self.G,self.E0,self.S)):
            raise ValueError('Nonfinite public operators')
        self.GH,self.SH=self.G.conj().T.copy(),self.S.conj().T.copy()
        self.identity=np.eye(self.n,dtype=complex)
        self.counts['cached_operator_imports']+=1

    @contextmanager
    def operation(self,kind,rhs=0,**cost):
        if self.guard is not None:
            self.guard()
        self.counts[kind]+=1
        self.counts['actual_complex_RHS']+=int(rhs)
        for name,value in cost.items():
            if isinstance(value,(float,int)):
                self.counts[name]+=value
        wall,cpu=time.perf_counter(),time.process_time()
        try:
            yield
        except Exception:
            self.counts['failed_operations']+=1
            raise
        finally:
            self.timings[kind]+=time.perf_counter()-wall
            self.cpu_timings[kind]+=time.process_time()-cpu

    def forward(self,m,label='trial',jac=False):
        m=np.asarray(m,float)
        if m.shape!=(self.n,) or not np.all(np.isfinite(m)):
            raise ValueError('Finite real material shape(n,) required')
        with self.operation('FULL_FORWARD',rhs=self.ns,LU=1,source_RHS=self.ns,full_factorization_dimension=self.n):
            A=self.identity-m[:,None]*self.G
            b=m[:,None]*self.E0
            factor=lu_factor(A,check_finite=False)
            j=lu_solve(factor,b,check_finite=False)
            E=self.E0+self.G@j
            y=self.S@j
            state=dict(m=m.copy(),j=j,E=E,Q=E,y=y,label=label,
                       state_relative_residual=float(np.linalg.norm(A@j-b)/max(np.linalg.norm(b),1e-300)))
        if not all(np.all(np.isfinite(x)) for x in (j,E,y)):
            raise FloatingPointError('Nonfinite nonlinear VIE state')
        self.state,self.factor=state,factor
        if jac:
            self.jacobian(state)
        return state

    def _matching(self,state):
        if state is not self.state or self.factor is None:
            raise RuntimeError('Derivative requires matching retained LU')

    def jacobian(self,state=None):
        state=self.state if state is None else state
        self._matching(state)
        if 'J' not in state:
            with self.operation('FULL_ANALYTIC_J',rhs=self.nr,receiver_transpose_RHS=self.nr,FullJ=1,material_J_builds=1):
                T=lu_solve(self.factor,self.S.T,trans=1,check_finite=False).T
                J=np.einsum('ri,is->sri',T,state['E']).reshape(self.ns*self.nr,self.n)
                state['Jcomplex']=J
                state['J']=np.stack((J.real,J.imag),axis=1).reshape(2*self.ns*self.nr,self.n)
        return state['J']

    def smallgram(self,state,target,gradient=True):
        self._matching(state)
        Q,R=state['E'],state['y']-np.asarray(target)
        fallback=False
        with self.operation('SMALL_GRAM',Gram=1,Gram_dimension=self.ns,Gram_factorizations=1,Gram_RHS=self.nr):
            H=Q.conj().T@Q
            H=(H+H.conj().T)*.5
            eigenvalues=np.linalg.eigvalsh(H)
            condition=float(eigenvalues[-1]/max(eigenvalues[0],1e-300))
            try:
                factor=cho_factor(H,lower=True,check_finite=False)
                X=cho_solve(factor,R.conj().T,check_finite=False).conj().T
                gram_res=float(np.linalg.norm(X@H-R)/max(np.linalg.norm(R),1e-300))
                fallback=not np.isfinite(condition) or condition>1e10 or gram_res>1e-8
            except np.linalg.LinAlgError:
                fallback=True
        if fallback:
            with self.operation('SMALL_GRAM_THIN_QR_FALLBACK',thin_QR=1,thin_QR_columns=self.ns):
                _,T=qr(Q,mode='economic',check_finite=False)
                # H=T^H T. Right solves avoid squaring conditioning during factorization.
                W=solve_triangular(T.T,R.T,lower=True,check_finite=False).T
                X=solve_triangular(T.conj(),W.T,lower=False,check_finite=False).T
                gram_res=float(np.linalg.norm(X@H-R)/max(np.linalg.norm(R),1e-300))
                if not np.all(np.isfinite(X)) or gram_res>1e-5:
                    raise FloatingPointError('Numerically rank-deficient source Gram; no regularizing jitter')
        psi=float(np.vdot(R,X).real)
        details=dict(gram_relative_residual=gram_res,gram_condition=condition,
                     gram_min_eigenvalue=float(eigenvalues[0]),gram_max_eigenvalue=float(eigenvalues[-1]),
                     thin_qr_fallback=bool(fallback))
        if not np.isfinite(psi) or psi < -1e-10:
            raise FloatingPointError('Invalid small-Gram objective')
        if not gradient:
            return max(psi,0.),None,details
        with self.operation('SMALL_GRAM_ADJOINT',rhs=self.ns,adjoint_RHS=self.ns,smallgram_gradients=1):
            C=X.conj().T@X
            U=self.SH@X-self.GH@(Q@C)
            adj=lu_solve(self.factor,U,trans=2,check_finite=False)
            g=2*np.real(np.sum(adj.conj()*Q,axis=1))
        if not np.all(np.isfinite(g)):
            raise FloatingPointError('Nonfinite small-Gram derivative')
        return max(psi,0.),g,details

    @staticmethod
    def pack(y):
        a=np.asarray(y).T
        return np.stack((a.real,a.imag),axis=-1).reshape(-1)


    def data_gradient(self,state,target,scale2):
        self._matching(state)
        with self.operation('DATA_ADJOINT',rhs=self.ns,adjoint_RHS=self.ns,data_adjoint_gradients=1):
            adj=lu_solve(self.factor,self.SH@(state['y']-target),trans=2,check_finite=False)
            return np.real(np.sum(adj.conj()*state['Q'],axis=1))/scale2

    def field_norm_squared(self,state):
        return float(np.vdot(state['Q'],state['Q']).real)



# R3-DIGITS material prolongation and common regularizer, unchanged formulas.
import numpy as np
from scipy.sparse import coo_matrix


def interpolation_1d(fine,coarse):
    x=np.linspace(0,1,int(fine))
    t=x*(int(coarse)-1)
    lo=np.minimum(np.floor(t).astype(int),int(coarse)-2)
    a=t-lo
    B=np.zeros((int(fine),int(coarse)))
    B[np.arange(int(fine)),lo]=1-a
    B[np.arange(int(fine)),lo+1]=a
    return B


def material_basis(grid_shape,stage):
    ny,nx=map(int,grid_shape)
    cy,cx=(int(stage),int(stage)) if np.isscalar(stage) else tuple(map(int,stage))
    if (ny,nx)==(cy,cx):
        return np.eye(ny*nx)
    return np.kron(interpolation_1d(ny,cy),interpolation_1d(nx,cx))


def resample_material(coarse,old_shape,new_shape):
    # New coefficients sample the preceding material field with linear interpolation.
    return material_basis(new_shape,old_shape)@np.asarray(coarse)


def difference_operator(grid_shape,width=19.5):
    ny,nx=map(int,grid_shape)
    rows,cols,values=[],[],[]
    edge=0
    for y in range(ny):
        for x in range(nx):
            i=y*nx+x
            for j in ([i+1] if x+1<nx else [])+([i+nx] if y+1<ny else []):
                rows.extend([edge,edge]);cols.extend([i,j]);values.extend([-1/width,1/width]);edge+=1
    return coo_matrix((values,(rows,cols)),shape=(edge,ny*nx)).tocsr()


class Regularizer:
    def __init__(self,grid_shape,smooth=1e-4,tv=1e-4,epsilon=.01):
        self.D=difference_operator(grid_shape)
        self.smooth,self.tv,self.epsilon=float(smooth),float(tv),float(epsilon)
        self.edges=self.D.shape[0]

    def evaluate(self,m):
        d=self.D@m
        q=np.sqrt(d*d+self.epsilon**2)
        value=self.smooth*np.mean(d*d)+self.tv*np.mean(q-self.epsilon)
        gd=(2*self.smooth*d+self.tv*d/q)/self.edges
        g=np.asarray(self.D.T@gd).reshape(-1)
        weights=(2*self.smooth+self.tv/q)/self.edges
        return float(value),g,weights
