"""Torch complex128 CUDA implementation of the frozen R3-DIGITS VIE formulas.

Scipy optimizers remain on the CPU. All physical dense solves, Full-J receiver
solves, and source-adjoint/Gram work use the retained CUDA LU and complex128.
"""
from contextlib import contextmanager
import numpy as np
import torch
from .physics import Model


class CudaModel(Model):
    def __init__(self,public=None,**geometry):
        if not torch.cuda.is_available():raise RuntimeError('CUDA backend requested but unavailable')
        super().__init__(public,**geometry)
        self.device=torch.device('cuda')
        with self.operation('CUDA_OPERATOR_IMPORT',cuda_operator_imports=1):
            self.tG=torch.as_tensor(self.G,dtype=torch.complex128,device=self.device)
            self.tE0=torch.as_tensor(self.E0,dtype=torch.complex128,device=self.device)
            self.tS=torch.as_tensor(self.S,dtype=torch.complex128,device=self.device)
            self.tGH=self.tG.mH.contiguous();self.tSH=self.tS.mH.contiguous()
            self.tidentity=torch.eye(self.n,dtype=torch.complex128,device=self.device)
        self.factor=None

    @staticmethod
    def synchronize():torch.cuda.synchronize()

    @contextmanager
    def operation(self,kind,rhs=0,**cost):
        # Synchronization stays inside the counted region so timings include kernels.
        with super().operation(kind,rhs=rhs,**cost):
            try:yield
            finally:torch.cuda.synchronize()

    def forward(self,m,label='trial',jac=False):
        m=np.asarray(m,float)
        if m.shape!=(self.n,) or not np.all(np.isfinite(m)):raise ValueError('Finite real material shape(n,) required')
        with self.operation('FULL_FORWARD',rhs=self.ns,LU=1,source_RHS=self.ns,full_factorization_dimension=self.n):
            tm=torch.as_tensor(m,dtype=torch.float64,device=self.device)
            A=self.tidentity-tm[:,None]*self.tG
            b=tm[:,None]*self.tE0
            factor=torch.linalg.lu_factor(A)
            j=torch.linalg.lu_solve(*factor,b)
            E=self.tE0+self.tG@j
            y=self.tS@j
            relative=float((torch.linalg.vector_norm(A@j-b)/torch.clamp(torch.linalg.vector_norm(b),min=1e-300)).item())
            state=dict(m=m.copy(),y=y.cpu().numpy(),label=label,state_relative_residual=relative,
                       _Q=E,_y=y)
        self.state,self.factor=state,factor
        if jac:self.jacobian(state)
        return state

    def jacobian(self,state=None):
        state=self.state if state is None else state
        self._matching(state)
        if 'J' not in state:
            with self.operation('FULL_ANALYTIC_J',rhs=self.nr,receiver_transpose_RHS=self.nr,FullJ=1,material_J_builds=1):
                # A.T X=S.T via conjugated A.H solve, reusing the same LU.
                T=torch.linalg.lu_solve(*self.factor,self.tS.T.conj(),adjoint=True).conj().T
                J=torch.einsum('ri,is->sri',T,state['_Q']).reshape(self.ns*self.nr,self.n)
                state['J']=torch.stack((J.real,J.imag),dim=1).reshape(2*self.ns*self.nr,self.n).cpu().numpy()
        return state['J']

    def data_gradient(self,state,target,scale2):
        self._matching(state)
        with self.operation('DATA_ADJOINT',rhs=self.ns,adjoint_RHS=self.ns,data_adjoint_gradients=1):
            target=torch.as_tensor(target,dtype=torch.complex128,device=self.device)
            adj=torch.linalg.lu_solve(*self.factor,self.tSH@(state['_y']-target),adjoint=True)
            return (torch.sum(adj.conj()*state['_Q'],dim=1).real/scale2).cpu().numpy()

    def field_norm_squared(self,state):
        return float(torch.sum(torch.abs(state['_Q'])**2).item())

    def smallgram(self,state,target,gradient=True):
        self._matching(state)
        target=torch.as_tensor(target,dtype=torch.complex128,device=self.device)
        Q,R=state['_Q'],state['_y']-target
        fallback=False
        with self.operation('SMALL_GRAM',Gram=1,Gram_dimension=self.ns,Gram_factorizations=1,Gram_RHS=self.nr):
            H=Q.mH@Q;H=(H+H.mH)*.5
            eigenvalues=torch.linalg.eigvalsh(H)
            condition=float((eigenvalues[-1]/torch.clamp(eigenvalues[0],min=1e-300)).item())
            L,info=torch.linalg.cholesky_ex(H,check_errors=False)
            if int(info.item())==0:
                X=torch.cholesky_solve(R.mH,L).mH
                gram_res=float((torch.linalg.vector_norm(X@H-R)/torch.clamp(torch.linalg.vector_norm(R),min=1e-300)).item())
                fallback=not np.isfinite(condition) or condition>1e10 or gram_res>1e-8
            else:fallback=True
        if fallback:
            with self.operation('SMALL_GRAM_THIN_QR_FALLBACK',thin_QR=1,thin_QR_columns=self.ns):
                _,T=torch.linalg.qr(Q,mode='reduced')
                W=torch.linalg.solve_triangular(T.T,R.T,upper=False).T
                X=torch.linalg.solve_triangular(T.conj(),W.T,upper=True).T
                gram_res=float((torch.linalg.vector_norm(X@H-R)/torch.clamp(torch.linalg.vector_norm(R),min=1e-300)).item())
                if not bool(torch.isfinite(X).all().item()) or gram_res>1e-5:
                    raise FloatingPointError('Numerically rank-deficient source Gram; no regularizing jitter')
        psi=float(torch.sum(R.conj()*X).real.item())
        details=dict(gram_relative_residual=gram_res,gram_condition=condition,
                     gram_min_eigenvalue=float(eigenvalues[0].item()),gram_max_eigenvalue=float(eigenvalues[-1].item()),
                     thin_qr_fallback=bool(fallback))
        if not np.isfinite(psi) or psi < -1e-10:raise FloatingPointError('Invalid small-Gram objective')
        if not gradient:return max(psi,0.),None,details
        with self.operation('SMALL_GRAM_ADJOINT',rhs=self.ns,adjoint_RHS=self.ns,smallgram_gradients=1):
            C=X.mH@X
            U=self.tSH@X-self.tGH@(Q@C)
            adj=torch.linalg.lu_solve(*self.factor,U,adjoint=True)
            g=(2*torch.sum(adj.conj()*Q,dim=1).real).cpu().numpy()
        if not np.all(np.isfinite(g)):raise FloatingPointError('Nonfinite small-Gram derivative')
        return max(psi,0.),g,details
