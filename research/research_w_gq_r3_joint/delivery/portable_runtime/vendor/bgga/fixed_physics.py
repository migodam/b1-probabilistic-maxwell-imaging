"""Paid full scalar-TM solves and exact fixed-space material/shape derivatives.

No truth/reference loader. W and X are static unknown parameters of each world;
changing parameters for an audit is not a motion model or a new observation.
"""
from dataclasses import dataclass
import time
import numpy as np
from scipy.sparse.linalg import splu
from skfem import asm
from wslam.fem_backend import FEM,_contrast_mass,_forcing,line_green,_sparse_bytes
from wslam.world import boundary_quadrature
from .fixed_mesh import public_fixed_mesh,TargetAssembler
from .operators import relative_backward

@dataclass
class FixedState:
    data: object
    u: object
    matrix: object
    rhs: object
    factor: object
    target: object
    world: object
    residual: dict
    timings: dict

class FixedEvaluator:
    def __init__(self,book,phase='G0',*,h=.17,wall_refinement=2,target_step=.125):
        self.book,self.phase=book,phase
        self.declared_completion=False
        tx=np.array([[-1.5,0.],[1.5,0.]])
        rx=np.column_stack((np.linspace(-2,2,16),np.zeros(16)))
        self.p=FEM(1.,tx,rx,book,h=h,order=3,phase=phase,
            mesh=public_fixed_mesh(h=h,grid_refinement=wall_refinement,target_step=target_step),
            residual_tolerance=1e-10)
        self.ass=TargetAssembler(self.p.basis,self.p._free_dofs)
        self.R=self.p._receiver_probe[:,self.p._free_dofs]
        self.direct=line_green(rx,tx,self.p.k).T
        self.incident=line_green(self.p.material_points,tx,self.p.k)
        self._key=None;self._last=None
    def background(self,world):
        p=self.p;eps=p._material(world.coefficient(p.freq,384,include_target=False))
        m=asm(_contrast_mass,p.material_basis,contrast=(eps-1).reshape(p._material_shape))
        L=(p._background_matrix-p.k**2*m[p._free_dofs][:,p._free_dofs]).tocsc()
        b=np.column_stack([asm(_forcing,p.material_basis,
            load=(p.k**2*(eps-1)*self.incident[:,j]).reshape(p._material_shape))[p._free_dofs] for j in range(len(p.tx))])
        return L,b
    def evaluate(self,world,*,label='fixed world'):
        p=self.p;t=time.process_time();wall=time.perf_counter()
        key=(world.support.tobytes(),np.asarray(world.epsinf).tobytes(),np.asarray(world.ell).tobytes(),
             *world.pose,world.radius,world.contrast,world.target_loss)
        if key==self._key and self._last is not None:return self._last
        self._last=None;self._key=None
        L,b=self.background(world)
        a=self.ass.assemble(world.pose,world.radius,p.k,world.contrast+1j*world.target_loss,
            lambda pts:line_green(pts,p.tx,p.k))
        matrix=(L-a.Vtarget).tocsc();rhs=b+a.bX
        self.book.reserve(self.phase,'fixed material-frequency LU',factorizations=1,world_label=label,declared_completion=self.declared_completion)
        factor=splu(matrix,permc_spec='COLAMD')
        storage=_sparse_bytes(factor.L)+_sparse_bytes(factor.U)
        if storage>p.max_factor_bytes:raise MemoryError('Actual fixed-factor storage guard exceeded')
        self.book.reserve(self.phase,'fixed full primal',rhs=len(p.tx),world_label=label,declared_completion=self.declared_completion)
        u=factor.solve(rhs);defect=relative_backward(matrix,u,rhs)
        if defect['backward']>1e-10 or not np.isfinite(u).all():raise RuntimeError('fixed full linear solve failed')
        data=self.direct+(self.R@u).T
        state=FixedState(data,u,matrix,rhs,factor,a,world,defect,
            dict(kernel_cpu_seconds=time.process_time()-t,wall_seconds=time.perf_counter()-wall,factor_bytes=storage))
        self.book.event(self.phase,'fixed_evaluation',world_label=label,backward=defect,
            ndofs=p.ndofs,actual_factor_bytes=storage,kernel_cpu_seconds=state.timings['kernel_cpu_seconds'])
        self._key=key;self._last=state
        return state
    def probe(self,points):
        """Local target-ROI FE point evaluation, no all-grid locator fallback."""
        points=np.asarray(points,float);ass=self.ass
        k=min(24,len(ass.roi_elements));_,near=ass._tree.query(points,k=k)
        if near.ndim==1:near=near[:,None]
        picked=np.full(len(points),-1,int)
        for j in range(near.shape[1]):
            todo=np.flatnonzero(picked<0)
            if not len(todo):break
            local=near[todo,j];tri=ass._triangles[local]
            x=tri[:,1]-tri[:,0];y=tri[:,2]-tri[:,0];d=points[todo]-tri[:,0]
            det=x[:,0]*y[:,1]-x[:,1]*y[:,0]
            a=(d[:,0]*y[:,1]-d[:,1]*y[:,0])/det
            b=(x[:,0]*d[:,1]-x[:,1]*d[:,0])/det
            good=(a>=-2e-11)&(b>=-2e-11)&(a+b<=1+2e-11)
            picked[todo[good]]=ass.roi_elements[local[good]]
        if np.any(picked<0):raise ValueError('target boundary local locator failed; no global fallback')
        return ass._known_element_probes(points,picked,np.arange(self.p.ndofs))[:,self.p._free_dofs]
    def material_field(self,state):
        p=self.p;dofs=np.zeros((p.ndofs,len(p.tx)),complex);dofs[p._free_dofs]=state.u
        scattered=np.column_stack([np.asarray(p.material_basis.interpolate(dofs[:,j])).reshape(-1) for j in range(len(p.tx))])
        return self.incident+scattered
    def derivatives(self,state,*,include_nuisance=True):
        """A=[shared wall epsinf,ell]; B=[px,py,contrast,ell,radius,gammaR,gammaI].

        Wall derivatives pay four native JVP RHS. Target shape/material terms
        share sixteen receiver-adjoint RHS. Real parameters are whitened by
        the caller, preserving all real/imag and source dimensions.
        """
        p=self.p;w=state.world;ns=len(p.tx)
        field=self.material_field(state);mask=w.occupied(p.material_points)
        mu=np.column_stack([asm(_forcing,p.material_basis,
            load=(p.k**2*mask*field[:,j]).reshape(p._material_shape))[p._free_dofs] for j in range(ns)])
        rhs=np.column_stack((mu,1j/p.freq*mu))
        self.book.reserve(self.phase,'native wall JVP',rhs=rhs.shape[1])
        sol=state.factor.solve(rhs)
        A=np.column_stack(((self.R@sol[:,:ns]).T.ravel(),(self.R@sol[:,ns:]).T.ravel()))
        if not include_nuisance:return A,None
        self.book.reserve(self.phase,'shared receiver adjoint',rhs=len(p.rx))
        q=state.factor.solve(self.R.conj().T.toarray(),trans='H')
        adj_res=relative_backward(state.matrix.conj().T,q,self.R.conj().T.toarray())
        if adj_res['backward']>1e-10:raise RuntimeError('receiver adjoint failed')
        via_adjoint=np.column_stack(((q.conj().T@mu).T.ravel(),(q.conj().T@(1j/p.freq*mu)).T.ravel()))
        duality=float(np.linalg.norm(A-via_adjoint)/max(np.linalg.norm(A),np.finfo(float).tiny))
        if duality>1e-9:raise RuntimeError('native JVP/receiver adjoint mismatch')
        state.derivative_audit=dict(wall_JVP_adjoint_relative=duality,adjoint_backward=adj_res,
            scope='same discrete equations, independent primal and conjugate-transpose solves')
        chi=w.contrast+1j*w.target_loss/p.freq
        if chi==0:raise ValueError('zero target derivative needs independent unit-mass assembly')
        prime=(state.target.Vtarget@state.u+state.target.bX)/chi
        contrast=(q.conj().T@prime).T.ravel();loss=1j/p.freq*contrast
        pts,normal,ds,radial=boundary_quadrature(w.pose,w.radius,384,degree=4)
        P=self.probe(pts);total=line_green(pts,p.tx,p.k)+P@state.u;adj=P@q
        shapes=[]
        for direction in (normal[:,0],normal[:,1],radial):
            shapes.append((p.k**2*chi*np.einsum('qr,q,qs->sr',adj.conj(),ds*direction,total,optimize=True)).ravel())
        B=np.column_stack((shapes[0],shapes[1],contrast,loss,shapes[2],state.data.ravel(),1j*state.data.ravel()))
        return A,B
