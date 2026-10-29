"""Exact algebraic Schur chart controls; no claim to new scattering theory."""
import numpy as np
def schur_step(hessian,gradient,nwall=12):
    """Exact nuisance elimination, equivalent to the same unconstrained GN step.

    This is parameter-space algebra only. It is deliberately not called SOM
    or a current/action compression and does not claim physical solve savings.
    """
    A=hessian[:nwall,:nwall];B=hessian[:nwall,nwall:];D=hessian[nwall:,nwall:];a=gradient[:nwall];d=gradient[nwall:]
    AB=np.linalg.solve(A,B);aa=np.linalg.solve(A,a);s=D-B.T@AB;t=-np.linalg.solve(s,d-B.T@aa);w=-aa-AB@t
    return np.r_[w,t]

class ForcingChart:
    """Mass-whitened weak-forcing chart, transported by current full operator.

    U lives in g=D^-1/2 q coordinates, never material coordinates. Reference
    S=R Lref^-1 D^1/2 uses 16 receiver adjoints. Each approximate Jacobian
    transports r independent loads through current L (r actual RHS).
    """
    def __init__(self,ev,state,rank=16,max_rank=16):
        from skfem import asm
        from wslam.fem_backend import _contrast_mass
        self.ev=ev;self.max_rank=max_rank;self.rank_requested=rank;self.last_rhs=None
        p=ev.p
        from bgga.fixed_mesh import TARGET_BOX,WALL_BOX
        points=p.material_points;mask=np.zeros(len(points),bool)
        for xl,xr,yb,yt in [WALL_BOX,TARGET_BOX]:mask|=(points[:,0]>=xl)&(points[:,0]<=xr)&(points[:,1]>=yb)&(points[:,1]<=yt)
        mass=asm(_contrast_mass,p.material_basis,contrast=mask.reshape(p._material_shape).astype(float))
        d=np.asarray(mass.diagonal()[p._free_dofs]).real
        self.dsqrt=np.sqrt(np.maximum(d,0.));self.dinvsqrt=np.zeros_like(d);np.divide(1.,self.dsqrt,out=self.dinvsqrt,where=self.dsqrt>1e-15)
        self.wall_mass=[];self.wall_force=[]
        for i in range(state.wall.epsinf.size):
            v=np.zeros(state.wall.shape);v.ravel()[i]=1.;M=ev.wall_mass_direction(state.wall,v);q=ev.wall_differential_rhs(state,v);self.wall_mass.append(M);self.wall_force.append(q-M@state.u)
        adj=ev._receiver_adjoint(state,'R3 reference forcing chart receiver adjoint')
        self.rebuild(adj,rank=rank,label='public uniform reference')
        q=self.rhs_columns(state);g=self.dinvsqrt[:,None,None]*q;coef=np.einsum('nr,nsp->rsp',self.U.conj(),g,optimize=True);S=adj.conj().T*self.dsqrt[None,:];approx=np.einsum('kr,rsp->skp',S@self.U,coef,optimize=True).reshape(state.data.size,q.shape[2]);full=(adj.conj().T@q.reshape(q.shape[0],-1)).reshape(len(ev.p.rx),len(ev.p.tx),q.shape[2]).transpose(1,0,2).reshape(state.data.size,q.shape[2]);self.reference_audit=dict(exact_jacobian=full,approximate_jacobian=approx,data=state.data,singular_values=self.singular_values,metric_diagonal=self.dsqrt**2)
    def rebuild(self,adj,rank=None,label='adaptive current receiver chart'):
        S=adj.conj().T*self.dsqrt[None,:]
        # S has only16rows; thin SVD is exact and stores at most16 forcing modes.
        _,s,Vh=np.linalg.svd(S,full_matrices=False);r=min(self.rank_requested if rank is None else rank,len(s));self.U=np.ascontiguousarray(Vh[:r].conj().T);self.singular_values=s
        self.ev.external_book.event('R3','forcing_chart_rebuild',freq=self.ev.p.freq,rank=r,label_of_reference=label,singular_values=s,metric='positive diagonal FE mass D on public wall union target ROI; g=D^-1/2 q')
    def rhs_columns(self,state):
        ev=self.ev;qs=[M@state.u+b for M,b in zip(self.wall_mass,self.wall_force)];qs+=[1j/ev.p.freq*q for q in qs.copy()];target=ev.target_differentials(state,names=('px','py','rx','ry','contrast','loss','angle'));qs+=[target[n] for n in ('px','py','rx','ry','contrast','loss','angle')];self.last_rhs=np.stack(qs,axis=2);return self.last_rhs
    def jacobian(self,state):
        q=self.rhs_columns(state);g=self.dinvsqrt[:,None,None]*q;coef=np.einsum('nr,nsp->rsp',self.U.conj(),g,optimize=True)
        fields=state.factor.solve(self.dsqrt[:,None]*self.U,label='R3 current operator chart transports',chart_rank=self.U.shape[1],forcing_metric='D^-1/2 q')
        recv=np.asarray(self.ev.R@fields);j=np.einsum('kr,rsp->skp',recv,coef,optimize=True).reshape(state.data.size,q.shape[2]);return j
    def audit(self,state,direction,approx_jacobian,tolerance=.1):
        q=np.einsum('nsp,p->ns',self.last_rhs,direction,optimize=True);du=state.factor.solve(q,label='R3 adaptive exact proposed-step tangent',source_columns=len(self.ev.p.tx));exact=np.asarray(self.ev.R@du).T.ravel();approx=approx_jacobian@direction;err=np.linalg.norm(exact-approx);rel=float(err/max(np.linalg.norm(exact),1e-14));bad=rel>tolerance
        event=dict(relative_defect=rel,absolute_defect=float(err),exact_norm=float(np.linalg.norm(exact)),tolerance=tolerance,failed=bool(bad),rank=self.U.shape[1],freq=self.ev.p.freq)
        self.ev.external_book.event('R3','adaptive_chart_step_audit',**event)
        if not bad:return event,None
        adj=self.ev._receiver_adjoint(state,'R3 adaptive fallback full receiver adjoint');cols=np.eye(state.wall.epsinf.size);A=self.ev._wall_jacobian(state,adj,cols,True);B=self.ev._target_jacobian(state,adj,('px','py','rx','ry','contrast','loss','angle'));full=np.column_stack((A,B))
        self.rebuild(adj,rank=self.rank_requested)
        # Enrich in source-resolved weak-forcing residuals of the failed step.
        g=self.dinvsqrt[:,None]*q;outside=g-self.U@(self.U.conj().T@g);Q,R=np.linalg.qr(outside,mode='reduced');keep=np.abs(np.diag(R))>1e-10*max(np.linalg.norm(g),1e-15);Q=Q[:,keep];space=self.max_rank-self.U.shape[1]
        if space>0 and Q.shape[1]:self.U=np.column_stack((self.U,Q[:,:space]));self.U=np.linalg.qr(self.U,mode='reduced')[0]
        self.ev.external_book.event('R3','adaptive_chart_fallback',freq=self.ev.p.freq,new_rank=self.U.shape[1],full_current_jacobian=True,enrichment='failed step residual weak forcing in mass metric')
        return event,full
