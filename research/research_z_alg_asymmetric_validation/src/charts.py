"""Two complete, bounded, real-material correctors for the same free coordinates."""
import time
import numpy as np
from scipy.linalg import qr

def partition(J):
    _,_,p=qr(J,mode='economic',pivoting=True)
    d=np.asarray(p[:J.shape[0]],dtype=int);f=np.asarray(sorted(set(range(J.shape[1]))-set(d)),dtype=int)
    return d,f

def solve_real(book,A,b,label):
    nrhs=1 if np.asarray(b).ndim==1 else np.asarray(b).shape[1]
    with book.op(label,algebraic_real_factorizations=1,algebraic_real_RHS=nrhs,algebraic_matrix_dimension=len(A)):
        try:return np.linalg.solve(A,b)
        except np.linalg.LinAlgError:
            with book.op(label+'_FALLBACK_SVD',algebraic_real_SVD=1,algebraic_real_RHS=nrhs):
                return np.linalg.lstsq(A,b,rcond=1e-12)[0]

def box_alpha(values,direction):
    alpha=1.
    for v,d in zip(values,direction):
        if d>0:alpha=min(alpha,.995*(4-v)/d)
        elif d<0:alpha=min(alpha,.995*(.5-v)/d)
    return max(0.,alpha)

class MaterialChart:
    def __init__(self,model,target,reference,dep,free):
        self.model=model;self.book=model.book;self.target=target.copy();self.reference=reference.copy();self.dep=np.asarray(dep);self.free=np.asarray(free)
        self.scale=max(np.linalg.norm(target),1e-12)
    def correct(self,initial,free_values,label,maxiter=24):
        m=np.asarray(initial).copy();m[self.free]=free_values;history=[];start=time.perf_counter();state=None
        if np.any(m<.5) or np.any(m>4):return dict(success=False,m=m,reason='INPUT_OUTSIDE_BOX',history=[],wall_seconds=0.)
        state=self.model.forward(m,label=label+'_initial')
        for it in range(maxiter):
            r=self.model.pack(state['y']-self.target)/self.scale
            history.append(dict(iteration=it,residual=float(np.linalg.norm(r)),m=m.tolist()))
            if np.linalg.norm(r)<=2e-13:
                return dict(success=True,m=m,state=state,history=history,wall_seconds=time.perf_counter()-start,iterations=it)
            J=self.model.jacobian(state)/self.scale
            step=solve_real(self.book,J[:,self.dep],-r,'MATERIAL_NEWTON_REAL_SOLVE')
            alpha=box_alpha(m[self.dep],step);accepted=False
            if alpha<1e-12:break
            for bt in range(12):
                trial=m.copy();trial[self.dep]+=alpha*step
                candidate=self.model.forward(trial,label=label+'_trial')
                rr=self.model.pack(candidate['y']-self.target)/self.scale
                accepted=np.linalg.norm(rr)<(1-1e-4*alpha)*np.linalg.norm(r)
                self.book.append(dict(event='CORRECTOR_TRIAL',method='material',iteration=it,backtrack=bt,alpha=alpha,accepted=bool(accepted),residual=float(np.linalg.norm(rr))))
                if accepted:m=trial;state=candidate;break
                alpha*=.5
            if not accepted:break
        return dict(success=False,m=m,state=state,history=history,reason='BOUNDARY_OR_NEWTON_STAGNATION',wall_seconds=time.perf_counter()-start,iterations=len(history))

class CurrentChart:
    def __init__(self,model,target,reference_state,dep,free):
        self.model=model;self.book=model.book;self.target=target.copy();self.dep=np.asarray(dep);self.free=np.asarray(free)
        n=model.n;ns=model.ns;self.nq=model.N.shape[1];self.nc=2*self.nq*ns
        with self.book.op('CURRENT_CHART_SETUP',algebraic_complex_pseudoinverse=1,algebraic_complex_factorizations=1,algebraic_complex_RHS=ns,algebraic_matrix_dimension=model.nr,Green_matvec_columns=self.nq+ns,current_real_variables=self.nc+len(dep)):
            self.N=model.N;self.GN=model.G@self.N
            self.j0=model.S.conj().T@np.linalg.solve(model.S@model.S.conj().T,target)
            self.e0=model.E0+model.G@self.j0
            self.scale=max(np.linalg.norm(reference_state['j'])/np.sqrt(n*ns),1e-6)
            self.rscale=max(np.linalg.norm(model.E0)/np.sqrt(n*ns),1e-6)
    def encode(self,j,m):
        c=self.N.conj().T@(j-self.j0)/self.scale
        return np.concatenate([c.T.real.ravel(),c.T.imag.ravel(),m[self.dep]])
    def decode(self,z,free_values):
        ns=self.model.ns;q=self.nq;cut=q*ns
        c=(z[:cut].reshape(ns,q)+1j*z[cut:2*cut].reshape(ns,q)).T
        j=self.j0+self.scale*self.N@c
        E=self.e0+self.scale*self.GN@c
        m=np.empty(self.model.n);m[self.dep]=z[self.nc:];m[self.free]=free_values
        return m,j,E
    def residual_jac(self,z,free_values,need_jac=True):
        tick=time.perf_counter()
        m,j,E=self.decode(z,free_values);n=self.model.n;ns=self.model.ns;q=self.nq;cut=q*ns
        r=(j-m[:,None]*E).T
        rr=np.concatenate([r.real.ravel(),r.imag.ravel()])/self.rscale
        self.book.append(dict(event='CURRENT_ALGEBRA_STATE',current_basis_matvec_columns=2*ns,shared_constitutive_cells=n,kernel_seconds=time.perf_counter()-tick))
        if not need_jac:return rr,None,(m,j,E)
        with self.book.op('CURRENT_ALGEBRA_JACOBIAN',algebraic_J_builds=1,Green_matvec_columns=0,constitutive_equations_real=2*n*ns):
            D=np.zeros((2*n*ns,self.nc+len(self.dep)))
            K=self.scale*(self.N-m[:,None]*self.GN)
            for s in range(ns):
                rows=slice(s*n,(s+1)*n);imrows=slice(n*ns+s*n,n*ns+(s+1)*n)
                cr=slice(s*q,(s+1)*q);ci=slice(cut+s*q,cut+(s+1)*q)
                D[rows,cr]=K.real;D[rows,ci]=-K.imag
                D[imrows,cr]=K.imag;D[imrows,ci]=K.real
                for h,i in enumerate(self.dep):D[s*n+i,self.nc+h]=-E[i,s].real;D[n*ns+s*n+i,self.nc+h]=-E[i,s].imag
        return rr,D/self.rscale,(m,j,E)
    def correct(self,initial_z,free_values,label,maxiter=24):
        z=initial_z.copy();history=[];start=time.perf_counter();decoded=None
        if np.any(np.asarray(free_values)<.5) or np.any(np.asarray(free_values)>4):return dict(success=False,z=z,reason='FREE_OUTSIDE_BOX',history=[],wall_seconds=0.)
        for it in range(maxiter):
            r,D,decoded=self.residual_jac(z,free_values);m,j,E=decoded
            history.append(dict(iteration=it,residual=float(np.linalg.norm(r)),m=m.tolist(),minimum_field_squared=float(np.min(np.sum(abs(E)**2,axis=1)))))
            if np.linalg.norm(r)<=2e-13:
                return dict(success=True,z=z,m=m,j=j,E=E,D=D,history=history,iterations=it,wall_seconds=time.perf_counter()-start)
            step=solve_real(self.book,D,-r,'CURRENT_NEWTON_REAL_SOLVE')
            alpha=box_alpha(z[self.nc:],step[self.nc:]);accepted=False
            if alpha<1e-12:break
            for bt in range(12):
                trial=z+alpha*step
                rr,_,_ = self.residual_jac(trial,free_values,False)
                accepted=np.linalg.norm(rr)<(1-1e-4*alpha)*np.linalg.norm(r)
                self.book.append(dict(event='CORRECTOR_TRIAL',method='current',iteration=it,backtrack=bt,alpha=alpha,accepted=bool(accepted),residual=float(np.linalg.norm(rr))))
                if accepted:z=trial;break
                alpha*=.5
            if not accepted:break
        m,j,E=self.decode(z,free_values)
        return dict(success=False,z=z,m=m,j=j,E=E,history=history,reason='BOUNDARY_OR_CURRENT_NEWTON_STAGNATION',wall_seconds=time.perf_counter()-start,iterations=len(history))
    def material_tangent(self,z,free_values):
        r,D,decoded=self.residual_jac(z,free_values)
        m,j,E=decoded;n=self.model.n;ns=self.model.ns
        R=np.zeros((2*n*ns,len(self.free)))
        for h,i in enumerate(self.free):
            for s in range(ns):R[s*n+i,h]=-E[i,s].real;R[n*ns+s*n+i,h]=-E[i,s].imag
        dz=solve_real(self.book,D,-R/self.rscale,'CURRENT_IMPLICIT_MATERIAL_TANGENT')
        tangent=np.zeros((n,len(self.free)));tangent[self.free]=np.eye(len(self.free));tangent[self.dep]=dz[self.nc:]
        return tangent
