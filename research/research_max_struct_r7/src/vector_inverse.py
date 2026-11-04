"""True 3D vector DDA inverse controls; scalar epsilon_r per voxel.

exp(-iwt); lambda0=Z0=eps_b=mu_b=1. Density j=chi E; dipole p=Vj.
Self s=-1/3+i k^3V/(6pi) yields alpha=V chi/(1-s chi), exactly
CM+RR. The optimized variable is real chi=eps_r-1, NOT alpha.
"""
import time
import numpy as np
from scipy import linalg as la
from scipy.optimize import minimize
from research.research_max_struct_r1.src.vector_backend import TensorDDA,plane_waves

class VectorInverse:
    def __init__(self, nside=5,h=.08):
        t=time.perf_counter();c=time.process_time()
        self.nside=nside;self.h=h;self.N=nside**3;self.M=3*self.N
        self.ijk=np.stack(np.meshgrid(*[np.arange(nside)]*3,indexing='ij'),axis=-1).reshape(-1,3)
        self.points=(self.ijk-(nside-1)/2)*h
        self.backend=TensorDDA(self.points,h**3)
        self.G=self.backend.G;self.s=self.backend.self_density
        self.directions=np.array([[0,0,1],[0,0,-1],[1,0,0],[-1,0,0],[0,1,0],[0,-1,0]],float)
        self.polarizations=np.array([[1,0,0],[0,1,0],[0,1,0],[0,0,1],[0,0,1],[1,0,0]],float)
        self.inc=plane_waves(self.points,self.backend.k,self.directions,self.polarizations)
        n=30;z=1-2*(np.arange(n)+.5)/n;ang=np.arange(n)*np.pi*(3-np.sqrt(5));r=np.sqrt(1-z*z)
        self.receivers=.8*np.column_stack([r*np.cos(ang),r*np.sin(ang),z])
        self.S=self.backend.external_electric_operator(self.receivers)
        block=self.ijk//2;same=np.all(block[:,None,:]==block[None,:,:],axis=-1)
        self.H=(self.G-self.s*np.eye(self.M))*np.repeat(np.repeat(same,3,axis=0),3,axis=1)
        self.edges=np.array([(i,j) for i in range(self.N) for j in range(i+1,self.N) if np.sum(abs(self.ijk[i]-self.ijk[j]))==1])
        self.E2=float(np.vdot(self.inc,self.inc).real);self.beta=1.;self.rho=10.
        self.lam=2e-4;self.l2=1e-6;self.box=(0.,4.)
        self.setup=dict(wall_s=time.perf_counter()-t,cpu_s=time.process_time()-c,N=self.N,vector_dimension=self.M,sources=6,receivers=30,complex_data=540,green_assemblies=1,external_green_assemblies=1)
        self.reset_counts()
    def reset_counts(self):
        self.counts=dict(physical_LU=0,physical_source_RHS=0,physical_adjoint_RHS=0,physical_calls=0,current_normal_factorizations=0,current_RHS=0,material_value_gradient_calls=0,failed_calls=[])
    def prior(self,x):
        a,b=self.edges.T;d=x[a]-x[b];f=self.lam*np.mean(d*d)+self.l2*np.mean(x*x)
        g=2*self.l2*x/self.N;v=2*self.lam*d/len(a);np.add.at(g,a,v);np.add.at(g,b,-v)
        return f,g
    def physical(self,x,Y=None,gradient=True):
        self.counts['physical_calls']+=1;self.counts['physical_LU']+=1
        d=np.repeat(x,3);K=np.eye(self.M)-d[:,None]*self.G
        fac=la.lu_factor(K,check_finite=False);J=la.lu_solve(fac,d[:,None]*self.inc,check_finite=False)
        self.counts['physical_source_RHS']+=6
        E=self.inc+self.G@J;pred=self.S@J
        lin=float(la.norm(K@J-d[:,None]*self.inc)/max(la.norm(d[:,None]*self.inc),1e-30))
        if Y is None:return pred,J,E,lin
        delta=pred-Y;Y2=float(np.vdot(Y,Y).real);f0=float(np.vdot(delta,delta).real/Y2);p,pg=self.prior(x)
        if not gradient:return f0+p,f0,J,E,lin
        adj=la.lu_solve(fac,self.S.conj().T@delta/Y2,trans=2,check_finite=False)
        self.counts['physical_adjoint_RHS']+=6
        g=2*np.real(np.sum(np.conj(adj)*E,axis=1)).reshape(self.N,3).sum(axis=1)+pg
        return f0+p,g
    def current_update(self,x,Y):
        d=np.repeat(x,3);K=np.eye(self.M)-d[:,None]*self.G;Y2=float(np.vdot(Y,Y).real)
        normal=self.S.conj().T@self.S/Y2+self.rho/self.E2*(K.conj().T@K)
        rhs=self.S.conj().T@Y/Y2+self.rho/self.E2*K.conj().T@(d[:,None]*self.inc)
        self.counts['current_normal_factorizations']+=1;self.counts['current_RHS']+=6
        # Strictly positive K*K in declared finite box; no artificial ridge.
        return la.cho_solve(la.cho_factor(normal,check_finite=False),rhs,check_finite=False)
    def material(self,x,J,Q,kind):
        self.counts['material_value_gradient_calls']+=1
        d=np.repeat(x,3);r=J-d[:,None]*Q
        if kind=='ordinary_joint':R=r;du=-Q;g=2*np.real(np.sum(np.conj(R)*du,axis=1))
        else:
            b=self.beta-self.s;den=1+b*d;a=d/den
            u=self.beta*r/den[:,None];du=-self.beta*(Q+b*J)/den[:,None]**2
            if kind=='dual_like_joint':R=u;g=2*np.real(np.sum(np.conj(R)*du,axis=1))
            elif kind=='A_material_joint':
                Hu=self.H@u;R=u+a[:,None]*Hu
                v=R+self.H.conj().T@(np.conj(a)[:,None]*R)
                g=2*np.real(np.sum(np.conj(v)*du+np.conj(R)*Hu/den[:,None]**2,axis=1))
            else:raise ValueError(kind)
        f=self.rho/self.E2*float(np.vdot(R,R).real);p,pg=self.prior(x)
        g=self.rho/self.E2*g.reshape(self.N,3).sum(axis=1)+pg
        return f+p,g
    def optimize(self,Y,kind,outer=15,inner=12,corrector=40):
        self.reset_counts();t=time.perf_counter();c=time.process_time();x=np.full(self.N,.25);x0=x.copy();history=[]
        best=x.copy();bestphi=self.physical(x,Y,False)[0]
        optimizer_logs=[]
        if kind=='direct_physical':
            res=minimize(lambda z:self.physical(z,Y),x,jac=True,method='L-BFGS-B',bounds=[self.box]*self.N,options=dict(maxiter=outer*inner,ftol=1e-12,gtol=1e-8,maxls=25))
            x=res.x;best=x.copy();optimizer_logs.append(dict(stage='direct',success=bool(res.success),message=str(res.message),nit=int(res.nit),nfev=int(res.nfev)))
        else:
            for it in range(outer):
                J=self.current_update(x,Y);Q=self.inc+self.G@J
                res=minimize(lambda z:self.material(z,J,Q,kind),x,jac=True,method='L-BFGS-B',bounds=[self.box]*self.N,options=dict(maxiter=inner,ftol=1e-12,gtol=1e-8,maxls=25))
                x=res.x
                phi,df,_,_,lin=self.physical(x,Y,False)
                history.append(dict(cycle=it,physical_phi=phi,data_residual_squared=df,material_loss=float(res.fun),material_nit=int(res.nit),material_nfev=int(res.nfev),linear_residual=lin))
                if phi<bestphi:bestphi=phi;best=x.copy()
                optimizer_logs.append(dict(stage='material',cycle=it,success=bool(res.success),message=str(res.message),nit=int(res.nit),nfev=int(res.nfev)))
        raw=best.copy();pre=self.physical(raw,Y,False)
        res=minimize(lambda z:self.physical(z,Y),raw,jac=True,method='L-BFGS-B',bounds=[self.box]*self.N,options=dict(maxiter=corrector,ftol=1e-12,gtol=1e-8,maxls=25))
        post=self.physical(res.x,Y,False);final=res.x if post[0]<=pre[0] else raw
        value,df,J,E,lin=self.physical(final,Y,False)
        optimizer_logs.append(dict(stage='common_physical_corrector',success=bool(res.success),message=str(res.message),nit=int(res.nit),nfev=int(res.nfev)))
        receipt=dict(method=kind,status='COMPLETED',selection='minimum public-data physical Phi before common corrector; truth unavailable',counts=dict(self.counts),wall_s=time.perf_counter()-t,cpu_s=time.process_time()-c,pre_corrector_phi=pre[0],pre_corrector_data_residual=float(np.sqrt(pre[1])),physical_phi=value,physical_data_residual=float(np.sqrt(df)),final_linear_residual=lin,optimizer=optimizer_logs,history=history)
        return dict(initial=x0,raw=raw,prediction=final,J=J,E=E),receipt
