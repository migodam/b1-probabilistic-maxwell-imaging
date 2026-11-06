"""Reduced full-vector DDA/FFT-TSOM adaptations, not a paper reproduction.
Dense Green propagation and exact reduced current LS; FFT constructs APCS bases.
"""
import time
import numpy as np
from scipy import linalg as la
from scipy.optimize import minimize
from scipy.sparse import csr_matrix
from research.research_max_struct_r1.src.vector_backend import TensorDDA
from research.research_max_struct_r1.src.a9_green import dyadic_green

class VectorCIE3D:
    def __init__(self,ns=6,width=.6,lossy=False):
        t=time.perf_counter();c=time.process_time();self.ns=ns;self.N=ns**3;self.M=3*self.N;self.h=width/ns;self.lossy=lossy
        self.ijk=np.indices((ns,)*3).reshape(3,-1).T;self.points=(self.ijk-(ns-1)/2)*self.h
        self.backend=TensorDDA(self.points,self.h**3);self.G=self.backend.G;self.s=self.backend.self_density
        a=2*np.pi*np.arange(20)/20;r=.6
        circles=[];pol=[]
        for fixed,center in [(2,[.04,0,-.02]),(0,[.02,0,-.03]),(1,[-.01,.02,0])]:
            q=np.zeros((20,3));axes=[j for j in range(3) if j!=fixed];q[:,axes[0]]=r*np.cos(a);q[:,axes[1]]=r*np.sin(a);q+=center;circles.append(q);e=np.zeros((20,3));e[:,fixed]=1;pol.append(e)
        self.receivers=np.concatenate(circles);self.polarizations=np.concatenate(pol);self.P=60
        ge=dyadic_green(self.points[:,None,:]-self.receivers[None,:,:],2*np.pi)
        self.inc=np.einsum('nspq,sq->nps',ge,self.polarizations).reshape(self.M,self.P)
        self.S=self.backend.external_electric_operator(self.receivers)
        block=self.ijk//2;same=np.all(block[:,None,:]==block[None,:,:],axis=-1)
        self.H=(self.G-self.s*np.eye(self.M))*np.repeat(np.repeat(same,3,0),3,1)
        self.Hop=csr_matrix(self.H)
        self.edges=np.array([(i,j) for i in range(self.N) for j in range(i+1,self.N) if np.sum(abs(self.ijk[i]-self.ijk[j]))==1])
        self.lam=2e-4;self.l2=1e-6;self.rho=10.;self.E2=np.sum(abs(self.inc)**2,axis=0);self.reset()
        self.setup=dict(wall_s=time.perf_counter()-t,cpu_s=time.process_time()-c,green_assemblies=1,external_assemblies=1,source_assemblies=1,N=self.N,vector_dimension=self.M)
    def reset(self):self.counts=dict(physical_LU=0,source_RHS=0,adjoint_RHS=0,reduced_current_Cholesky=0,current_RHS=0,material_calls=0,receiver_SVD=0,fft_basis_columns=0)
    def decode(self,z):return z[:self.N]+1j*z[self.N:] if self.lossy else np.asarray(z,complex)
    def encode(self,x):return np.r_[x.real,x.imag] if self.lossy else x.real
    def grad(self,w):return np.r_[2*w.real,-2*w.imag] if self.lossy else 2*w.real
    def prior(self,z):
        x=self.decode(z);a,b=self.edges.T;d=x[a]-x[b];g=self.l2*x/self.N;v=self.lam*d/len(a);np.add.at(g,a,v);np.add.at(g,b,-v)
        return self.lam*np.mean(abs(d)**2)+self.l2*np.mean(abs(x)**2),self.grad(np.conj(g))
    def physical(self,z,Y=None,gradient=True):
        x=self.decode(z);d=np.repeat(x,3);K=np.eye(self.M)-d[:,None]*self.G;fac=la.lu_factor(K,check_finite=False);J=la.lu_solve(fac,d[:,None]*self.inc,check_finite=False)
        self.counts['physical_LU']+=1;self.counts['source_RHS']+=self.P
        E=self.inc+self.G@J;pred=self.S@J
        if Y is None:return pred,J,E
        Y2=np.sum(abs(Y)**2,axis=0);D=pred-Y;f0=float(np.mean(np.sum(abs(D)**2,axis=0)/Y2));pr,pg=self.prior(z)
        if not gradient:return f0+pr,f0,J,E
        adj=la.lu_solve(fac,self.S.conj().T@(D/Y2[None,:]/self.P),trans=2,check_finite=False);self.counts['adjoint_RHS']+=self.P
        w=np.sum(np.conj(adj)*E,axis=1).reshape(self.N,3).sum(1)
        return f0+pr,self.grad(w)+pg
    def prepare(self,Y):
        t=time.perf_counter();u,s,vh=la.svd(self.S,full_matrices=False,check_finite=False);rank=int(np.sum(s>=.1*s[0]));self.V=vh[:rank].conj().T;self.Jd=self.V@((u[:,:rank].conj().T@Y)/s[:rank,None]);self.bases={}
        for mf in [1,2]:
            freq=np.fft.fftfreq(self.ns)*self.ns;ids=np.argwhere(np.all(np.stack(np.meshgrid(*[np.isin(freq,np.r_[np.arange(mf),np.arange(-mf,0)])]*3,indexing='ij'),-1),-1));cols=[]
            for idx in ids:
                coeff=np.zeros((self.ns,)*3,complex);coeff[tuple(idx)]=1;f=np.fft.ifftn(coeff,norm='ortho').reshape(-1)
                for component in range(3):
                    q=np.zeros((self.N,3),complex);q[:,component]=f;cols.append(q.reshape(-1))
            B=np.array(cols).T;B-=self.V@(self.V.conj().T@B)
            # SVD drops numerical zero columns only; same exact affine family.
            ub,sb,_=la.svd(B,full_matrices=False,check_finite=False);B=ub[:,sb>1e-10];self.bases[mf]=B
        self.counts['receiver_SVD']+=1;self.counts['fft_basis_columns']+=sum(b.shape[1] for b in self.bases.values())
        return dict(rank=rank,svd_threshold_relative=.1,singular_values=s.tolist(),basis_complex_dimensions={str(k):v.shape[1] for k,v in self.bases.items()},wall_s=time.perf_counter()-t,implementation='actual 3D inverse FFT basis construction; dense Maxwell propagation; exact reduced quadratic current solve; no FFT complexity claim')
    def material(self,z,J,kind,beta):
        self.counts['material_calls']+=1;x=self.decode(z);d=np.repeat(x,3);Q=self.inc+self.G@J;r=J-d[:,None]*Q
        if kind=='LSIE_TSOM':R=r;du=-Q;w=np.sum(np.conj(R)*du/self.E2,axis=1)
        else:
            b=(beta-self.s) if kind=='A_TSOM' else beta;den=1+b*d;u=beta*r/den[:,None];du=-beta*(Q+b*J)/den[:,None]**2
            if kind=='A_TSOM':
                aa=d/den;Hu=self.Hop@u;R=u+aa[:,None]*Hu;v=R+self.Hop.conj().T@(np.conj(aa)[:,None]*R)
                w=np.sum((np.conj(v)*du+np.conj(R)*Hu/den[:,None]**2)/self.E2,axis=1)
            else:R=u;w=np.sum(np.conj(R)*du/self.E2,axis=1)
        pr,pg=self.prior(z);return self.rho/self.P*float(np.sum(abs(R)**2/self.E2))+pr,self.rho/self.P*self.grad(w.reshape(self.N,3).sum(1))+pg
    def current(self,z,Y,kind,beta,mf):
        x=self.decode(z);d=np.repeat(x,3);K=np.eye(self.M)-d[:,None]*self.G;rhs=d[:,None]*self.inc
        if kind!='LSIE_TSOM':
            b=(beta-self.s) if kind=='A_TSOM' else beta;t=beta/(1+b*d);K=t[:,None]*K;rhs=t[:,None]*rhs
            if kind=='A_TSOM':
                aa=d/(1+b*d);K=K+aa[:,None]*(self.Hop@K);rhs=rhs+aa[:,None]*(self.Hop@rhs)
        B=self.bases[mf];SB=self.S@B;KB=K@B;D=Y-self.S@self.Jd;R=rhs-K@self.Jd;Y2=np.sum(abs(Y)**2,axis=0)
        # Different per-source normalization entails genuinely distinct reduced normals.
        J=self.Jd.copy();Gs=SB.conj().T@SB;Gk=KB.conj().T@KB;rs=SB.conj().T@D;rk=KB.conj().T@R
        for p in range(self.P):
            N=Gs/Y2[p]+self.rho*Gk/self.E2[p];rr=rs[:,p]/Y2[p]+self.rho*rk[:,p]/self.E2[p]
            gam=la.cho_solve(la.cho_factor(N,check_finite=False),rr,check_finite=False);J[:,p]+=B@gam
        self.counts['reduced_current_Cholesky']+=self.P;self.counts['current_RHS']+=self.P;return J
    def optimize(self,Y,kind,outer=16,inner=12,corrector=300):
        self.reset();t=time.perf_counter();c=time.process_time();z=self.encode(np.full(self.N,.25+(.05j if self.lossy else 0)));initial=z.copy();best=z.copy();bestphi=self.physical(z,Y,False)[0];logs=[];curve=[];trials=[];bounds=[(0,4)]*self.N+([(0,1.2)]*self.N if self.lossy else [])
        def callback(x):
            nonlocal best,bestphi
            phi=self.physical(x,Y,False)[0];trials.append(x.copy());curve.append(dict(phi=phi,physical_LU=self.counts['physical_LU'],cpu_s=time.process_time()-c))
            if phi<bestphi:bestphi=phi;best=x.copy()
        prep=None
        if kind=='direct_physical':
            res=minimize(lambda q:self.physical(q,Y),z,jac=True,bounds=bounds,method='L-BFGS-B',callback=callback,options=dict(maxiter=outer*inner,ftol=1e-12,gtol=1e-8,maxls=25));logs.append(dict(stage='direct',nit=res.nit,nfev=res.nfev,message=str(res.message),success=bool(res.success)));callback(res.x)
        else:
            prep=self.prepare(Y)
            for it in range(outer):
                mf=1 if it<outer//2 else 2;beta=6 if it<outer//2 else 1;beta_current=(3 if kind=='CIE_dualbeta_TSOM' and it>=outer//2 else beta)
                J=self.current(z,Y,kind,beta_current,mf);res=minimize(lambda q:self.material(q,J,kind,beta),z,jac=True,bounds=bounds,method='L-BFGS-B',options=dict(maxiter=inner,ftol=1e-12,gtol=1e-8,maxls=25));z=res.x;callback(z);logs.append(dict(stage='joint',cycle=it,mf=mf,beta_material=beta,beta_current=beta_current,nit=res.nit,nfev=res.nfev,message=str(res.message),success=bool(res.success)))
        raw=best.copy();pre=self.physical(raw,Y,False);res=minimize(lambda q:self.physical(q,Y),raw,jac=True,bounds=bounds,method='L-BFGS-B',callback=callback,options=dict(maxiter=corrector,ftol=1e-13,gtol=1e-9,maxls=30));callback(res.x);final=best;f,df,J,E=self.physical(final,Y,False);logs.append(dict(stage='common_long_corrector',nit=res.nit,nfev=res.nfev,message=str(res.message),success=bool(res.success)))
        return dict(initial=self.decode(initial),raw=self.decode(raw),prediction=self.decode(final),J=J,E=E,trials=np.array([self.decode(q) for q in trials])),dict(method=kind,wall_s=time.perf_counter()-t,cpu_s=time.process_time()-c,counts=self.counts.copy(),physical_phi=f,physical_data_residual=np.sqrt(df),pre_corrector_phi=pre[0],pre_corrector_residual=np.sqrt(pre[1]),optimizer=logs,curve=curve,subspace=prep,truth_access='NONE during all optimization and selection')

def box_fraction(points,h,center,half):
    lo=np.maximum(points-h/2,np.array(center)-np.array(half));hi=np.minimum(points+h/2,np.array(center)+np.array(half));return np.prod(np.maximum(hi-lo,0),axis=1)/h**3

def objects(model):
    p=model.points;h=model.h;four=np.zeros(model.N)
    for center,chi in [([-.14,.14,-.14],1),([.14,-.14,-.14],2),([-.14,-.14,.14],3),([.14,.14,.14],4)]:four+=chi*box_fraction(p,h,center,[.08]*3)
    outer=box_fraction(p,h,[0,0,0],[.2]*3);inner=box_fraction(p,h,[0,0,0],[.1]*3);coat=(.5+.3j)*outer+(.5+.5j)*inner
    # Explicit adapted Austria dimensions, as numerical dimensions absent source prose.
    hollow=box_fraction(p,h,[0,.05,0],[.15,.13,.15])-box_fraction(p,h,[0,.05,0],[.09,.08,.09])
    for x in [-.2,.2]:hollow+=box_fraction(p,h,[x,-.2,0],[.035,.035,.18])
    return dict(four_cubes=four,hollow_cube_rods=1.5*hollow,lossy_coated_cube=coat)
