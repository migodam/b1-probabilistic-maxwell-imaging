"""R1 CM+RR DDA on uniform 3D grid using A18_C padded FFT propagation form.
No A18 state/hash machinery imported; geometry/scalar material newly specified.
"""
import time
import numpy as np
from scipy import fft,linalg as la
from scipy.sparse.linalg import LinearOperator,gmres
from research.research_max_struct_r1.src.a9_green import dyadic_green
from research.research_max_struct_r8.src.vector_cie3d import box_fraction
class VectorFFT3D:
    def __init__(self,ns=12,width=3.):
        tick=time.perf_counter();self.ns=ns;self.N=ns**3;self.M=3*self.N;self.h=width/ns;self.shape=(ns,)*3;self.pad=(2*ns,)*3;self.ijk=np.indices(self.shape).reshape(3,-1).T;self.points=(self.ijk-(ns-1)/2)*self.h;self.s=-1/3+1j*(2*np.pi)**3*self.h**3/(6*np.pi)
        idx=np.arange(2*ns);idx=np.where(idx<=ns,idx,idx-2*ns);dv=np.stack(np.meshgrid(idx,idx,idx,indexing='ij'),-1)*self.h;zero=np.all(dv==0,axis=-1);dv[zero]=[1,0,0];g=self.h**3*dyadic_green(dv,2*np.pi);g[zero]=self.s*np.eye(3);self.kernel=fft.fftn(g,axes=(0,1,2));self.counts=dict(fft_matvec=0,primal_RHS=0,adjoint_RHS=0,gmres_iterations=0,failed_rhs=0,physical_evaluations=0)
        a=2*np.pi*np.arange(20)/20;circles=[];pol=[]
        for fixed,center in [(2,[.2,0,-.1]),(0,[.1,0,-.15]),(1,[-.05,.1,0])]:
            q=np.zeros((20,3));axes=[j for j in range(3) if j!=fixed];q[:,axes[0]]=3*np.cos(a);q[:,axes[1]]=3*np.sin(a);q+=center;circles.append(q);e=np.zeros((20,3));e[:,fixed]=1;pol.append(e)
        self.receivers=np.concatenate(circles);self.polarizations=np.concatenate(pol);self.P=60;gg=dyadic_green(self.points[:,None,:]-self.receivers[None,:,:],2*np.pi);self.inc=np.einsum('nspq,sq->nps',gg,self.polarizations).reshape(self.M,self.P);gg=dyadic_green(self.receivers[:,None,:]-self.points[None,:,:],2*np.pi);self.S=self.h**3*gg.transpose(0,2,1,3).reshape(180,self.M);self.setup_wall=time.perf_counter()-tick;self.max_residual=0.;self.failures=[]
    def G(self,v,adj=False):
        self.counts['fft_matvec']+=1;z=v.reshape(self.shape+(3,));hat=fft.fftn(z,s=self.pad,axes=(0,1,2));kg=self.kernel.conj().swapaxes(-1,-2) if adj else self.kernel;out=fft.ifftn(np.einsum('...ij,...j->...i',kg,hat),axes=(0,1,2))[:self.ns,:self.ns,:self.ns].reshape(-1);return out
    def solve(self,x,rhs,adj=False):
        d=np.repeat(x,3);op=LinearOperator((self.M,self.M),matvec=(lambda v:v-self.G(d.conj()*v,True)) if adj else (lambda v:v-d*self.G(v)),dtype=complex);diag=1-d*self.s;pre=LinearOperator((self.M,self.M),matvec=lambda v:v/(diag.conj() if adj else diag),dtype=complex);sol=np.empty_like(rhs,dtype=complex)
        for p in range(rhs.shape[1]):
            it=[0]
            def cb(r):it[0]+=1
            sol[:,p],info=gmres(op,rhs[:,p],M=pre,rtol=2e-6,atol=0,restart=60,maxiter=8,callback=cb,callback_type='pr_norm');res=la.norm(op@sol[:,p]-rhs[:,p])/max(la.norm(rhs[:,p]),1e-30);self.max_residual=max(self.max_residual,res);self.counts['gmres_iterations']+=it[0];self.counts['adjoint_RHS' if adj else 'primal_RHS']+=1
            if info:self.counts['failed_rhs']+=1;self.failures.append(dict(adjoint=adj,source=p,info=info,residual=res));raise RuntimeError('GMRES cap failed: '+str(self.failures[-1]))
        return sol
    def physical(self,x,Y=None,gradient=True):
        self.counts['physical_evaluations']+=1;d=np.repeat(x,3);J=self.solve(x,d[:,None]*self.inc);pred=self.S@J
        if Y is None:return pred,J
        Y2=np.sum(abs(Y)**2,axis=0);D=pred-Y;f=float(np.mean(np.sum(abs(D)**2,axis=0)/Y2));lam=2e-5;l2=1e-6;cube=x.reshape(self.shape);prior=l2*np.mean(x*x);pg=2*l2*x/self.N
        for ax in range(3):
            delta=np.diff(cube,axis=ax);prior+=lam*np.mean(delta*delta);gg=2*lam*delta/delta.size;lo=[slice(None)]*3;hi=lo.copy();lo[ax]=slice(None,-1);hi[ax]=slice(1,None);pv=pg.reshape(self.shape);pv[tuple(lo)]-=gg;pv[tuple(hi)]+=gg
        if not gradient:return f+prior,f,J
        adj=self.solve(x,self.S.conj().T@(D/Y2[None,:]/self.P),True);E=self.inc+np.column_stack([self.G(J[:,p]) for p in range(self.P)]);w=2*np.real(np.sum(np.conj(adj)*E,axis=1)).reshape(self.N,3).sum(1)+pg;return f+prior,w
    def four_cubes(self):
        x=np.zeros(self.N)
        for center,chi in [([-.7,.7,-.7],1),([.7,-.7,-.7],2),([-.7,-.7,.7],3),([.7,.7,.7],4)]:x+=chi*box_fraction(self.points,self.h,center,[.4]*3)
        return x
