from physics import *
from scipy.optimize import brentq
def normcols(a):return (a.abs()**2).sum(0)
def eigspace(model,V):
    """Full restricted Hermitian spectrum, using identical spatial dyadic kernel."""
    n=model.n;N=model.N
    coord=torch.stack(torch.meshgrid(*(torch.arange(n,device=DEV),)*3,indexing='ij')).reshape(3,-1).T
    ker=torch.fft.ifftn(model.kernel,dim=(-3,-2,-1)).permute(2,3,4,0,1)
    dense=torch.empty((model.M,model.M),device=DEV,dtype=CD)
    for i in range(0,N,48):
        d=(coord[i:i+48,None,:]-coord[None,:,:])%(2*n)
        block=ker[d[:,:,0],d[:,:,1],d[:,:,2]].permute(0,2,1,3)
        dense[3*i:3*(i+len(block))]=block.reshape(3*len(block),model.M)
    # Largest eigenvalue of G*G estimated with converged deterministic Lanczos.
    from scipy.sparse.linalg import LinearOperator,eigsh
    def gram(a):
        v=t(a);return npy(dense.mH@(dense@v))
    rng=np.random.default_rng(5168)
    lmax=float(eigsh(LinearOperator((model.M,model.M),matvec=gram,dtype=np.complex64),k=1,which='LA',v0=rng.normal(size=model.M).astype(np.complex64),tol=1e-5,return_eigenvectors=False)[0])
    dense-= (dense@V)@V.mH
    herm=dense.mH@dense;del dense
    vals,U=torch.linalg.eigh(herm);sel=vals>.01*lmax
    model.count['route_dense_Gram']+=1;model.count['route_full_eigh']+=1
    return U[:,sel].contiguous(),dict(route_rank=int(sel.sum()),lmax=lmax,smallest_kept=float(vals[sel][0]),max_discard=float(vals[~sel][-1]))
def cg(mod,rhs,tol=3e-5):
    def A(a):return a+36*mod.G(mod.G(a),True)
    x=torch.zeros_like(rhs);r=rhs.clone();p=r.clone();rr=normcols(r);base=rr.clamp_min(1e-30)
    for i in range(120):
        v=A(p);a=rr/(p.conj()*v).sum(0).real.clamp_min(1e-30)
        x=x+p*a;r=r-v*a;new=normcols(r)
        if (new/base).max().item()<tol**2:break
        p=r+p*(new/rr.clamp_min(1e-30));rr=new
    err=(normcols(A(x)-rhs)/base).max().sqrt().item()
    mod.count['H_cg_iterations']+=i+1;mod.count['H_cg_rhs']+=rhs.shape[1]
    if err>3*tol:raise ArithmeticError(('CG',err))
    return x,err
def material(x,q):
    n=x.shape[0]//3;X=x.reshape(n,-1);Q=q.reshape(n,-1)
    a=normcols(X.T);c=normcols(Q.T);b=(X.conj()*Q).real.sum(1);d=c-a
    stable=-.5*(d+torch.copysign(torch.hypot(d,2*b),d));z=torch.zeros_like(d)
    candidates=torch.stack([z,z+4/6,stable/b,-b/stable],1)
    good=torch.isfinite(candidates)&(candidates>=0)&(candidates<=4/6)
    vals=torch.where(good,(a[:,None]-2*b[:,None]*candidates+c[:,None]*candidates**2)/(1+candidates**2),torch.inf)
    alpha=candidates.gather(1,vals.argmin(1)[:,None]);cq=(Q+alpha*X)/(1+alpha**2)
    return (alpha*cq).reshape_as(x),cq.reshape_as(q),6*alpha.flatten()
class Projector:
    def __init__(self,mod,V,Y,delta):
        start=time.perf_counter();self.mod=mod;self.Y=Y;self.delta=delta;self.T=6*mod.S;self.logs=[]
        self.U,self.info=eigspace(mod,V);self.Z,setup_err=cg(mod,self.T.mH)
        M=self.T@self.Z;ev,W=torch.linalg.eigh((M+M.mH)/2);self.ev=npy(ev).astype(float);self.W=W
        self.info.update(setup_cg=setup_err,schur_min=float(ev.min()),setup_s=time.perf_counter()-start)
    def route(self,v):return self.U@(self.U.mH@v)
    def affine(self,x,q):
        m=self.mod;xu,cgres=cg(m,x+6*m.G(q-m.E,True));r=self.T@xu-self.Y
        z=self.W.mH@r;powers=npy(normcols(z.T)).astype(float);ev=self.ev
        floor=ev<=max(ev.max()*1e-8,1e-30);target=self.delta**2;mu=0.
        if powers.sum()>target:
            if powers[floor].sum()>target:raise ArithmeticError('Unreachable data ball')
            solve=lambda s:sum(powers[floor])+np.sum(powers[~floor]/(1+s*ev[~floor])**2)-target
            hi=1/ev.max()
            while solve(hi)>0:hi*=4
            mu=brentq(solve,0,hi);weight=np.where(floor,0,1/(ev+1/mu))
            xu=xu-self.Z@(self.W@(t(weight)[:,None]*z))
        aq=m.E+6*m.G(xu);res=float(torch.linalg.norm(self.T@xu-self.Y))
        self.logs.append(dict(cg_residual=cgres,mu=mu,data_over_delta=res/self.delta))
        return xu,aq
