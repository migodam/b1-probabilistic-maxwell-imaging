"""Full nonlinear pulse-collocation TM VIE with counted source/adjoint actions."""
import itertools,json,time
import numpy as np
from scipy.linalg import lu_factor,lu_solve,null_space
from scipy.special import hankel1,jv
from accounting import P,write

class Model:
    def __init__(self,book,centers,sources,receivers,radius=.04,k=2*np.pi):
        self.book=book;self.n=len(centers);self.ns=len(sources);self.nr=len(receivers);self.k=k;self.radius=radius
        self.centers=np.asarray(centers,dtype=np.float64);self.sources=np.asarray(sources,dtype=np.float64);self.receivers=np.asarray(receivers,dtype=np.float64)
        self.state=None;self.factor=None
        with book.op('HANKEL_GREEN_SETUP',Green_kernel_evaluations=self.n**2+self.n*self.ns+self.nr*self.n,receiver_null_SVD=1,algebraic_setup=1):
            d=np.linalg.norm(self.centers[:,None]-self.centers[None,:],axis=2)
            assert np.min(d+np.eye(self.n)*100)>2*radius
            np.fill_diagonal(d,1.0)
            self.G=self.external(d)
            np.fill_diagonal(self.G,.5j*np.pi*k*radius*hankel1(1,k*radius)-1)
            self.E0=.25j*hankel1(0,k*np.linalg.norm(self.centers[:,None]-self.sources[None,:],axis=2))
            self.S=self.external(np.linalg.norm(self.receivers[:,None]-self.centers[None,:],axis=2))
            self.N=null_space(self.S)
    def external(self,d):return .5j*np.pi*self.k*self.radius*jv(1,self.k*self.radius)*hankel1(0,self.k*d)
    def forward(self,m,label='state',jac=False,independent=False):
        m=np.asarray(m,dtype=np.float64);assert m.shape==(self.n,)
        # No unbounded live-LU dictionary. Every attempted new state is billed.
        self.factor=None;self.state=None
        with self.book.op('INDEPENDENT_FORWARD' if independent else 'FULL_FORWARD',rhs=self.ns,LU=1,new_states=1,Green_matvec_columns=self.ns,label=label,material_legal=bool(np.all((m>=.5)&(m<=4)))) as action:
            tick=time.perf_counter()
            L=np.eye(self.n,dtype=np.complex128)-m[:,None]*self.G
            rhs=m[:,None]*self.E0
            if independent:j=np.linalg.solve(L,rhs);fac=None
            else:fac=lu_factor(L,check_finite=False);j=lu_solve(fac,rhs,check_finite=False)
            E=self.E0+self.G@j;y=self.S@j
            state_res=np.linalg.norm(L@j-rhs)/max(np.linalg.norm(rhs),1e-300)
            sigma_min=np.linalg.svd(L,compute_uv=False)[-1]
            state=dict(m=m.copy(),j=j,E=E,y=y,L=L,state_relative_residual=state_res,resolvent_sigma_min=sigma_min,action=action)
            self.book.append(dict(event='PHYSICS_KERNEL',action=action,kernel_seconds=time.perf_counter()-tick,serialization_excluded=True))
            self.book.save(action,m=m,j=j,E=E,y=y,state_relative_residual=state_res,resolvent_singular_values=np.linalg.svd(L,compute_uv=False))
        self.state=state;self.factor=fac
        if jac:self.jacobian(state)
        return state
    def jacobian(self,state=None):
        state=state or self.state
        if 'J' in state:return state['J']
        if state is not self.state or self.factor is None:raise RuntimeError('J needs the one live LU at the same state')
        with self.book.op('FULL_MATERIAL_J_ADJOINT',rhs=self.nr,adjoint_RHS=self.nr,material_J_builds=1) as action:
            tick=time.perf_counter()
            # Transpose solve (not conjugate transpose): analytic complex derivative.
            T=lu_solve(self.factor,self.S.T,trans=1,check_finite=False).T
            J=np.einsum('ri,is->sri',T,state['E']).reshape(self.ns*self.nr,self.n)
            state['Jcomplex']=J;state['J']=self.pack_jac(J)
            self.book.append(dict(event='PHYSICS_KERNEL',action=action,kernel_seconds=time.perf_counter()-tick,serialization_excluded=True))
            self.book.save(action,m=state['m'],Jcomplex=J,Jreal=state['J'],receiver_transfer=T)
        return state['J']
    def pack(self,y):
        a=np.asarray(y).T
        return np.stack([a.real,a.imag],axis=-1).reshape(-1)
    def pack_jac(self,J):return np.stack([J.real,J.imag],axis=1).reshape(2*self.ns*self.nr,self.n)
    def permutation_audit(self):
        gaps=[]
        for p in itertools.permutations(range(self.n)):
            if p==tuple(range(self.n)):continue
            p=np.asarray(p)
            gap=max(np.linalg.norm(self.G[np.ix_(p,p)]-self.G)/np.linalg.norm(self.G),
                    np.linalg.norm(self.E0[p]-self.E0)/np.linalg.norm(self.E0),
                    np.linalg.norm(self.S[:,p]-self.S)/np.linalg.norm(self.S))
            gaps.append(gap)
        return dict(nonidentity_permutations=len(gaps),minimum_full_labeled_operator_defect=float(min(gaps)),exact_nontrivial_permutations_at_1e_12=int(np.sum(np.asarray(gaps)<1e-12)),incident_singular_values=np.linalg.svd(self.E0,compute_uv=False).tolist())

def model_for(book,stage='E1'):
    contract=json.loads((P/'EXPERIMENT_CONTRACT.json').read_text());c=contract['E1']
    centers=c['centers']
    if stage=='E2':centers=[centers[i] for i in contract['E2']['centers_from_E1']]
    return Model(book,centers,c['sources'],c['receivers'])

def quality(model,state,target,mref):
    m=state['m'];pair=np.linalg.norm(state['y']-target)
    norm=np.linalg.norm(target)
    return dict(material_RMS=float(np.linalg.norm(m-mref)/np.sqrt(model.n)),pair_relative=float(pair/max(norm,1e-300)),
        pair_whitened=float(pair/.001),Phi=float(pair**2/(2*.001**2)),state_relative=float(state['state_relative_residual']),
        sigma_min=float(state['resolvent_sigma_min']),legal=bool(np.all((m>=.5-1e-12)&(m<=4+1e-12))),
        qualified=bool(pair/max(norm,1e-300)<=1e-10 and state['state_relative_residual']<=1e-10 and np.all((m>=.5-1e-12)&(m<=4+1e-12)) and state['resolvent_sigma_min']>=1e-6))
