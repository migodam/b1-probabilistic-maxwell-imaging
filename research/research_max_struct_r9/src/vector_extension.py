"""R9 Phase2 strict real-vector adapter over R1/R8 physics and readouts.

No new Green/factorization backend. Current is ordinary CIE for every CIE/A
arm; only the material readout changes for A/A+. No truth loading here.
"""
import time
from collections import Counter
import numpy as np
from scipy import linalg as la
from scipy.optimize import minimize
from scipy.sparse import csr_matrix
from scipy.stats import chi2
from research.research_max_struct_r1.src.vector_backend import TensorDDA
from research.research_max_struct_r1.src.a9_green import dyadic_green
from research.research_max_struct_r8.src.vector_cie3d import box_fraction

ARMS=('cie_single','dual_source','global_dual','a_material','a_plus','direct',
      'global_dual_qualified','a_plus_qualified')

def paid_density_state(G,inc,x,S,cost):
    """Unchanged R8 density convention, retained LU/all-source RHS."""
    d=np.repeat(np.asarray(x,float),3);M=len(d);P=inc.shape[1]
    K=np.eye(M,dtype=complex)-d[:,None]*G
    factor=la.lu_factor(K,check_finite=False)
    J=la.lu_solve(factor,d[:,None]*inc,check_finite=False)
    cost['physical_LU']+=1;cost['physical_source_RHS']+=P
    E=inc+G@J;cost['physical_G_columns']+=P
    pred=None if S is None else S@J
    if S is not None:cost['physical_S_columns']+=P
    lin=float(la.norm(J-d[:,None]*E)/max(la.norm(d[:,None]*inc),1e-30))
    return dict(m=np.asarray(x,float).copy(),factor=factor,J=J,E=E,Y=pred,
                linear_relative_residual=lin)

class StrictVectorInverse:
    def __init__(self,nside,width,config,auxiliary=True):
        wall,cpu=time.perf_counter(),time.process_time();self.nside=nside;self.N=nside**3;self.M=3*self.N
        self.width=width;self.h=width/nside;self.config=config;self.solver=config['solver']
        self.lo,self.hi=map(float,config['physics']['material_box']);self.P=60
        self.ijk=np.indices((nside,)*3).reshape(3,-1).T;self.points=(self.ijk-(nside-1)/2)*self.h
        self.backend=TensorDDA(self.points,self.h**3);self.G=self.backend.G;self.s=self.backend.self_density
        angle=2*np.pi*np.arange(20)/20;circles=[];polarizations=[]
        for fixed,center in zip((2,0,1),config['acquisition']['original_3lambda_circle_centers']):
            q=np.zeros((20,3));axes=[j for j in range(3) if j!=fixed]
            q[:,axes[0]]=width*np.cos(angle);q[:,axes[1]]=width*np.sin(angle)
            q+=np.asarray(center)*width/3;circles.append(q)
            p=np.zeros((20,3));p[:,fixed]=1;polarizations.append(p)
        self.receivers=np.concatenate(circles);self.polarizations=np.concatenate(polarizations)
        if np.any(np.max(abs(self.receivers),axis=1)<=width/2):raise ValueError('Source/rx inside DoI')
        electric=dyadic_green(self.points[:,None,:]-self.receivers[None,:,:],2*np.pi)
        self.inc=np.einsum('nspq,sq->nps',electric,self.polarizations).reshape(self.M,self.P)
        self.S=self.backend.external_electric_operator(self.receivers);self.E2=np.sum(abs(self.inc)**2,axis=0)
        ids=np.arange(self.N).reshape((nside,)*3);edges=[]
        self.neighbors=np.full((self.N,6),self.N,int)
        for axis in range(3):
            a=[slice(None)]*3;b=a.copy();a[axis]=slice(None,-1);b[axis]=slice(1,None)
            left,right=ids[tuple(a)].ravel(),ids[tuple(b)].ravel();edges.append(np.c_[left,right])
            self.neighbors[left,2*axis]=right;self.neighbors[right,2*axis+1]=left
        self.edges=np.concatenate(edges);self.degree=np.sum(self.neighbors<self.N,axis=1)
        self.lam=self.solver['prior_smooth'];self.l2=self.solver['prior_l2'];self.rho=self.solver['state_penalty_rho']
        self.Hop=self.Hblocks=self.block_cells=self.block_rows=None
        if auxiliary:
            if nside%2:raise ValueError('Frozen2³localblocks require even nside')
            bases=np.indices((nside//2,)*3).reshape(3,-1).T*2
            colors=np.indices((2,)*3).reshape(3,-1).T
            coords=bases[:,None,:]+colors[None,:,:]
            self.block_cells=np.ravel_multi_index(coords.transpose(2,0,1),(nside,)*3)
            self.block_rows=np.repeat(self.block_cells,3,axis=1)*3+np.tile(np.arange(3),8)
            rr=self.block_rows;self.Hblocks=self.G[rr[:,:,None],rr[:,None,:]]-self.s*np.eye(24)[None]
            rows=np.broadcast_to(rr[:,:,None],self.Hblocks.shape).ravel()
            cols=np.broadcast_to(rr[:,None,:],self.Hblocks.shape).ravel()
            self.Hop=csr_matrix((self.Hblocks.ravel(),(rows,cols)),shape=(self.M,self.M));self.Hop.eliminate_zeros()
        self.counts=Counter();self.state_cache=None;self.Y=None;self.best=None
        self.setup=dict(wall_seconds=time.perf_counter()-wall,CPU_seconds=time.process_time()-cpu,
            material_dof=self.N,vector_dimension=self.M,sources=self.P,receiver_positions=60,complex_data=10800,
            Green_assemblies=1,source_operator_assemblies=1,receiver_operator_assemblies=1,
            sparse_local_H=auxiliary,source_radius=width,source_centers_scale=width/3)

    def prior(self,x):
        a,b=self.edges.T;d=x[a]-x[b];f=self.lam*np.mean(d*d)+self.l2*np.mean(x*x)
        g=2*self.l2*x/self.N;v=2*self.lam*d/len(a);np.add.at(g,a,v);np.add.at(g,b,-v)
        return float(f),g

    def state(self,x):
        x=np.asarray(x,float)
        if x.shape!=(self.N,) or not np.isfinite(x).all():raise ValueError('Invalid full real material')
        if np.min(x)<self.lo-1e-12 or np.max(x)>self.hi+1e-12:raise ValueError('Material outside frozen shared box')
        if self.state_cache is not None and np.array_equal(x,self.state_cache['m']):
            self.counts['physical_state_cache_hits']+=1;return self.state_cache
        self.state_cache=paid_density_state(self.G,self.inc,x,self.S,self.counts)
        return self.state_cache

    def evaluate(self,x,gradient=False):
        self.counts['physical_gradient_requests' if gradient else 'physical_value_requests']+=1
        st=self.state(x)
        if 'Phi' not in st:
            D=st['Y']-self.Y;pr,pg=self.prior(st['m'])
            data=float(np.mean(np.sum(abs(D)**2,axis=0)*self.Yweights))
            st.update(Phi=data+pr,data_squared=data,prior=pr,prior_gradient=pg)
        if gradient and 'gradient' not in st:
            delta=st['Y']-self.Y
            adj=la.lu_solve(st['factor'],self.S.conj().T@(delta*self.Yweights[None,:]/self.P),trans=2,check_finite=False)
            self.counts['physical_adjoint_RHS']+=self.P
            st['gradient']=2*np.real(np.sum(np.conj(adj)*st['E'],axis=1)).reshape(self.N,3).sum(1)+st['prior_gradient']
        elif gradient:self.counts['physical_gradient_cache_hits']+=1
        self.observe(st)
        return (st['Phi'],st['gradient']) if gradient else st

    def observe(self,st):
        def snapshot(s):
            return {k:(v.copy() if isinstance(v,np.ndarray) else v) for k,v in s.items()
                    if k in ('m','J','E','Y','Phi','data_squared','gradient','linear_relative_residual')}
        if self.best is None or st['Phi']<self.best['Phi']:
            self.best=snapshot(st)
        elif np.array_equal(st['m'],self.best['m']) and 'gradient' in st:
            self.best['gradient']=st['gradient'].copy()
        signature=(self.counts['physical_LU'],self.counts['physical_adjoint_RHS'])
        if signature!=self.last_budget_signature:
            self.budget.append(dict(stage=self.phase,CPU_seconds=time.process_time()-self.online_cpu,
                wall_seconds=time.perf_counter()-self.online_wall,LU=signature[0],source_RHS=self.counts['physical_source_RHS'],
                adjoint_RHS=signature[1],current_Phi=st['Phi'],original_Phi=self.best['Phi']))
            self.budget_material.append(self.best['m'].copy());self.last_budget_signature=signature

    def prepare(self,noise_meta,qualified):
        wall,cpu=time.perf_counter(),time.process_time()
        U,s,Vh=la.svd(self.S,full_matrices=False,check_finite=False);self.counts['receiver_SVD']+=1
        threshold=np.finfo(float).eps*max(self.S.shape)*s[0];supported=int(np.sum(s>threshold))
        for key,value in noise_meta.items():
            if any(word in key.lower() for word in ('whiten','covariance','precision')) and value not in (None,False):
                raise ValueError('Unsupported noise whitening/covariance flag')
        decision=dict(rule='fixed_010',threshold_relative=.1,supported_rank=supported,singular_values=s.tolist())
        if not qualified:rank=int(np.sum(s>=.1*s[0]))
        else:
            if noise_meta['noise_model']=='legacy_fixed_energy_bounded':
                radius2=float((.1/.9)**2*np.vdot(self.Y,self.Y).real)
                decision.update(rule='legacy_deterministic_noise_bound',radius_formula='(.1/.9)^2*||Y||F^2; no chi-square')
            elif noise_meta['noise_model']=='iid_circular_complex':
                sigma=np.asarray(noise_meta['sigma'])
                if sigma.size!=1 or np.iscomplexobj(sigma) or float(sigma)<0:raise ValueError('Unsupported non-scalar sigma')
                radius2=float(float(sigma)**2*chi2.ppf(.99,2*self.Y.size)/2)
                decision.update(rule='iid_CN_full_noise_ball_99',sigma=float(sigma),real_noise_dof=2*self.Y.size)
            else:raise ValueError('Unsupported noise/whitening model')
            coord=U.conj().T@self.Y
            outside=0. if U.shape[1]==self.Y.shape[0] else float(la.norm(self.Y-U@coord)**2)
            tails=outside+np.r_[np.cumsum(np.sum(abs(coord)**2,axis=1)[::-1])[::-1],0.]
            eligible=np.flatnonzero(tails[:supported+1]<=radius2)
            decision.update(full_radius_squared=radius2,discarded_energy_squared_by_rank=tails.tolist(),
                scope='noise-compatible hardlock only; discretization mismatch NOT bounded; no GT/per-r testing')
            if not len(eligible):raise ValueError('No numerically supported DPCS prefix meets frozen noise-only radius')
            rank=int(eligible[0])
        if rank>supported:raise ValueError('Fixed DPCS includes numerically unsupported mode')
        self.V=Vh[:rank].conj().T
        self.Jd=self.V@((U[:,:rank].conj().T@self.Y)/s[:rank,None])
        self.GJd=self.G@self.Jd;self.SJd=self.S@self.Jd
        self.counts['prepare_G_columns']+=self.P;self.counts['prepare_S_columns']+=self.P
        self.bases={};freq=np.fft.fftfreq(self.nside)*self.nside
        for mf in self.solver['MF']:
            mask=np.isin(freq,np.r_[np.arange(mf),np.arange(-mf,0)])
            indices=np.argwhere(np.all(np.stack(np.meshgrid(mask,mask,mask,indexing='ij'),axis=-1),axis=-1))
            columns=[];keys=[]
            for idx in indices:
                coefficient=np.zeros((self.nside,)*3,complex);coefficient[tuple(idx)]=1
                field=np.fft.ifftn(coefficient,norm='ortho').ravel();self.counts['basis_IFFT_columns']+=1
                for component in range(3):
                    col=np.zeros((self.N,3),complex);col[:,component]=field;columns.append(col.ravel())
                    keys.append(tuple(int(round(freq[i])) for i in idx)+(component,))
            B=np.array(columns).T;B-=self.V@(self.V.conj().T@B)
            self.bases[mf]=dict(B=B,GB=self.G@B,SB=self.S@B,keys=keys)
            self.counts['prepare_G_columns']+=B.shape[1];self.counts['prepare_S_columns']+=B.shape[1]
        decision.update(rank=rank,whiten_apcs=False,basis_complex_dimensions={str(k):v['B'].shape[1] for k,v in self.bases.items()},
            CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall)
        return decision

    def current_step(self,x,gamma,old_gradient,old_direction,beta,stage,sourcewise):
        B,GB,SB=(stage[k] for k in ('B','GB','SB'));d=np.repeat(x,3);p=beta/(1+beta*d)
        KB=p[:,None]*(B-d[:,None]*GB)
        T=p[:,None]*(self.Jd-d[:,None]*(self.inc+self.GJd))
        data=self.SJd+SB@gamma-self.Y;state=T+KB@gamma
        grad=SB.conj().T@(data*self.Yweights)+self.rho*KB.conj().T@(state/self.E2)
        axis=0 if sourcewise else None
        if old_gradient is None:direction=-grad
        else:
            pr=np.maximum(0.,np.real(np.sum(np.conj(grad)*(grad-old_gradient),axis=axis))/np.maximum(np.sum(abs(old_gradient)**2,axis=axis),1e-300))
            direction=-grad+old_direction*pr
        slope=np.real(np.sum(np.conj(grad)*direction,axis=axis))
        if sourcewise:
            bad=slope>=0;direction[:,bad]=-grad[:,bad];self.counts['source_PR_restarts']+=int(bad.sum())
        elif slope>=0:direction=-grad;self.counts['global_PR_restarts']+=1
        slope=np.real(np.sum(np.conj(grad)*direction,axis=axis));ds=SB@direction;dk=KB@direction
        denom=np.sum(abs(ds)**2,axis=0)*self.Yweights+self.rho*np.sum(abs(dk)**2,axis=0)/self.E2
        if not sourcewise:denom=float(denom.sum())
        alpha=np.divide(-slope,denom,out=np.zeros_like(slope),where=denom>1e-300)
        change=direction*alpha;ratio=np.sum(abs(change)**2,axis=0)/np.maximum(np.sum(abs(gamma)**2,axis=0),1e-300)
        delta=float(np.sqrt(np.mean(ratio))) if np.isfinite(ratio).all() else None
        gamma=gamma+change;J=self.Jd+B@gamma;Q=self.inc+self.GJd+GB@gamma
        self.counts['current_PR_steps']+=1;self.counts['current_cached_GB_reuses']+=1;self.counts['current_cached_SB_reuses']+=1
        return gamma,grad,direction,J,Q,delta

    def material(self,x,J,Q,metric,beta):
        self.counts['material_value_gradient_calls']+=1;d=np.repeat(x,3)
        b=beta-self.s if metric=='A' else beta;den=1+b*d
        u=beta*(J-d[:,None]*Q)/den[:,None];du=-beta*(Q+b*J)/den[:,None]**2
        if metric=='A':
            Hu=self.Hop@u;R=u+(d/den)[:,None]*Hu
            v=R+self.Hop.conj().T@(np.conj(d/den)[:,None]*R)
            w=np.sum((np.conj(v)*du+np.conj(R)*Hu/den[:,None]**2)/self.E2,axis=1)
            self.counts['material_sparse_H_actions']+=2
        else:R=u;w=np.sum(np.conj(R)*du/self.E2,axis=1)
        pr,pg=self.prior(x);value=self.rho/self.P*float(np.sum(abs(R)**2/self.E2))+pr
        grad=self.rho/self.P*2*np.real(w.reshape(self.N,3).sum(1))+pg
        return value,grad

    def aplus(self,x,J,Q,beta):
        """Eight colors: disjoint H blocks and no same-color prior edge.
        Each batch is the same exact conditional-arc proposal as R8, followed
        by complete unchanged-prior acceptance; no tensor/complex extension.
        """
        x=x.copy();cells,rows=self.block_cells,self.block_rows;b=beta-self.s
        jb,qb=J[rows],Q[rows];v=-beta*(qb+b*jb);z=x/(1+b*x)
        zr=np.repeat(z[cells],3,axis=1);u=beta*jb+zr[:,:,None]*v
        Hu=self.Hblocks@u;R=u+zr[:,:,None]*Hu
        for sweep in range(self.solver['aplus_sweeps']):
            for color in range(8):
                sl=slice(3*color,3*color+3);ids=cells[:,color]
                HiVi=self.Hblocks[:,:,sl]@v[:,sl,:]
                D=zr[:,:,None]*HiVi;D[:,sl,:]=v[:,sl,:]+Hu[:,sl,:]
                C=R-z[ids,None,None]*D
                q=self.rho/self.P*np.sum(abs(D)**2/self.E2,axis=(1,2))
                h=self.rho/self.P*np.sum(np.conj(C)*D/self.E2,axis=(1,2))
                aa=q*b.real+np.real(h*np.conj(b)**2);bb=q+2*np.real(h*np.conj(b));cc=h.real
                disc=bb*bb-4*aa*cc;ok=disc>=-1e-14*np.maximum.reduce((bb*bb,abs(4*aa*cc),np.full_like(bb,1e-300)))
                square=np.sqrt(np.maximum(disc,0));stable=-.5*(bb+np.copysign(square,bb))
                r1=np.divide(stable,aa,out=np.full_like(aa,np.nan),where=abs(aa)>1e-30)
                r2=np.divide(cc,stable,out=np.full_like(cc,np.nan),where=abs(stable)>1e-30)
                linear=np.divide(-cc,bb,out=np.full_like(cc,np.nan),where=abs(bb)>1e-30)
                r1=np.where(abs(aa)<=1e-30,linear,r1);r2=np.where(abs(aa)<=1e-30,np.nan,r2)
                roots=np.c_[r1,r2];legal=ok[:,None]&np.isfinite(roots)&(roots>=self.lo)&(roots<=self.hi)
                candidates=np.c_[x[ids],np.full(len(ids),self.lo),np.full(len(ids),self.hi),np.where(legal,roots,x[ids,None])]
                neighbors=np.r_[x,0.][self.neighbors[ids]]
                a=self.lam*self.degree[ids]/len(self.edges)+self.l2/self.N
                c=-self.lam*np.sum(neighbors,axis=1)/len(self.edges)
                zz=candidates/(1+b*candidates)
                loss=q[:,None]*abs(zz)**2+2*np.real(h[:,None]*zz)+a[:,None]*candidates**2+2*c[:,None]*candidates
                new=candidates[np.arange(len(ids)),np.argmin(loss,axis=1)];newz=new/(1+b*new);dz=newz-z[ids]
                R+=dz[:,None,None]*D;u[:,sl,:]+=dz[:,None,None]*v[:,sl,:];Hu+=dz[:,None,None]*HiVi
                self.counts['aplus_coordinate_changes']+=int(np.sum(new!=x[ids]));x[ids]=new;z[ids]=newz
                zr[:,sl]=newz[:,None];self.counts['aplus_color_batches']+=1
                self.counts['aplus_arc_candidates']+=candidates.size;self.counts['aplus_feasible_arc_roots']+=int(legal.sum())
            self.counts['aplus_sweeps']+=1
        return x

    @staticmethod
    def copy_state(st):
        return {k:(v.copy() if isinstance(v,np.ndarray) else v) for k,v in st.items()
                if k in ('m','J','E','Y','Phi','data_squared','gradient','linear_relative_residual')}

    def run(self,Y,noise_meta,arm,progress=None):
        if arm not in ARMS:raise ValueError('Unknown frozen arm')
        self.counts=Counter();self.state_cache=None;self.Y=np.asarray(Y,complex)
        if self.Y.shape!=(180,self.P) or not np.isfinite(self.Y).all():raise ValueError('Invalid public Y')
        energy=np.sum(abs(self.Y)**2,axis=0)
        if np.any(energy<=0):raise ValueError('Zero public source energy')
        self.Yweights=1/energy;self.best=None;self.budget=[];self.budget_material=[];self.last_budget_signature=None
        self.online_wall,self.online_cpu=time.perf_counter(),time.process_time();self.phase='initial'
        x=np.full(self.N,self.config['physics']['initial_chi']);initial=self.copy_state(self.evaluate(x))
        logs=[];joint_history=[];stages=[];snapshots={};dpcs=None
        if arm!='direct':
            dpcs=self.prepare(noise_meta,arm.endswith('qualified'))
            gamma=None;previous_stage=None;last_J=None
            sourcewise=arm in ('cie_single','dual_source');material_A=arm in ('a_material','a_plus','a_plus_qualified')
            for k,(mf,cycles) in enumerate(zip(self.solver['MF'],self.solver['cycles_per_stage'])):
                stage=self.bases[mf];newgamma=np.zeros((len(stage['keys']),self.P),complex)
                if gamma is not None:
                    oldmap={key:i for i,key in enumerate(previous_stage['keys'])}
                    for j,key in enumerate(stage['keys']):
                        if key in oldmap:newgamma[j]=gamma[oldmap[key]]
                gamma=newgamma;grad=direction=None;transport=0.
                if last_J is not None:transport=float(la.norm(self.Jd+stage['B']@gamma-last_J)/max(la.norm(last_J),1e-30))
                beta_m=self.solver['beta_material'][k]
                beta_c=(self.solver['beta_current_single'] if arm=='cie_single' else self.solver['beta_current_dual'])[k]
                self.phase=f'joint_stage_{k}'
                for it in range(cycles):
                    gamma,grad,direction,J,Q,delta=self.current_step(x,gamma,grad,direction,beta_c,stage,sourcewise)
                    metric='A' if material_A else 'CIE'
                    if arm in ('a_plus','a_plus_qualified'):
                        before=self.material(x,J,Q,metric,beta_m)[0];candidate=self.aplus(x,J,Q,beta_m)
                        after=self.material(candidate,J,Q,metric,beta_m)[0];accepted=bool(after<=before+1e-13)
                        if accepted:x=candidate
                        else:self.counts['aplus_full_loss_rollbacks']+=1
                        info=dict(nit=2,nfev=2,success=accepted,message='two8-colorrealMobiusarcproposals; wholeoriginalquadraticprioraccept',before=before,after=after)
                    else:
                        res=minimize(lambda z:self.material(z,J,Q,metric,beta_m),x,jac=True,method='L-BFGS-B',bounds=[(self.lo,self.hi)]*self.N,
                            options=dict(maxiter=self.solver['material_inner_maxiter'],ftol=self.solver['material_ftol'],gtol=self.solver['material_gtol'],maxls=self.solver['material_maxls']))
                        x=res.x;info=dict(nit=int(res.nit),nfev=int(res.nfev),success=bool(res.success),message=str(res.message))
                    st=self.evaluate(x);joint_history.append(dict(stage=k,iteration=it+1,beta_current=beta_c,beta_material=beta_m,
                        paper_delta_sqrt_mean=delta,physical_Phi=st['Phi'],physical_residual=np.sqrt(st['data_squared']),**info))
                    if progress and (it+1)%20==0:progress(dict(stage=self.phase,iteration=it+1,physical_LU=self.counts['physical_LU'],Phi=self.best['Phi']))
                stages.append(dict(stage=k,m=x.copy(),J=J.copy(),Q=Q.copy(),gamma=gamma.copy(),mf=mf,
                    beta_current=beta_c,beta_material=beta_m,transport_relative_error=transport))
                previous_stage=stage;last_J=J.copy()
        raw=self.copy_state(self.best);self.phase='physical_corrector' if arm!='direct' else 'direct_physical'
        maximum=self.solver['direct_maxiter'] if arm=='direct' else self.solver['corrector_maxiter']
        evaluations=accepted_steps=0
        def objective(z):
            nonlocal evaluations
            value,g=self.evaluate(z,True);evaluations+=1
            if evaluations==self.solver['checkpoint_evaluations']:snapshots['checkpoint80_evaluation']=self.copy_state(self.best)
            return value,g
        def callback(z):
            nonlocal accepted_steps
            accepted_steps+=1
            if not np.array_equal(z,self.state_cache['m']):raise RuntimeError('Accepted point not last paid physical state')
            if accepted_steps==self.solver['checkpoint_accepted']:snapshots['checkpoint80_accepted']=self.copy_state(self.best)
            if progress and accepted_steps%40==0:progress(dict(stage=self.phase,accepted=accepted_steps,physical_LU=self.counts['physical_LU'],Phi=self.best['Phi']))
        res=minimize(objective,raw['m'],jac=True,method='L-BFGS-B',bounds=[(self.lo,self.hi)]*self.N,callback=callback,
            options=dict(maxiter=maximum,ftol=self.solver['physical_ftol'],gtol=self.solver['physical_gtol'],maxls=self.solver['physical_maxls']))
        logs.append(dict(stage=self.phase,nit=int(res.nit),nfev=int(res.nfev),success=bool(res.success),message=str(res.message),allowance=maximum))
        final=self.copy_state(self.best)
        receipt=dict(status='SEALED_PUBLIC_PREDICTION',arm=arm,truth_access='NONE',current_metric='ordinary_CIE_only' if arm!='direct' else 'NONE',
            current_aggregation='sourcewise' if arm in ('cie_single','dual_source') else 'global',
            material_metric='A' if arm in ('a_material','a_plus','a_plus_qualified') else 'CIE' if arm!='direct' else 'original_physics',
            actual_material_dof=self.N,vector_dimension=self.M,observations=int(Y.size),material_box=[self.lo,self.hi],dpcs=dpcs,
            online_CPU_seconds=time.process_time()-self.online_cpu,online_wall_seconds=time.perf_counter()-self.online_wall,counts=dict(self.counts),
            original_Phi=final['Phi'],physical_source_RMS_relative_residual=float(np.sqrt(final['data_squared'])),linear_relative_residual=final['linear_relative_residual'],
            optimizer=logs,joint_history=joint_history,checkpoint80_accepted_reached='checkpoint80_accepted' in snapshots,
            checkpoint80_evaluation_reached='checkpoint80_evaluation' in snapshots,selection='minimum original Phi among every paid host/joint/physical evaluation; no GT',
            Fourier='native3Dprojected-IDFTgamma; exactnestedmasktransport; no coefficientwhitening; denseGpropagation',
            physics='inherited R1/R8 full dyadic CM+RR DDA, densityj=chiE, epsilon=1+chi; not a new backend or original-paper physical reproduction',
            finalization_new_LU=0,finalization_new_RHS=0)
        return dict(initial=initial,raw=raw,final=final,stages=stages,snapshots=snapshots,
                    budget=self.budget,budget_material=np.array(self.budget_material)),receipt


def target_material(case,points,h,contract):
    if case=='hollow_cube_rods':
        c=contract[case];x=box_fraction(points,h,c['center'],[c['outer_side']/2]*3)-box_fraction(points,h,c['center'],[c['inner_side']/2]*3)
        for center in c['rods_centers']:x+=box_fraction(points,h,center,c['rod_halfwidths'])
        return (c['epsilon_r']-1)*x
    if case=='two_cubes_eps9':
        c=contract[case];return (c['epsilon_r']-1)*sum(box_fraction(points,h,center,[c['side']/2]*3) for center in c['centers'])
    raise ValueError('Legacy material must only be copied by generator')
