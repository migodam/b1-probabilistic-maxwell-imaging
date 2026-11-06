"""Real-isotropic 3D fixed-current Green denominator on complex3N currents.

The original CM+RR forward is never changed. Material has N real coordinates;
each coordinate multiplies its three interleaved components for every source.
All readout operations remain on the input device, including the small LU.
"""
from collections import Counter
import torch


class VectorLowRankDenominator:
    def __init__(self,U,GU,projection='galerkin',cost=None):
        if U.ndim!=2 or U.shape!=GU.shape or U.shape[0]%3:
            raise ValueError('Matching complex3N by rank U/GU required')
        self.U=U;self.B=U@(U.mH@GU) if projection=='galerkin' else GU
        if projection not in ('galerkin','right'):raise ValueError(projection)
        self.M,self.rank=U.shape;self.N=self.M//3;self.projection=projection
        self.cost=Counter() if cost is None else cost

    def action(self,X):
        self.cost['denominator_Green_actions']+=1
        self.cost['denominator_Green_complex_RHS']+=X.shape[1]
        return self.B@(self.U.mH@X)

    def factor(self,m,beta):
        if m.shape!=(self.N,):raise ValueError('N real scalar materials required')
        d=torch.repeat_interleave(m,3);a=1+beta*d;F=(d/a)[:,None]*self.B
        T=torch.eye(self.rank,device=m.device,dtype=self.U.dtype)-self.U.mH@F
        LU,pivots,info=torch.linalg.lu_factor_ex(T,check_errors=False)
        self.cost['denominator_small_LU']+=1
        self.cost['denominator_small_LU_dimension_sum']+=self.rank
        return a,F,LU,pivots,info

    def solve(self,factor,rhs,adjoint=False):
        a,F,LU,pivots,info=factor
        self.cost['denominator_small_triangular_complex_RHS']+=rhs.shape[1]
        if adjoint:
            v=torch.linalg.lu_solve(LU,pivots,F.mH@rhs,adjoint=True)
            return (rhs+self.U@v)/a[:,None]
        b=rhs/a[:,None]
        return b+F@torch.linalg.lu_solve(LU,pivots,self.U.mH@b)

    @staticmethod
    def weighted(X,omega):
        return X*omega if omega.ndim==1 else X@omega

    def value_gradient(self,J,Q,m,beta,omega,rho=1.,reference_factor=None,gradient=True):
        factor=self.factor(m,beta) if reference_factor is None else reference_factor
        self.latest_factor=factor
        w=self.solve(factor,J-torch.repeat_interleave(m,3)[:,None]*Q)
        weighted=self.weighted(w,omega)
        value=.5*rho*beta**2*torch.sum((w.conj()*weighted).real)
        self.cost['material_inner_evaluations']+=1
        if gradient:
            a=self.solve(factor,weighted,adjoint=True)
            h=Q+beta*w-self.action(w) if reference_factor is None else Q
            g=-rho*beta**2*torch.sum((a.conj()*h).real.reshape(self.N,3,-1),dim=(1,2))
        else:g=None
        return value,g,w,factor[-1]

    def readout(self,J,Q,m,beta,omega,mode,upper=8.,steps=8,rho=1.):
        """R0 projected descent, eight paid steps; old/ordinary seed accepted safely.

        Lagged uses A(m_old) throughout this readout only. Its exact line search
        is along the feasible projected direction with alpha constrained [0,1].
        Full Armijo makes at most eight trial evaluations per paid step.
        """
        if mode=='contractedLagged':return self.contracted_readout(J,Q,m,beta,omega,upper,steps,rho)
        if mode not in ('full','lagged'):raise ValueError(mode)
        if steps!=8:raise ValueError('Frozen eight-step material-readout budget')
        before=Counter(self.cost);reference=self.factor(m,beta) if mode=='lagged' else None
        n=self.N;V=(Q+beta*J).reshape(n,3,-1);T=(beta*J).reshape(n,3,-1)
        if omega.ndim!=1:raise ValueError('Readout ordinary seed currently expects diagonal source metric')
        den=torch.sum(abs(V)**2*omega,dim=(1,2));num=torch.sum((V.conj()*T).real*omega,dim=(1,2))
        r=torch.where(den>0,num/den,torch.zeros_like(num)).clamp(0,beta*upper/(1+beta*upper))
        seed=(r/(beta*(1-r))).clamp(0,upper)
        old_v,old_g,old_w,info=self.value_gradient(J,Q,m,beta,omega,rho,reference)
        old_factor=self.latest_factor
        seed_v,seed_g,seed_w,seed_info=self.value_gradient(J,Q,seed,beta,omega,rho,reference)
        seed_factor=self.latest_factor
        use_seed=torch.isfinite(seed_v)&(seed_v<old_v)&(seed_info==0)
        x=torch.where(use_seed,seed,m).clone();value=torch.where(use_seed,seed_v,old_v)
        w=torch.where(use_seed,seed_w,old_w)
        current_factor=tuple(torch.where(use_seed,s,o) for s,o in zip(seed_factor,old_factor))
        g=torch.where(use_seed,seed_g,old_g);records=[];lu_infos=[info,seed_info];trial_counts=[]
        for step in range(steps):
            # Positive scalar scaling only selects a feasible direction; the
            # fixed-A line minimizer below is recomputed after projection.
            scale=torch.clamp(g.abs().max(),min=torch.finfo(g.dtype).tiny)
            direction=torch.clamp(x-g/scale,0,upper)-x
            slope=torch.dot(g,direction)
            if mode=='lagged':
                dw=-self.solve(reference,torch.repeat_interleave(direction,3)[:,None]*Q)
                curvature=rho*beta**2*torch.sum((dw.conj()*self.weighted(dw,omega)).real)
                alpha=torch.where(curvature>0,-slope/curvature,torch.zeros_like(slope)).clamp(0,1)
                candidate=x+alpha*direction
                trial,tg,tw,ti=self.value_gradient(J,Q,candidate,beta,omega,rho,reference)
                ok=torch.isfinite(trial)&(trial<=value)&(ti==0)
                x=torch.where(ok,candidate,x);g=torch.where(ok,tg,g);value=torch.where(ok,trial,value);w=torch.where(ok,tw,w)
                lu_infos.append(ti);trials=1
            else:
                h=Q+beta*w-self.action(w)
                dw=-self.solve(current_factor,torch.repeat_interleave(direction,3)[:,None]*h)
                curvature=rho*beta**2*torch.sum((dw.conj()*self.weighted(dw,omega)).real)
                alpha=torch.where(curvature>0,-slope/curvature,torch.zeros_like(slope)).clamp(0,1)
                self.cost['full_direction_GN_initializations']+=1
                # Only one acceptance scalar crosses to the host per trial.
                # U/B, material, currents, gradients and all LU stay resident.
                for trial_index in range(8):
                    candidate=x+alpha*direction
                    trial,tg,tw,ti=self.value_gradient(J,Q,candidate,beta,omega,rho)
                    ok=torch.isfinite(trial)&(ti==0)&(trial<=value+1e-4*alpha*slope)
                    lu_infos.append(ti)
                    if bool(ok.item()):
                        x=candidate;g=tg;value=trial;w=tw;current_factor=self.latest_factor;break
                    alpha=alpha*.5
                trials=trial_index+1
            records.append(torch.stack((value,alpha,slope,g.abs().max())))
            trial_counts.append(trials)
        self.cost['matrix_material_readouts']+=1;self.cost['matrix_material_projected_steps']+=steps
        packed=torch.stack([old_v,seed_v,value,use_seed.to(x.dtype),torch.stack(lu_infos).abs().max().to(x.dtype)])
        return x.clamp(0,upper),dict(mode=mode,projection=self.projection,rank=self.rank,steps=steps,
            packed_metrics=packed,device_step_trace=torch.stack(records),
            per_step_trials=trial_counts,source_metric='diagonal all-source supplied omega',prior='R0',
            counters={k:int(v-before.get(k,0)) for k,v in self.cost.items() if v-before.get(k,0)})

    def factor_contracted(self,m_reference,beta):
        factor=self.factor(m_reference,beta)
        # A=(I+beta Dm) T, with T=I-Dz K. F=z B and the small
        # Woodbury LU are identical, but a T-solve has no diagonal divide.
        return (torch.ones_like(factor[0]),)+factor[1:]

    def contracted_value_gradient(self,J,V,z,beta,omega,rho,reference_factor,gradient=True):
        w=self.solve(reference_factor,J-torch.repeat_interleave(z,3)[:,None]*V)
        weighted=self.weighted(w,omega)
        value=.5*rho*beta**2*torch.sum((w.conj()*weighted).real)
        if gradient:
            b=self.solve(reference_factor,weighted,adjoint=True)
            g=-rho*beta**2*torch.sum((b.conj()*V).real.reshape(self.N,3,-1),dim=(1,2))
        else:g=None
        self.cost['material_inner_evaluations']+=1
        return value,g,w,reference_factor[-1]

    def contracted_readout(self,J,Q,m,beta,omega,upper=8.,steps=8,rho=1.):
        """Root-contracted z-coordinate fixed-T quadratic, R0; not an MM claim."""
        if steps!=8 or omega.ndim!=1:raise ValueError('Frozen eight steps and diagonal source metric')
        before=Counter(self.cost);reference=self.factor_contracted(m,beta)
        V=Q+beta*J;vv=V.reshape(self.N,3,-1);jj=J.reshape(self.N,3,-1)
        den=torch.sum(abs(vv)**2*omega,dim=(1,2));num=torch.sum((vv.conj()*jj).real*omega,dim=(1,2))
        zmax=upper/(1+beta*upper);seed=torch.where(den>0,num/den,torch.zeros_like(num)).clamp(0,zmax)
        old=m/(1+beta*m)
        old_v,old_g,_,info=self.contracted_value_gradient(J,V,old,beta,omega,rho,reference)
        seed_v,seed_g,_,_=self.contracted_value_gradient(J,V,seed,beta,omega,rho,reference)
        use_seed=torch.isfinite(seed_v)&(seed_v<old_v)&(info==0)
        z=torch.where(use_seed,seed,old);value=torch.where(use_seed,seed_v,old_v);g=torch.where(use_seed,seed_g,old_g)
        records=[]
        for step in range(steps):
            scale=torch.clamp(g.abs().max(),min=torch.finfo(g.dtype).tiny)
            direction=torch.clamp(z-g/scale,0,zmax)-z;slope=torch.dot(g,direction)
            dw=-self.solve(reference,torch.repeat_interleave(direction,3)[:,None]*V)
            curvature=rho*beta**2*torch.sum((dw.conj()*self.weighted(dw,omega)).real)
            alpha=torch.where(curvature>0,-slope/curvature,torch.zeros_like(slope)).clamp(0,1)
            candidate=z+alpha*direction
            trial,tg,_,_=self.contracted_value_gradient(J,V,candidate,beta,omega,rho,reference)
            ok=torch.isfinite(trial)&(trial<=value)&(info==0)
            z=torch.where(ok,candidate,z);g=torch.where(ok,tg,g);value=torch.where(ok,trial,value)
            records.append(torch.stack((value,alpha,slope,g.abs().max())))
        self.cost['matrix_material_readouts']+=1;self.cost['matrix_material_projected_steps']+=steps
        self.cost['contracted_lagged_material_readouts']+=1
        # Guard only the last ulp of the exact z-box inverse at m=upper.
        result=(z/(1-beta*z)).clamp(0,upper)
        return result,dict(mode='contractedLagged',projection=self.projection,rank=self.rank,steps=steps,
            packed_metrics=torch.stack((old_v,seed_v,value,use_seed.to(m.dtype),info.abs().to(m.dtype))),
            device_step_trace=torch.stack(records),per_step_trials=[1]*steps,
            source_metric='diagonal all-source supplied omega',prior='R0',
            coordinate='z=m/(1+beta*m), fixed T=I-Dz_ref K, exact feasible-direction line search',
            counters={k:int(v-before.get(k,0)) for k,v in self.cost.items() if v-before.get(k,0)})
