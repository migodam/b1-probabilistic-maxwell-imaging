"""R10 parameterized R9 real-isotropic Fourier-CIE readouts.
Known-algorithm adaptation: upper8 throughout, m0.25, native MF2/3/4,
min20/max500 paper coefficient stopping, checkpointed cooperative GPU gate.
GlobalDual and AplusGlobal have identical ordinary CIE current updates.
SourceCIE differs in current beta AND sourcewise PR; not a matched-A ablation.
No NIE, tensor/complex-material claim, or new Green backend.
"""
import time
import numpy as np
import torch
from research.research_max_struct_r10.src.vector_basis import ImplicitFourier

def cuda_material_readout(J,Q,weights,beta,upper=8.):
    n=J.shape[0]//3;V=(Q+beta*J).reshape(n,3,-1);T=(beta*J).reshape(n,3,-1)
    den=torch.sum(abs(V)**2*weights,dim=(1,2));num=torch.sum(torch.real(V.conj()*T)*weights,dim=(1,2))
    r=torch.where(den>0,num/den,torch.zeros_like(num));r=torch.clamp(r,0,beta*upper/(1+beta*upper))
    return torch.clamp(r/(beta*(1-r)),0,upper)

class LocalArcReadout:
    def __init__(self,model,upper=8.):
        self.upper=float(upper)
        from research.research_max_struct_r9.src.vector_original_extension import LocalArcReadout
        local=LocalArcReadout(model.cpu);self.ids=torch.as_tensor(local.ids,device=model.device);self.H=model.array(local.H);self.s=local.self_term
    def __call__(self,J,Q,old,weights,beta):
        ids=self.ids;nb=len(ids);b=beta-self.s
        jb=J.reshape(-1,3,J.shape[1])[ids].reshape(nb,24,-1);qb=Q.reshape(-1,3,Q.shape[1])[ids].reshape(nb,24,-1)
        h=qb+b*jb
        def state(m):
            z=m/(1+b*m);zr=torch.repeat_interleave(z,3,dim=1)
            v=jb-zr[:,:,None]*h;w=torch.einsum('ij,bjs->bis',self.H,v)
            return z,v,w,beta*(v+zr[:,:,None]*w)
        def loss(r):return torch.sum(abs(r)**2*weights,dim=(1,2))
        oldb=old[ids];m=oldb.clone();z,v,w,r=state(m);initial=loss(r)
        seed=cuda_material_readout(J,Q,weights,beta,self.upper)[ids];zs,vs,ws,rs=state(seed);use=loss(rs)<initial
        m=torch.where(use[:,None],seed,m);z=torch.where(use[:,None],zs,z)
        v=torch.where(use[:,None,None],vs,v);w=torch.where(use[:,None,None],ws,w);r=torch.where(use[:,None,None],rs,r)
        for sweep in range(2):
            for slot in range(8):
                sl=slice(3*slot,3*slot+3);zi=z[:,slot].clone();hi=h[:,sl]
                B=-beta*torch.repeat_interleave(z,3,dim=1)[:,:,None]*torch.einsum('ij,bjs->bis',self.H[:,sl],hi)
                B[:,sl]=beta*(-hi+w[:,sl]);C=r-zi[:,None,None]*B
                q=torch.sum(abs(B)**2*weights,dim=(1,2));hh=torch.sum(C.conj()*B*weights,dim=(1,2))
                aa=q*b.real+torch.real(hh*b.conjugate()**2);bb=q+2*torch.real(hh*b.conjugate());cc=hh.real
                disc=bb*bb-4*aa*cc;root=torch.sqrt(torch.clamp(disc,min=0))
                stable=-.5*(bb+torch.where(bb>=0,1.,-1.)*root)
                nan=torch.full_like(aa,float('nan'))
                r1=torch.where(aa!=0,stable/aa,nan);r2=torch.where(stable!=0,cc/stable,nan)
                linear=torch.where(bb!=0,-cc/bb,nan);r1=torch.where(aa==0,linear,r1)
                candidates=torch.stack([m[:,slot],torch.zeros_like(aa),torch.full_like(aa,self.upper),r1,r2],dim=1)
                valid=torch.isfinite(candidates)&(candidates>=0)&(candidates<=self.upper)
                valid[:,3:]&=(disc[:,None]>=0)|(aa[:,None]==0)
                zc=candidates/(1+b*candidates);values=q[:,None]*abs(zc)**2+2*torch.real(hh[:,None]*zc)
                values=torch.where(valid,values,torch.full_like(values,float('inf')))
                best=torch.argmin(values,dim=1);mn=torch.gather(candidates,1,best[:,None])[:,0];zn=mn/(1+b*mn);dz=zn-zi
                r=r+dz[:,None,None]*B;dv=-dz[:,None,None]*hi
                w=w+torch.einsum('ij,bjs->bis',self.H[:,sl],dv);v[:,sl]=v[:,sl]+dv;m[:,slot]=mn;z[:,slot]=zn
        _,_,_,final=state(m);final_loss=loss(final)
        ok=torch.isfinite(final_loss)&(final_loss<=initial*(1+1e-12)+1e-25);m=torch.where(ok[:,None],m,oldb)
        out=old.clone();out[ids.flatten()]=m.flatten()
        return out,dict(blocks=nb,rejected_blocks=int((~ok).sum().item()),coordinate_sweeps=2,
            fixed_current_loss_before=float(initial.sum().item()),fixed_current_loss_after=float(torch.where(ok,final_loss,initial).sum().item()),
            variant='R9vectorAplus_with_full_A_loss_accepted_ordinary_CIE_seed_proposal',seed_proposals=nb,
            seed_accepted_blocks=int(use.sum().item()),seed_and_final_full_block_evaluations=3)

def auxiliary_solve(model,receiver,Y,arm,config,callback=None,endpoint=None):
    if arm not in ('sourceCIE','sourceDual','GlobalDual','AplusGlobal'):raise ValueError(arm)
    wd=receiver.wd;ws=receiver.ws
    m=torch.full((model.N,),config['initial_chi'],device=model.device,dtype=torch.float64);J=receiver.Jd.clone();Q=receiver.Qd.clone();D=receiver.Dd.clone()
    arc=LocalArcReadout(model,config['material_upper']) if arm=='AplusGlobal' else None
    history=[];stages=[];old_basis=None;gamma=None
    if endpoint:endpoint('initial_auxiliary',m,dict(J=J,Q=Q,gamma=None,stage=None))
    for stage,mf in enumerate(config['mf']):
        bm=config['beta_material'][stage];bj=(config['beta_material'] if arm=='sourceCIE' else config['beta_current_dual'])[stage]
        basis=ImplicitFourier(model,receiver,mf)
        gamma=torch.zeros((basis.K,60),device=model.device,dtype=torch.complex128) if gamma is None else basis.transport(old_basis,gamma)
        previous_gradient=None;previous_direction=None;stage_start=time.perf_counter()
        for iteration in range(config['max_iterations_per_stage']):
            if iteration%20==0:model.gate.wait('PR_'+arm+'_stage'+str(stage)+'_iteration'+str(iteration),torch.cuda.memory_reserved() if model.device.type=='cuda' else 0)
            wall,cput=time.perf_counter(),time.process_time();old_gamma=gamma
            chi=torch.repeat_interleave(m,3)[:,None];p=bj/(1+bj*chi)
            residual=p*(J-chi*Q);v=p.conj()*residual*ws
            gradient=basis.adjoint(model.S.conj().T@(D*wd)+v-model.green(chi*v,True))
            dim=0 if arm in ('sourceCIE','sourceDual') else None
            if previous_gradient is None:direction=-gradient
            else:
                den=torch.sum(abs(previous_gradient)**2,dim=dim)
                numerator=torch.real(torch.sum(gradient.conj()*(gradient-previous_gradient),dim=dim))
                coefficient=torch.clamp(torch.where(den>0,numerator/den,torch.zeros_like(den)),min=0)
                direction=-gradient+previous_direction*coefficient
                descent=torch.real(torch.sum(gradient.conj()*direction,dim=dim))
                direction=torch.where((descent>=0)[None,:],-gradient,direction) if dim==0 else torch.where(descent>=0,-gradient,direction)
            u=basis.apply(direction);gu=model.green(u);su=model.S@u;hu=p*(u-chi*gu)
            numerator=-torch.real(torch.sum(gradient.conj()*direction,dim=dim))
            denominator=torch.sum(abs(su)**2*wd,dim=dim)+torch.sum(abs(hu)**2*ws,dim=dim)
            alpha=torch.where(denominator>0,numerator/denominator,torch.zeros_like(numerator))
            gamma=gamma+direction*alpha;J=J+u*alpha;Q=Q+gu*alpha;D=D+su*alpha
            previous_gradient=gradient;previous_direction=direction
            if arc:m,material_info=arc(J,Q,m,ws,bm)
            else:m=cuda_material_readout(J,Q,ws,bm,config['material_upper']);material_info={}
            den=torch.sum(abs(old_gamma)**2,dim=0);num=torch.sum(abs(gamma-old_gamma)**2,dim=0)
            valid=bool(torch.all((den>0)|(num==0)).item())
            delta_sq=float(torch.mean(torch.where(den>0,num/den,torch.zeros_like(num))).item()) if valid else None
            delta=float(np.sqrt(delta_sq)) if valid else None
            chi=torch.repeat_interleave(m,3)[:,None];p=bj/(1+bj*chi)
            data=float(torch.mean(torch.sum(abs(D)**2,dim=0)*wd).item())
            state=float(torch.mean(torch.sum(abs(p*(J-chi*Q))**2,dim=0)*ws).item())
            if model.device.type=='cuda':torch.cuda.synchronize()
            record=dict(stage=stage,mf=mf,iteration=iteration+1,beta_material=bm,beta_current=bj,beta_active=True,
                current_representation='ordinary_real_isotropic_CIE_on_full3N_current',
                coefficient_delta_RMS=delta,coefficient_delta_squared_mean=delta_sq,auxiliary_data=data,auxiliary_state=state,
                auxiliary_sum=data+state,peak_contrast=float(m.max().item()),material_readout=material_info,
                CPU_seconds=time.process_time()-cput,wall_seconds=time.perf_counter()-wall)
            history.append(record)
            if callback:callback(record,m,J,Q,gamma,basis)
            if iteration+1>=config['minimum_iterations_per_stage'] and delta is not None and delta<config['coefficient_delta_RMS_tolerance']:stop='PAPER_COEFFICIENT_STOP';break
        else:stop='CAP'
        stages.append(dict(stage=stage,mf=mf,iterations=iteration+1,stop=stop,wall_seconds=time.perf_counter()-stage_start))
        if endpoint:endpoint('stage'+str(stage+1),m,dict(stages[-1],J=J,Q=Q,gamma=gamma,basis=basis))
        old_basis=basis
    return dict(material=m,history=history,stages=stages,receiver=receiver.info)
