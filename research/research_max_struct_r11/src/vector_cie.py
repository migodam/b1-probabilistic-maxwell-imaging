"""R11 material-readout-only extension of R10 true3D ordinary GlobalDual.
The inherited global PR current step and Fourier transport are copied without
mathematical change. The original physical forward is imported unchanged.
"""
import time
import numpy as np
import torch
from research.research_max_struct_r10.src.vector_basis import ImplicitFourier
from research.research_max_struct_r10.src.vector_cie import cuda_material_readout

def auxiliary_solve(model,receiver,Y,arm,config,callback=None,endpoint=None,material_operator=None):
    if arm!='GlobalDual':raise ValueError('R11 readout ablation uses matched ordinary GlobalDual current')
    wd=receiver.wd;ws=receiver.ws
    m=torch.full((model.N,),config['initial_chi'],device=model.device,dtype=torch.float64);J=receiver.Jd.clone();Q=receiver.Qd.clone();D=receiver.Dd.clone()
    mode=config['material_mode']
    if mode not in ('ordinary','full','lagged','contractedLagged'):raise ValueError(mode)
    if mode!='ordinary' and material_operator is None:raise ValueError('New readout requires explicit public operator')
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
            if mode=='ordinary':
                m=cuda_material_readout(J,Q,ws,bm,config['material_upper']);material_info={}
            else:
                m,material_info=material_operator.readout(J,Q,m,bm,ws,mode,upper=config['material_upper'],rho=1/model.P)
                packed=material_info.pop('packed_metrics').cpu().tolist()
                material_info['steps_trace']=material_info.pop('device_step_trace').cpu().tolist()
                material_info.update(old_own_loss=packed[0],ordinary_seed_own_loss=packed[1],final_own_loss=packed[2],seed_selected=bool(packed[3]),max_small_LU_info=packed[4])
                if packed[4]!=0:raise RuntimeError('Readout small LU failed; preserved trace records')
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
