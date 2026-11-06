"""Resident complex128 Torch/cuFFT wrapper of the existing CM+RR DDA arrays.

Convolution layout follows A18's cached Torch FFT primitive. Krylov vectors stay
resident; restarted left-preconditioned GMRES uses CGS2 and small Torch QR.
This is a numerical-backend adaptation, not a new scientific inverse algorithm.
"""
from collections import Counter
from pathlib import Path
import time
import numpy as np
import torch


class CudaOriginalModel:
    def __init__(self,cpu,*,fft_layout='strided',block_rhs=None):
        if not torch.cuda.is_available():raise RuntimeError('CUDA required')
        wall,cput=time.perf_counter(),time.process_time()
        torch.backends.cuda.matmul.allow_tf32=False
        self.cpu=cpu;self.device=torch.device('cuda');self.N=cpu.N;self.M=cpu.M;self.P=60
        self.side=cpu.side;self.full_N=cpu.full_N;self.s=cpu.s;self.rtol=cpu.rtol
        self.restart=cpu.restart;self.cycles=cpu.cycles;self.block_rhs=cpu.block_rhs;self.fft_layout='strided'
        self.shape=cpu.convolution_shape;self.pad=cpu.pad;self.convolution_N=cpu.convolution_N
        self.rows=torch.as_tensor(cpu.green_rows,device='cuda',dtype=torch.long)
        self.kernel=torch.as_tensor(cpu.kernel,device='cuda',dtype=torch.complex128)
        self.inc=torch.as_tensor(cpu.inc,device='cuda',dtype=torch.complex128)
        self.S=torch.as_tensor(cpu.S,device='cuda',dtype=torch.complex128)
        self.counts=Counter();self.solve_records=[];self.max_residual=0.;self._physical_cache=None
        self.configure_fft(fft_layout,block_rhs or cpu.block_rhs)
        torch.cuda.synchronize()
        self.setup=dict(cpu.setup,CUDA_upload_CPU_seconds=time.process_time()-cput,
            CUDA_upload_wall_seconds=time.perf_counter()-wall,backend='Torchcomplex128residentcuFFT_CGS2_GMRES60',
            GMRES_check_interval=8,GMRES_original_residual_tolerance=self.rtol,
            backend_change='CGS2 instead of SciPy MGS; QR/check every8 steps and cycle end; up to7 extra Arnoldi steps',
            device=torch.cuda.get_device_name(0),torch_version=torch.__version__)
        self.setup.update(FFT_layout=self.fft_layout,FFT_block_rhs=self.block_rhs,external_assembly_chunk=cpu.block_rhs)

    def configure_fft(self,layout,block_rhs):
        if layout not in ('strided','packed') or block_rhs not in (1,2,4,16):raise ValueError('Frozen engineering candidate only')
        if layout!=self.fft_layout:
            order=(3,4,0,1,2) if layout=='packed' else (2,3,4,0,1)
            self.kernel=self.kernel.permute(order).contiguous();self.fft_layout=layout
        self.block_rhs=block_rhs
        if hasattr(self,'setup'):self.setup.update(FFT_layout=layout,FFT_block_rhs=block_rhs)

    def array(self,x):return torch.as_tensor(x,device='cuda',dtype=torch.complex128)

    def green(self,X,adjoint=False):
        one=X.ndim==1;X=X[:,None] if one else X;out=torch.empty_like(X)
        kg=self.kernel.transpose(0,1).conj() if adjoint and self.fft_layout=='packed' else (self.kernel.transpose(-1,-2).conj() if adjoint else self.kernel)
        for i in range(0,X.shape[1],self.block_rhs):
            j=min(i+self.block_rhs,X.shape[1]);b=j-i
            full=torch.zeros((3*self.convolution_N,b),device='cuda',dtype=torch.complex128)
            full[self.rows]=X[:,i:j];cube=full.reshape(self.shape+(3,b))
            if self.fft_layout=='packed':
                cube=cube.permute(4,3,0,1,2).contiguous()
                hat=torch.fft.fftn(cube,s=self.pad,dim=(-3,-2,-1))
                prod=torch.einsum('ijxyz,bjxyz->bixyz',kg,hat)
                cropped=torch.fft.ifftn(prod,dim=(-3,-2,-1))[..., :self.shape[0],:self.shape[1],:self.shape[2]]
                value=cropped.permute(2,3,4,1,0).reshape(3*self.convolution_N,b)
            else:
                hat=torch.fft.fftn(cube,s=self.pad,dim=(0,1,2))
                prod=torch.einsum('...ij,...jb->...ib',kg,hat)
                value=torch.fft.ifftn(prod,dim=(0,1,2))[tuple(slice(0,d) for d in self.shape)].reshape(3*self.convolution_N,b)
            out[:,i:j]=value[self.rows]
            self.counts['Green_adjoint_RHS' if adjoint else 'Green_forward_RHS']+=b
            self.counts['Green_block_calls']+=1;self.counts['scalar_FFT3_transforms']+=6*b
        return out[:,0] if one else out

    def solve_one(self,chi,b,adjoint=False):
        """Authoritative stopping uses original unpreconditioned residual."""
        d=torch.repeat_interleave(chi,3);diag=1-d*self.s
        if adjoint:diag=diag.conj()
        def action(v):return v-self.green(d*v,True) if adjoint else v-d*self.green(v)
        bnorm=torch.linalg.vector_norm(b);bn=float(bnorm.item())
        if bn==0:return torch.zeros_like(b),dict(iterations=0,cycles=0,relative_residual=0.,info=0,projected_checks=0,original_residual_checks=0)
        maxdiag=torch.max(torch.abs(diag));x=torch.zeros_like(b);r=b.clone();iterations=0;checks=0;actual_checks=0
        last_rr=None
        for cycle in range(self.cycles):
            z=r/diag;beta=torch.linalg.vector_norm(z)
            V=torch.zeros((self.M,self.restart+1),device='cuda',dtype=torch.complex128)
            H=torch.zeros((self.restart+1,self.restart),device='cuda',dtype=torch.complex128)
            V[:,0]=z/beta
            for j in range(self.restart):
                w=action(V[:,j])/diag;vv=V[:,:j+1]
                h=vv.conj().T@w;w=w-vv@h
                correction=vv.conj().T@w;w=w-vv@correction
                H[:j+1,j]=h+correction
                hn=torch.linalg.vector_norm(w);H[j+1,j]=hn
                V[:,j+1]=w/torch.clamp(hn,min=1e-300)
                iterations+=1;self.counts['gmres_Arnoldi_steps']+=1
                if (j+1)%8 and j+1!=self.restart:continue
                checks+=1;hh=H[:j+2,:j+1]
                q,rr=torch.linalg.qr(hh,mode='reduced')
                y=torch.linalg.solve_triangular(rr,(q[0,:].conj()*beta)[:,None],upper=True)[:,0]
                e=torch.zeros(j+2,device='cuda',dtype=torch.complex128);e[0]=beta
                predicted=torch.linalg.vector_norm(e-hh@y)
                conservative=float((maxdiag*predicted/bnorm).item())
                self.counts['small_QR_solves']+=1;self.counts['GMRES_scalar_convergence_checks']+=1
                if not np.isfinite(conservative):raise RuntimeError('Nonfinite resident GMRES projected residual')
                if conservative<=self.rtol or j+1==self.restart:
                    candidate=x+V[:,:j+1]@y;true_r=b-action(candidate)
                    last_rr=float((torch.linalg.vector_norm(true_r)/bnorm).item());actual_checks+=1
                    self.counts['GMRES_original_residual_checks']+=1
                    if last_rr<=self.rtol:
                        return candidate,dict(iterations=iterations,cycles=cycle+1,relative_residual=last_rr,info=0,
                            projected_checks=checks,original_residual_checks=actual_checks)
                    if j+1==self.restart:x=candidate;r=true_r
        return x,dict(iterations=iterations,cycles=self.cycles,relative_residual=last_rr,info=self.cycles,
            projected_checks=checks,original_residual_checks=actual_checks)

    def solve(self,chi,rhs,*,adjoint=False,label='',progress=None,checkpoint_dir=None):
        chi=torch.as_tensor(chi,device='cuda',dtype=torch.float64)
        rhs=self.array(rhs) if not isinstance(rhs,torch.Tensor) else rhs
        out=torch.empty_like(rhs);checkpoint_dir=Path(checkpoint_dir) if checkpoint_dir is not None else None
        if checkpoint_dir is not None:
            if not self.cpu.generation_only or adjoint:raise ValueError('Generation-only checkpoints')
            checkpoint_dir.mkdir(parents=True,exist_ok=True);material_file=checkpoint_dir/'material.npz'
            chi_host=chi.cpu().numpy()
            if material_file.exists():
                saved=np.load(material_file)
                if not np.allclose(saved['chi'],chi_host,rtol=0,atol=1e-14) or not np.array_equal(saved['active'],self.cpu.active):raise ValueError('Checkpoint material mismatch')
            else:np.savez_compressed(material_file,chi=chi_host,active=self.cpu.active)
        for source in range(rhs.shape[1]):
            file=checkpoint_dir/('rhs%03d.npz'%source) if checkpoint_dir else None
            if file is not None and file.exists():
                saved=np.load(file)
                rhs_host=rhs[:,source].cpu().numpy()
                rhs_difference=float(np.linalg.norm(saved['rhs']-rhs_host)/max(np.linalg.norm(rhs_host),1e-300))
                if rhs_difference>1e-12:raise ValueError('Checkpoint RHS differs beyond cross-platform roundoff')
                candidate=self.array(saved['solution']);d=torch.repeat_interleave(chi,3)
                residual=float((torch.linalg.vector_norm(candidate-d*self.green(candidate)-rhs[:,source])/torch.linalg.vector_norm(rhs[:,source])).item())
                self.counts['resumed_RHS_original_residual_verifications']+=1
                if residual<=self.rtol:
                    out[:,source]=candidate;self.counts['resumed_qualified_primal_RHS']+=1
                    self.max_residual=max(self.max_residual,residual)
                    if progress:progress(dict(label=label,source=source,status='RESUMED_QUALIFIED_RHS',relative_RHS_platform_difference=rhs_difference,
                        original_relative_residual=residual,additional_Arnoldi_steps=0))
                    continue
                self.counts['checkpoint_RHS_requiring_repair_solve']+=1
            wall,cput=time.perf_counter(),time.process_time();before=self.counts['gmres_Arnoldi_steps']
            try:out[:,source],record=self.solve_one(chi,rhs[:,source],adjoint)
            except BaseException as exc:
                self.counts['interrupted_or_failed_RHS']+=1
                record=dict(iterations=self.counts['gmres_Arnoldi_steps']-before,status='INTERRUPTED_OR_FAILED',error_type=type(exc).__name__)
                record.update(label=label,source=source,adjoint=adjoint,CPU_seconds=time.process_time()-cput,wall_seconds=time.perf_counter()-wall)
                self.solve_records.append(record)
                if progress:progress(record)
                raise
            torch.cuda.synchronize();record.update(label=label,source=source,adjoint=adjoint,
                CPU_seconds=time.process_time()-cput,wall_seconds=time.perf_counter()-wall)
            self.solve_records.append(record);self.counts['adjoint_RHS' if adjoint else 'primal_RHS']+=1
            self.max_residual=max(self.max_residual,record['relative_residual'])
            if progress:progress(record)
            if record['info'] or record['relative_residual']>self.rtol:
                self.counts['failed_qualified_RHS']+=1;raise RuntimeError('QUALIFIED_RESIDENT_GMRES_FAILED '+str(record))
            if file is not None:
                solution=out[:,source].cpu().numpy();rhs_host=rhs[:,source].cpu().numpy()
                prediction=(self.S@out[:,source]).cpu().numpy()
                np.savez_compressed(file,rhs=rhs_host,solution=solution,Y_clean_column=prediction,
                    relative_residual=record['relative_residual'],iterations=record['iterations'],rtol=self.rtol,
                    backend='TorchCUDAcomplex128CGS2')
                self.counts['qualified_source_checkpoints_written']+=1
        return out

    def physical(self,chi,Y=None,*,label='',progress=None,checkpoint_dir=None):
        chi=torch.as_tensor(chi,device='cuda',dtype=torch.float64)
        self.counts['physical_evaluations']+=1
        J=self.solve(chi,torch.repeat_interleave(chi,3)[:,None]*self.inc,label=label,progress=progress,checkpoint_dir=checkpoint_dir)
        prediction=self.S@J
        if Y is None:return prediction.cpu().numpy()
        Y=self.array(Y) if not isinstance(Y,torch.Tensor) else Y
        loss=float(torch.mean(torch.sum(abs(prediction-Y)**2,dim=0)/torch.sum(abs(Y)**2,dim=0)).item())
        return dict(prediction=prediction.cpu().numpy(),data_loss_mean=loss,relative_source_RMS=float(np.sqrt(loss)))


class CudaReceiverSubspace:
    def __init__(self,model,Y,rtol=.1):
        wall,cput=time.perf_counter(),time.process_time();S=model.S
        ev,U=torch.linalg.eigh(S@S.conj().T);ev=torch.clamp(ev.flip((0,)),min=0);U=U.flip((1,))
        singular=torch.sqrt(ev);keep=singular>=rtol*singular[0];self.rank=int(keep.sum().item())
        U=U[:,keep];sig=singular[keep];self.V=(S.conj().T@U)/sig[None,:]
        self.Jd=self.V@((U.conj().T@Y)/sig[:,None])
        self.Qd=model.inc+model.green(self.Jd);self.Dd=S@self.Jd-Y;torch.cuda.synchronize()
        self.info=dict(rank=self.rank,relative_cutoff=rtol,singular_values=singular.cpu().tolist(),
            right_vector_orthogonality=float(torch.linalg.vector_norm(self.V.conj().T@self.V-torch.eye(self.rank,device='cuda',dtype=torch.complex128)).item()),
            CPU_seconds=time.process_time()-cput,wall_seconds=time.perf_counter()-wall,
            rule_scope='Explicit0.1noise-scale adaptation; not author exactrank',shared_across_five_arms=True)
        model.counts.update(receiver_Gram_builds=1,receiver_small_eigendecompositions=1)
    def project(self,X):return X-self.V@(self.V.conj().T@X)


class CudaImplicitFourier:
    def __init__(self,model,receiver,mf):
        self.model=model;self.receiver=receiver;self.mf=mf;self.n=model.side
        freq=np.rint(np.fft.fftfreq(self.n)*self.n).astype(int)
        mask=np.isin(freq,np.r_[np.arange(mf),np.arange(-mf,0)])
        self.ids_np=np.flatnonzero((mask[:,None,None]&mask[None,:,None]&mask[None,None,:]).ravel())
        self.ids=torch.as_tensor(self.ids_np,device='cuda');self.K=3*len(self.ids_np)
    def apply(self,gamma):
        b=gamma.shape[1];value=torch.empty((self.model.M,b),device='cuda',dtype=torch.complex128)
        for i in range(0,b,self.model.block_rhs):
            j=min(i+self.model.block_rhs,b);grid=torch.zeros((self.n**3,3,j-i),device='cuda',dtype=torch.complex128)
            grid[self.ids]=gamma[:,i:j].reshape(len(self.ids),3,j-i)
            cube=grid.reshape((self.n,)*3+(3,j-i))
            if self.model.fft_layout=='packed':
                cube=cube.permute(4,3,0,1,2).contiguous()
                value[:,i:j]=torch.fft.ifftn(cube,dim=(-3,-2,-1),norm='ortho').permute(2,3,4,1,0).reshape(self.model.M,j-i)
            else:value[:,i:j]=torch.fft.ifftn(cube,dim=(0,1,2),norm='ortho').reshape(self.model.M,j-i)
        self.model.counts['coefficient_IFFT_source_RHS']+=b
        return self.receiver.project(value)
    def adjoint(self,value):
        value=self.receiver.project(value);b=value.shape[1];gamma=torch.empty((self.K,b),device='cuda',dtype=torch.complex128)
        for i in range(0,b,self.model.block_rhs):
            j=min(i+self.model.block_rhs,b)
            spatial=value[:,i:j].reshape((self.n,)*3+(3,j-i))
            if self.model.fft_layout=='packed':
                packed=spatial.permute(4,3,0,1,2).contiguous()
                cube=torch.fft.fftn(packed,dim=(-3,-2,-1),norm='ortho').permute(2,3,4,1,0).reshape(self.n**3,3,j-i)
            else:cube=torch.fft.fftn(spatial,dim=(0,1,2),norm='ortho').reshape(self.n**3,3,j-i)
            gamma[:,i:j]=cube[self.ids].reshape(self.K,j-i)
        self.model.counts['coefficient_FFT_source_RHS']+=b
        return gamma
    def transport(self,old,gamma):
        positions=np.searchsorted(self.ids_np,old.ids_np)
        if not np.array_equal(self.ids_np[positions],old.ids_np):raise ValueError('Non-nested MF')
        new=torch.zeros((len(self.ids),3,gamma.shape[1]),device='cuda',dtype=torch.complex128)
        new[torch.as_tensor(positions,device='cuda')]=gamma.reshape(len(old.ids),3,gamma.shape[1])
        return new.reshape(self.K,gamma.shape[1])


def cuda_material_readout(J,Q,weights,beta):
    n=J.shape[0]//3;V=(Q+beta*J).reshape(n,3,-1);T=(beta*J).reshape(n,3,-1)
    den=torch.sum(abs(V)**2*weights,dim=(1,2));num=torch.sum(torch.real(V.conj()*T)*weights,dim=(1,2))
    r=torch.where(den>0,num/den,torch.zeros_like(num));r=torch.clamp(r,0,beta*4/(1+beta*4))
    return torch.clamp(r/(beta*(1-r)),0,4)


class CudaLocalArcReadout:
    def __init__(self,model):
        from .vector_original_extension import LocalArcReadout
        local=LocalArcReadout(model.cpu);self.ids=torch.as_tensor(local.ids,device='cuda');self.H=model.array(local.H);self.s=local.self_term
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
        seed=cuda_material_readout(J,Q,weights,beta)[ids];zs,vs,ws,rs=state(seed);use=loss(rs)<initial
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
                candidates=torch.stack([m[:,slot],torch.zeros_like(aa),torch.full_like(aa,4.),r1,r2],dim=1)
                valid=torch.isfinite(candidates)&(candidates>=0)&(candidates<=4)
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


def cuda_auxiliary_solve(model,receiver,Y,arm,config,callback=None,endpoint=None):
    if arm not in ('sourceCIE','sourceDual','GlobalDual','AplusGlobal','LSIEsource'):raise ValueError(arm)
    wd=1/torch.sum(abs(Y)**2,dim=0);ws=1/torch.sum(abs(model.inc)**2,dim=0)
    m=torch.zeros(model.N,device='cuda',dtype=torch.float64);J=receiver.Jd.clone();Q=receiver.Qd.clone();D=receiver.Dd.clone()
    arc=CudaLocalArcReadout(model) if arm=='AplusGlobal' else None
    history=[];stages=[];old_basis=None;gamma=None
    if endpoint:endpoint('air',m,None)
    for stage,mf in enumerate(config['mf']):
        bm=config['beta_material'][stage];bj=(config['beta_material'] if arm=='sourceCIE' else config['beta_current_dual'])[stage]
        basis=CudaImplicitFourier(model,receiver,mf)
        gamma=torch.zeros((basis.K,60),device='cuda',dtype=torch.complex128) if gamma is None else basis.transport(old_basis,gamma)
        previous_gradient=None;previous_direction=None;stage_start=time.perf_counter()
        for iteration in range(config['max_iterations_per_stage']):
            wall,cput=time.perf_counter(),time.process_time();old_gamma=gamma
            chi=torch.repeat_interleave(m,3)[:,None];p=torch.ones_like(chi) if arm=='LSIEsource' else bj/(1+bj*chi)
            residual=p*(J-chi*Q);v=p.conj()*residual*ws
            gradient=basis.adjoint(model.S.conj().T@(D*wd)+v-model.green(chi*v,True))
            dim=0 if arm in ('sourceCIE','sourceDual','LSIEsource') else None
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
            if arm=='LSIEsource':
                qr=Q.reshape(model.N,3,60);jr=J.reshape(model.N,3,60)
                den=torch.sum(abs(qr)**2*ws,dim=(1,2));num=torch.sum(torch.real(qr.conj()*jr)*ws,dim=(1,2))
                m=torch.clamp(torch.where(den>0,num/den,torch.zeros_like(num)),0,4);material_info={}
            elif arc:m,material_info=arc(J,Q,m,ws,bm)
            else:m=cuda_material_readout(J,Q,ws,bm);material_info={}
            den=torch.sum(abs(old_gamma)**2,dim=0);num=torch.sum(abs(gamma-old_gamma)**2,dim=0)
            valid=bool(torch.all((den>0)|(num==0)).item())
            delta_sq=float(torch.mean(torch.where(den>0,num/den,torch.zeros_like(num))).item()) if valid else None
            delta=float(np.sqrt(delta_sq)) if valid else None
            chi=torch.repeat_interleave(m,3)[:,None];p=torch.ones_like(chi) if arm=='LSIEsource' else bj/(1+bj*chi)
            data=float(torch.mean(torch.sum(abs(D)**2,dim=0)*wd).item())
            state=float(torch.mean(torch.sum(abs(p*(J-chi*Q))**2,dim=0)*ws).item());torch.cuda.synchronize()
            record=dict(stage=stage,mf=mf,iteration=iteration+1,beta_material=bm,beta_current=bj,beta_active=arm!='LSIEsource',
                current_representation='ordinary_LSIE_p1' if arm=='LSIEsource' else 'ordinary_scalar_CIE',
                coefficient_delta_RMS=delta,coefficient_delta_squared_mean=delta_sq,auxiliary_data=data,auxiliary_state=state,
                auxiliary_sum=data+state,peak_contrast=float(m.max().item()),material_readout=material_info,
                CPU_seconds=time.process_time()-cput,wall_seconds=time.perf_counter()-wall)
            history.append(record)
            if callback:callback(record,m)
            if delta is not None and delta<config['coefficient_delta_RMS_tolerance']:stop='PAPER_COEFFICIENT_STOP';break
        else:stop='CAP'
        stages.append(dict(stage=stage,mf=mf,iterations=iteration+1,stop=stop,wall_seconds=time.perf_counter()-stage_start))
        if endpoint:endpoint('stage'+str(stage+1),m,stages[-1])
        old_basis=basis
    return dict(material=m,history=history,stages=stages,receiver=receiver.info)
