"""Original3lambda geometry, R1 CM+RR DDA and implicit native FFT-TSOM.

Independent DDA adaptation, not author integration/software. No work on import.
Fine-data-only active-current restriction is exact when the prescribed chi=0;
the inverse model always exposes all27000 real cells. No material prior.
"""
from collections import Counter
from pathlib import Path
import time
import numpy as np
from scipy import fft,linalg as la
from scipy.sparse.linalg import LinearOperator,gmres
from research.research_max_struct_r1.src.a9_green import dyadic_green
from research.research_max_struct_r8.src.vector_cie3d import box_fraction


def original_geometry():
    angle=2*np.pi*np.arange(20)/20;points=[];pol=[]
    for fixed,center in [(2,[.2,0,-.1]),(0,[.1,0,-.15]),(1,[-.05,.1,0])]:
        q=np.zeros((20,3));axes=[i for i in range(3) if i!=fixed]
        q[:,axes[0]]=3*np.cos(angle);q[:,axes[1]]=3*np.sin(angle);q+=center
        p=np.zeros((20,3));p[:,fixed]=1;points.append(q);pol.append(p)
    return np.concatenate(points),np.concatenate(pol)


def grid_points(side):
    ijk=np.indices((side,)*3).reshape(3,-1).T
    return (ijk-(side-1)/2)*(3./side)


def four_cubes(side):
    points=grid_points(side);h=3./side;chi=np.zeros(side**3)
    for center,value in [([-.7,.7,-.7],1.),([.7,-.7,-.7],2.),
                          ([-.7,-.7,.7],3.),([.7,.7,.7],4.)]:
        chi+=value*box_fraction(points,h,center,[.4]*3)
    return chi


class OriginalVectorModel:
    """Cached dyadic FFT, block external fields and qualified library GMRES."""
    def __init__(self,side,*,block_rhs=4,active_indices=None,generation_only=False,crop_generation_fft=False,
                 gmres_rtol=2e-6,gmres_restart=60,gmres_cycles=20):
        wall,cpu=time.perf_counter(),time.process_time()
        if side%2 or block_rhs not in (1,2,4):raise ValueError('Even grid and RHS block1/2/4 required')
        if active_indices is not None and not generation_only:
            raise ValueError('Known support is permitted only for synthetic data generation')
        self.side=side;self.width=3.;self.h=3./side;self.k=2*np.pi
        self.shape=(side,)*3;self.pad=(2*side,)*3;self.full_N=side**3
        self.active=np.arange(self.full_N) if active_indices is None else np.asarray(active_indices,int)
        if len(np.unique(self.active))!=len(self.active) or np.any(self.active<0) or np.any(self.active>=self.full_N):raise ValueError('Invalid active indices')
        self.generation_only=generation_only;self.N=len(self.active);self.M=3*self.N
        self.row_indices=(3*self.active[:,None]+np.arange(3)).ravel()
        self.points=grid_points(side)[self.active];self.P=60;self.nr=180
        self.block_rhs=block_rhs;self.rtol=gmres_rtol;self.restart=gmres_restart;self.cycles=gmres_cycles
        self.s=-1/3+1j*self.k**3*self.h**3/(6*np.pi)
        self.counts=Counter();self.solve_records=[];self.max_residual=0.;self._physical_cache=None
        self.sources,self.polarizations=original_geometry();self.receivers=self.sources.copy()
        self.configure_embedding(crop_generation_fft)
        self.inc=np.empty((self.M,60),complex);self.S=np.empty((180,self.M),complex)
        for i in range(0,60,block_rhs):
            j=min(i+block_rhs,60)
            g=dyadic_green(self.points[:,None,:]-self.sources[None,i:j,:],self.k)
            self.inc[:,i:j]=np.einsum('nspq,sq->nps',g,self.polarizations[i:j]).reshape(self.M,j-i)
            g=dyadic_green(self.receivers[i:j,None,:]-self.points[None,:,:],self.k)
            self.S[3*i:3*j]=self.h**3*g.transpose(0,2,1,3).reshape(3*(j-i),self.M)
        self.counts.update(incident_dipole_fields=60,receiver_vector_rows=180)
        self.setup=dict(CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,
            full_cells=self.full_N,stored_current_cells=self.N,complex_current_components=self.M,
            known_support_restriction=bool(active_indices is not None),block_rhs=block_rhs,
            kernel_bytes=self.kernel.nbytes,incident_bytes=self.inc.nbytes,receiver_bytes=self.S.nbytes,
            source_radius=3.,geometry_width=3.,self_real=float(self.s.real),self_imag=float(self.s.imag),
            convolution_shape=list(self.convolution_shape),convolution_padding=list(self.pad),
            generation_fft_crop=bool(crop_generation_fft),embedding_origin_indices=self.embedding_origin.tolist(),
            physical_model='R1CM+RRpointdipoleDDA; notauthorcellintegration')

    def configure_embedding(self,crop=False):
        """Exact zero-exterior elimination; h/absolute points/S/inc stay fixed."""
        if crop and not self.generation_only:raise ValueError('FFT bounding-box crop is generation-only')
        coordinates=np.array(np.unravel_index(self.active,self.shape)).T
        self.embedding_origin=coordinates.min(axis=0) if crop else np.zeros(3,int)
        self.convolution_shape=tuple((coordinates.max(axis=0)-self.embedding_origin+1).tolist()) if crop else self.shape
        self.convolution_N=int(np.prod(self.convolution_shape))
        embedding=np.ravel_multi_index((coordinates-self.embedding_origin).T,self.convolution_shape)
        self.green_rows=(3*embedding[:,None]+np.arange(3)).ravel()
        self.pad=tuple(fft.next_fast_len(2*d) for d in self.convolution_shape) if crop else tuple(2*d for d in self.shape)
        lags=[np.where(np.arange(p)<=p//2,np.arange(p),np.arange(p)-p) for p in self.pad]
        dv=np.stack(np.meshgrid(*lags,indexing='ij'),axis=-1)*self.h
        zero=np.all(dv==0,axis=-1);dv[zero]=[1,0,0]
        kernel=self.h**3*dyadic_green(dv,self.k);kernel[zero]=self.s*np.eye(3)
        self.kernel=fft.fftn(kernel,axes=(0,1,2));self.counts['Green_kernel_builds']+=1
        if hasattr(self,'setup'):
            self.setup.update(convolution_shape=list(self.convolution_shape),convolution_padding=list(self.pad),
                generation_fft_crop=bool(crop),embedding_origin_indices=self.embedding_origin.tolist(),kernel_bytes=self.kernel.nbytes)

    def green(self,X,adjoint=False):
        X=np.asarray(X,complex);one=X.ndim==1;X=X[:,None] if one else X
        if X.shape[0]!=self.M:raise ValueError('Current dimension mismatch')
        out=np.empty_like(X)
        # Conjugate+component transpose is the Hermitian Fourier multiplier.
        kg=self.kernel.conj().swapaxes(-1,-2) if adjoint else self.kernel
        for i in range(0,X.shape[1],self.block_rhs):
            j=min(i+self.block_rhs,X.shape[1]);b=j-i
            full=np.zeros((3*self.convolution_N,b),complex);full[self.green_rows]=X[:,i:j]
            cube=full.reshape(self.convolution_shape+(3,b))
            hat=fft.fftn(cube,s=self.pad,axes=(0,1,2))
            prod=np.einsum('...ij,...jb->...ib',kg,hat,optimize=True)
            value=fft.ifftn(prod,axes=(0,1,2))[tuple(slice(0,d) for d in self.convolution_shape)].reshape(3*self.convolution_N,b)
            out[:,i:j]=value[self.green_rows]
            self.counts['Green_adjoint_RHS' if adjoint else 'Green_forward_RHS']+=b
            self.counts['Green_block_calls']+=1
            self.counts['scalar_FFT3_transforms']+=6*b
        return out[:,0] if one else out

    def solve(self,chi,rhs,*,adjoint=False,label='',progress=None,checkpoint_dir=None):
        chi=np.asarray(chi,float);d=np.repeat(chi,3)
        if chi.shape!=(self.N,):raise ValueError('Material dimension mismatch')
        rhs=np.asarray(rhs,complex);out=np.empty_like(rhs)
        checkpoint_dir=Path(checkpoint_dir) if checkpoint_dir is not None else None
        if checkpoint_dir is not None:
            if not self.generation_only or adjoint:raise ValueError('Checkpoint support is fine-generation-only')
            checkpoint_dir.mkdir(parents=True,exist_ok=True)
            material_file=checkpoint_dir/'material.npz'
            if material_file.exists():
                saved=np.load(material_file)
                if not np.array_equal(saved['chi'],chi) or not np.array_equal(saved['active'],self.active):raise ValueError('Checkpoint material mismatch')
            else:np.savez_compressed(material_file,chi=chi,active=self.active)
        def action(v):return v-self.green(d*v,True) if adjoint else v-d*self.green(v)
        op=LinearOperator((self.M,self.M),matvec=action,dtype=complex)
        diag=1-d*self.s
        pre=LinearOperator((self.M,self.M),matvec=lambda v:v/(diag.conj() if adjoint else diag),dtype=complex)
        for source in range(rhs.shape[1]):
            checkpoint=checkpoint_dir/('rhs%03d.npz'%source) if checkpoint_dir else None
            if checkpoint is not None and checkpoint.exists():
                saved=np.load(checkpoint)
                if not np.array_equal(saved['rhs'],rhs[:,source]):raise ValueError('Checkpoint RHS mismatch')
                out[:,source]=saved['solution'];self.counts['resumed_qualified_primal_RHS']+=1
                if progress:progress(dict(label=label,source=source,status='RESUMED_QUALIFIED_RHS',additional_solve_CPU_seconds=0))
                continue
            wall,cpu=time.perf_counter(),time.process_time();iterations=[0]
            def callback(_):iterations[0]+=1
            try:
                out[:,source],info=gmres(op,rhs[:,source],M=pre,rtol=self.rtol,atol=0,
                    restart=self.restart,maxiter=self.cycles,callback=callback,callback_type='pr_norm')
            except BaseException as exc:
                self.counts['interrupted_RHS']+=1;self.counts['gmres_iterations']+=iterations[0]
                record=dict(label=label,source=source,adjoint=adjoint,iterations=iterations[0],status='INTERRUPTED',
                    error_type=type(exc).__name__,CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall)
                self.solve_records.append(record)
                if progress:progress(record)
                raise
            residual=float(la.norm(action(out[:,source])-rhs[:,source])/max(la.norm(rhs[:,source]),1e-300))
            self.max_residual=max(self.max_residual,residual)
            self.counts['adjoint_RHS' if adjoint else 'primal_RHS']+=1
            self.counts['gmres_iterations']+=iterations[0]
            record=dict(label=label,source=source,adjoint=adjoint,iterations=iterations[0],info=int(info),
                relative_residual=residual,CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall)
            self.solve_records.append(record)
            if progress:progress(record)
            if info or residual>max(2*self.rtol,1e-12):
                self.counts['failed_RHS']+=1
                raise RuntimeError('QUALIFIED_GMRES_FAILED: '+str(record))
            if checkpoint is not None:
                np.savez_compressed(checkpoint,rhs=rhs[:,source],solution=out[:,source],relative_residual=residual,
                    iterations=iterations[0],rtol=self.rtol)
                self.counts['qualified_source_checkpoints_written']+=1
        return out

    def physical(self,chi,Y=None,*,label='',progress=None,checkpoint_dir=None):
        if self._physical_cache is not None and np.array_equal(chi,self._physical_cache[0]):
            prediction=self._physical_cache[1];self.counts['cached_physical_hits']+=1
        else:
            self.counts['physical_evaluations']+=1
            J=self.solve(chi,np.repeat(chi,3)[:,None]*self.inc,label=label,progress=progress,checkpoint_dir=checkpoint_dir)
            prediction=self.S@J
            self._physical_cache=(np.array(chi).copy(),prediction.copy())
        if Y is None:return prediction.copy()
        # Paper Eq15 data portion averaged over sources. No added prior.
        loss=float(np.mean(np.sum(abs(prediction-Y)**2,axis=0)/np.sum(abs(Y)**2,axis=0)))
        return dict(prediction=prediction.copy(),data_loss_mean=loss,relative_source_RMS=float(np.sqrt(loss)))


class ReceiverSubspace:
    def __init__(self,model,Y,rtol=.1):
        wall,cpu=time.perf_counter(),time.process_time();S=model.S
        gram=S@S.conj().T
        ev,U=la.eigh(gram,check_finite=False);order=np.argsort(ev)[::-1]
        ev=np.maximum(ev[order],0.);U=U[:,order];singular=np.sqrt(ev)
        keep=singular>=rtol*singular[0];self.singular=singular;self.rank=int(keep.sum())
        if not self.rank:raise ValueError('No retained receiver modes')
        U=U[:,keep];sig=singular[keep]
        self.V=(S.conj().T@U)/sig[None,:]
        self.Jd=self.V@((U.conj().T@Y)/sig[:,None])
        self.info=dict(rank=self.rank,relative_cutoff=rtol,cutoff_scope='Explicitadaptation; notauthorrank',
            singular_values=singular.tolist(),right_vector_orthogonality=float(la.norm(self.V.conj().T@self.V-np.eye(self.rank))),
            method='180x180receiverGramHermitianEVD; nofullB/SB/KB/APCSSVD',
            CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall)
        model.counts.update(receiver_Gram_builds=1,receiver_small_eigendecompositions=1)

    def project(self,X):return X-self.V@(self.V.conj().T@X)


class ImplicitFourier:
    """Native coefficient Euclidean metric; orthonormalFFT common to all arms."""
    def __init__(self,model,receiver,mf):
        if model.N!=model.full_N:raise ValueError('Inverse basis may not use known support')
        self.model=model;self.receiver=receiver;self.mf=mf;self.n=model.side
        freq=np.fft.fftfreq(self.n)*self.n
        mask=np.isin(np.rint(freq).astype(int),np.r_[np.arange(mf),np.arange(-mf,0)])
        self.ids=np.flatnonzero((mask[:,None,None]&mask[None,:,None]&mask[None,None,:]).ravel())
        self.K=3*len(self.ids)

    def apply(self,gamma):
        b=gamma.shape[1];value=np.empty((self.model.M,b),complex)
        for i in range(0,b,self.model.block_rhs):
            j=min(i+self.model.block_rhs,b);grid=np.zeros((self.n**3,3,j-i),complex)
            grid[self.ids]=gamma[:,i:j].reshape(len(self.ids),3,j-i)
            value[:,i:j]=fft.ifftn(grid.reshape((self.n,)*3+(3,j-i)),axes=(0,1,2),norm='ortho').reshape(self.model.M,j-i)
        self.model.counts['coefficient_IFFT_source_RHS']+=b
        return self.receiver.project(value)

    def adjoint(self,value):
        value=self.receiver.project(value);b=value.shape[1];gamma=np.empty((self.K,b),complex)
        for i in range(0,b,self.model.block_rhs):
            j=min(i+self.model.block_rhs,b)
            cube=fft.fftn(value[:,i:j].reshape((self.n,)*3+(3,j-i)),axes=(0,1,2),norm='ortho').reshape(self.n**3,3,j-i)
            gamma[:,i:j]=cube[self.ids].reshape(self.K,j-i)
        self.model.counts['coefficient_FFT_source_RHS']+=b
        return gamma

    def transport(self,old,gamma):
        new=np.zeros((len(self.ids),3,gamma.shape[1]),complex)
        positions=np.searchsorted(self.ids,old.ids)
        if not np.array_equal(self.ids[positions],old.ids):raise ValueError('Fourier masks are not nested')
        new[positions]=gamma.reshape(len(old.ids),3,gamma.shape[1])
        return new.reshape(self.K,gamma.shape[1])


def material_readout(J,Q,weights,beta,upper=4.):
    """Exact shared-real R=0 CIE material minimizer (no material floor)."""
    n=J.shape[0]//3
    V=(Q+beta*J).reshape(n,3,-1);T=(beta*J).reshape(n,3,-1)
    den=np.sum(abs(V)**2*weights,axis=(1,2))
    num=np.sum(np.real(V.conj()*T)*weights,axis=(1,2))
    ratio=np.divide(num,den,out=np.zeros(n),where=den>0)
    ratio=np.clip(ratio,0,beta*upper/(1+beta*upper))
    return np.clip(ratio/(beta*(1-ratio)),0,upper)


class LocalArcReadout:
    """R=0 real-material disjoint 2^3-cell one-hop A+ coordinate readout.

    All blocks at a given local slot are independent and updated together.
    The endpoint is accepted only against its full fixed-current block metric.
    No replacement of the full Maxwell operator or full-physics rank-one claim.
    """
    def __init__(self,model):
        n=model.side
        self.ids=np.arange(n**3).reshape(n//2,2,n//2,2,n//2,2).transpose(0,2,4,1,3,5).reshape(-1,8)
        p=np.indices((2,2,2)).reshape(3,-1).T*model.h
        dv=p[:,None,:]-p[None,:,:];zero=np.all(dv==0,axis=-1);dv[zero]=[1,0,0]
        H=model.h**3*dyadic_green(dv,model.k);H[zero]=0
        self.H=H.transpose(0,2,1,3).reshape(24,24);self.self_term=model.s

    def __call__(self,J,Q,old,weights,beta):
        ids=self.ids;nb=len(ids);b=beta-self.self_term
        jb=J.reshape(-1,3,J.shape[1])[ids].reshape(nb,24,-1)
        qb=Q.reshape(-1,3,Q.shape[1])[ids].reshape(nb,24,-1)
        h=qb+b*jb
        def state(m):
            z=m/(1+b*m);zr=np.repeat(z,3,axis=1)
            v=jb-zr[:,:,None]*h;w=np.einsum('ij,bjs->bis',self.H,v,optimize=True)
            r=beta*(v+zr[:,:,None]*w)
            return z,v,w,r
        def loss(r):return np.sum(abs(r)**2*weights,axis=(1,2))
        m=old[ids].copy();z,v,w,r=state(m);initial=loss(r)
        seed=material_readout(J,Q,weights,beta)[ids];zs,vs,ws,rs=state(seed)
        use=loss(rs)<initial;seed_accepted=int(use.sum());m[use]=seed[use];z[use]=zs[use];v[use]=vs[use];w[use]=ws[use];r[use]=rs[use]
        for sweep in range(2):
            for slot in range(8):
                sl=slice(3*slot,3*slot+3);zi=z[:,slot].copy();hi=h[:,sl]
                B=-beta*np.repeat(z,3,axis=1)[:,:,None]*np.einsum('ij,bjs->bis',self.H[:,sl],hi,optimize=True)
                B[:,sl]=beta*(-hi+w[:,sl])
                C=r-zi[:,None,None]*B
                q=np.sum(abs(B)**2*weights,axis=(1,2));hh=np.sum(C.conj()*B*weights,axis=(1,2))
                aa=q*b.real+np.real(hh*b.conjugate()**2)
                bb=q+2*np.real(hh*b.conjugate());cc=hh.real
                # Stable quadratic roots; degenerate linear/constant cases retained.
                disc=bb*bb-4*aa*cc;root=np.sqrt(np.maximum(disc,0))
                stable=-.5*(bb+np.where(bb>=0,1.,-1.)*root)
                r1=np.divide(stable,aa,out=np.full(nb,np.nan),where=aa!=0)
                r2=np.divide(cc,stable,out=np.full(nb,np.nan),where=stable!=0)
                linear=np.divide(-cc,bb,out=np.full(nb,np.nan),where=bb!=0)
                r1=np.where(aa==0,linear,r1)
                candidates=np.stack([m[:,slot],np.zeros(nb),np.full(nb,4.),r1,r2],axis=1)
                valid=np.isfinite(candidates)&(candidates>=0)&(candidates<=4)
                valid[:,3:]&=(disc[:,None]>=0)|(aa[:,None]==0)
                zc=candidates/(1+b*candidates)
                values=q[:,None]*abs(zc)**2+2*np.real(hh[:,None]*zc)
                values[~valid]=np.inf
                best=np.argmin(values,axis=1);mn=candidates[np.arange(nb),best];zn=mn/(1+b*mn);dz=zn-zi
                r+=dz[:,None,None]*B
                dv=-dz[:,None,None]*hi
                w+=np.einsum('ij,bjs->bis',self.H[:,sl],dv,optimize=True)
                v[:,sl]+=dv;m[:,slot]=mn;z[:,slot]=zn
        # Full recomputation prevents accumulated coordinate roundoff from escaping acceptance.
        _,_,_,final=state(m);ok=np.isfinite(loss(final))&(loss(final)<=initial*(1+1e-12)+1e-25)
        m[~ok]=old[ids][~ok]
        out=old.copy();out[ids.ravel()]=m.ravel()
        return out,dict(blocks=nb,rejected_blocks=int((~ok).sum()),coordinate_sweeps=2,
            fixed_current_loss_before=float(initial.sum()),fixed_current_loss_after=float(np.where(ok,loss(final),initial).sum()),
            variant='R9vectorAplus_with_full_A_loss_accepted_ordinary_CIE_seed_proposal',
            seed_proposals=nb,seed_accepted_blocks=seed_accepted,seed_and_final_full_block_evaluations=3)


def auxiliary_solve(model,Y,arm,config,callback=None,endpoint=None):
    """Paper-style two-stage PR+ loop, no extra physical corrector or prior."""
    if arm not in ('sourceCIE','sourceDual','GlobalDual','AplusGlobal','LSIEsource'):raise ValueError(arm)
    receiver=ReceiverSubspace(model,Y,config['dpcs_relative_cutoff'])
    wd=1/np.sum(abs(Y)**2,axis=0);ws=1/np.sum(abs(model.inc)**2,axis=0)
    m=np.zeros(model.N);J=receiver.Jd.copy();Q=model.inc+model.green(J);D=model.S@J-Y
    arc=LocalArcReadout(model) if arm=='AplusGlobal' else None
    history=[];stages=[];old_basis=None;gamma=None
    if endpoint:endpoint('air',m,None)
    for stage,mf in enumerate(config['mf']):
        beta_m=config['beta_material'][stage]
        beta_j=(config['beta_material'] if arm=='sourceCIE' else config['beta_current_dual'])[stage]
        basis=ImplicitFourier(model,receiver,mf)
        gamma=np.zeros((basis.K,60),complex) if gamma is None else basis.transport(old_basis,gamma)
        previous_gradient=None;previous_direction=None;stage_start=time.perf_counter()
        for iteration in range(config['max_iterations_per_stage']):
            wall,cpu=time.perf_counter(),time.process_time();old_gamma=gamma.copy()
            chi=np.repeat(m,3)[:,None]
            p=np.ones_like(chi) if arm=='LSIEsource' else beta_j/(1+beta_j*chi)
            residual=p*(J-chi*Q);v=p.conj()*residual*ws
            pull=model.S.conj().T@(D*wd)+v-model.green(chi*v,True)
            gradient=basis.adjoint(pull)
            sourcewise=arm in ('sourceCIE','sourceDual','LSIEsource');axes=0 if sourcewise else None
            if previous_gradient is None:direction=-gradient
            else:
                denominator=np.sum(abs(previous_gradient)**2,axis=axes)
                coefficient=np.maximum(0,np.divide(np.real(np.sum(gradient.conj()*(gradient-previous_gradient),axis=axes)),denominator,
                    out=np.zeros_like(denominator,dtype=float),where=denominator>0))
                direction=-gradient+previous_direction*coefficient
                descent=np.real(np.sum(gradient.conj()*direction,axis=axes))
                if sourcewise:direction[:,descent>=0]=-gradient[:,descent>=0]
                elif descent>=0:direction=-gradient
            u=basis.apply(direction);gu=model.green(u);su=model.S@u
            hu=p*(u-chi*gu)
            numerator=-np.real(np.sum(gradient.conj()*direction,axis=axes))
            denominator=np.sum(abs(su)**2*wd,axis=axes)+np.sum(abs(hu)**2*ws,axis=axes)
            alpha=np.divide(numerator,denominator,out=np.zeros_like(numerator,dtype=float),where=denominator>0)
            gamma+=direction*alpha;J+=u*alpha;Q+=gu*alpha;D+=su*alpha
            previous_gradient=gradient;previous_direction=direction
            if arm=='LSIEsource':
                qr=Q.reshape(model.N,3,60);jr=J.reshape(model.N,3,60)
                den=np.sum(abs(qr)**2*ws,axis=(1,2));num=np.sum(np.real(qr.conj()*jr)*ws,axis=(1,2))
                m=np.clip(np.divide(num,den,out=np.zeros_like(num),where=den>0),0,4);material_info={}
            elif arc:m,material_info=arc(J,Q,m,ws,beta_m)
            else:m=material_readout(J,Q,ws,beta_m);material_info={}
            den=np.sum(abs(old_gamma)**2,axis=0);num=np.sum(abs(gamma-old_gamma)**2,axis=0)
            valid=np.all((den>0)|(num==0))
            delta_sq=float(np.mean(np.divide(num,den,out=np.zeros(60),where=den>0))) if valid else None
            delta=float(np.sqrt(delta_sq)) if valid else None
            chi=np.repeat(m,3)[:,None];p=np.ones_like(chi) if arm=='LSIEsource' else beta_j/(1+beta_j*chi)
            data=float(np.mean(np.sum(abs(D)**2,axis=0)*wd));state=float(np.mean(np.sum(abs(p*(J-chi*Q))**2,axis=0)*ws))
            record=dict(stage=stage,mf=mf,iteration=iteration+1,beta_material=beta_m,beta_current=beta_j,
                beta_active=arm!='LSIEsource',current_representation='ordinary_LSIE_p1' if arm=='LSIEsource' else 'ordinary_scalar_CIE',
                coefficient_delta_RMS=delta,coefficient_delta_squared_mean=delta_sq,auxiliary_data=data,
                auxiliary_state=state,auxiliary_sum=data+state,peak_contrast=float(m.max()),
                CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,material_readout=material_info)
            history.append(record)
            if callback:callback(record)
            if delta is not None and delta<config['coefficient_delta_RMS_tolerance']:
                stop='PAPER_COEFFICIENT_STOP';break
        else:stop='CAP'
        stages.append(dict(stage=stage,mf=mf,iterations=iteration+1,stop=stop,wall_seconds=time.perf_counter()-stage_start))
        if endpoint:endpoint('stage'+str(stage+1),m,stages[-1])
        old_basis=basis
    return dict(material=m,history=history,stages=stages,receiver=receiver.info)
