"""Unchanged native R9 projected-IDFT and nested frequency transport.
Device generalized to model.device; no coefficient whitening or alternate basis.
"""
import numpy as np
import torch

class ImplicitFourier:
    def __init__(self,model,receiver,mf):
        self.model=model;self.receiver=receiver;self.mf=mf;self.n=model.side
        freq=np.rint(np.fft.fftfreq(self.n)*self.n).astype(int)
        mask=np.isin(freq,np.r_[np.arange(mf),np.arange(-mf,0)])
        self.ids_np=np.flatnonzero((mask[:,None,None]&mask[None,:,None]&mask[None,None,:]).ravel())
        self.ids=torch.as_tensor(self.ids_np,device=self.model.device);self.K=3*len(self.ids_np)
    def apply(self,gamma):
        b=gamma.shape[1];value=torch.empty((self.model.M,b),device=self.model.device,dtype=torch.complex128)
        for i in range(0,b,self.model.block_rhs):
            j=min(i+self.model.block_rhs,b);grid=torch.zeros((self.n**3,3,j-i),device=self.model.device,dtype=torch.complex128)
            grid[self.ids]=gamma[:,i:j].reshape(len(self.ids),3,j-i)
            cube=grid.reshape((self.n,)*3+(3,j-i))
            if self.model.fft_layout=='packed':
                cube=cube.permute(4,3,0,1,2).contiguous()
                value[:,i:j]=torch.fft.ifftn(cube,dim=(-3,-2,-1),norm='ortho').permute(2,3,4,1,0).reshape(self.model.M,j-i)
            else:value[:,i:j]=torch.fft.ifftn(cube,dim=(0,1,2),norm='ortho').reshape(self.model.M,j-i)
        self.model.counts['coefficient_IFFT_source_RHS']+=b
        return self.receiver.project(value)
    def adjoint(self,value):
        value=self.receiver.project(value);b=value.shape[1];gamma=torch.empty((self.K,b),device=self.model.device,dtype=torch.complex128)
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
        new=torch.zeros((len(self.ids),3,gamma.shape[1]),device=self.model.device,dtype=torch.complex128)
        new[torch.as_tensor(positions,device=self.model.device)]=gamma.reshape(len(old.ids),3,gamma.shape[1])
        return new.reshape(self.K,gamma.shape[1])
