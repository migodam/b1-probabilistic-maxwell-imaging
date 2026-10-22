"""Conventional fixed surrogate forces; supplied arrays only, no physics IO.

Both forces include the unchanged standard Gaussian66 prior in physical
probit coordinates x. They are proposal dynamics only. A separate caller
must perform the exact Full endpoint Hamiltonian MH decision.
"""
import numpy as np
from scipy.special import ndtr


def _fixed(a,shape,name):
    a=np.array(a,copy=True)
    if a.shape!=shape or not np.isfinite(a).all():raise ValueError(f'invalid {name}{shape}')
    a.setflags(write=False);return a


def native(x):
    x=np.asarray(x,float)
    if x.shape!=(66,):raise ValueError('full physical probit66 required')
    q=x.copy();q[64:]=ndtr(q[64:]);return q


def probit_derivative(x):
    return np.r_[np.ones(64),np.exp(-.5*np.asarray(x)[64:]**2)/np.sqrt(2*np.pi)]


def white(data,sigma):
    data=np.asarray(data)
    return np.concatenate((data.real.T,data.imag.T),axis=1).ravel()*np.sqrt(2)/sigma


def unwhite(vector,sigma):
    blocks=np.asarray(vector).reshape(8,64)
    return (blocks[:,:32].T+1j*blocks[:,32:].T)*sigma/np.sqrt(2)


class LinearResponseForce:
    kind='FROZEN_LINEAR_WHITE_RESPONSE'
    def __init__(self,x0,F0,Jx,y,sigma):
        self.x0=_fixed(x0,(66,),'x0');self.F0=_fixed(F0,(32,8),'F0')
        self.Jx=_fixed(Jx,(512,66),'Jx');self.y=_fixed(y,(32,8),'y')
        if not np.isfinite(sigma) or sigma<=0:raise ValueError('positive sigma required')
        self.sigma=float(sigma);self.r0=white(self.F0-self.y,self.sigma)
    def residual(self,x):return self.r0+self.Jx@(np.asarray(x)-self.x0)
    def prediction(self,x):return self.y+unwhite(self.residual(x),self.sigma)
    def energy(self,x):
        r=self.residual(x);return .5*float(r@r+np.asarray(x)@np.asarray(x))
    def gradient(self,x):return -self.Jx.T@self.residual(x)-np.asarray(x)


class MaterialResponseForce:
    kind='FROZEN_FINITE_MATERIAL_RESPONSE'
    def __init__(self,prior,x0,chi0,F0,Cchi,y,sigma):
        self.prior=prior;self.x0=_fixed(x0,(66,),'x0')
        self.chi0=_fixed(chi0,(576,),'chi0');self.F0=_fixed(F0,(32,8),'F0')
        self.Cchi=_fixed(Cchi,(256,576),'Cchi');self.y=_fixed(y,(32,8),'y')
        if not np.isfinite(sigma) or sigma<=0:raise ValueError('positive sigma required')
        self.sigma=float(sigma)
    def prediction(self,x):
        chi=self.prior.decode(native(x))[0]
        # Cchi rows are receiver-major/source-minor, exactly core's reshape.
        return self.F0+(self.Cchi@(chi-self.chi0)).reshape(32,8)
    def energy(self,x):
        r=self.prediction(x)-self.y
        return float(np.vdot(r,r).real/self.sigma**2+.5*np.asarray(x)@np.asarray(x))
    def gradient(self,x):
        q=native(x);chi=self.prior.decode(q)[0]
        residual=(self.F0+(self.Cchi@(chi-self.chi0)).reshape(32,8)-self.y).ravel()
        adjoint=self.Cchi.conj().T@residual
        likelihood=2*np.real(self.prior.derivative(q).conj().T@adjoint)/self.sigma**2
        return -likelihood*probit_derivative(x)-np.asarray(x)


def response_error(force,x,true_F,true_Phi):
    predicted=force.prediction(x);difference=predicted-true_F
    r=predicted-force.y;proxy_Phi=float(np.vdot(r,r).real/force.sigma**2)
    return dict(proxy_Phi=proxy_Phi,Phi_error=proxy_Phi-float(true_Phi),
        complex_response_RMS=float(np.sqrt(np.mean(np.abs(difference)**2))),
        white_response_RMS=float(np.sqrt(np.mean(white(difference,force.sigma)**2))),
        response_relative_error=float(np.linalg.norm(difference)/max(1e-30,np.linalg.norm(true_F))))
