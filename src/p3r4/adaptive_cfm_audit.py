"""Isolated CPU float64 continuous-CFM integration diagnostic, no deployment.

For dx/dt=v, dl/dt=-div(v). Backward1->0 with l(1)=0 yields
l_back=+integral_0^1 div(v)dt. Thus logq1=logp0(x0)-l_back.
The caller adds the frozen raw->network affine-map log determinant.
"""
from __future__ import annotations
import time
import numpy as np
import torch
from scipy.integrate import solve_ivp


class FullTraceODE:
    """Frozen real-valued network with an exact66-diagonal autograd trace."""
    def __init__(self,model,condition,guard=None):
        if model.kind!='cfm' or model.dimension!=66:raise ValueError('full66 CFM required')
        self.model=model.double().eval()
        self.model.requires_grad_(False)
        self.condition=torch.as_tensor(condition,dtype=torch.float64,device='cpu')
        if self.condition.shape!=(1,model.context_dim):raise ValueError('one frozen normalized condition required')
        self.guard=guard or (lambda:None)
        self.velocity_calls=0;self.trace_component_VJPs=0

    def rhs(self,t,state,trace):
        self.guard();self.velocity_calls+=1
        x=torch.tensor(np.asarray(state[:66]),dtype=torch.float64).reshape(1,66)
        if trace:
            x.requires_grad_(True)
            with torch.enable_grad():
                velocity=self.model.velocity(x,float(t),self.condition)
                divergence=x.new_zeros(())
                for d in range(66):
                    derivative=torch.autograd.grad(velocity[0,d],x,retain_graph=d<65)[0][0,d]
                    divergence=divergence+derivative
                    self.trace_component_VJPs+=1
            result=np.r_[velocity.detach().numpy()[0],-float(divergence.detach())]
        else:
            with torch.no_grad():result=self.model.velocity(x,float(t),self.condition).numpy()[0]
        if not np.isfinite(result).all():raise FloatingPointError('nonfinite full64 ODE RHS')
        return result

    def backward_density(self,standardized_x,base_mean,base_chol,outer_logdet,*,atol,rtol):
        start,cpu=time.monotonic(),time.process_time()
        before_calls,before_trace=self.velocity_calls,self.trace_component_VJPs
        # Every solve starts independently at the same true raw target point.
        solution=solve_ivp(lambda t,z:self.rhs(t,z,True),(1.,0.),np.r_[standardized_x,0.],
            method='DOP853',atol=atol,rtol=rtol)
        terminal=solution.y[:,-1]
        x0=torch.tensor(terminal[:66],dtype=torch.float64)[None]
        kwargs={}
        if base_mean is not None:
            kwargs=dict(base_mean=torch.tensor(base_mean,dtype=torch.float64),base_chol=torch.tensor(base_chol,dtype=torch.float64))
        with torch.no_grad():base_logp=float(self.model.base_log_prob(x0,**kwargs)[0])
        logq=base_logp-terminal[66]+float(outer_logdet)
        return dict(success=bool(solution.success),message=solution.message,
            raw_logq_x66=logq,raw_NLL_x66=-logq,base_logp=base_logp,
            backward_accumulated_negative_divergence=float(terminal[66]),x0_standardized=terminal[:66],
            accepted_time_points=len(solution.t),last_t=float(solution.t[-1]),
            velocity_calls=self.velocity_calls-before_calls,trace_component_VJPs=self.trace_component_VJPs-before_trace,
            scipy_nfev=solution.nfev,wall_seconds=time.monotonic()-start,cpu_seconds=time.process_time()-cpu)

    def forward_cycle(self,x0,*,atol,rtol):
        before=self.velocity_calls;start=time.monotonic()
        # Coordinate reconstruction needs only v; divergence was fully traced
        # in the backward augmented density solve, not stochastically sampled.
        solution=solve_ivp(lambda t,z:self.rhs(t,z,False),(0.,1.),np.asarray(x0),method='DOP853',atol=atol,rtol=rtol)
        return dict(success=bool(solution.success),message=solution.message,x1_standardized=solution.y[:,-1],
            last_t=float(solution.t[-1]),velocity_calls=self.velocity_calls-before,scipy_nfev=solution.nfev,
            wall_seconds=time.monotonic()-start)


def fixed32_double(model,condition,x,base_mean,base_chol,outer_logdet):
    kwargs={}
    if base_mean is not None:kwargs=dict(base_mean=torch.tensor(base_mean,dtype=torch.float64),base_chol=torch.tensor(base_chol,dtype=torch.float64))
    with torch.no_grad():
        value=model.log_prob(torch.tensor(x,dtype=torch.float64)[None],torch.tensor(condition,dtype=torch.float64),
            allow_approximate=True,nfe=32,**kwargs)
    return float(value[0])+float(outer_logdet)
