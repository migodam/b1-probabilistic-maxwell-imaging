"""R10 full original-Phi value/adjoint adapter of approved R9 derivative.
Real chi box [0,8], R0, mean-source normalized Phi without half.
Exact-material paid J retained for value-to-gradient upgrade; no extra forward.
No truth access. Both auxiliary and physical states are saved as 3N vectors.
"""
import time
import numpy as np
import torch
from scipy.optimize import minimize

class PhysicalObjective:
    def __init__(self,model,Y,progress=None,upper=8.):
        self.upper=float(upper)
        self.model=model;self.Y=model.array(Y);self.P=self.Y.shape[1]
        self.wd=1/torch.sum(abs(self.Y)**2,dim=0);self.records={};self.pending_current={};self.events=[];self.progress=progress
    @staticmethod
    def key(m):return np.asarray(m,dtype=np.float64).tobytes()
    def evaluate(self,m,owner,gradient=True):
        m=np.asarray(m,dtype=np.float64)
        if m.shape!=(self.model.N,) or not np.isfinite(m).all() or m.min()<0 or m.max()>self.upper:raise ValueError('Outside frozen fullmaterialbox')
        key=self.key(m);wall,cpu=time.perf_counter(),time.process_time();before=dict(self.model.counts);first=len(self.model.solve_records)
        hit=key in self.records;upgrade=False
        if hit and (not gradient or self.records[key].get('gradient') is not None):
            self.model.counts['physical_refinement_cache_hits']+=1;r=self.records[key]
        else:
            mt=torch.as_tensor(m,device=self.model.device,dtype=torch.float64);d=torch.repeat_interleave(mt,3)[:,None]
            air=bool(np.all(m==0));upgrade=hit
            if hit:
                r=self.records[key];J=self.pending_current.get(key);pred=self.model.array(r['prediction'])
                self.model.counts['paid_forward_reused_for_gradient']+=1
            else:
                if air:
                    J=None;pred=torch.zeros_like(self.Y);self.model.counts['analytic_air_value_calls']+=1
                else:
                    J=self.model.solve(mt,d*self.model.inc,label=owner+'_primal',progress=self.progress)
                    pred=self.model.S@J
                delta=pred-self.Y;phi=float(torch.mean(torch.sum(abs(delta)**2,dim=0)*self.wd).item())
                r=dict(id=f'state{len(self.records):03d}',m=m.copy(),prediction=pred.cpu().numpy(),Phi=phi,gradient=None,
                    first_owner=owner,air_analytic=air,uses=[],initial_paid_event=len(self.events))
                self.records[key]=r
            if gradient:
                delta=pred-self.Y;rhs=self.model.S.conj().T@(delta*self.wd)
                if air:
                    Q=self.model.inc;P=rhs;self.model.counts['analytic_air_gradient_calls']+=1
                    defect=0.
                else:
                    Q=self.model.inc+self.model.green(J)
                    P=self.model.solve(mt,rhs,adjoint=True,label=owner+'_adjoint',progress=self.progress)
                    defect=float((torch.linalg.vector_norm(J-d*Q)/torch.clamp(torch.linalg.vector_norm(d*self.model.inc),min=1e-300)).item())
                g=(2/self.P)*torch.real(torch.sum(P.conj()*Q,dim=1)).reshape(self.model.N,3).sum(dim=1)
                r.update(gradient=g.cpu().numpy(),constitutive_relative_residual=defect)
                r['J']=np.zeros((self.model.M,self.P),complex) if air else J.cpu().numpy()
                r['Q']=Q.cpu().numpy()
                r['projected_gradient_inf']=float(np.max(abs(m-np.clip(m-r['gradient'],0,self.upper))))
                self.pending_current.pop(key,None)
            elif J is not None:
                self.pending_current[key]=J
                r['J']=J.cpu().numpy();r['Q']=(self.model.inc+self.model.green(J)).cpu().numpy()
            self.model.counts['physical_refinement_value_gradient_evaluations' if gradient else 'physical_refinement_value_only_evaluations']+=1
        r['uses'].append(owner)
        records=self.model.solve_records[first:]
        event=dict(index=len(self.events),owner=owner,state_id=r['id'],cache_hit=hit,gradient_upgrade_from_paid_J=upgrade,
            gradient_requested=gradient,Phi=r['Phi'],projected_gradient_inf=r.get('projected_gradient_inf'),
            CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,
            counts_delta={k:int(v-before.get(k,0)) for k,v in self.model.counts.items() if v-before.get(k,0)},
            actual_primal_RHS=sum(not x['adjoint'] for x in records),actual_adjoint_RHS=sum(x['adjoint'] for x in records),
            maximum_relative_linear_residual=max([x.get('relative_residual',0.) for x in records] or[0.]))
        self.events.append(event);return r


class PhysicalRequestCap(RuntimeError): pass


def physical_corrector(objective, initial, config, owner, save, emit):
    """Common 80 accepted / strict200 objective requests, exact cache hits free.

    Callback copies already-paid states, and never asks for another solve.
    Endpoint selection uses only public Phi, never private material/error.
    """
    initial=np.asarray(initial,float);requests=accepted=0;best=None;result=None
    begin=len(objective.events);wall,cpu=time.perf_counter(),time.process_time()
    def consider(record):
        nonlocal best
        if best is None or record['Phi']<best['Phi']:best=record
    def fun(m):
        nonlocal requests
        if requests>=config['strict_maxfun']:raise PhysicalRequestCap('STRICT_OBJECTIVE_REQUEST_CAP')
        requests+=1;record=objective.evaluate(m,owner+'_request'+str(requests),True);consider(record)
        emit(dict(kind='physical_evaluation',request=requests,Phi=record['Phi'],PG_inf=record['projected_gradient_inf']))
        return record['Phi'],record['gradient']
    def callback(m):
        nonlocal accepted
        accepted+=1;record=objective.records[objective.key(m)];consider(record)
        if accepted in config['checkpoints_accepted']:save('physical_accepted'+str(accepted),record)
        emit(dict(kind='physical_accepted',accepted=accepted,Phi=record['Phi'],PG_inf=record['projected_gradient_inf']))
    try:
        result=minimize(fun,initial,method='L-BFGS-B',jac=True,
            bounds=[(0.,objective.upper)]*objective.model.N,callback=callback,
            options=dict(maxiter=config['accepted_maxiter'],maxfun=config['strict_maxfun'],
                maxls=config['maxls'],ftol=config['ftol'],gtol=config['gtol']))
        stop=str(result.message)
    except PhysicalRequestCap as exc:stop=str(exc)
    if best is None:raise RuntimeError('No paid physical-corrector state')
    save('physical_selected',best)
    return best,dict(requests=requests,accepted=accepted,stop=stop,
        scipy_success=bool(result.success) if result is not None else False,
        selected_state_id=best['id'],Phi=best['Phi'],PG_inf=best['projected_gradient_inf'],
        CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,
        events=objective.events[begin:],scope='Nested corrector interval, not another additive whole process',
        selection='Lowest public Phi among paid corrector states, no GT; rejected trials paid')

