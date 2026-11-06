"""Independent short original-Phi controls over frozen CM CUDA operators.
No truth/data-family geometry enters evaluation or selection. No new backend.
"""
import time
import numpy as np
import torch
from scipy.optimize import minimize

class RequestBudgetStop(RuntimeError):pass

class PhysicalObjective:
    def __init__(self,model,Y,progress=None):
        self.model=model;self.Y=model.array(Y);self.P=self.Y.shape[1]
        self.wd=1/torch.sum(abs(self.Y)**2,dim=0);self.records={};self.pending_current={};self.events=[];self.progress=progress
    @staticmethod
    def key(m):return np.asarray(m,dtype=np.float64).tobytes()
    def evaluate(self,m,owner,gradient=True):
        m=np.asarray(m,dtype=np.float64)
        if m.shape!=(self.model.N,) or not np.isfinite(m).all() or m.min()<0 or m.max()>4:raise ValueError('Outside frozen fullmaterialbox')
        key=self.key(m);wall,cpu=time.perf_counter(),time.process_time();before=dict(self.model.counts);first=len(self.model.solve_records)
        hit=key in self.records;upgrade=False
        if hit and (not gradient or self.records[key].get('gradient') is not None):
            self.model.counts['physical_refinement_cache_hits']+=1;r=self.records[key]
        else:
            mt=torch.as_tensor(m,device='cuda',dtype=torch.float64);d=torch.repeat_interleave(mt,3)[:,None]
            air=bool(np.all(m==0));upgrade=hit
            if hit:
                r=self.records[key];J=self.pending_current[key];pred=self.model.array(r['prediction'])
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
                r['projected_gradient_inf']=float(np.max(abs(m-np.clip(m-r['gradient'],0,4))))
                self.pending_current.pop(key,None)
            elif J is not None:self.pending_current[key]=J
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

    def wiring(self,config):
        air=self.evaluate(np.zeros(self.model.N),'gradient_wiring_air',True);g=air['gradient'];scale=float(np.max(abs(g)))
        direction=np.maximum(-g,0)/scale if scale>0 else np.ones_like(g)
        fallback=not np.any(direction>0)
        if fallback:direction=np.ones_like(g)
        predicted=float(np.dot(g,direction));rows=[]
        for h in (config['h'],config['second_h']):
            value=self.evaluate(h*direction,f'gradient_wiring_FD_{h}',False)
            fd=(value['Phi']-air['Phi'])/h;absolute=abs(fd-predicted);relative=absolute/max(abs(fd),abs(predicted),1e-12)
            rows.append(dict(h=h,predicted_derivative=predicted,finite_difference=fd,absolute_difference=absolute,relative_difference=relative,
                state_id=value['id'],air_state_id=air['id']))
            if not (relative>1e-3 and absolute>1e-8):break
        last=rows[-1];good=np.isfinite(last['finite_difference']) and not(last['relative_difference']>1e-2 and last['absolute_difference']>1e-7)
        return dict(status='CONSISTENT_FIXED_AIR_WIRING' if good else 'REQUIRES_ROOT_REVIEW',checks=rows,
            all_positive_fallback=fallback,direction_max=float(direction.max()),air_projected_gradient_inf=air['projected_gradient_inf'],
            scope='one orconditionallytwo predeclaredone-sidedair probes; noGT/nootheridentityreplay')

def run_short(objective,start,name,config,save,emit):
    wall,cpu=time.perf_counter(),time.process_time();first_event=len(objective.events);first_count=dict(objective.model.counts)
    requests=0;accepted=0;best=None;rows=[];stop=None;res=None
    def observe(r):
        nonlocal best
        if best is None or r['Phi']<best['Phi']:best=r
    def fun(m):
        nonlocal requests
        if requests>=config['strict_request_cap']:raise RequestBudgetStop('STRICT_MAXFUN_24_BEFORE_NEW_REQUEST')
        requests+=1;r=objective.evaluate(m,f'{name}_request{requests:02d}',True);observe(r)
        row=dict(request=requests,state_id=r['id'],Phi=r['Phi'],PG_inf=r['projected_gradient_inf'],event_index=len(objective.events)-1)
        rows.append(row);emit(dict(kind='evaluation',start=name,**row));return r['Phi'],r['gradient']
    def callback(m):
        nonlocal accepted
        accepted+=1;key=objective.key(m)
        if key not in objective.records or objective.records[key]['gradient'] is None:raise RuntimeError('Accepted state notalready paid')
        r=objective.records[key];observe(r);save(f'accepted_{accepted:02d}',r)
        emit(dict(kind='accepted',start=name,accepted=accepted,Phi=r['Phi'],PG_inf=r['projected_gradient_inf'],state_id=r['id']))
    initial=objective.evaluate(start,name+'_initial',True);observe(initial);save('initial',initial)
    try:
        res=minimize(fun,np.array(start,float),method='L-BFGS-B',jac=True,bounds=[(0.,4.)]*objective.model.N,callback=callback,
            options={k:config[k] for k in ('maxiter','maxfun','maxls','ftol','gtol')})
        stop=str(res.message)
    except RequestBudgetStop as exc:stop=str(exc)
    save('prediction',best)
    return dict(status='TERMINAL_PUBLIC_PHYSICAL_REFINEMENT',start=name,requests=requests,accepted_steps=accepted,stop=stop,
        scipy_success=bool(res.success) if res is not None else False,scipy_nit=int(res.nit) if res is not None else None,
        initial_state_id=initial['id'],selected_state_id=best['id'],initial_Phi=initial['Phi'],Phi=best['Phi'],
        initial_projected_gradient_inf=initial['projected_gradient_inf'],projected_gradient_inf=best['projected_gradient_inf'],
        selection='MinimumpublicPhiamongallpaidstatesofthisstart; GTnotused; selectedtrialneednotbeaccepted',
        CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,events=objective.events[first_event:],requests_trace=rows,
        counts={k:int(v-first_count.get(k,0)) for k,v in objective.model.counts.items() if v-first_count.get(k,0)},
        additive=False,included_in='RUN_PROCESS_COST.json')
