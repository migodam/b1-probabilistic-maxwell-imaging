"""Small full-wave continuation controls; existing Maxwell operators are reused.

No claim of global convexity. LSIE full-current penalty is known WRI/CSI
variable projection, distinct from the inherited fixed-beta CIE profile.
"""
from collections import Counter
import time
import numpy as np
from scipy.linalg import cho_factor, cho_solve, lu_solve
from scipy.optimize import minimize


class LandscapeRun:
    def __init__(self, model, Y, output, snapshot_stride=10):
        self.model, self.Y, self.output = model, Y, output
        self.yd = float(np.vdot(Y,Y).real)
        self.ed = float(np.vdot(model.E,model.E).real)
        self.wall, self.cpu = time.perf_counter(), time.process_time()
        self.cost = Counter(); self.history=[]; self.stages=[]
        self.budget=[]; self.budget_m=[]; self.best=None; self.calls=0
        self.snapshot_stride=snapshot_stride

    def physical_record(self, m, label, force=False):
        state=self.model.state(m)
        reg=self.model.reg.evaluate(m)[0]
        phi=float(.5*np.vdot(state['Y']-self.Y,state['Y']-self.Y).real/self.yd+reg)
        if self.best is None or phi<self.best['Phi']:
            self.best=dict(m=np.array(m).copy(),Phi=phi,label=label)
        self.calls+=1
        if force or self.calls==1 or self.calls%self.snapshot_stride==0:
            self.budget.append(dict(CPU_seconds=time.process_time()-self.cpu,
                wall_seconds=time.perf_counter()-self.wall,Phi=self.best['Phi'],label=label,
                physical_cost=dict(self.model.cost),auxiliary_cost=dict(self.cost)))
            self.budget_m.append(self.best['m'].astype(np.float32))
        return phi,state

    def physical_value_gradient(self,m,target,label):
        # Reuse the existing complete forward/adjoint; keep original denominator
        # and regularizer while the target moves along the data homotopy.
        val,g=self.model.data_value_gradient(m,target)
        rv,rg,_=self.model.reg.evaluate(m)
        ratio=float(np.vdot(target,target).real)/self.yd
        value=(val-rv)*ratio+rv; gradient=(g-rg)*ratio+rg
        phi,_=self.physical_record(m,label)
        self.history.append(dict(kind=label,objective=float(value),original_Phi=phi,
                                 gradient_inf=float(np.linalg.norm(gradient,np.inf))))
        return 1e6*float(value),1e6*gradient

    def minimize(self,fun,m,maximum,label,callback=None):
        wall,cpu=time.perf_counter(),time.process_time()
        result=minimize(fun,m,jac=True,method='L-BFGS-B',bounds=[(.02,20.)]*self.model.n,
            callback=callback,options=dict(maxiter=int(maximum),maxfun=4*int(maximum),
                                           ftol=1e-12,gtol=.01,maxls=30))
        phi,state=self.physical_record(result.x,label+'_endpoint',force=True)
        pg=result.x-np.clip(result.x-result.jac/1e6,.02,20.)
        row=dict(label=label,iterations=int(result.nit),evaluations=int(result.nfev),
                 stop=str(result.message),success=bool(result.success),objective=float(result.fun/1e6),
                 original_Phi=phi,projected_gradient_inf=float(np.linalg.norm(pg,np.inf)),
                 CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall)
        self.stages.append(row)
        np.savez_compressed(self.output/(label+'_physical.npz'),m=result.x,j=state['j'],
                            Q=state['Q'],Y_prediction=state['Y'])
        return result.x.copy()

    def profile_scale(self,m0):
        state=self.model.state(m0)
        with self.model._physics.operation('PROFILE_SCALE_RECEIVER_SOLVE',rhs=self.model.nr,
                                            receiver_RHS=self.model.nr):
            O=lu_solve(state['lu'],self.model.S.conj().T,trans=2,check_finite=False).conj().T
        gram=O@O.conj().T
        self.cost['port_Gram_builds']+=1;self.cost['port_eigendecompositions']+=1
        return float(np.linalg.eigvalsh(gram)[-1]*self.ed/self.yd)

    def profile_setup(self):
        S=self.model.S
        self.SS=S.conj().T@S/self.yd
        self.SY=S.conj().T@self.Y/self.yd
        self.cost['receiver_Gram_builds']+=1

    def profile_value_gradient(self,m,penalty,label,details=False):
        G,E,S=self.model.G,self.model.E,self.model.S
        L=self.model._physics.identity-m[:,None]*G
        b=m[:,None]*E
        K=self.SS+(penalty/self.ed)*(L.conj().T@L)
        rhs=self.SY+(penalty/self.ed)*(L.conj().T@b)
        fac=cho_factor(K,lower=True,check_finite=False)
        J=cho_solve(fac,rhs,check_finite=False)
        Q=E+G@J; e=J-m[:,None]*Q; d=S@J-self.Y
        rv,rg,_=self.model.reg.evaluate(m)
        data=float(.5*np.vdot(d,d).real/self.yd)
        state=float(.5*penalty*np.vdot(e,e).real/self.ed)
        grad=-(penalty/self.ed)*np.real(np.sum(e.conj()*Q,axis=1))+rg
        self.cost.update(current_Cholesky=1,current_RHS=self.model.ns,dense_state_Gram_builds=1,
                         Green_columns=self.model.ns,receiver_columns=self.model.ns)
        value=float(data+state+rv)
        self.history.append(dict(kind=label,objective=value,data=data,state=state,penalty=penalty,
                                 gradient_inf=float(np.linalg.norm(grad,np.inf))))
        if details:
            normal=float(np.linalg.norm(K@J-rhs)/np.linalg.norm(rhs))
            return value,grad,J,normal
        return 1e6*value,1e6*grad

    def run(self,m0,algorithm,config):
        initial_phi,initial_state=self.physical_record(m0,'initial',force=True)
        m=m0.copy(); extra={}
        if algorithm=='data_homotopy':
            Y0=initial_state['Y'].copy()
            for i,t in enumerate(config['data_homotopy_t']):
                target=(1-float(t))*Y0+float(t)*self.Y
                label=f'data_t{i}_{t}'
                m=self.minimize(lambda x:self.physical_value_gradient(x,target,label),m,
                                config['data_stage_maxiter'],label)
            extra['path']='Y_t=(1-t)F(m0)+tY; intermediate data are auxiliary, not newly measured'
        elif algorithm in ('profile_fixed','profile_homotopy','profile_restarted'):
            scale=self.profile_scale(m0);self.profile_setup()
            factors=(config['profile_factors'] if algorithm=='profile_homotopy' else
                     [1.]*len(config['profile_factors']) if algorithm=='profile_restarted' else [1.])
            maximum=config['profile_stage_maxiter'] if algorithm!='profile_fixed' else len(config['profile_factors'])*config['profile_stage_maxiter']
            extra.update(penalty_scale=scale,penalty_scale_definition='||S L(m0)^-1||2² * ||E||F² / ||Y||F²')
            for i,factor in enumerate(factors):
                penalty=scale*factor;label=f'profile_{i}_{factor}'; accepted=[0]
                def callback(x):
                    accepted[0]+=1
                    if accepted[0]%25==0:self.physical_record(x,label+f'_accepted{accepted[0]}',force=True)
                m=self.minimize(lambda x:self.profile_value_gradient(x,penalty,label),m,maximum,label,callback)
                value,g,J,normal=self.profile_value_gradient(m,penalty,label,True)
                self.stages[-1]['profile_normal_relative_residual']=normal
                np.savez_compressed(self.output/(label+'_auxiliary.npz'),m=m,J=J,penalty=penalty)
        elif algorithm not in ('physical_direct','weak_profile_handoff'):raise ValueError(algorithm)
        # Continuation transports its final iterate. No intermediate original-Phi
        # veto; the original-Phi incumbent is retained separately for final output.
        m=self.minimize(lambda x:self.physical_value_gradient(x,self.Y,'original_corrector'),m,
                         config['direct_maxiter'] if algorithm=='physical_direct' else config['corrector_maxiter'],
                         'original_corrector')
        phi,state=self.physical_record(self.best['m'],'selected_prediction',force=True)
        np.savez_compressed(self.output/'prediction.npz',m=self.best['m'],Y_prediction=state['Y'],
                            j=state['j'],Q=state['Q'])
        np.savez_compressed(self.output/'BUDGET_MATERIALS.npz',m=np.asarray(self.budget_m))
        return dict(status='COMPLETED_PUBLIC_PREDICTION',algorithm=algorithm,initial_Phi=initial_phi,
            Phi=phi,relative_residual=float(np.linalg.norm(state['Y']-self.Y)/np.sqrt(self.yd)),
            constitutive_defect=state['state_relative_residual'],selected_label=self.best['label'],
            online_selection_rule='Minimum original Phi among paid physical states; transported final auxiliary iterate starts original corrector',
            stages=self.stages,history=self.history,budget=self.budget,extra=extra,
            cost=dict(CPU_seconds=time.process_time()-self.cpu,wall_seconds=time.perf_counter()-self.wall,
                      physical=dict(self.model.cost),auxiliary=dict(self.cost)),truth_access='NONE')
