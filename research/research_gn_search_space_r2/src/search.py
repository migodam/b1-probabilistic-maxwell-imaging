"""Public-only full-VIE continuation and finite directional searches."""
import importlib.util
import time
from pathlib import Path
import numpy as np
from scipy.linalg import svd
from scipy.optimize import least_squares
from threadpoolctl import threadpool_limits

HIST=Path(__file__).resolve().parents[2]/'research_z_alg_r4_diversity'
spec=importlib.util.spec_from_file_location('r2_original_physics',HIST/'src/physics.py')
physics=importlib.util.module_from_spec(spec);spec.loader.exec_module(physics)
Model,Regularizer=physics.Model,physics.Regularizer

class Stop(RuntimeError):pass
class PhaseStop(RuntimeError):pass

def safe(v):
    if isinstance(v,dict):return {str(k):safe(x)for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [safe(x)for x in v]
    if isinstance(v,np.ndarray):return safe(v.tolist())
    if isinstance(v,np.generic):return safe(v.item())
    if isinstance(v,float)and not np.isfinite(v):return None
    return v

def load_public(observation_id):
    p=HIST/'public_data';out={}
    for path in (p/'operators.npz',p/(observation_id+'.npz')):
        with np.load(path,allow_pickle=False)as f:out.update({k:f[k]for k in f.files})
    return out

def detail_bands(m):
    a=np.asarray(m).reshape(24,24)
    def means(n):return np.repeat(np.repeat(a.reshape(24//n,n,24//n,n).mean((1,3)),n,0),n,1)
    return [(a-means(2)).ravel(),(means(2)-means(4)).ravel()]

def ray_shrink(m,h):
    positive=h>0;negative=h<0;t=1.
    if np.any(positive):t=min(t,float(np.min((20.-m[positive])/h[positive])))
    if np.any(negative):t=min(t,float(np.min((.02-m[negative])/h[negative])))
    t=max(0.,t)
    if 0.<t<1.:t=float(np.nextafter(t,0.))
    while t>0 and (np.any(m+t*h<.02)or np.any(m+t*h>20.)):t=float(np.nextafter(t,0.))
    return t*h,t

class Search:
    def __init__(self,public,anchor,task):
        self.started=time.perf_counter();self.cpu=time.process_time();self.task=task
        self.model=Model(public);self.target=np.asarray(public['Y']);self.scale2=float(np.vdot(self.target,self.target).real)
        self.reg=Regularizer((24,24),smooth=1e-4*task.get('reg_scale',1.),tv=1e-4*task.get('reg_scale',1.))
        self.stopping_profile=task.get('stopping_profile','frozen')
        self.anchor=np.asarray(anchor,float).copy();self.best=None;self.trace=[];self.attempts=[];self.arrays={};self.phase_end=np.inf
        self.model.guard=self.guard
    def elapsed(self):return time.perf_counter()-self.started
    def guard(self):
        if self.elapsed()>=45.:raise Stop('WALL_SECONDS')
        if self.model.counts['LU']>=500:raise Stop('MAX_LU')
        if self.model.counts['actual_complex_RHS']>=30000:raise Stop('MAX_RHS')
        if self.elapsed()>=self.phase_end:raise PhaseStop('PHASE_SECONDS')
    def evaluate(self,m,label):
        state=self.model.forward(m,label=label)
        r=float(np.linalg.norm(state['y']-self.target)/np.sqrt(self.scale2));reg=self.reg.evaluate(m)[0]
        row=dict(objective=.5*r*r+reg,relative_residual=r,regularization=reg,label=label,seconds=self.elapsed(),counts=dict(self.model.counts))
        self.trace.append(row)
        item={**row,'m':np.asarray(m).copy()}
        if self.best is None or row['objective']<self.best['objective']:self.best=item
        return state,item
    def correct(self,m,label,seconds=45.,max_nfev=2000):
        start=self.elapsed();counts=dict(self.model.counts);self.phase_end=min(45.,start+seconds)
        cached_x=cached_state=cached_row=None;local_best=None;reason=None;status='FINISHED';message=None;nfev=None;optimality=None;active_count=None;projected_gradient_inf=None;gradient_inf=None
        def at(x):
            nonlocal cached_x,cached_state,cached_row,local_best
            self.guard()
            if cached_x is None or not np.array_equal(x,cached_x):
                cached_state,cached_row=self.evaluate(x,label);cached_x=np.asarray(x).copy()
                if local_best is None or cached_row['objective']<local_best['objective']:local_best=cached_row
            return cached_state
        def residual(x):
            state=at(x);d=self.reg.D@x;q=np.sqrt(d*d+self.reg.epsilon**2)
            return np.concatenate((self.model.pack(state['y']-self.target)/np.sqrt(self.scale2),np.sqrt(2*self.reg.smooth/self.reg.edges)*d,np.sqrt(2*self.reg.tv/self.reg.edges)*d/np.sqrt(q+self.reg.epsilon)))
        def jac(x):
            state=at(x);J=self.model.jacobian(state)/np.sqrt(self.scale2)
            d=self.reg.D@x;q=np.sqrt(d*d+self.reg.epsilon**2);D=self.reg.D.toarray()
            td=np.sqrt(2*self.reg.tv/self.reg.edges)*(1/np.sqrt(q+self.reg.epsilon)-d*d/(2*q*(q+self.reg.epsilon)**1.5))
            with self.model.operation('TRF_ANALYTIC_RESIDUAL_J',trf_analytic_jacobians=1,parameter_dimension=576):return np.vstack((J,np.sqrt(2*self.reg.smooth/self.reg.edges)*D,td[:,None]*D))
        try:
            fit=least_squares(residual,m,jac=jac,bounds=(.02,20.),method='trf',x_scale='jac',tr_solver='exact',max_nfev=max_nfev,ftol=1e-12 if self.stopping_profile=='tight' else 1e-10,xtol=1e-12 if self.stopping_profile=='tight' else 1e-10,gtol=1e-12 if self.stopping_profile=='tight' else 1e-8)
            message=str(fit.message);nfev=int(fit.nfev);optimality=float(fit.optimality);active_count=int(np.count_nonzero(fit.active_mask));gradient_inf=float(np.linalg.norm(fit.grad,np.inf))
            projected=fit.grad.copy();projected[(fit.x<=.02+1e-10)&(projected>0)]=0.;projected[(fit.x>=20.-1e-10)&(projected<0)]=0.;projected_gradient_inf=float(np.linalg.norm(projected,np.inf))
        except PhaseStop as exc:status='PHASE_STOP';reason=str(exc)
        except Stop as exc:status='BUDGET_STOP';reason=str(exc)
        except Exception as exc:status='FAILED';reason=type(exc).__name__+': '+str(exc)
        finally:
            self.phase_end=np.inf
            self.attempts.append(dict(label=label,status=status,reason=reason,message=message,nfev=nfev,optimality=optimality,active_bound_count=active_count,gradient_inf=gradient_inf,projected_gradient_inf=projected_gradient_inf,start_seconds=start,wall_seconds=self.elapsed()-start,start_counts=counts,end_counts=dict(self.model.counts),best=None if local_best is None else {k:v for k,v in local_best.items()if k!='m'}))
        return local_best
    def anchor_spectrum(self,state):
        self.guard();J=self.model.jacobian(state)/np.sqrt(self.scale2)
        with self.model.operation('ANCHOR_J_SVD',dense_SVD=1):_,s,Vt=svd(J,full_matrices=False,check_finite=False,lapack_driver='gesdd')
        self.arrays.update(anchor_J_singular_values=s,anchor_J_right_vectors=Vt.T)
        return s,Vt.T
    def run(self):
        status='COMPLETED';reason=None;geometry=None
        try:
            state,base=self.evaluate(self.anchor,'anchor');self.arrays['anchor']=self.anchor
            arm=self.task['arm'];needs_spectrum=self.task['experiment']=='e3' or self.task.get('reg_scale')==1.
            if arm=='green_pullback':
                from .space import prepare_geometry
                geometry=prepare_geometry(self.model,state,compute_j_svd=True)
                s=np.asarray(geometry['J_singular_values'])/np.sqrt(self.scale2);V=np.asarray(geometry['J_right_vectors'])
                self.arrays.update(anchor_J_singular_values=s,anchor_J_right_vectors=V,green_dark_basis=geometry['approx_basis'],green_eigenvalues=geometry['eigenvalues'])
            elif needs_spectrum:s,V=self.anchor_spectrum(state)
            if arm=='continue_trf':self.correct(self.anchor,'continue_trf')
            else:
                proposal_rows=[];self.proposals=proposal_rows;raw=[];filtered=[];seeds=[];afters=[]
                for band_index,band in enumerate(detail_bands(self.anchor)):
                    rms=float(np.sqrt(np.mean(band*band)))
                    for radius in (.03,.1):
                        for sign in (-1,1):
                            self.guard();p=sign*radius*band/max(rms,1e-300)
                            if arm=='green_pullback':
                                Q=geometry['approx_basis'];h=Q@(Q.T@p)
                            elif arm=='soft_jacobian':
                                epsilon=1e-4 if self.task['noise']=='clean' else 1e-3
                                lam=(epsilon/(radius*np.sqrt(576)))**2
                                h=V@((lam/(s*s+lam))*(V.T@p))
                            else:h=p
                            h,alpha=ray_shrink(self.anchor,h);seed=self.anchor+h
                            i=len(raw);raw.append(p);filtered.append(h);seeds.append(seed)
                            self.arrays.update(proposal_raw=np.asarray(raw),proposal_filtered=np.asarray(filtered),proposal_before=np.asarray(seeds),proposal_after=np.asarray(afters))
                            before_state,before=self.evaluate(seed,f'proposal_{i}_before')
                            corrected=self.correct(seed,f'proposal_{i}_correct',seconds=3.,max_nfev=8)
                            after=corrected if corrected is not None else before;afters.append(after['m']);self.arrays['proposal_after']=np.asarray(afters)
                            strong=s*.03*np.sqrt(576)>(1e-4 if self.task['noise']=='clean' else 1e-3)
                            shift=after['m']-seed
                            proposal_rows.append(dict(index=i,band=band_index+1,radius=radius,sign=sign,ray_scale=alpha,raw_RMS=float(np.sqrt(np.mean(p*p))),filtered_RMS=float(np.sqrt(np.mean(h*h))),before={k:v for k,v in before.items()if k!='m'},after={k:v for k,v in after.items()if k!='m'},correction_strong_norm=float(np.linalg.norm(V[:,strong].T@shift)),correction_weak_norm=float(np.linalg.norm(V[:,~strong].T@shift))))
                self.arrays.update(proposal_raw=np.asarray(raw),proposal_filtered=np.asarray(filtered),proposal_before=np.asarray(seeds),proposal_after=np.asarray(afters))
                self.proposals=proposal_rows
                if self.best is not None:self.correct(self.best['m'],'incumbent_refine')
        except Stop as exc:status='BUDGET_STOP';reason=str(exc)
        except Exception as exc:status='FAILED';reason=type(exc).__name__+': '+str(exc)
        stops=[a for a in self.attempts if a['status']=='BUDGET_STOP']
        if stops:status='BUDGET_STOP';reason=stops[-1]['reason']
        self.arrays['top1']=self.best['m']if self.best is not None else self.anchor
        result=dict(schema='gn-search-r2-v1',task=self.task,status=status,stop_reason=reason,wall_seconds=self.elapsed(),cpu_seconds=time.process_time()-self.cpu,costcounts=dict(self.model.counts),operation_seconds=dict(self.model.timings),top1=None if self.best is None else {k:v for k,v in self.best.items()if k!='m'},trace=self.trace,attempts=self.attempts,proposals=getattr(self,'proposals',[]),geometry=None if geometry is None else {'diagnostics':geometry['diagnostics'],'costs':geometry['costs']},truth_access='NONE_PUBLIC_INPUT_ONLY',ranking='original nonlinear .5||F(m)-Y||^2/||Y||^2 + R',budget=dict(wall_seconds=45,max_LU=500,max_complex_RHS=30000),wall_overshoot_seconds=max(0.,self.elapsed()-45.))
        return safe(result),self.arrays

def solve(public,anchor,task):
    with threadpool_limits(limits=1):return Search(public,anchor,task).run()
