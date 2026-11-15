"""Four public-data solvers with a fixed physical grid and shared budgets."""
import json
import time
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares, minimize
from threadpoolctl import threadpool_limits, threadpool_info
from .physics import Model, Regularizer, material_basis, resample_material

DEFAULTS = dict(wall_seconds=45., max_lu=2000, max_rhs=60000,
                stages=[8,12,24], stage_fractions=[.20,.25,.55],
                lower=.02, upper=20., initial_material=.1,
                smooth_weight=1e-4, tv_weight=1e-4, tv_epsilon=.01,
                dual_fraction=.45, maxiter=2000, maxfun=4000,
                ftol=1e-12, gtol=1e-8, trf_ftol=1e-10,
                trf_xtol=1e-10, trf_gtol=1e-8,
                candidate_limit=16, retain_material_trajectory=False)
METHODS = ('native_trf','reciprocal_trf','data_lbfgs','globaldual')

class BudgetStop(RuntimeError): pass
class PhaseStop(RuntimeError): pass


def json_safe(value):
    if isinstance(value, dict): return {str(k):json_safe(v) for k,v in value.items()}
    if isinstance(value, (list,tuple)): return [json_safe(v) for v in value]
    if isinstance(value, np.ndarray): return json_safe(value.tolist())
    if isinstance(value, np.generic): return json_safe(value.item())
    if isinstance(value, float) and not np.isfinite(value): return None
    return value


def load_public(source):
    if isinstance(source, (str,Path)):
        path=Path(source)
        with np.load(path,allow_pickle=False) as f: public={k:f[k] for k in f.files}
        if 'operators_path' in public:
            with np.load(path.parent/str(np.asarray(public['operators_path']).item()),allow_pickle=False) as f:
                public={**{k:f[k] for k in f.files},**public}
        return public
    return source


class Run:
    def __init__(self, public, method, config, backend):
        self.started,self.cpu_started=time.perf_counter(),time.process_time()
        self.cfg={**DEFAULTS,**(config or {})}
        if 'max_LU' in self.cfg:self.cfg['max_lu']=self.cfg['max_LU']
        if method not in METHODS:raise ValueError('Unknown method '+str(method))
        self.method,self.backend=method,backend
        if backend=='cuda':
            from .cuda_physics import CudaModel
            self.model=CudaModel(public)
        elif backend=='cpu':self.model=Model(public)
        else:raise ValueError('backend must be cpu or cuda')
        self.target=np.asarray(public.get('Y',public.get('y')),complex)
        self.shape=tuple(map(int,np.asarray(public.get('grid_shape',public.get('shape',[24,24])))))
        if self.target.shape!=(self.model.nr,self.model.ns):raise ValueError('Y shape must be (receivers,sources)')
        if self.shape!=(24,24) or (self.model.n,self.model.ns,self.model.nr)!=(576,12,48):
            raise ValueError('Comparison requires complete 24x24,12-source,48-receiver physics')
        for key,val in [('lower',.02),('upper',20.),('initial_material',.1),('smooth_weight',1e-4),('tv_weight',1e-4)]:
            if self.cfg[key]!=val:raise ValueError('Frozen common setting '+key)
        if list(self.cfg['stages'])!=[8,12,24] or list(self.cfg['stage_fractions'])!=[.2,.25,.55]:
            raise ValueError('Frozen 8->12->24 stages with .2/.25/.55 shares required')
        self.scale2=max(float(np.vdot(self.target,self.target).real),1e-300)
        self.reg=Regularizer(self.shape,self.cfg['smooth_weight'],self.cfg['tv_weight'],self.cfg['tv_epsilon'])
        self.trace=[];self.events=[];self.candidates=[];self.phase_rows=[]
        self.materials=[];self.material_seconds=[]
        self.best=None;self.stage=None;self.phase_end=np.inf
        self.phase_lu_end=self.phase_rhs_end=np.inf
        self.status='COMPLETED_NUMERICAL';self.stop_reason=None
        self.uniform_reference=None
        self.model.guard=self.guard

    def elapsed(self):return time.perf_counter()-self.started

    def guard(self):
        if self.elapsed()>=float(self.cfg['wall_seconds']):raise BudgetStop('WALL_SECONDS')
        if self.model.counts['LU']>=int(self.cfg['max_lu']):raise BudgetStop('MAX_LU')
        if self.model.counts['actual_complex_RHS']>=int(self.cfg['max_rhs']):raise BudgetStop('MAX_COMPLEX_RHS')
        if self.elapsed()>=self.phase_end:raise PhaseStop('PHASE_SECONDS')
        if self.model.counts['LU']>=self.phase_lu_end:raise PhaseStop('PHASE_LU')
        if self.model.counts['actual_complex_RHS']>=self.phase_rhs_end:raise PhaseStop('PHASE_COMPLEX_RHS')

    def evaluate(self,c,B,origin):
        state=self.model.forward(B@c,label=origin)
        residual=float(np.linalg.norm(state['y']-self.target)/np.sqrt(self.scale2))
        regularization=self.reg.evaluate(state['m'])[0]
        scalar=dict(objective=.5*residual**2+regularization,ranking_score=.5*residual**2+regularization,
                    data_objective=.5*residual**2,regularization=regularization,
                    relative_residual=residual,residual_norm=residual*np.sqrt(self.scale2),
                    origin=origin,stage=self.stage,seconds=self.elapsed(),counts=dict(self.model.counts),
                    state_relative_residual=state['state_relative_residual'])
        row={**scalar,'m':state['m'].copy(),'coarse_material':np.asarray(c).copy()}
        self.trace.append(scalar)
        if self.best is None or row['objective']<self.best['objective']:self.best=row
        if self.cfg['retain_material_trajectory']:
            self.materials.append(state['m'].copy());self.material_seconds.append(scalar['seconds'])
        if self.uniform_reference is None:
            self.uniform_reference=dict(mean_source_field_norm_squared=self.model.field_norm_squared(state)/self.model.ns,
                                        material=.1,seconds=self.elapsed(),counts=dict(self.model.counts))
        return state,row

    def correct(self,c,B,origin,wall_end,lu_end,rhs_end,dual=False,trf=False,reciprocal=False):
        wall_share=max(0.,float(wall_end)-self.elapsed())
        lu_share=max(0,int(lu_end)-self.model.counts['LU'])
        rhs_share=max(0,int(rhs_end)-self.model.counts['actual_complex_RHS'])
        self.phase_end=float(wall_end)
        self.phase_lu_end=int(lu_end)
        self.phase_rhs_end=int(rhs_end)
        start=time.perf_counter();cpu=time.process_time();start_counts=dict(self.model.counts)
        best=None;cached_x=None;cached_state=None;cached_row=None
        lo,hi=self.cfg['lower'],self.cfg['upper'];width=hi-lo
        c=np.clip(np.asarray(c,float),lo,hi)
        def unpack(x):return 1./x if reciprocal else x
        def state_at(x):
            nonlocal cached_x,cached_state,cached_row,best
            self.guard()
            if cached_x is None or not np.array_equal(x,cached_x):
                cached_state,cached_row=self.evaluate(unpack(x),B,origin)
                cached_x=np.asarray(x).copy()
                if best is None or cached_row['objective']<best['objective']:best=cached_row
            return cached_state,cached_row
        def objective(c):
            state,row=state_at(c)
            reg,greg,_=self.reg.evaluate(state['m'])
            if dual:
                value,g,details=self.model.smallgram(state,self.target,gradient=True)
                factor=.5*self.uniform_reference['mean_source_field_norm_squared']/self.scale2
                row['proposal_objective']=value*factor+reg
                self.trace[-1]['proposal_objective']=row['proposal_objective']
                self.trace[-1]['gram']=details
                return value*factor+reg,B.T@(factor*g+greg)
            return row['objective'],B.T@(self.model.data_gradient(state,self.target,self.scale2)+greg)
        def residual(x):
            state,_=state_at(x)
            d=self.reg.D@state['m'];q=np.sqrt(d*d+self.reg.epsilon**2)
            return np.concatenate((self.model.pack(state['y']-self.target)/np.sqrt(self.scale2),
                np.sqrt(2*self.reg.smooth/self.reg.edges)*d,
                np.sqrt(2*self.reg.tv/self.reg.edges)*d/np.sqrt(q+self.reg.epsilon)))
        def jacobian(x):
            state,_=state_at(x)
            chain=B*(-unpack(x)**2)[None,:] if reciprocal else B
            J=self.model.jacobian(state)@chain/np.sqrt(self.scale2)
            d=self.reg.D@state['m'];q=np.sqrt(d*d+self.reg.epsilon**2)
            DB=self.reg.D@chain
            td=np.sqrt(2*self.reg.tv/self.reg.edges)*(1/np.sqrt(q+self.reg.epsilon)-d*d/(2*q*(q+self.reg.epsilon)**1.5))
            with self.model.operation('TRF_ANALYTIC_RESIDUAL_J',trf_analytic_jacobians=1,parameter_dimension=B.shape[1]):
                return np.vstack((J,np.sqrt(2*self.reg.smooth/self.reg.edges)*DB,td[:,None]*DB))
        status='FINISHED';message='';reason=None;optimizer_status=None;success=None
        try:
            if trf:
                x=1./c if reciprocal else c
                bounds=(1./hi,1./lo) if reciprocal else (lo,hi)
                fit=least_squares(residual,x,jac=jacobian,method='trf',bounds=bounds,x_scale='jac',
                    tr_solver='exact',max_nfev=int(self.cfg['maxfun']),ftol=self.cfg['trf_ftol'],
                    xtol=self.cfg['trf_xtol'],gtol=self.cfg['trf_gtol'])
                final_x=fit.x
            else:
                def scaled(z):
                    value,gradient=objective(lo+width*z)
                    return value,width*gradient
                fit=minimize(scaled,(c-lo)/width,jac=True,method='L-BFGS-B',bounds=[(0.,1.)]*len(c),
                    options=dict(maxiter=int(self.cfg['maxiter']),maxfun=int(self.cfg['maxfun']),
                        ftol=self.cfg['ftol'],gtol=self.cfg['gtol'],maxls=25))
                final_x=lo+width*fit.x
            message=str(fit.message);optimizer_status=int(fit.status);success=bool(fit.success)
            state_at(final_x)
        except PhaseStop as exc:status='PHASE_BUDGET_STOP';reason=str(exc)
        except BudgetStop as exc:
            status='GLOBAL_BUDGET_STOP';reason=str(exc)
            self.status='BUDGET_STOP';self.stop_reason=reason
        finally:
            self.events.append(dict(kind='CORRECTION',origin=origin,stage=self.stage,dual=dual,
                trf=trf,reciprocal=reciprocal,status=status,stop_reason=reason,message=message,
                optimizer_status=optimizer_status,optimizer_success=success,
                wall_seconds=time.perf_counter()-start,cpu_seconds=time.process_time()-cpu,
                allocated_wall_seconds=float(wall_share),allocated_lu=int(lu_share),allocated_rhs=int(rhs_share),
                start_counts=start_counts,end_counts=dict(self.model.counts),seconds=self.elapsed(),
                optimizer_linear_algebra='scipy TRF dense exact SVD' if trf else 'scipy L-BFGS-B'))
            if best is not None:
                best={**best,'corrector_status':status,'corrector_message':message,'phase_stop_reason':reason}
                self.candidates.append(best)
                self.candidates=sorted(self.candidates,key=lambda r:r['objective'])[:int(self.cfg['candidate_limit'])]
        return best

    def run(self):
        c=np.full(64,.1);old_shape=(8,8)
        cumulative=np.cumsum(self.cfg['stage_fractions'])
        for index,(stage,fraction) in enumerate(zip(self.cfg['stages'],cumulative)):
            self.stage=stage
            if index:c=resample_material(c,old_shape,(stage,stage))
            B=material_basis(self.shape,stage)
            wall_end=float(self.cfg['wall_seconds'])*fraction
            lu_end=int(round(self.cfg['max_lu']*fraction))
            rhs_end=int(round(self.cfg['max_rhs']*fraction))
            first=len(self.trace)
            if stage==24 and self.method in ('data_lbfgs','globaldual'):
                fraction1=float(self.cfg['dual_fraction'])
                now=self.elapsed()
                phase1_wall=now+max(0.,wall_end-now)*fraction1
                current_lu=self.model.counts['LU'];current_rhs=self.model.counts['actual_complex_RHS']
                phase1_lu=current_lu+int(max(0,lu_end-current_lu)*fraction1)
                phase1_rhs=current_rhs+int(max(0,rhs_end-current_rhs)*fraction1)
                row=self.correct(c,B,'globaldual_proposal' if self.method=='globaldual' else 'data_matched_phase1',
                    phase1_wall,phase1_lu,phase1_rhs,dual=self.method=='globaldual')
                if row is not None:c=row['coarse_material']
                if self.stop_reason is None:
                    row=self.correct(c,B,'globaldual_data_correct' if self.method=='globaldual' else 'data_matched_restart',
                        wall_end,lu_end,rhs_end)
                    if row is not None:c=row['coarse_material']
            else:
                trf=self.method in ('native_trf','reciprocal_trf')
                row=self.correct(c,B,self.method+'_stage'+str(stage),wall_end,lu_end,rhs_end,trf=trf,
                    reciprocal=self.method=='reciprocal_trf')
                if row is not None:c=row['coarse_material']
            self.events.append(dict(kind='MULTISCALE_STAGE',stage=stage,parameter_dimension=B.shape[1],
                physical_dimension=self.model.n,has_evaluated_candidate=len(self.trace)>first,
                cumulative_wall_end=wall_end,cumulative_lu_end=lu_end,cumulative_rhs_end=rhs_end,
                seconds=self.elapsed(),counts=dict(self.model.counts)))
            old_shape=(stage,stage)
            if self.stop_reason is not None:break
        return self.finish()

    def finish(self):
        if self.backend=='cuda':self.model.synchronize()
        wall=self.elapsed();counts=dict(self.model.counts)
        top=self.best
        def metadata(row):return {k:v for k,v in row.items() if k not in ('m','coarse_material')}
        # Only correction endpoints/best corrected proposals are candidates. Scalar trace is an evaluation log.
        candidates=[{**metadata(r),'m':r['m']} for r in self.candidates[:16]]
        result=json_safe(dict(schema='r4-diversity-solver-v1',method=self.method,backend=self.backend,
            status=self.status,stop_reason=self.stop_reason,config=self.cfg,
            wall_seconds=wall,wall_s=wall,cpu_seconds=time.process_time()-self.cpu_started,
            m=None if top is None else top['m'],top1=None if top is None else metadata(top),
            candidates=candidates,trace=self.trace,events=self.events,costcounts=counts,
            budget_overshoot_seconds=max(0.,wall-self.cfg['wall_seconds']),
            lu_budget_overshoot=max(0,counts.get('LU',0)-self.cfg['max_lu']),
            complex_rhs_budget_overshoot=max(0,counts.get('actual_complex_RHS',0)-self.cfg['max_rhs']),
            accounting=dict(full_factorizations=counts.get('LU',0),full_complex_rhs=counts.get('actual_complex_RHS',0),
                source_rhs=counts.get('source_RHS',0),adjoint_rhs=counts.get('adjoint_RHS',0),
                fullj_receiver_transpose_rhs=counts.get('receiver_transpose_RHS',0),
                full_jacobians=counts.get('FullJ',0),small_gram=counts.get('Gram',0),
                wall_by_operation=dict(self.model.timings),cpu_by_operation=dict(self.model.cpu_timings)),
            all_unknowns_released=any(t['stage']==24 for t in self.trace),full_material_dimension=576,
            uniform_reference=self.uniform_reference,threadpools=threadpool_info(),
            ranking='original nonlinear .5||F(m)-Y||^2/||Y||^2 + common R3-DIGITS regularizer',
            truth_access='NONE_PUBLIC_INPUT_ONLY',candidate_count_is_branch_count=False,
            candidate_semantics='correction-phase best physical candidates; not distinct roots',
            optimizer_archive_is_multiroot_evidence=False))
        if self.cfg['retain_material_trajectory']:
            result['trajectory_arrays']={'m':np.asarray(self.materials,float),'seconds':np.asarray(self.material_seconds,float)}
        return result


def solve(public,method='native_trf',config=None,backend='cpu'):
    """Solve from public Y/operators. Optional trajectory_arrays must be saved as NPZ, not JSON."""
    if backend=='cuda':
        import torch
        torch.cuda.synchronize()
    with threadpool_limits(limits=1):
        run=Run(load_public(public),method,config,backend)
        return run.run()

solve_case=solve
