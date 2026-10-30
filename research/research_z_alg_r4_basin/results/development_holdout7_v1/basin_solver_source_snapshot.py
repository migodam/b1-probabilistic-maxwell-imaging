"""Public-only bounded multiscale solvers and separated rational-valley correction."""
import time
from collections import Counter
import numpy as np
from scipy.linalg import lu_solve, cho_factor, cho_solve
from scipy.optimize import minimize, least_squares
from threadpoolctl import threadpool_limits, threadpool_info
from research.research_z_alg_r3_digits.src.solver import load_public, json_safe
from .physics import Model, Regularizer, material_basis, resample_material
from .rational_rays import scan_ray

DEFAULTS=dict(wall_seconds=75.,max_lu=2000,stages=[8,12,24],stage_fractions=[.20,.25,.55],
              lower=.5,upper=20.,initial_material=1.,smooth_weight=1e-4,tv_weight=1e-4,
              tv_epsilon=.01,maxiter=2000,maxfun=4000,ftol=1e-12,gtol=1e-8,
              trf_ftol=1e-10,trf_xtol=1e-10,trf_gtol=1e-8,
              dual_fraction=.45,anchor_fraction=.25,ray_fraction=.25,
              ray_points=33,ray_refine_iterations=10,ray_refine_tolerance=1e-5,
              ray_directions=['amplitude','gn'],candidate_rms=.25,
              gn_damping=.01)


class BudgetStop(RuntimeError): pass
class PhaseStop(RuntimeError): pass


class BasinRun:
    def __init__(self,public,method,config,initial):
        self.started,self.cpu_started=time.perf_counter(),time.process_time()
        self.cfg=dict(DEFAULTS);self.cfg.update(config or {})
        if 'max_LU' in self.cfg:self.cfg['max_lu']=self.cfg['max_LU']
        if 'maxLU' in self.cfg:self.cfg['max_lu']=self.cfg['maxLU']
        if method not in ('data_lbfgs','globaldual','trf','basin_bridge','ray_best_only'):
            raise ValueError('Unknown method '+method)
        self.method=method;self.model=Model(public)
        self.target=np.asarray(public.get('Y',public.get('y')),complex)
        self.shape=tuple(map(int,np.asarray(public.get('grid_shape',public.get('shape',[24,24])))))
        self.scale2=float(np.vdot(self.target,self.target).real)
        self.reg=Regularizer(self.shape,self.cfg['smooth_weight'],self.cfg['tv_weight'],self.cfg['tv_epsilon'])
        self.initial=initial
        self.trajectory=[];self.candidates=[];self.raw_proposals=[];self.rays=[];self.events=[]
        self.status='COMPLETED_NUMERICAL';self.stop_reason=None
        self.phase_end=np.inf;self.stage=None;self.uniform_reference=None
        self.model.guard=self.guard

    def elapsed(self):return time.perf_counter()-self.started

    def guard(self):
        if self.elapsed()>=float(self.cfg['wall_seconds']):raise BudgetStop('WALL_SECONDS')
        if self.model.counts['LU']>=int(self.cfg['max_lu']):raise BudgetStop('MAX_LU')
        if self.elapsed()>=self.phase_end:raise PhaseStop('PHASE_SECONDS')

    def row(self,state,origin,coarse=None):
        residual=float(np.linalg.norm(state['y']-self.target)/np.sqrt(self.scale2))
        regularization=self.reg.evaluate(state['m'])[0]
        return dict(m=state['m'].copy(),coarse_material=None if coarse is None else np.asarray(coarse).copy(),
                    objective=.5*residual*residual+regularization,ranking_score=.5*residual*residual+regularization,
                    data_objective=.5*residual*residual,regularization=regularization,
                    relative_residual=residual,residual_norm=residual*np.sqrt(self.scale2),origin=origin,
                    stage=self.stage,seconds=self.elapsed(),counts=dict(self.model.counts),
                    state_relative_residual=state['state_relative_residual'])

    def evaluate(self,c,B,origin):
        m=B@c
        state=self.model.forward(m,label=origin)
        row=self.row(state,origin,c)
        self.trajectory.append(row)
        if self.uniform_reference is None and np.max(np.abs(m-1))<1e-12:
            self.uniform_reference=dict(mean_source_field_norm_squared=float(np.vdot(state['Q'],state['Q']).real)/self.model.ns,
                                        seconds=self.elapsed(),counts=dict(self.model.counts))
        return state,row

    def data_gradient(self,state):
        with self.model.operation('DATA_ADJOINT',rhs=self.model.ns,adjoint_RHS=self.model.ns,data_adjoint_gradients=1):
            adj=lu_solve(self.model.factor,self.model.SH@(state['y']-self.target),trans=2,check_finite=False)
            return np.real(np.sum(adj.conj()*state['Q'],axis=1))/self.scale2

    def correct(self,start,B,end,origin,dual=False,trf=False):
        """Every correction records every paid physical trial and returns its best original score."""
        self.phase_end=end
        c=np.clip(np.asarray(start,float),self.cfg['lower'],self.cfg['upper'])
        best=None;cache_c=None;cache_state=None;cache_row=None
        def state_at(x):
            nonlocal cache_c,cache_state,cache_row,best
            self.guard()
            if cache_c is None or not np.array_equal(x,cache_c):
                cache_state,cache_row=self.evaluate(x,B,origin)
                cache_c=np.asarray(x).copy()
                if best is None or cache_row['objective']<best['objective']:best=cache_row
            # An optimizer Jacobian follows the same cached forward and owns its LU.
            return cache_state,cache_row
        def objective(x):
            state,row=state_at(x)
            if dual:
                value,g,details=self.model.smallgram(state,self.target,gradient=True)
                factor=.5*self.uniform_reference['mean_source_field_norm_squared']/self.scale2
                reg,greg,_=self.reg.evaluate(state['m'])
                return value*factor+reg,B.T@(g*factor+greg)
            g=self.data_gradient(state)+self.reg.evaluate(state['m'])[1]
            return row['objective'],B.T@g
        def residual(x):
            state,_=state_at(x)
            d=self.reg.D@state['m'];q=np.sqrt(d*d+self.reg.epsilon**2)
            return np.concatenate((self.model.pack(state['y']-self.target)/np.sqrt(self.scale2),
                                   np.sqrt(2*self.reg.smooth/self.reg.edges)*d,
                                   np.sqrt(2*self.reg.tv/self.reg.edges)*d/np.sqrt(q+self.reg.epsilon)))
        def jacobian(x):
            state,_=state_at(x)
            J=self.model.jacobian(state)@B/np.sqrt(self.scale2)
            d=self.reg.D@state['m'];q=np.sqrt(d*d+self.reg.epsilon**2)
            DB=self.reg.D@B
            tv_derivative=np.sqrt(2*self.reg.tv/self.reg.edges)*(1/np.sqrt(q+self.reg.epsilon)-d*d/(2*q*(q+self.reg.epsilon)**1.5))
            with self.model.operation('TRF_ANALYTIC_RESIDUAL_J',trf_analytic_jacobians=1,parameter_dimension=B.shape[1]):
                return np.vstack((J,np.sqrt(2*self.reg.smooth/self.reg.edges)*DB,tv_derivative[:,None]*DB))
        status='FINISHED';message=''
        correction_started,correction_cpu=time.perf_counter(),time.process_time()
        try:
            if trf:
                fit=least_squares(residual,c,jac=jacobian,method='trf',bounds=(self.cfg['lower'],self.cfg['upper']),
                                  ftol=self.cfg['trf_ftol'],xtol=self.cfg['trf_xtol'],gtol=self.cfg['trf_gtol'],
                                  max_nfev=int(self.cfg['maxfun']),tr_solver='exact',x_scale='jac')
            else:
                fit=minimize(objective,c,method='L-BFGS-B',jac=True,
                             bounds=[(self.cfg['lower'],self.cfg['upper'])]*len(c),
                             options=dict(maxiter=int(self.cfg['maxiter']),maxfun=int(self.cfg['maxfun']),
                                          ftol=self.cfg['ftol'],gtol=self.cfg['gtol'],maxls=25))
            message=str(fit.message)
            state_at(fit.x)
        except PhaseStop:status='PHASE_BUDGET_STOP'
        except BudgetStop:
            status='GLOBAL_BUDGET_STOP'
            raise
        finally:
            self.events.append(dict(kind='CORRECTION',origin=origin,dual=dual,trf=trf,status=status,
                                    message=message,stage=self.stage,seconds=self.elapsed(),
                                    has_physical_candidate=best is not None,
                                    wall_seconds=time.perf_counter()-correction_started,
                                    cpu_seconds=time.process_time()-correction_cpu,
                                    optimizer_linear_algebra=('scipy TRF exact dense internal SVD, inclusive time; opaque call count' if trf else 'L-BFGS-B internal limited-memory work, inclusive time')))
        return best

    def retain(self,row,**metadata):
        if row is None:return
        row=dict(row);row.update(metadata)
        # Corrected candidates are deduplicated by material RMS, after correction.
        for i,old in enumerate(self.candidates):
            if np.sqrt(np.mean((old['m']-row['m'])**2))<self.cfg['candidate_rms']:
                if row['objective']<old['objective']:self.candidates[i]=row
                return
        self.candidates.append(row)

    def direction(self,state,name):
        m=state['m']
        if name=='amplitude':return m-1.
        if name=='shape':return m-np.mean(m)
        if name=='offset':return np.ones_like(m)
        J=self.model.jacobian(state)/np.sqrt(self.scale2)
        r=self.model.pack(state['y']-self.target)/np.sqrt(self.scale2)
        _,greg,weights=self.reg.evaluate(m)
        g=J.T@r+greg
        if name=='data':return -g/max(np.sqrt(np.mean(g*g)),1e-30)
        if name=='dual':
            _,g,_=self.model.smallgram(state,self.target)
            return -g/max(np.sqrt(np.mean(g*g)),1e-30)
        if name!='gn':raise ValueError('Unknown ray direction '+name)
        with self.model.operation('RAY_GN_DIRECTION',normal_factorizations=1,normal_dimension=self.model.n):
            Hess=J.T@J+(self.reg.D.T@self.reg.D.multiply(weights[:,None])).toarray()
            damping=self.cfg['gn_damping']*max(float(np.max(np.diag(Hess))),1e-12)
            Hess.flat[::self.model.n+1]+=damping
            p=cho_solve(cho_factor(Hess,lower=True,check_finite=False),-g,check_finite=False)
        return np.clip(m+p,self.cfg['lower'],self.cfg['upper'])-m

    def bridge(self,anchor,B,end):
        self.phase_end=end
        anchor_c=anchor['coarse_material']
        # Final stage is the identity physical basis. Each ray uses the common current anchor.
        proposals=[]
        scan_end=self.elapsed()+max(0.,end-self.elapsed())*self.cfg['ray_fraction']
        self.phase_end=scan_end
        for name in self.cfg['ray_directions']:
            record=dict(direction=name,status='STARTED',nodes=[],minima=[])
            self.rays.append(record)
            try:
                state,_=self.evaluate(anchor_c,B,'ray_anchor_'+name)
                p=self.direction(state,name)
                found,record=scan_ray(self.model,state,p,self.target,self.scale2,self.reg,self.cfg,name,record)
            except PhaseStop:
                record['status']='SCAN_PHASE_BUDGET_STOP'
                break
            for proposal in found:
                proposal['proposal_id']=len(self.raw_proposals)
                self.raw_proposals.append(proposal)
                # Pre-correction dedup only removes materially identical ray minima.
                if not any(np.sqrt(np.mean((proposal['m']-old['m'])**2))<self.cfg['candidate_rms'] for old in proposals):
                    proposals.append(proposal)
        self.phase_end=end
        if self.method=='ray_best_only' and proposals:
            selected=[min(proposals,key=lambda q:q['raw_score'])]
            for q in proposals:q['selection']='RAW_TOP1_ABLATION' if q is selected[0] else 'DISCARDED_RAW_SCORE_ABLATION'
        else:
            # Enumeration order, never raw-score ordering: all distinct valleys receive a correction slot.
            selected=proposals
            for q in selected:q['selection']='CORRECT_ALL_DISTINCT_VALLEYS'
        for i,proposal in enumerate(selected):
            remaining=len(selected)-i
            slot_end=self.elapsed()+max(0.,end-self.elapsed())/remaining
            row=self.correct(proposal['m'],B,slot_end,'ray_correct_'+str(proposal['proposal_id']))
            proposal['corrected']=row is not None
            proposal['corrected_score']=None if row is None else row['objective']
            proposal['correction_seconds_end']=self.elapsed()
            self.retain(row,proposal_id=proposal['proposal_id'],ray=proposal['ray'],t=proposal['t'],
                        candidate_kind='CORRECTED_RATIONAL_VALLEY')

    def run(self):
        stages=self.cfg['stages'];fractions=np.asarray(self.cfg['stage_fractions'],float)
        ends=np.cumsum(fractions/np.sum(fractions))*self.cfg['wall_seconds']
        coarse=np.full(int(stages[0])**2,float(self.cfg['initial_material']))
        if self.initial is not None:
            initial=np.asarray(self.initial,float).reshape(self.shape)
            yi=np.rint(np.linspace(0,self.shape[0]-1,int(stages[0]))).astype(int)
            xi=np.rint(np.linspace(0,self.shape[1]-1,int(stages[0]))).astype(int)
            coarse=initial[np.ix_(yi,xi)].reshape(-1)
        old_shape=(int(stages[0]),int(stages[0]))
        for index,(stage,end) in enumerate(zip(stages,ends)):
            self.stage=int(stage);self.phase_end=float(end)
            if index:coarse=resample_material(coarse,old_shape,(stage,stage))
            B=material_basis(self.shape,stage)
            final=index==len(stages)-1
            stage_start=len(self.trajectory)
            try:
                if final and self.method in ('basin_bridge','ray_best_only'):
                    anchor_end=self.elapsed()+max(0.,end-self.elapsed())*self.cfg['anchor_fraction']
                    anchor=self.correct(coarse,B,anchor_end,'common_data_anchor')
                    self.retain(anchor,candidate_kind='COMMON_LOCAL_ANCHOR')
                    if anchor is not None:self.bridge(anchor,B,float(end))
                elif final and self.method=='globaldual':
                    dual_end=self.elapsed()+max(0.,end-self.elapsed())*self.cfg['dual_fraction']
                    # Use the dual optimizer's physically best visited state as its proposal seed.
                    dualrow=self.correct(coarse,B,dual_end,'globaldual_proposal',dual=True)
                    if dualrow is not None:coarse=dualrow['coarse_material']
                    row=self.correct(coarse,B,float(end),'globaldual_data_correct')
                    self.retain(row,candidate_kind='GLOBALDUAL_DATA_CORRECTED')
                else:
                    row=self.correct(coarse,B,float(end),'trf' if self.method=='trf' else 'data_lbfgs',trf=self.method=='trf')
                    if final:self.retain(row,candidate_kind='LOCAL_BASELINE')
            except PhaseStop:
                self.events.append(dict(kind='RAY_OR_STAGE_PHASE_STOP',stage=self.stage,seconds=self.elapsed()))
            finally:
                rows=self.trajectory[stage_start:]
                if rows:coarse=min(rows,key=lambda a:a['objective'])['coarse_material'].copy()
                self.events.append(dict(kind='MULTISCALE_STAGE',stage=stage,parameter_dimension=B.shape[1],
                                        fixed_physics_dimension=self.model.n,seconds=self.elapsed(),
                                        has_evaluated_candidate=bool(rows),counts=dict(self.model.counts)))
            old_shape=(stage,stage)
        return self.finish()

    def finish(self):
        # A stopped correction still exposes all completed full-physics evaluations.
        finalrows=[r for r in self.trajectory if r['stage']==self.cfg['stages'][-1]]
        if finalrows:self.retain(min(finalrows,key=lambda r:r['objective']),candidate_kind='BEST_COMPLETED_FULL_PHYSICS')
        candidates=sorted(self.candidates,key=lambda r:r['objective'])
        top=candidates[0] if candidates else (min(self.trajectory,key=lambda r:r['objective']) if self.trajectory else None)
        counts=dict(self.model.counts);wall=self.elapsed()
        return json_safe(dict(schema='r4-basin-solver-v1',method=self.method,status=self.status,stop_reason=self.stop_reason,
                    config=self.cfg,wall_seconds=wall,cpu_seconds=time.process_time()-self.cpu_started,
                    budget_seconds=self.cfg['wall_seconds'],budget_overshoot_seconds=max(0.,wall-self.cfg['wall_seconds']),
                    m=None if top is None else top['m'],top1=top,candidates=candidates,archive=candidates,
                    trajectory=self.trajectory,material_trajectory=self.trajectory,raw_proposals=self.raw_proposals,
                    rays=self.rays,events=self.events,costcounts=counts,
                    accounting=dict(full_factorizations=counts.get('LU',0),full_complex_rhs=counts.get('actual_complex_RHS',0),
                                    source_rhs=counts.get('source_RHS',0),adjoint_rhs=counts.get('adjoint_RHS',0),
                                    fullj_receiver_transpose_rhs=counts.get('receiver_transpose_RHS',0),
                                    ray_kernel_rhs=counts.get('ray_kernel_RHS',0),ray_illumination_rhs=counts.get('ray_illumination_RHS',0),
                                    dense_complex_schur=counts.get('dense_complex_schur',0),counts=counts,
                                    wall_by_operation=dict(self.model.timings),cpu_by_operation=dict(self.model.cpu_timings)),
                    all_unknowns_released=bool(finalrows),full_material_dimension=self.model.n,
                    uniform_m1_reference=self.uniform_reference,threadpools=threadpool_info(),
                    ranking='full nonlinear .5||F(m)-Y||²/||Y||² + identical generic smoothing/TV',
                    truth_access='NONE_PUBLIC_INPUT_ONLY',candidate_count_is_branch_count=False))


def solve(public,method='basin_bridge',config=None,initial=None):
    with threadpool_limits(limits=1):
        run=BasinRun(load_public(public),method,config,initial)
        try:return run.run()
        except BudgetStop as exc:
            run.status='BUDGET_STOP';run.stop_reason=str(exc)
            return run.finish()


solve_case=solve
