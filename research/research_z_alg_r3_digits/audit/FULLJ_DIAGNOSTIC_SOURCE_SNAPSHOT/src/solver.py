"""Public-input-only multiscale bounded Full-J and small-Gram proposal solvers."""
import json
import time
from pathlib import Path
from collections import Counter
import numpy as np
from scipy.linalg import cho_factor,cho_solve
from scipy.optimize import minimize
from threadpoolctl import threadpool_limits,threadpool_info
from .physics import Model
from .multiscale import material_basis,resample_material,Regularizer


DEFAULTS=dict(wall_seconds=75.,stages=[8,12,24],stage_fractions=[.23,.27,.50],
              smooth_weight=1e-4,tv_weight=1e-4,tv_epsilon=.01,
              max_gn_iterations=1000,initial_damping=.01,damping_min=1e-8,
              dual_fraction=.45,dual_maxiter=1000,dual_maxfun=2000,
              dual_ftol=1e-12,dual_gtol=1e-8,archive_size=32,
              min_material_rms=1e-4,lower=.5,upper=20.,initial_material=1.,
              relative_step_tolerance=1e-6,projected_gradient_tolerance=1e-8,
              common_coarse_native=True,coarse_gn_max_iterations=320)


class BudgetStop(RuntimeError):
    pass


class StageStop(RuntimeError):
    pass


def json_safe(x):
    if isinstance(x,dict):return {str(k):json_safe(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [json_safe(v) for v in x]
    if isinstance(x,np.ndarray):return json_safe(x.tolist())
    if isinstance(x,np.generic):return json_safe(x.item())
    if isinstance(x,float) and not np.isfinite(x):return None
    return x


def load_public(source,operators=None):
    allowed=('G','H','U0','S','E0','Y','y','centers','sources','receivers','radius','k',
             'grid_shape','shape','sigma','id','case_id','operators_path','box')
    if isinstance(source,(str,Path)):
        path=Path(source)
        with np.load(path,allow_pickle=False) as f:public={k:f[k] for k in allowed if k in f.files}
        if operators is None and 'operators_path' in public:
            operators=path.parent/str(public['operators_path'].item())
    else:public={k:source[k] for k in allowed if k in source}
    if operators is not None:
        if isinstance(operators,(str,Path)):
            with np.load(operators,allow_pickle=False) as f:ops={k:f[k] for k in allowed if k in f.files}
        else:ops={k:operators[k] for k in allowed if k in operators}
        ops.update(public);public=ops
    return public


class Run:
    def __init__(self,public,method,config):
        self.started,self.cpu_started=time.perf_counter(),time.process_time()
        self.cfg=dict(DEFAULTS);self.cfg.update(config or {})
        aliases={'fullj':'fullj_native','native':'fullj_native','reciprocal':'fullj_reciprocal',
                 'globaldual_smallgram':'globaldual'}
        self.method=aliases.get(method,method)
        if self.method not in ('fullj_native','fullj_reciprocal','globaldual'):
            raise ValueError('Unknown method '+method)
        self.model=Model(public)
        self.target=np.asarray(public.get('Y',public.get('y')),complex)
        if self.target.shape!=(self.model.nr,self.model.ns) or not np.all(np.isfinite(self.target)):
            raise ValueError('Finite Y shape(nr,ns) required')
        self.grid_shape=tuple(map(int,np.asarray(public.get('grid_shape',public.get('shape',[24,24])))))
        if np.prod(self.grid_shape)!=self.model.n:raise ValueError('Physical grid size mismatch')
        if self.cfg['stages'][-1]!=self.grid_shape[0] or len(set(self.grid_shape))!=1:
            raise ValueError('Final stage must release every physical material cell')
        if self.cfg['lower']!=.5 or self.cfg['upper']!=20. or self.cfg['initial_material']!=1.:
            raise ValueError('Frozen common bounds .5..20 and initial material1 required')
        self.scale2=max(float(np.vdot(self.target,self.target).real),1e-300)
        self.reg=Regularizer(self.grid_shape,self.cfg['smooth_weight'],self.cfg['tv_weight'],self.cfg['tv_epsilon'])
        self.archive=[];self.stage_best={};self.trace=[];self.events=[];self.failures=Counter()
        self.stage_deadline=np.inf;self.phase_deadline=np.inf;self.dual_factor=None;self.uniform_reference=None
        self.active_stage=None;self.status='COMPLETED_NUMERICAL';self.stop_reason=None
        self.model.guard=self.guard

    def elapsed(self):return time.perf_counter()-self.started

    def guard(self):
        elapsed=self.elapsed()
        if elapsed>=float(self.cfg['wall_seconds']):raise BudgetStop('WALL_SECONDS')
        if elapsed>=min(self.stage_deadline,self.phase_deadline):raise StageStop('STAGE_OR_PHASE_SECONDS')

    def score(self,state,origin,coarse):
        if self.uniform_reference is None and np.max(np.abs(state['m']-self.cfg['initial_material']))<1e-12:
            self.uniform_reference=dict(material=1.,mean_source_field_norm_squared=float(np.vdot(state['E'],state['E']).real)/self.model.ns,
                                        counts=dict(self.model.counts),seconds=self.elapsed())
            self.dual_factor=.5*self.uniform_reference['mean_source_field_norm_squared']/self.scale2
        R=state['y']-self.target
        rr=float(np.linalg.norm(R)/np.sqrt(self.scale2))
        reg,_,_=self.reg.evaluate(state['m'])
        score=.5*rr*rr+reg
        row=dict(m=state['m'].copy(),coarse_material=np.asarray(coarse).copy(),objective=score,ranking_score=score,
                 data_objective=.5*rr*rr,regularization=reg,relative_residual=rr,
                 residual_norm=float(np.linalg.norm(R)),origin=origin,stage=self.active_stage,
                 discovery_seconds=self.elapsed(),discovery_counts=dict(self.model.counts),
                 state_relative_residual=state['state_relative_residual'])
        if self.active_stage not in self.stage_best or score<self.stage_best[self.active_stage]['objective']:
            self.stage_best[self.active_stage]=row
        duplicate=any(np.sqrt(np.mean((row['m']-old['m'])**2))<=self.cfg['min_material_rms'] for old in self.archive)
        if not duplicate:
            self.archive.append(row)
            self.archive.sort(key=lambda a:a['objective'])
            self.archive=self.archive[:int(self.cfg['archive_size'])]
            if not self.trace or score<self.trace[-1]['objective']:
                self.trace.append(dict(objective=score,relative_residual=rr,seconds=self.elapsed(),
                                       counts=dict(self.model.counts),origin=origin,stage=self.active_stage))
        return row

    def evaluate(self,m,origin,coarse):
        state=self.model.forward(m,label=origin)
        row=self.score(state,origin,coarse)
        return state,row

    def gn(self,start,B,seconds_end,reciprocal=False,origin='GN',max_iterations=None):
        self.phase_deadline=seconds_end
        lo,hi=float(self.cfg['lower']),float(self.cfg['upper'])
        c=np.clip(np.asarray(start,float),lo,hi)
        x=1/c if reciprocal else c.copy()
        xlo,xhi=(1/hi,1/lo) if reciprocal else (lo,hi)
        identity_basis=B.shape[1]==self.model.n
        DB=None if identity_basis else self.reg.D@B
        state,row=self.evaluate(c if identity_basis else B@c,origin+'_initial',c)
        damping=float(self.cfg['initial_damping'])
        accepted=0
        iteration_limit=int(self.cfg['max_gn_iterations'] if max_iterations is None else max_iterations)
        for iteration in range(iteration_limit):
            self.guard()
            c=1/x if reciprocal else x
            deriv=-c*c if reciprocal else np.ones_like(c)
            Jfull=self.model.jacobian(state)
            with self.model.operation('GN_NORMAL_EQUATIONS',normal_equations=1,parameter_dimension=len(x)):
                J=(Jfull if identity_basis else Jfull@B)*deriv[None,:]/np.sqrt(self.scale2)
                r=self.model.pack(state['y']-self.target)/np.sqrt(self.scale2)
                _,greg,weights=self.reg.evaluate(state['m'])
                g=J.T@r+(greg if identity_basis else B.T@greg)*deriv
                Hess=J.T@J
                if identity_basis:
                    Dcsr=self.reg.D.multiply(deriv[None,:])
                    Hess+=(Dcsr.T@Dcsr.multiply(weights[:,None])).toarray()
                else:
                    D=DB*deriv[None,:]
                    Hess+=D.T@(weights[:,None]*D)
                diag=np.maximum(np.diag(Hess),max(float(np.max(np.diag(Hess)))*1e-8,1e-12))
            projected=x-np.clip(x-g,xlo,xhi)
            if np.max(np.abs(projected))<self.cfg['projected_gradient_tolerance']:
                self.events.append(dict(kind='GN_CONVERGED_PROJECTED_GRADIENT',stage=self.active_stage,iteration=iteration))
                break
            moved=False
            for attempt in range(8):
                self.guard()
                with self.model.operation('LM_NORMAL_SOLVE',normal_factorizations=1,normal_dimension=len(x)):
                    A=Hess.copy();A.flat[::len(x)+1]+=damping*diag
                    try:p=cho_solve(cho_factor(A,lower=True,check_finite=False),-g,check_finite=False)
                    except np.linalg.LinAlgError:
                        damping*=10;self.failures['LM_CHOLESKY']+=1;continue
                trialx=np.clip(x+p,xlo,xhi)
                step=trialx-x
                if np.linalg.norm(step)<=self.cfg['relative_step_tolerance']*(1+np.linalg.norm(x)):
                    damping*=10
                    if attempt==7:break
                    continue
                trialc=1/trialx if reciprocal else trialx
                try:trial,trialrow=self.evaluate(trialc if identity_basis else B@trialc,origin+'_trial',trialc)
                except (FloatingPointError,np.linalg.LinAlgError) as exc:
                    self.failures[type(exc).__name__]+=1;damping*=10;continue
                predicted=float(-g@step-.5*step@Hess@step)
                gain=(row['objective']-trialrow['objective'])/max(predicted,1e-300)
                if trialrow['objective']<row['objective']:
                    x,state,row=trialx,trial,trialrow
                    accepted+=1;moved=True
                    damping=max(float(self.cfg['damping_min']),damping*(.35 if gain>.75 else .8))
                    break
                damping=min(damping*5,1e18)
                # Failed trial changed the retained physical LU; restore only if another derivative is needed.
            self.events.append(dict(kind='GN_ITERATION',stage=self.active_stage,iteration=iteration,
                                    accepted=moved,objective=row['objective'],damping=damping,seconds=self.elapsed()))
            if not moved:
                self.events.append(dict(kind='GN_NO_ACCEPTED_STEP',stage=self.active_stage,iteration=iteration))
                break
            # Last accepted trial owns the retained LU, except rejected attempts followed by no move (break above).
        return (1/x if reciprocal else x),dict(iterations=iteration+1,accepted=accepted,iteration_limit=iteration_limit)

    def dual(self,start,B,seconds_end):
        self.phase_deadline=seconds_end
        if self.dual_factor is None:
            raise RuntimeError('GlobalDual requires paid uniform-m1 reference from common coarse initialization')
        normalization=1/self.dual_factor
        self.events.append(dict(kind='DUAL_NORMALIZATION',stage=self.active_stage,normalization=normalization,
                                reference='fixed uniform-m1 field from first common coarse forward',
                                mean_source_field_norm_squared=self.uniform_reference['mean_source_field_norm_squared']))
        calls=0
        identity_basis=B.shape[1]==self.model.n
        accepted_iterations=0
        def accepted_callback(c):
            nonlocal accepted_iterations
            accepted_iterations+=1
            self.model.counts['dual_optimizer_accepted_iterations']+=1
            self.events.append(dict(kind='DUAL_ACCEPTED_ITERATION',stage=self.active_stage,
                                    accepted_iterations=accepted_iterations,seconds=self.elapsed(),additional_RHS=0))
        def objective(c):
            nonlocal calls
            self.guard();calls+=1
            state=self.model.forward(c if identity_basis else B@c,label='dual_proposal')
            row=self.score(state,'dual_proposal',c)
            value,g,details=self.model.smallgram(state,self.target,gradient=True)
            reg,greg,_=self.reg.evaluate(state['m'])
            self.events.append(dict(kind='DUAL_EVALUATION',stage=self.active_stage,seconds=self.elapsed(),
                                    proposal_objective=value/normalization+reg,ranking_score=row['objective'],**details))
            gradient=g/normalization+greg
            return value/normalization+reg,(gradient if identity_basis else B.T@gradient)
        fit=minimize(objective,np.asarray(start),jac=True,method='L-BFGS-B',callback=accepted_callback,
                     bounds=[(self.cfg['lower'],self.cfg['upper'])]*B.shape[1],
                     options=dict(maxiter=int(self.cfg['dual_maxiter']),maxfun=int(self.cfg['dual_maxfun']),
                                  ftol=float(self.cfg['dual_ftol']),gtol=float(self.cfg['dual_gtol']),maxls=25))
        self.events.append(dict(kind='DUAL_OPTIMIZER',stage=self.active_stage,status=int(fit.status),message=str(fit.message),
                                nfev=int(fit.nfev),iterations=int(fit.nit),proposal_objective=float(fit.fun)))
        return fit.x

    def run(self):
        stages=self.cfg['stages'];fractions=np.asarray(self.cfg['stage_fractions'],float)
        fractions=fractions/np.sum(fractions)
        ends=np.cumsum(fractions)*float(self.cfg['wall_seconds'])
        coarse=np.full(int(stages[0])**2,float(self.cfg['initial_material']))
        old_shape=(int(stages[0]),int(stages[0]))
        for index,(stage,stage_end) in enumerate(zip(stages,ends)):
            self.active_stage=int(stage)
            self.stage_deadline=float(stage_end)
            self.phase_deadline=np.inf
            shape=(int(stage),int(stage))
            if index:coarse=resample_material(coarse,old_shape,shape)
            B=material_basis(self.grid_shape,stage)
            if stage==self.grid_shape[0] and not np.array_equal(B,np.eye(self.model.n)):
                raise RuntimeError('Final material basis must be identity')
            begun=self.elapsed();before=dict(self.model.counts)
            record=dict(kind='MULTISCALE_STAGE',stage=int(stage),parameter_dimension=B.shape[1],
                        fixed_physics_dimension=self.model.n,started_seconds=begun,deadline_seconds=float(stage_end),
                        initial_material=(coarse.copy() if B.shape[1]==self.model.n else B@coarse),initial_coarse_material=coarse.copy(),
                        coarse_seed_counts=dict(self.model.counts),
                        optimizer=('common_native_gn' if self.cfg['common_coarse_native'] and index<len(stages)-1 else self.method))
            # Keep a feasible stage state even when a preceding phase consumed most of its allocation.
            stage_rows_before=len(self.archive)
            try:
                if self.cfg['common_coarse_native'] and index<len(stages)-1:
                    coarse,summary=self.gn(coarse,B,float(stage_end),False,origin='common_coarse_native',
                                           max_iterations=self.cfg['coarse_gn_max_iterations'])
                elif self.method=='globaldual':
                    dual_end=begun+max(0.,float(stage_end)-begun)*float(self.cfg['dual_fraction'])
                    try:self.dual(coarse,B,dual_end)
                    except StageStop:record['dual_status']='PHASE_BUDGET_STOP'
                    self.phase_deadline=np.inf
                    # Use best original-data+common-reg proposal representable in this stage.
                    candidates=[self.stage_best[int(stage)]] if int(stage) in self.stage_best else []
                    if candidates:
                        best=min(candidates,key=lambda a:a['objective'])
                        # B is a tensor interpolation; least-squares map is tiny (p<=576), identity at final stage.
                        coarse=best['coarse_material'].copy()
                        coarse=np.clip(coarse,self.cfg['lower'],self.cfg['upper'])
                    coarse,summary=self.gn(coarse,B,float(stage_end),False,origin='dual_fullj_polish')
                else:
                    coarse,summary=self.gn(coarse,B,float(stage_end),self.method=='fullj_reciprocal',origin=self.method)
                record.update(status='OPTIMIZER_FINISHED',**summary)
            except StageStop:
                record['status']='STAGE_BUDGET_STOP'
            except BudgetStop:
                record['status']='WALL_BUDGET_STOP';self.status='BUDGET_STOP';self.stop_reason='WALL_SECONDS'
            except (ValueError,FloatingPointError,np.linalg.LinAlgError) as exc:
                self.failures[type(exc).__name__]+=1;record.update(status='FAILED_PAID',error=repr(exc))
            finally:
                candidates=[self.stage_best[int(stage)]] if int(stage) in self.stage_best else []
                if candidates:
                    best=min(candidates,key=lambda a:a['objective'])
                    coarse=best['coarse_material'].copy()
                    coarse=np.clip(coarse,self.cfg['lower'],self.cfg['upper'])
                record.update(wall_seconds=self.elapsed()-begun,
                              coarse_iteration_cap_truncated=bool(record.get('optimizer')=='common_native_gn' and record.get('status') in ('STAGE_BUDGET_STOP','WALL_BUDGET_STOP')),
                              completed_gn_iterations=sum(e['kind']=='GN_ITERATION' and e.get('stage')==int(stage) for e in self.events),
                              counts={k:v-before.get(k,0) for k,v in self.model.counts.items()},
                              has_evaluated_candidate=bool(candidates))
                self.events.append(record)
            old_shape=shape
            if self.status=='BUDGET_STOP':break
        return self.finish()

    def finish(self):
        wall=self.elapsed();counts=dict(self.model.counts)
        archive=sorted(self.archive,key=lambda a:a['objective'])
        top=archive[0] if archive else None
        stages=[e for e in self.events if e['kind']=='MULTISCALE_STAGE']
        fullstage=any(e['parameter_dimension']==self.model.n and e['has_evaluated_candidate'] for e in stages)
        return json_safe(dict(schema='r3-digits-solver-v1',method=self.method,status=self.status,
                    stop_reason=self.stop_reason,wall_seconds=wall,cpu_seconds=time.process_time()-self.cpu_started,
                    budget_seconds=self.cfg['wall_seconds'],budget_overshoot_seconds=max(0.,wall-self.cfg['wall_seconds']),
                    config=self.cfg,m=top['m'] if top else None,top1=top,candidates=archive,archive=archive,
                    ranking='full nonlinear .5||F(m)-Y||²/||Y||² + identical generic smoothing/TV',
                    full_material_dimension=self.model.n,all_unknowns_released=fullstage,
                    accounting=dict(full_factorizations=counts.get('LU',0),full_complex_rhs=counts.get('actual_complex_RHS',0),
                                    source_rhs=counts.get('source_RHS',0),adjoint_rhs=counts.get('adjoint_RHS',0),
                                    fullj_receiver_transpose_rhs=counts.get('receiver_transpose_RHS',0),
                                    fullj_builds=counts.get('FullJ',0),gram_factorizations=counts.get('Gram',0),counts=counts,
                                    wall_by_operation=dict(self.model.timings),cpu_by_operation=dict(self.model.cpu_timings)),
                    failures=dict(self.failures),uniform_m1_reference=self.uniform_reference,best_trace=self.trace,events=self.events,
                    stages=stages,threadpools=threadpool_info(),truth_access='NONE_PUBLIC_INPUT_ONLY'))


def solve_case(public,method='fullj_native',config=None,operators=None):
    with threadpool_limits(limits=1):
        run=Run(load_public(public,operators),method,config)
        try:
            return run.run()
        except BaseException as exc:
            run.status='FAILED_PAID_UNEXPECTED'
            run.stop_reason=type(exc).__name__
            run.failures[type(exc).__name__]+=1
            run.events.append(dict(kind='UNEXPECTED_FAILURE_PAID',error=repr(exc),seconds=run.elapsed()))
            return run.finish()
