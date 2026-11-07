"""Cold, public-input CIE reference plus data-only analytic TRF correction.

The vendored run_cie is unchanged. Only its parent-model interface, common
data-only prior choice, finite budget and output transport are adapted here.
"""
from contextlib import contextmanager
import time
import numpy as np
from scipy.optimize import least_squares
from threadpoolctl import threadpool_limits
from research.research_z_alg_r4_diversity.src.physics import Model
from research.research_gn_nn_preservation_r4.vendor.cie_tsom_baseline import run_cie

BUDGET=dict(wall_seconds=55.,max_LU=600,max_complex_RHS=36000)
CURRENT_FRACTION=.55

class BudgetStop(RuntimeError):pass
class DataGate(RuntimeError):pass

class BudgetModel(Model):
    def __init__(self,public,owner):self.owner=owner;super().__init__(public)
    @contextmanager
    def operation(self,kind,rhs=0,**cost):
        self.owner.guard()
        if self.counts['LU']+int(cost.get('LU',0))>BUDGET['max_LU']:raise BudgetStop('MAX_LU_PROSPECTIVE')
        if self.counts['actual_complex_RHS']+int(rhs)>BUDGET['max_complex_RHS']:raise BudgetStop('MAX_RHS_PROSPECTIVE')
        with super().operation(kind,rhs=rhs,**cost):yield

class CIEModelAdapter:
    """run_cie interface over the unchanged frozen R4 physical Model."""
    def __init__(self,owner):
        self.owner=owner;self.model=owner.model
        self.G,self.E,self.S=self.model.G,self.model.E0,self.model.S
        self.shape=(24,24);self.cost=self.model.counts
    def state(self,m):
        state=self.model.forward(np.asarray(m,float),label='CIE_physical_endpoint')
        self.owner.add_candidate(state,'anchor'if not self.owner.candidates else f'cie_endpoint_{len(self.owner.candidates)}','anchor'if not self.owner.candidates else'raw')
        return dict(state,Y=state['y'],Q=state['Q'])

def frozen_config(method,Y,wall_seconds):
    if method not in ('cie_single','cie_dual'):raise ValueError('Only explicit single/dual CIE reference arms are supported')
    stages=[dict(mf=3,beta=1.,maxiter=500),dict(mf=6,beta=.2,maxiter=500),dict(mf=10,beta=.2,maxiter=500)]
    if method=='cie_dual':
        for stage in stages[1:]:stage['beta_current']=1.
    return dict(method='cie',stages=stages,dpcs_rtol=1e-3,coefficient_tol=1e-8,miniter=50,
        rho=1.,material_box=(.02,20.),regularizer=None,physical_sigma=float(np.linalg.norm(Y)),
        max_wall_seconds=max(0.,wall_seconds),sourcewise_cg=True)

class CIEReference:
    def __init__(self,public,*,method,noise,preparation_seconds=0.):
        self.started=time.perf_counter();self.cpu=time.process_time();self.preparation_seconds=float(preparation_seconds)
        self.method=method;self.Y=np.asarray(public['Y'],complex);self.scale2=squared(self.Y)
        if not np.isfinite(self.scale2)or self.scale2<=0:raise ValueError('Finite nonzero public Y required')
        self.data_gate=.001 if float(noise)==0 else .003
        self.model=BudgetModel(public,self);self.model.guard=self.guard
        self.anchor=np.full(576,.1);self.arrays=dict(anchor=self.anchor.copy())
        self.candidates=[];self.trace=[];self.attempts=[];self.current=None;self.gate_active=False
    def elapsed(self):return time.perf_counter()-self.started+self.preparation_seconds
    def guard(self):
        if self.elapsed()>=BUDGET['wall_seconds']:raise BudgetStop('WALL_SECONDS')
        if self.model.counts['LU']>BUDGET['max_LU']:raise BudgetStop('MAX_LU')
        if self.model.counts['actual_complex_RHS']>BUDGET['max_complex_RHS']:raise BudgetStop('MAX_RHS')
    def add_candidate(self,state,label,kind):
        residual=float(np.linalg.norm(state['y']-self.Y)/np.sqrt(self.scale2))
        row=dict(index=len(self.candidates),label=label,kind=kind,objective=.5*residual**2,
            relative_residual=residual,regularization=0.,old_data_TV_score=.5*residual**2,
            physical_data_feasible=bool(residual<=self.data_gate),full_VIE_equation_defect=state['state_relative_residual'],
            seconds=self.elapsed(),counts=dict(self.model.counts))
        self.candidates.append(row);self.arrays[f'candidate_{row["index"]}_m']=np.asarray(state['m']).copy()
        if label=='anchor':self.arrays['anchor_j0']=state['j'].copy()
        return row
    def correct(self,seed):
        """Same bounded data-only analytic TRF and gate as the NN corrector.

        Wrapper uses the remaining cold 55 s and 600/36000 physical caps;
        it does not inherit the R2 helper's hard-coded 45 s phase endpoint.
        """
        start=self.elapsed();begin=dict(self.model.counts);cached_x=cached_state=None
        status='FINISHED';reason=None;fit=None;best=None
        def at(x):
            nonlocal cached_x,cached_state,best
            self.guard()
            if cached_x is None or not np.array_equal(x,cached_x):
                cached_state=self.model.forward(x,label='CIE_data_TRF_trial');cached_x=np.asarray(x).copy()
                r=float(np.linalg.norm(cached_state['y']-self.Y)/np.sqrt(self.scale2))
                row=dict(objective=.5*r*r,relative_residual=r,label='CIE_data_TRF',seconds=self.elapsed(),counts=dict(self.model.counts))
                self.trace.append(row)
                if best is None or row['objective']<best[0]['objective']:best=(row,dict(cached_state))
                if r<=self.data_gate:raise DataGate('PHYSICAL_DATA_GATE')
            return cached_state
        def residual(x):return self.model.pack(at(x)['y']-self.Y)/np.sqrt(self.scale2)
        def jac(x):
            state=at(x)
            return self.model.jacobian(state)/np.sqrt(self.scale2)
        try:
            fit=least_squares(residual,np.asarray(seed,float),jac=jac,bounds=(.02,20.),method='trf',
                x_scale='jac',tr_solver='exact',max_nfev=1000,ftol=1e-12,xtol=1e-12,gtol=1e-12)
        except DataGate as e:status='PHASE_STOP';reason=str(e)
        except BudgetStop as e:status='BUDGET_STOP';reason=str(e)
        except Exception as e:status='FAILED';reason=type(e).__name__+': '+str(e)
        if best is not None:self.add_candidate(best[1],'corrected_1','corrected')
        self.attempts.append(dict(label='CIE_data_TRF',status=status,reason=reason,
            start_seconds=start,wall_seconds=self.elapsed()-start,start_counts=begin,end_counts=dict(self.model.counts),
            nfev=None if fit is None else int(fit.nfev),message=None if fit is None else str(fit.message)))
    def run(self):
        status='COMPLETED';reason=None;current_status='NOT_RUN';config=None
        try:
            config=frozen_config(self.method,self.Y,BUDGET['wall_seconds']*CURRENT_FRACTION-self.elapsed())
            config['budget_check']=self.guard
            self.current=run_cie(CIEModelAdapter(self),self.Y,self.anchor,config)
            current_status=self.current['status']
        except BudgetStop as e:status='BUDGET_STOP';reason=str(e);current_status='CAP_INTERRUPTED'
        except Exception as e:status='HOLD_CIE_IMPLEMENTATION_FAILURE';reason=type(e).__name__+': '+str(e);current_status='FAILED'
        current_wall=self.elapsed()
        if self.candidates and status!='HOLD_CIE_IMPLEMENTATION_FAILURE':
            best=min(self.candidates,key=lambda c:c['objective']);seed=self.arrays[f'candidate_{best["index"]}_m'].copy()
            # Explicit duplicate slot transports the current-phase physical incumbent;
            # no additional solve or independent exploration is claimed.
            raw=dict(best,index=len(self.candidates),label='raw_1',kind='raw',duplicate_of=best['index'])
            self.candidates.append(raw);self.arrays[f'candidate_{raw["index"]}_m']=seed.copy()
            self.arrays.update(m_proposal_unclipped=seed.copy(),m_proposal_clipped=seed.copy())
            if best['physical_data_feasible']:
                self.attempts.append(dict(label='CIE_data_TRF',status='NOT_NEEDED_DATA_GATE',reason='PHYSICAL_DATA_GATE',wall_seconds=0.))
            elif self.elapsed()<BUDGET['wall_seconds']:
                self.correct(seed)
        feasible=[c for c in self.candidates if c['physical_data_feasible']]
        chosen=min(feasible,key=lambda c:c['objective'])if feasible else None
        old=min(self.candidates,key=lambda c:c['objective'])if self.candidates else None
        if chosen is not None:self.arrays['top1']=self.arrays[f'candidate_{chosen["index"]}_m'].copy()
        if old is not None:self.arrays['old_score_top1']=self.arrays[f'candidate_{old["index"]}_m'].copy()
        self.arrays['candidate_materials']=np.asarray([self.arrays[f'candidate_{c["index"]}_m']for c in self.candidates]).reshape(-1,576)
        if 'm_proposal_unclipped'not in self.arrays:self.arrays.update(m_proposal_unclipped=self.anchor.copy(),m_proposal_clipped=self.anchor.copy())
        stops=[r for r in self.attempts if r['status']=='BUDGET_STOP']
        if stops:status='BUDGET_STOP';reason=stops[-1]['reason']
        config_receipt=None if config is None else {k:v for k,v in config.items()if k!='budget_check'}
        receipt=dict(schema='gn-nn-r4-cie-reference-v1',arm=self.method,status=status,stop_reason=reason,
            selection_status='PASS'if chosen else'FAIL_NO_FEASIBLE_CANDIDATE',accounting_status='HOLD_INCOMPLETE_CURRENT_INTERNAL_COUNTERS'if self.current is None else'COMPLETE',
            noise='clean'if self.data_gate==.001 else'noisy',data_gate=self.data_gate,truth_access='NONE_PUBLIC_INPUT_ONLY',
            wall_seconds=self.elapsed(),cpu_seconds=time.process_time()-self.cpu,budget=BUDGET,
            wall_overshoot_seconds=max(0.,self.elapsed()-BUDGET['wall_seconds']),
            preparation_cost=dict(wall_seconds=self.preparation_seconds,LU=0,actual_complex_RHS=0,complete=True),
            cost_scope='COLD_PUBLIC_UNIFORM_START_CURRENT_PLUS_DATA_TRF_55_SECONDS',
            costcounts=dict(self.model.counts),operation_seconds=dict(self.model.timings),
            current_phase=dict(status=current_status,wall_seconds=current_wall,allocation_fraction=CURRENT_FRACTION,
                config=config_receipt,cost=None if self.current is None else self.current['cost'],
                history=[]if self.current is None else self.current['history'],
                maturity='schedule completion/MAXITER/WALL_BUDGET are not convergence certificates'),
            candidates=self.candidates,top1=chosen,old_score_top1=old,trace=self.trace,attempts=self.attempts,
            auxiliary=None,control=self.method,online_selection='minimum original data-only objective among physically feasible frozen candidates; otherwise FAIL',
            clipping=dict(bounds=[.02,20.]),
            adaptations=['2D scalar-TM real material','explicit dense Fourier columns, no FFT complexity claim',
                'zero spatial prior to match NN data-only correction','same physical .001/.003 gate; remaining-budget analytic TRF',
                'cold public uniform m=.1; no GN anchor or label input'],
            provenance='unchanged vendored run_cie, HEAD 4baa5b0f88db1a49e9880c8e159e8ba336cc7da1',
            implementation_validation='EXECUTED_WITH_FROZEN_PHYSICAL_ENDPOINTS'if self.current is not None else'HOLD_RUNTIME_FAILURE')
        return receipt,self.arrays

def squared(x):return float(np.vdot(x,x).real)

def run_reference(public,*,method='cie_single',noise=0.,preparation_seconds=0.):
    """No truth arguments. Calling this function performs physical optimization."""
    with threadpool_limits(limits=1):return CIEReference(public,method=method,noise=noise,preparation_seconds=preparation_seconds).run()
