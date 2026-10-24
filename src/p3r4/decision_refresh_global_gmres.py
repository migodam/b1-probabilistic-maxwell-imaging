"""Same-paid-cache global GMRES for eight-source matrix currents.

One vector in C^(N*8), Frobenius complex inner product, max4 Arnoldi steps.
No columnwise solves, automatic Full, unpaid inverse, gradient, or GT access.
"""
from collections import OrderedDict
import time
import numpy as np
from .decision_refresh_fast import DecisionRefreshFast
from .interval_chain import interval_decision


class DecisionRefreshGlobalGMRES(DecisionRefreshFast):
    def __post_init__(self):
        super().__post_init__()
        if self.incident.shape[1]!=8:raise ValueError('all8 native source columns required')
        self._basis_cache=OrderedDict();self._reduced_solver=None

    def configure_reduced_solver(self,callback):
        if not callable(callback):raise ValueError('explicit charged reduced-LS callback required')
        self._reduced_solver=callback

    def retain_endpoint_pair(self,current,candidate):
        keys={self._key(current.chi),self._key(candidate.chi)}
        for key in list(self._basis_cache):
            if key not in keys:del self._basis_cache[key]
        if len(self._basis_cache)>2:raise AssertionError('bounded two-endpoint basis cache')

    @staticmethod
    def _key(chi):
        a=np.asarray(chi);return a.dtype.str,a.tobytes()

    def _new_basis(self,state):
        if state.refinements!=0:raise AssertionError('nonzero refinement basis must be retained, never silently restart')
        residual=state.point['residual'];z=self.R0@residual;self.counts['global_gmres_R0_source_columns']+=z.shape[1]
        beta=float(np.linalg.norm(z));V=[] if beta==0 else [z/beta]
        for v in V:v.flags.writeable=False
        initial=np.array(state.point['current'],copy=True);initial.flags.writeable=False
        self.counts['global_gmres_basis_copy_bytes']+=initial.nbytes
        return dict(initial=initial,V=V,H=np.zeros((5,4),complex),beta=beta,solutions=[],profiles=[],done=beta==0)

    def refine_state(self,state,observation):
        if state.exact or state.refinements>=4:raise ValueError('approximate endpoint with fewer than4iterations required')
        if self._reduced_solver is None:raise AssertionError('small LS must have explicit charged solver')
        key=self._key(state.chi);entry=self._basis_cache.get(key)
        if entry is None:
            if len(self._basis_cache)>=2:raise AssertionError('caller must retain onlycurrent/candidate endpoint basis')
            entry=self._new_basis(state);self._basis_cache[key]=entry
        wanted=state.refinements+1;start=time.perf_counter();hit=wanted<=len(entry['solutions'])
        if not hit and not entry['done']:
            m=len(entry['solutions']);v=entry['V'][m]
            gv=self.G@v;lv=v-state.chi[:,None]*gv;w=self.R0@lv
            self.counts['global_gmres_G_source_columns']+=v.shape[1];self.counts['global_gmres_R0_source_columns']+=v.shape[1]
            self.counts['global_gmres_operator_result_bytes']+=gv.nbytes+lv.nbytes+w.nbytes
            original_norm=float(np.linalg.norm(w))
            # Two-pass modified Gram-Schmidt uses Frobenius conjugate products.
            for _ in range(2):
                for j,basis in enumerate(entry['V'][:m+1]):
                    h=np.vdot(basis,w);entry['H'][j,m]+=h;w-=h*basis
                    self.counts['global_gmres_Frobenius_inner_products']+=1;self.counts['global_gmres_inner_product_entries']+=w.size
            next_norm=float(np.linalg.norm(w));entry['H'][m+1,m]=next_norm
            breakdown=next_norm<=1e-14*max(1.,original_norm)
            if not breakdown:
                new=w/next_norm;new.flags.writeable=False;entry['V'].append(new);self.counts['global_gmres_basis_copy_bytes']+=new.nbytes
            rows,cols=m+2,m+1;H=entry['H'][:rows,:cols].copy();rhs=np.zeros((rows,1),complex);rhs[0,0]=entry['beta']
            tick=time.perf_counter();coef,charged=self._reduced_solver(H,rhs);solve_wall=time.perf_counter()-tick
            coef=np.asarray(coef)
            if int(charged)!=1 or coef.shape!=(cols,1) or not np.isfinite(coef).all():raise AssertionError('one actual complex reducedLS RHS and finite coefficients required')
            current=entry['initial'].copy()
            for j in range(cols):current+=coef[j,0]*entry['V'][j]
            current.flags.writeable=False;entry['solutions'].append(current)
            profile=dict(iteration=cols,reduced_rows=rows,reduced_columns=cols,reduced_complex_RHS=1,reduced_solver_wall_seconds=solve_wall,reduced_solver_wall_scope='callback inclfee-booking andLS',reduced_numeric_solver_wall_seconds=getattr(self,'_last_reduced_numeric_seconds',None),
                Hessenberg_real=H.real.tolist(),Hessenberg_imag=H.imag.tolist(),coeff_real=coef[:,0].real.tolist(),coeff_imag=coef[:,0].imag.tolist(),
                preconditioned_initial_residual_F=entry['beta'],Arnoldi_next_norm=next_norm,happy_breakdown=breakdown)
            entry['profiles'].append(profile);entry['done']=breakdown or cols==4
            self.counts['global_gmres_Arnoldi_iterations']+=1;self.counts['global_gmres_reduced_LS_complex_RHS']+=1
            self.counts['global_gmres_small_solver_nanoseconds']+=int(solve_wall*1e9)
            self.counts['global_gmres_solution_copy_bytes']+=current.nbytes+H.nbytes+rhs.nbytes
        elif hit:self.counts['global_gmres_solution_cache_hits']+=1
        if entry['solutions']:
            index=min(wanted,len(entry['solutions']))-1;current=entry['solutions'][index];profile=dict(entry['profiles'][index])
        else:current=entry['initial'];profile=dict(iteration=0,reduced_complex_RHS=0,happy_breakdown=True,preconditioned_initial_residual_F=0.)
        point=self.inspect(state.chi,current,observation,self.bound_constants(state.chi));profile.update(requested_refinement=wanted,solution_cache_hit=hit)
        point['global_gmres_profile']=profile
        self.counts['global_gmres_refinement_wall_nanoseconds']+=int((time.perf_counter()-start)*1e9)
        self.counts['global_gmres_basis_cache_peak_bytes']=max(self.counts['global_gmres_basis_cache_peak_bytes'],self.basis_cache_bytes())
        return point

    def basis_cache_bytes(self):
        return sum(e['initial'].nbytes+e['H'].nbytes+sum(v.nbytes for v in e['V'])+sum(v.nbytes for v in e['solutions']) for e in self._basis_cache.values())


def resolve_global(cache,current,candidate,observation,logu,paid_full,max_refinements=4,record=None,*,fixed_accuracy=False):
    if not isinstance(cache,DecisionRefreshGlobalGMRES) or max_refinements!=4:raise ValueError('samecacheglobalGMRES max4')
    cache.retain_endpoint_pair(current,candidate);history=[];full_calls=0
    def accurate(state):
        error=state.point.get('error_white_upper')
        return state.exact or (error is not None and error<=.05 and state.upper-state.lower<=.05)
    for attempt in range(12):
        rawdecision=interval_decision(current,candidate,logu);accuracy_ok=accurate(current) and accurate(candidate)
        decision=rawdecision if not fixed_accuracy or accuracy_ok else 'UNKNOWN'
        history.append(dict(attempt=attempt,current_lower=current.lower,current_upper=current.upper,candidate_lower=candidate.lower,candidate_upper=candidate.upper,
            decision=decision,raw_MH_decision=rawdecision,fixed_accuracy=fixed_accuracy,current_accuracy=accurate(current),candidate_accuracy=accurate(candidate),same_logu=logu))
        if decision!='UNKNOWN':return dict(current=current,candidate=candidate,accepted=decision=='ACCEPT',history=history,Full_calls=full_calls)
        available=[(candidate.upper-candidate.lower,1,'candidate',candidate),(current.upper-current.lower,0,'current',current)]
        available=sorted([v for v in available if not v[3].exact],key=lambda v:(v[0],v[1]),reverse=True)
        if not available:raise AssertionError('exactendpoints must resolve')
        if fixed_accuracy and not accuracy_ok:available=[v for v in available if not accurate(v[3])]
        refinable=[v for v in available if v[3].refinements<4 and v[3].point.get('inverse_upper') is not None]
        _,_,role,state=refinable[0] if refinable else available[0]
        if refinable:
            state.point=cache.refine_state(state,observation);state.refinements+=1
            if record:record(state,role,'global_GMRES_refined')
        else:
            exact,charged=paid_full(state,role)
            if int(charged)!=cache.incident.shape[1] or not exact.exact or exact.lower!=exact.upper or not np.array_equal(exact.x,state.x) or not np.array_equal(exact.chi,state.chi):raise AssertionError('sameworld allsourceFullfallback required')
            full_calls+=1
            if role=='current':current=exact
            else:candidate=exact
            if record:record(exact,role,'paid_Full')
    raise RuntimeError('bounded iterations exhausted')
