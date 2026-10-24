"""Two-endpoint interval pCN decisions with explicit paid Full callbacks.

No random draws, physics backends, GT access, or automatic exact-cache claims.
The fixed original DecisionRefresh R0 remains unchanged after every acceptance.
"""
from dataclasses import dataclass
import numpy as np


@dataclass
class IntervalState:
    x:np.ndarray
    q:np.ndarray
    chi:np.ndarray
    point:dict
    exact:bool=False
    refinements:int=0
    endpoint:str=''

    @property
    def lower(self):return float(self.point['phi_lower'])
    @property
    def upper(self):return float(self.point['phi_upper'])


def interval_decision(current,candidate,logu):
    if not np.isfinite(logu) or logu>0:raise ValueError('fixed same logU required')
    lower=min(0.,current.lower-candidate.upper)
    upper=min(0.,current.upper-candidate.lower)
    if logu<lower:return 'ACCEPT'
    if logu>=upper:return 'REJECT'
    return 'UNKNOWN'


def approximate_state(cache,x,q,chi,observation):
    current=cache.j0+cache.R0@((chi-cache.chi0)[:,None]*cache.E0)
    cache.counts['R0_matvec_source_columns']+=current.shape[1]
    point=cache.inspect(chi,current,observation,cache.bound_constants(chi))
    return IntervalState(np.array(x,copy=True),np.array(q,copy=True),np.array(chi,copy=True),point)


def resolve_step(cache,current,candidate,observation,logu,paid_full,max_refinements=4,record=None):
    """Refine both potentials, retaining a single U; Full only via callback.

    paid_full(state,role)->(exact_state, actual_RHS) must return8 source RHS
    (or cache source count in small algebra checks). record persists each new
    approximation immediately in caller's own arm; other-arm output is unused.
    """
    if not callable(paid_full) or not 0<=max_refinements<=4:raise ValueError('explicit Full callback,max4')
    history=[];full_calls=0;sources=cache.incident.shape[1]
    for attempt in range(2*(max_refinements+1)+2):
        decision=interval_decision(current,candidate,logu)
        history.append(dict(attempt=attempt,current_lower=current.lower,current_upper=current.upper,
                            candidate_lower=candidate.lower,candidate_upper=candidate.upper,decision=decision,logu=logu))
        if decision!='UNKNOWN':return dict(current=current,candidate=candidate,accepted=decision=='ACCEPT',history=history,Full_calls=full_calls)
        # Widest uncertain endpoint first; deterministic candidate tie break.
        available=[(candidate.upper-candidate.lower,1,'candidate',candidate),
                   (current.upper-current.lower,0,'current',current)]
        available=sorted((item for item in available if not item[3].exact),reverse=True,key=lambda item:(item[0],item[1]))
        if not available:raise RuntimeError('exact endpoint decision cannot remain unknown')
        refinable=[item for item in available if item[3].refinements<max_refinements and item[3].point.get('inverse_upper') is not None]
        _,_,role,state=refinable[0] if refinable else available[0]
        if refinable:
            approx=state.point['current']+cache.R0@state.point['residual'];cache.counts['R0_matvec_source_columns']+=sources
            state.point=cache.inspect(state.chi,approx,observation,cache.bound_constants(state.chi));state.refinements+=1
            if record is not None:record(state,role,'refined')
        else:
            exact,charged=paid_full(state,role)
            if int(charged)!=sources or not exact.exact or exact.lower!=exact.upper:raise AssertionError('Full refresh must charge all source RHS and return exact potential')
            if not np.array_equal(exact.x,state.x) or not np.array_equal(exact.chi,state.chi):raise AssertionError('Full callback changed proposed state')
            full_calls+=1
            if role=='current':current=exact
            else:candidate=exact
            if record is not None:record(exact,role,'paid_Full')
    raise RuntimeError('frozen refinement limit exhausted without exact fallback resolution')
