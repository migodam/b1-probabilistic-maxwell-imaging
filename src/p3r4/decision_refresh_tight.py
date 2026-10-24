"""Optional development revision; the original DecisionRefresh stays unchanged.

Exact-arithmetic upper bounds: ||B||2 <= min(||B||F, sqrt(||B||1||B||inf)).
No factorization, inverse acquisition, solve, randomized probe or GT is used.
Computed norms use ordinary float64 reductions; NOT directed interval certified.
All cached B acquisition cost remains owned by the caller's physical ledger.
"""
import time
import numpy as np
from .decision_refresh import DecisionRefresh


def inverse_norm_upper_bounds(B,*,known_frobenius=None):
    """Two safe mathematical bounds, with explicit floating-point limitation.

    Added setup over an existing Frobenius cache is one entrywise absolute
    array and row/column sum reductions. No matrix product or LU is added.
    """
    tick=time.perf_counter();B=np.asarray(B)
    if B.ndim!=2 or not B.size or not np.isfinite(B).all():raise ValueError('finite nonempty inverse/preconditioner array required')
    absolute=np.abs(B)
    one=float(np.max(np.sum(absolute,axis=0)))
    infinity=float(np.max(np.sum(absolute,axis=1)))
    frobenius=float(np.linalg.norm(B)) if known_frobenius is None else float(known_frobenius)
    if not np.isfinite(frobenius) or frobenius<=0:raise ValueError('positive finite Frobenius norm required')
    one_infinity=float(np.sqrt(one)*np.sqrt(infinity))
    if not np.isfinite(one_infinity):one_infinity=np.inf
    upper=min(frobenius,one_infinity)
    return dict(frobenius=frobenius,one_norm=one,infinity_norm=infinity,one_infinity_upper=one_infinity,
        selected_upper=upper,selected_kind='sqrt_1_infinity' if one_infinity<frobenius else 'Frobenius',
        tightening_factor=frobenius/upper,
        scope='CONDITIONAL_EXACT_ARITHMETIC_FLOAT64_NOT_DIRECTED_INTERVAL_CERTIFIED',
        profile=dict(wall_seconds=time.perf_counter()-tick,absolute_entries=B.size,
          row_column_sum_entries=2*B.size,added_matrix_products=0,new_factorizations=0,
          new_inverse_columns=0,new_complex_RHS=0))


def neumann_constants(eta0,b_upper,delta_operator_frobenius):
    if not all(np.isfinite(v) for v in (eta0,b_upper,delta_operator_frobenius)) or eta0<0 or b_upper<=0 or delta_operator_frobenius<0:
        raise ValueError('finite nonnegative defect/perturbation and positive inverse norm required')
    eta=float(eta0+b_upper*delta_operator_frobenius)
    return dict(left_inverse_defect_F=float(eta0),inverse_operator_norm_upper=float(b_upper),
        perturbation_operator_F=float(delta_operator_frobenius),eta_upper=eta,
        inverse_upper=float(b_upper/(1-eta)) if eta<1 else None,
        sufficient_condition_valid=bool(eta<1),
        scope='CONDITIONAL_EXACT_ARITHMETIC_FLOAT64_NOT_DIRECTED_INTERVAL_CERTIFIED')


class DecisionRefreshTight(DecisionRefresh):
    """Opt-in norm-only change; identical iteration/proposal/decision mechanics.

    No automatic replacement of existing caches or production chain is made.
    Construct a new instance after a rebase or changed B; arrays must stay fixed.
    """
    def __post_init__(self):
        super().__post_init__()
        self.inverse_norm_bounds=inverse_norm_upper_bounds(self.R0,known_frobenius=self.R0F)
        self.R0_operator_upper=self.inverse_norm_bounds['selected_upper']

    def bound_constants(self,chi):
        chi=np.asarray(chi)
        if chi.shape!=self.chi0.shape or not np.isfinite(chi).all():raise ValueError('candidate material shape')
        self.counts['bound_constant_evaluations']+=1
        self.counts['bound_norm_matrix_entries']+=self.G.size
        delta=float(np.linalg.norm((chi-self.chi0)[:,None]*self.G))
        result=neumann_constants(self.left_inverse_defect_F,self.R0_operator_upper,delta)
        result.update(norm_upper_kind=self.inverse_norm_bounds['selected_kind'],legacy_R0_F=self.R0F)
        return result


def phi_interval_from_saved_residual(residual_white_norm,residual_F,sigma,S_F,constants):
    """Replay a norm bound at the SAME saved approximate state; no new state.

    This audits numerical enclosures only on saved proposals. Current-potential
    intervals and both-endpoint MH ratios remain caller responsibilities.
    """
    r=float(residual_white_norm);rho=float(residual_F)
    if not all(np.isfinite(v) for v in (r,rho,sigma,S_F)) or min(r,rho,S_F)<0 or sigma<=0:raise ValueError('invalid saved residual/norm parameters')
    inverse=constants['inverse_upper']
    if inverse is None:return dict(phi_lower=0.,phi_upper=None,error_white_upper=None,available=False)
    eps=float(np.sqrt(2)/sigma*S_F*inverse*rho)
    return dict(phi_lower=.5*max(0.,r-eps)**2,phi_upper=.5*(r+eps)**2,
        error_white_upper=eps,available=True,
        scope='CONDITIONAL_EXACT_ARITHMETIC_FLOAT64_NOT_DIRECTED_INTERVAL_CERTIFIED')
