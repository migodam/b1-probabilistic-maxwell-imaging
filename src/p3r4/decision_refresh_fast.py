"""Opt-in tight-bound algebra/cache optimization, unchanged exact inequalities.

||diag(dchi)G||F²=sum_i |dchi_i|² sum_j |G_ij|². No dense candidate
perturbation array is formed. Float64 reductions are not interval-certified.
"""
from collections import OrderedDict
from types import MappingProxyType
import time
import numpy as np
from .decision_refresh_tight import DecisionRefreshTight,neumann_constants


class DecisionRefreshFast(DecisionRefreshTight):
    def __post_init__(self):
        start=time.perf_counter();copied=0
        # Immutable owned arrays make cached bounds invariant to caller mutation.
        for key in ('chi0','j0','E0','G','S','incident','R0'):
            owned=np.array(getattr(self,key),copy=True);copied+=owned.nbytes;setattr(self,key,owned)
        super().__post_init__()
        self.row_g2=np.einsum('ij,ij->i',self.G.real,self.G.real)+np.einsum('ij,ij->i',self.G.imag,self.G.imag)
        for key in ('chi0','j0','E0','G','S','incident','R0','row_g2'):getattr(self,key).flags.writeable=False
        self._bound_cache=OrderedDict()
        self.fast_setup_profile=dict(wall_seconds=time.perf_counter()-start,owned_array_copy_bytes=copied,
            row_g2_bytes=self.row_g2.nbytes,row_g2_reduction_entries=2*self.G.size,bound_cache_entries_max=2,
            added_complex_RHS=0,added_factorizations=0,added_inverse_columns=0,
            inherited_inverse_audit_dense_products=2)

    def bound_constants(self,chi):
        tick=time.perf_counter();chi=np.asarray(chi)
        if chi.shape!=self.chi0.shape or not np.isfinite(chi).all():raise ValueError('candidate material shape')
        payload=chi.tobytes();key=(chi.dtype.str,payload);self.counts['bound_cache_key_copy_bytes']+=len(payload)
        self.counts['bound_constant_requests']+=1
        if key in self._bound_cache:
            self.counts['bound_constant_cache_hits']+=1
            value=self._bound_cache.pop(key);self._bound_cache[key]=value
            self.counts['bound_constant_wall_nanoseconds']+=int((time.perf_counter()-tick)*1e9)
            return value
        delta=chi-self.chi0
        magnitude2=delta.real*delta.real+delta.imag*delta.imag
        perturbation=float(np.sqrt(np.dot(magnitude2,self.row_g2)))
        result=neumann_constants(self.left_inverse_defect_F,self.R0_operator_upper,perturbation)
        result.update(norm_upper_kind=self.inverse_norm_bounds['selected_kind'],legacy_R0_F=self.R0F,
            perturbation_norm_method='exact rowG2 algebra, ordinaryfloat64 reduction',candidate_dense_matrix_allocated=False)
        frozen=MappingProxyType(result);self._bound_cache[key]=frozen
        if len(self._bound_cache)>2:self._bound_cache.popitem(last=False)
        self.counts['bound_constant_evaluations']+=1
        self.counts['bound_norm_row_entries']+=len(delta)
        self.counts['bound_delta_vector_copy_bytes']+=delta.nbytes+magnitude2.nbytes
        self.counts['bound_constant_wall_nanoseconds']+=int((time.perf_counter()-tick)*1e9)
        return frozen

    def inspect(self,*args,**kwargs):
        point=super().inspect(*args,**kwargs)
        self.counts['inspection_saved_array_copy_bytes']+=point['current'].nbytes
        return point
