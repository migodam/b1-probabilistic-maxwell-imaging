"""Factor ownership and exact finite-grid likelihood refresh correction."""
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class Factor:
    factor_id: str
    kind: str
    version: int
    observation_id: str | None = None
    payload: object = None


class FactorRegistry:
    def __init__(self):
        self._factors = {}
        self._observations = {}

    def register(self,factor_id,kind,version=1,observation_id=None,payload=None):
        if kind not in ('prior','transition','likelihood'):
            raise ValueError('factor type must be prior, transition, or likelihood')
        if not isinstance(version,int) or isinstance(version,bool) or version<1:
            raise ValueError('positive integer version required')
        if not factor_id:
            raise ValueError('factor ID required')
        if kind == 'likelihood' and (observation_id is None or observation_id == ''):
            raise ValueError('likelihood requires a unique observation ID')
        if kind != 'likelihood' and observation_id is not None:
            raise ValueError('only likelihood factors own observation IDs')
        old = self._factors.get(factor_id)
        if old and (version <= old.version or kind != old.kind or observation_id != old.observation_id):
            raise ValueError('replacement must preserve type/observation and strictly increase version')
        if observation_id is not None and observation_id in self._observations and self._observations[observation_id] != factor_id:
            raise ValueError('duplicate likelihood observation ID')
        f = Factor(factor_id,kind,version,observation_id,payload)
        self._factors[factor_id] = f
        if observation_id is not None:
            self._observations[observation_id] = factor_id
        return f

    def get(self,factor_id):
        return self._factors[factor_id]

    def factors(self):
        return tuple(self._factors.values())

    def __len__(self):
        return len(self._factors)


def refresh_logweights(logweights,old_loglikelihood,new_loglikelihood,normalize=True):
    """Apply log w_new = log w_old + log L_new - log L_old exactly.

    Finite likelihood values are required: zero-support changes require explicit
    reconstruction from the base measure and cannot be recovered from weights.
    """
    w,old,new = (np.asarray(v,dtype=float) for v in (logweights,old_loglikelihood,new_loglikelihood))
    if w.ndim != 1 or w.size == 0 or old.shape != w.shape or new.shape != w.shape:
        raise ValueError('matching nonempty one-dimensional finite-grid arrays required')
    if not np.all(np.isfinite(old)) or not np.all(np.isfinite(new)) or np.any(np.isnan(w)) or np.any(np.isposinf(w)):
        raise ValueError('likelihood logs must be finite; weights may contain negative infinity')
    result = w + new - old
    maximum = np.max(result)
    if not np.isfinite(maximum):
        raise ValueError('weights have empty support')
    if normalize:
        result = result-(maximum+np.log(np.exp(result-maximum).sum()))
    return result


likelihood_refresh = refresh_logweights
