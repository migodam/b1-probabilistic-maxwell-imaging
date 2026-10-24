"""Explicit additive best-validated-query candidates, never accepted-step claims."""
import numpy as np

def ordinary_best_query(queries):
 eligible=[q for q in queries if q['prior_legal'] and q['forward_valid'] and np.isfinite(q['Phi'])]
 if not eligible:raise ValueError('no valid paid query or anchor')
 return min(eligible,key=lambda q:(q['Phi'],q['UID']))

def material_best_query(queries,mu,*,phi_threshold=282.873792):
 mu=np.asarray(mu,complex)
 eligible=[q for q in queries if q['prior_legal'] and q['forward_valid'] and np.isfinite(q['Phi']) and q['Phi']<=phi_threshold]
 if not eligible:raise ValueError('no originalnoiseballvalidquery')
 def objective(q):
  diff=np.asarray(q['chi'])-mu
  if diff.shape!=mu.shape or not np.isfinite(diff).all():raise ValueError('finiteoriginalmutargetrequired')
  return float(np.vdot(diff,diff).real/diff.size)
 return min(eligible,key=lambda q:(objective(q),q['UID']))
