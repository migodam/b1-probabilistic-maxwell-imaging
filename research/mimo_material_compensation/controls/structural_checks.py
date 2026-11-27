"""Reproducible small algebraic checks; none uses the user's array or truth data.

Run: python structural_checks.py
Only NumPy is required. Assertions check identities, not a Maxwell performance claim.
"""
import itertools, json
from pathlib import Path
import numpy as np

rng = np.random.default_rng(20261127)
out = {}

# Exact reciprocal full-port acquisition comparison, with the same total energy.
m = 16
Z = rng.normal(size=(m,m)) + 1j*rng.normal(size=(m,m))
Z = (Z+Z.T)/2
T = np.arange(0,m,2); R = np.arange(1,m,2)
B = Z[np.ix_(R,T)]/np.sqrt(8)
C = Z/np.sqrt(16)
simulated_B = (C[np.ix_(R,T)] + C[np.ix_(T,R)].T)/np.sqrt(2)
out['full_vs_cut'] = {
    'distance_C_over_B': float(np.linalg.norm(C)/np.linalg.norm(B)),
    'deterministic_ball_ratio_C_over_B_when_rC_is_2rB': float(np.linalg.norm(C)/(2*np.linalg.norm(B))),
    'noiseless_B_simulation_error': float(np.linalg.norm(simulated_B-B)),
    'noise_variance_after_pair_combination': 1.0,
}
assert np.linalg.norm(C) >= np.linalg.norm(B)-1e-12
assert np.linalg.norm(simulated_B-B) < 1e-12

# All 6,435 unordered balanced partitions and the optimal 15-round simplex code.
edges = list(itertools.combinations(range(m),2))
parts = [(0,)+s for s in itertools.combinations(range(1,m),7)]
counts = np.zeros(len(edges), dtype=int)
for p in parts:
    s = set(p)
    counts += [(i in s)!=(j in s) for i,j in edges]
codes = np.array([[(-1)**((a & i).bit_count()%2) for i in range(m)] for a in range(1,16)])
cover = np.array([np.sum(codes[:,i]!=codes[:,j]) for i,j in edges])
four = codes[[0,1,3,7],:]
four_cover = np.array([np.sum(four[:,i]!=four[:,j]) for i,j in edges])
Acut = np.array([i<8 for i in range(m)])
Bcut = np.array([i%2==0 for i in range(m)])
ua = {(i,j) for i,j in edges if Acut[i]!=Acut[j]}
ub = {(i,j) for i,j in edges if Bcut[i]!=Bcut[j]}
out['coverage'] = {
    'unordered_balanced_partitions': len(parts),
    'uniform_per_edge_probability': sorted(set((counts/len(parts)).tolist())),
    'simplex_rounds': len(codes),
    'simplex_transmit_counts': np.sum(codes>0,axis=1).tolist(),
    'simplex_edge_cover_counts': sorted(set(cover.tolist())),
    'simplex_gram_rank': int(np.linalg.matrix_rank(codes.T@codes)),
    'four_round_worst_fraction': float(np.min(four_cover)/4),
    'A_B_unique_complex_offdiagonal_entries': len(ua|ub),
    'A_B_intersection_edges': len(ua&ub),
}
assert len(parts)==6435 and np.all(cover==8) and np.all(np.sum(codes>0,axis=1)==8)

# Feedback LFT finite-difference identity for total scattering matrices.
def sym_contraction(n):
    z=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)); z=(z+z.T)/2
    return .55*z/np.linalg.norm(z,2)
S0=sym_contraction(4); S1=sym_contraction(4); K=sym_contraction(4)*.5
I=np.eye(4)
def feedback(S): return np.linalg.solve(I-S@K,S)
lhs=feedback(S1)-feedback(S0)
rhs=np.linalg.solve(I-S1@K,S1-S0)@np.linalg.inv(I-K@S0)
out['feedback_identity_relative_error']=float(np.linalg.norm(lhs-rhs)/np.linalg.norm(lhs))
assert out['feedback_identity_relative_error'] < 1e-12

# Rank-one perturbations can span ALL reciprocal port directions without resonance.
m=4; c=2.0; eye=np.eye(m)
vectors=[eye[:,i] for i in range(m)]+[eye[:,i]+eye[:,j] for i,j in itertools.combinations(range(m),2)]
tri=np.triu_indices(m)
cols=[]; max_rank_tail=0.; max_inverse=0.
for b in vectors:
    for t in [.3,1.1]:
        imp=c*eye+1j*t*np.outer(b,b)
        inv=np.linalg.inv(imp)
        dt=2*(eye/c-inv)
        ss=np.linalg.svd(dt,compute_uv=False)
        max_rank_tail=max(max_rank_tail,ss[1])
        max_inverse=max(max_inverse,np.linalg.norm(inv,2))
        z=dt[tri]
        cols.append(np.r_[z.real,z.imag])
rank=np.linalg.matrix_rank(np.stack(cols,axis=1))
out['rank_one_family']={
    'm':m,'real_symmetric_complex_space_dimension':m*(m+1),
    'secant_real_span_rank':int(rank),
    'largest_second_singular_value':float(max_rank_tail),
    'largest_impedance_inverse_norm':float(max_inverse),
    'uniform_inverse_bound':1/c,
    'model':'passive reciprocal lumped reactance network; not a free-space material array',
}
assert rank==m*(m+1) and max_rank_tail<1e-12 and max_inverse<=1/c+1e-12

Path('structural_checks.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
