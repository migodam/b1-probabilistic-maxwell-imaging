"""GT-free real material pullbacks of a fixed complex current subspace.

No module global truth, labels, observations or endpoint arrays are loaded here.
The caller must provide a matching retained-LU state from the existing R4 model.
"""
from pathlib import Path
import time
import numpy as np
from scipy.linalg import svd, eigh, eigvalsh, lu_solve, cholesky, solve_triangular

DEFAULT_SPECTRA = Path(__file__).resolve().parents[2]/'research_gn_green_error_20261106/results/fixed_spectra.npz'

def euclidean_orth(v):
    v=np.asarray(v)
    if not v.shape[1]: return np.empty((v.shape[0],0),dtype=v.dtype)
    # Full-rank generalized eigenvectors have positive Gram; no jitter/ridge.
    g=v.T@v
    c=cholesky((g+g.T)*.5,lower=True,check_finite=False)
    return solve_triangular(c,v.T,lower=True,check_finite=False).T

def current_complement_from_cache(path=None,tau=.1,null_rtol=1e-12):
    with np.load(path or DEFAULT_SPECTRA,allow_pickle=False) as z:
        a=z['in_V'][:,z['in_s']>=tau*z['in_s'][0]]
        b=z['ob_V'][:,z['ob_s']>null_rtol*z['ob_s'][0]]
    joined=np.column_stack((a,b))
    u,s,_=svd(joined,full_matrices=False,check_finite=False,lapack_driver='gesvd')
    rank=int(np.count_nonzero(s>null_rtol*s[0]))
    return u[:,:rank]

def real_j_svd(j):
    _,s,vh=svd(j,full_matrices=j.shape[0]<j.shape[1],check_finite=False,lapack_driver='gesvd')
    ss=np.zeros(j.shape[1]);ss[:len(s)]=s
    return ss,vh.T

def prepare_geometry(model,state,*,fixed_spectra_path=None,current_complement=None,tau=.1,
                     leakage_threshold=.01,null_rtol=1e-12,compute_j_svd=True):
    """Return paid exact tangent geometry without using observations or GT.

    approx_basis/U_dark: Euclidean orthonormal columns with generalized relative
    current leakage <= threshold. leakage_vectors are HT-normalized generalized
    eigenvectors, ascending eigenvalues. J_right_vectors columns descend by J
    singular value; all n columns and zero padding are retained. numerical_kernel
    is computed from the UNSQUARED realified leakage map (relative rank 1e-12).
    """
    model._matching(state)
    begin=time.perf_counter();before=dict(model.counts)
    u=(current_complement_from_cache(fixed_spectra_path,tau,null_rtol)
       if current_complement is None else np.asarray(current_complement))
    q=state['Q']
    with model.operation('EXACT_CURRENT_TANGENT_INVERSE',rhs=model.n,current_tangent_RHS=model.n):
        w=lu_solve(model.factor,model.identity,check_finite=False)
    with model.operation('REAL_MATERIAL_CURRENT_GEOMETRY',current_geometry_builds=1):
        uw=u.conj().T@w
        # Per-source leakage rows, so shared real material is enforced directly.
        ac=np.concatenate([uw*q[:,s][None,:] for s in range(q.shape[1])],axis=0)
        ar=np.concatenate((ac.real,ac.imag),axis=0)
        ht=np.real((w.conj().T@w)*(q@q.conj().T).T)
        ht=(ht+ht.T)*.5
        hc=ar.T@ar;hc=(hc+hc.T)*.5
        vals,vec=eigh(hc,ht,check_finite=False,driver='gv')
        _,sa,vha=svd(ar,full_matrices=ar.shape[0]<model.n,check_finite=False,lapack_driver='gesvd')
        rank=int(np.count_nonzero(sa>null_rtol*sa[0])) if len(sa) and sa[0] else 0
        kernel=vha[rank:].T.copy()
        basis=euclidean_orth(vec[:,vals<=leakage_threshold])
        ht_eig=eigvalsh(ht,check_finite=False,driver='evr')
        diag=dict(current_complement_rank=u.shape[1],current_dark_dimension=model.n-u.shape[1],
                  material_approx_rank=basis.shape[1],material_sensitivity_rank=int(np.count_nonzero(vals<=.1)),
                  numerical_leakage_rank=rank,numerical_material_kernel_dimension=model.n-rank,
                  null_rtol=null_rtol,leakage_threshold=leakage_threshold,
                  leakage_sigma_max=float(sa[0]),leakage_sigma_min=float(sa[-1]),leakage_sigma_min_relative=float(sa[-1]/sa[0]),
                  min_relative_leakage=float(vals[0]),max_relative_leakage=float(vals[-1]),
                  ht_min_eigenvalue=float(ht_eig[0]),ht_condition=float(ht_eig[-1]/ht_eig[0]),
                  kernel_leakage_norm=float(np.linalg.norm(ar@kernel)),
                  approx_orthogonality_error=float(np.linalg.norm(basis.T@basis-np.eye(basis.shape[1]))),
                  generalized_eigen_residual=float(np.linalg.norm(hc@vec-(ht@vec)*vals)/max(np.linalg.norm(hc@vec),1e-300)))
    j=model.jacobian(state)
    with model.operation('MATERIAL_J_SVD',material_svd=int(compute_j_svd)):
        sj,vj=real_j_svd(j) if compute_j_svd else (None,None)
    diag['kernel_J_norm']=float(np.linalg.norm(j@kernel))
    diag['geometry_wall_seconds']=time.perf_counter()-begin
    costs={k:int(v-before.get(k,0)) for k,v in model.counts.items() if v-before.get(k,0)}
    return dict(approx_basis=basis,U_dark=basis,eigenvalues=vals,leakage_vectors=vec,
                numerical_kernel=kernel,leakage_singular_values=sa,current_complement=u,
                J=j,J_singular_values=sj,J_right_vectors=vj,HT=ht,HC=hc,
                W=w,Q=q,leakage_operator=ar,diagnostics=diag,costs=costs)
