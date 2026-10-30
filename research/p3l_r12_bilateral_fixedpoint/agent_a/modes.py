"""Decomposition-only bilateral sensitivity modes. No Maxwell solves or truth IO.

All parameter vectors are REAL.  Dmaterial is d(d_i)/d(x_j), complex (q,p).
Current fields: L=O(I-DC)^-1, R=(I-CD)^-1 ER.  Host fields are O,ER.
Whitening is frozen sqrt(2)/sigma source-major real32,imag32 per source.
The default output is material-RMS orthonormal, not latent-Euclidean orthonormal.
"""
from __future__ import annotations
import numpy as np
from scipy.linalg import eigh, svd

SIGMA = .00119211798514519
MATERIAL_RMS_DENOMINATOR = 576


def _real_symmetric(A):
    A = np.asarray(A).real
    return (A + A.T) * .5


def _sign_columns(B):
    B = np.array(B, dtype=float, copy=True)
    if B.size:
        piv = np.argmax(np.abs(B), axis=0)
        B *= np.where(B[piv, np.arange(B.shape[1])] < 0, -1., 1.)
    return B


def material_metric(Dmaterial, rms_denominator=MATERIAL_RMS_DENOMINATOR):
    """Pullback of squared full-grid complex-material RMS to real latent x."""
    D = np.asarray(Dmaterial, complex)
    if D.ndim != 2 or rms_denominator <= 0:
        raise ValueError('Dmaterial must be complex (cells,real parameters)')
    return _real_symmetric(D.conj().T @ D) / float(rms_denominator)


def material_canonical_basis(Dmaterial, rms_denominator=MATERIAL_RMS_DENOMINATOR,
                             gauge_rtol=1e-10, return_info=False):
    """T spans the differential material image, T.T G T=I.

    Relative cutoff applies to singular values, NOT eigenvalues.  Removed directions
    have locally zero or numerically negligible material derivative; this is only a
    local decoder quotient, not proof of a global parameterization gauge.
    """
    D = np.asarray(Dmaterial, complex)
    if D.ndim != 2 or not np.isfinite(D).all():
        raise ValueError('finite complex 2D Dmaterial required')
    B = np.concatenate((D.real, D.imag), axis=0) / np.sqrt(rms_denominator)
    _, s, vh = svd(B, full_matrices=False, check_finite=False)
    tol = max(float(gauge_rtol), np.finfo(float).eps * max(B.shape)) * (s[0] if len(s) else 0.)
    keep = s > tol
    T = _sign_columns(vh[keep].T / s[keep])
    info = dict(material_rank=int(keep.sum()), latent_dimension=int(D.shape[1]),
                removed_local_decoder_nullity=int(D.shape[1]-keep.sum()),
                material_singular_values=s.tolist(), singular_value_cutoff=float(tol),
                gauge_rtol=float(gauge_rtol), rms_denominator=float(rms_denominator),
                quotient_scope='local differential material image; no global gauge claim')
    return (T, info) if return_info else T


def bilateral_gram(L, R):
    """Complex pixel-material J*J via exact Hadamard product; no inverse/solve."""
    L, R = np.asarray(L, complex), np.asarray(R, complex)
    if L.ndim != 2 or R.ndim != 2 or L.shape[1] != R.shape[0]:
        raise ValueError('L(receiver,cells), R(cells,source) required')
    G = (L.conj().T @ L) * (R.conj() @ R.T)
    return (G + G.conj().T) * .5


def pair_modes(L, R, n_pairs=None):
    """Weighted pair row operator K for diagonal material changes.

    L=U s V*, R=P t Q*.  K[(a,b),i]=s[a]t[b]V*[a,i]P[i,b].
    Retain top products s[a]t[b].  Full K is related to vec(L diag(d) R)
    by output-unitary rotations, hence K* K is the bilateral Gram.
    """
    L, R = np.asarray(L, complex), np.asarray(R, complex)
    _, sl, vl = svd(L, full_matrices=False, check_finite=False)
    pr, sr, _ = svd(R, full_matrices=False, check_finite=False)
    strength = (sl[:, None] * sr[None, :]).ravel()
    # Stable ordering makes ties reproducible without assigning scientific meaning.
    order = np.argsort(-strength, kind='stable')
    if n_pairs is not None:
        order = order[:min(max(int(n_pairs), 0), len(order))]
    ia, ib = np.unravel_index(order, (len(sl), len(sr)))
    K = strength[order, None] * vl[ia, :] * pr[:, ib].T
    return K, dict(pair_receiver_indices=ia.tolist(), pair_source_indices=ib.tolist(),
                   pair_strengths=strength[order].tolist(), retained_pairs=len(order),
                   available_pairs=len(strength), receiver_singular_values=sl.tolist(),
                   source_singular_values=sr.tolist())


def pullback_gram(Gpixel, Dmaterial, sigma=SIGMA):
    """Real-parameter whitened Hessian: 2/sigma^2 Re(D* Gpixel D)."""
    D = np.asarray(Dmaterial, complex)
    return _real_symmetric(D.conj().T @ np.asarray(Gpixel, complex) @ D) * (2./sigma**2)


def sensitivity_hessian(Jwhite, Dmaterial, kind='full', L=None, R=None,
                        sigma=SIGMA, n_pairs=32):
    """Local data metric for requested physical ablation; no prior Hessian added."""
    k = str(kind).lower().replace('-', '_')
    prior_aware = k.endswith('_prior')
    if prior_aware:
        H, info = sensitivity_hessian(Jwhite,Dmaterial,k[:-6],L,R,sigma,n_pairs)
        return H+np.eye(H.shape[0]),dict(info,prior_hessian_added=True,prior_curvature='original Gaussianized I')
    if k in ('full', 'fullj', 'full_j', 'full_latent'):
        J = np.asarray(Jwhite)
        if np.iscomplexobj(J) or J.ndim != 2:
            raise ValueError('Jwhite must already be a REAL whitened Jacobian')
        return _real_symmetric(J.T @ J), dict(kind='FullJ', source='provided current Jwhite')
    if L is None or R is None:
        raise ValueError('cached L and R are required; this module never obtains them by solves')
    L, R = np.asarray(L, complex), np.asarray(R, complex)
    if k in ('fullgram', 'full_gram', 'bilateral', 'host', 'host_bilateral'):
        G, info = bilateral_gram(L, R), dict(kind='host_bilateral' if k in ('host', 'host_bilateral') else 'fullGram', field_scope='caller-supplied host' if k in ('host', 'host_bilateral') else 'caller-supplied current')
    elif k in ('receiver', 'receiver_only'):
        G, info = L.conj().T @ L, dict(kind='receiver_only', discarded='source illumination')
    elif k in ('source', 'source_only'):
        G, info = R.conj() @ R.T, dict(kind='source_only', discarded='receiver observability')
    elif k in ('pair', 'pair_truncated'):
        K, info = pair_modes(L, R, n_pairs=n_pairs)
        G, info = K.conj().T @ K, dict(kind='pair_truncated', **info)
    else:
        raise ValueError('unknown basis kind: '+str(kind))
    return pullback_gram(G, Dmaterial, sigma=sigma), info


def proposal_basis(Jwhite, Dmaterial, kind='full', rank=16, L=None, R=None, *,
                   sigma=SIGMA, n_pairs=32, rms_denominator=MATERIAL_RMS_DENOMINATOR,
                   gauge_rtol=1e-10, return_info=False):
    """Return weakest local directions after removing decoder material-null space.

    Default: W.T Gmaterial W=I. Coordinates therefore measure full-grid material
    RMS to first order.  W need not be Euclidean-orthogonal in latent coordinates.
    Suffix '_prior' augments the factor with T, adding the original Gaussian
    curvature I in latent coordinates; this supplies no new observational data.
    kind='full_latent' instead returns latent-Euclidean-orthonormal directions
    within the same differential material-image subspace; prior geometry differs.
    """
    k = str(kind).lower().replace('-', '_')
    prior_aware = k.endswith('_prior')
    if prior_aware:k=k[:-6]
    T, mi = material_canonical_basis(Dmaterial, rms_denominator, gauge_rtol, True)
    if k == 'full_latent' and T.shape[1]:
        T, _ = np.linalg.qr(T, mode='reduced')
    if k in ('full', 'fullj', 'full_j', 'full_latent'):
        J = np.asarray(Jwhite)
        if np.iscomplexobj(J) or J.ndim != 2 or J.shape[1] != T.shape[0]:
            raise ValueError('Jwhite must be REAL whitened Jacobian matching Dmaterial')
        A = J @ T
        info = dict(kind='FullJ', source='provided current Jwhite')
    else:
        if L is None or R is None:
            raise ValueError('cached L,R required; this module never obtains fields by solves')
        L, R = np.asarray(L, complex), np.asarray(R, complex)
        if k in ('fullgram', 'full_gram', 'bilateral', 'host', 'host_bilateral'):
            # Explicit exact factor avoids squaring condition before pullback.
            K = np.einsum('ri,is->sri', L, R, optimize=True).reshape(-1,L.shape[1])
            info = dict(kind='host_bilateral' if k in ('host','host_bilateral') else 'fullGram',
                        field_scope='caller-supplied host' if k in ('host','host_bilateral') else 'caller-supplied current')
        elif k in ('receiver','receiver_only'):
            K, info = L, dict(kind='receiver_only',discarded='source illumination')
        elif k in ('source','source_only'):
            K, info = R.T, dict(kind='source_only',discarded='receiver observability')
        elif k in ('pair','pair_truncated'):
            K, info = pair_modes(L,R,n_pairs=n_pairs)
            info = dict(kind='pair_truncated',**info)
        else:
            raise ValueError('unknown basis kind: '+str(kind))
        Acomplex = K @ (np.asarray(Dmaterial,complex) @ T)
        A = np.concatenate((Acomplex.real,Acomplex.imag),axis=0)*(np.sqrt(2)/sigma)
    if prior_aware:
        A = np.vstack((A,T))
        info = dict(info,geometry='original Gaussianized prior-aware; data information unchanged')
    if T.shape[1]:
        # Full Vh supplies true null directions when the ablation has fewer rows
        # than material coordinates. Normal equations would lose weak accuracy.
        _, singular, Vh = svd(A, full_matrices=A.shape[0]<A.shape[1], check_finite=False)
        descending = np.r_[singular**2,np.zeros(T.shape[1]-len(singular))]
        order = np.argsort(descending,kind='stable')
        eigenvalues = descending[order]
        take = min(max(int(rank),0),T.shape[1])
        W = _sign_columns(T @ Vh[order[:take]].T)
    else:
        eigenvalues = np.empty(0)
        W = np.empty((T.shape[0],0))
    info = dict(**info, **mi, selected_rank=int(W.shape[1]),
                eigenvalues=eigenvalues.tolist(), selected_eigenvalues=eigenvalues[:W.shape[1]].tolist(),
                extraction='SVD of explicit real factor after material canonicalization; no normal equations',
                output_metric='latent Euclidean' if k == 'full_latent' else 'full-grid complex-material RMS',
                physical_solve_RHS=0, physical_LU=0, prior_hessian_added=prior_aware,
                prior_curvature='original Gaussianized I' if prior_aware else 'none')
    return (W, info) if return_info else W


def weak_basis(Jwhite, Dmaterial, rank=16, **kwargs):
    return proposal_basis(Jwhite, Dmaterial, kind='full', rank=rank, **kwargs)


def coordinate_preserving_complement(W, Dmaterial, *,
                                    rms_denominator=MATERIAL_RMS_DENOMINATOR,
                                    rtol=1e-10):
    """Euclidean orthonormal Z with (W.T G_anchor) Z=0.

    Anchor G and W must remain fixed within one predictor/corrector. Recomputing
    either at each corrector would define a moving horizontal rule, not preserve
    the same finite weak coordinate. Local decoder-null directions are included;
    callers may separately restrict to their declared differential material image.
    """
    from scipy.linalg import null_space
    W = np.asarray(W, float)
    Cw = W.T @ material_metric(Dmaterial, rms_denominator)
    Z = _sign_columns(null_space(Cw, rcond=rtol))
    return Z, Cw


def proposal_frame(Jwhite, Dmaterial, kind='full', rank=16, L=None, R=None, **kwargs):
    """Convenience dict containing weak W, fixed coordinate Cw and complement Z."""
    W, info = proposal_basis(Jwhite, Dmaterial, kind, rank, L, R,
                             return_info=True, **kwargs)
    Z, Cw = coordinate_preserving_complement(
        W, Dmaterial,
        rms_denominator=kwargs.get('rms_denominator', MATERIAL_RMS_DENOMINATOR))
    return dict(W=W, Z=Z, Cw=Cw, info=info)
