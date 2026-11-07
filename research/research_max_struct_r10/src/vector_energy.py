"""Separate vector-current and scalar-material energies; no truth/file I/O."""
import numpy as np
import torch


def scalar_material_energy(material, N, h):
    if tuple(material.shape) != (N,): raise ValueError('N real material vector required')
    return dict(domain='realN_scalar_material', integrated_squared_energy=float(h**3*np.sum(abs(material)**2)))


def current_receiver_energy(J, receiver, model):
    if J.ndim != 2 or J.shape[0] != 3*model.N:
        raise ValueError('Receiver singular basis acts only on complex3N current columns')
    coordinates=receiver.V.conj().T@J; retained=receiver.V@coordinates; complement=J-retained
    total=float(torch.sum(abs(J)**2).item()); strong=float(torch.sum(abs(retained)**2).item())
    other=float(torch.sum(abs(complement)**2).item())
    return dict(domain='complex3N_vector_current_all_sources', h3=model.cpu.h**3,
        integrated_squared_energy=model.cpu.h**3*total,
        retained_receiver_energy_fraction=strong/max(total,1e-300),
        receiver_complement_energy_fraction=other/max(total,1e-300),
        energy_closure_relative=(strong+other-total)/max(total,1e-300),
        interpretation='Finite retained receiver prefix vs its full orthogonal complement; complement is NOT a strict Green kernel')


def affine_representation_energy(J, gamma, basis, receiver, model):
    """Measure actual projected-IDFT current, not the coefficient norm proxy."""
    contribution=basis.apply(gamma); reconstruction=receiver.Jd+contribution
    removed=torch.sum(abs(gamma)**2)-torch.sum(abs(contribution)**2)
    total=torch.sum(abs(J)**2)
    return dict(domain='complex3N_vector_current',
        actual_affine_reconstruction_relative_error=float((torch.linalg.vector_norm(J-reconstruction)/torch.clamp(torch.linalg.vector_norm(J),min=1e-300)).item()),
        actual_projected_APCS_squared_norm=float(torch.sum(abs(contribution)**2).item()),
        raw_orthonormal_IDFT_squared_norm=float(torch.sum(abs(gamma)**2).item()),
        projection_removed_squared_norm=float(removed.item()),
        gamma_norm_is_projected_current_norm=False,
        affine_DPCS_APCS_cross_term=float((2*torch.real(torch.sum(receiver.Jd.conj()*contribution))).item()),
        actual_current_squared_norm=float(total.item()))


def lift_material_difference(delta_material, Q):
    """Explicit state/source-dependent diagnostic lift, not material SVD energy."""
    delta_material=np.asarray(delta_material)
    if delta_material.ndim!=1 or Q.ndim!=2 or Q.shape[0]!=3*len(delta_material):
        raise ValueError('N material difference and complex3N-by-source Q required')
    return np.repeat(delta_material,3)[:,None]*Q
