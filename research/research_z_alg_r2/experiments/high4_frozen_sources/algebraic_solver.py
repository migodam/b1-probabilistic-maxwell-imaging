"""Observation-blind geometry setup and exact identifiable inverse."""
import numpy as np

try:
    from .physics import ReducedKernel
except ImportError:
    from physics import ReducedKernel


def reconstruct_identifiable(model, y, rank_rtol=1e-12, node_threshold=1e-20,
                             constitutive_rtol=1e-8, box=(.5, 4.)):
    """Invert receiver currents, then recover their shared real material.

    This is a deterministic algebraic inverse, not a material optimization.
    It applies only when the fixed receiver matrix has full column rank.
    A returned candidate always requires a fresh independent Maxwell forward.
    """
    y = np.asarray(y, dtype=complex)
    if y.shape != (model.nr, model.ns):
        raise ValueError('Expected y shape (receivers,sources)')
    singular = model.S_s
    rank = int(np.sum(singular > singular[0]*rank_rtol)) if singular.size else 0
    result = dict(method='RECEIVER_SVD_SHARED_REAL_CONSTITUTIVE',
                  receiver_rank=rank, receiver_columns=model.n,
                  receiver_singular_values=singular.copy(),
                  rank_relative_cutoff=rank_rtol,
                  completeness='NOT_APPLICABLE', candidate=None)
    if rank < model.n:
        result.update(status='RANK_DEFICIENT_RECEIVER_REQUIRES_BRANCH_ENUMERATION')
        return result
    with model.operation('ALGEBRAIC_IDENTIFIABLE_INVERSE',
                         receiver_pseudoinverse_applications=model.ns,
                         constitutive_material_reconstructions=1):
        j = model.S_vh[:model.n].conj().T@((model.S_u[:, :model.n].conj().T@y)/singular[:model.n, None])
        E = model.E0+model.G@j
        field_power = np.sum(np.abs(E)**2, axis=1)
        if np.any(field_power <= node_threshold):
            result.update(status='UNRESOLVED_FIELD_NODE', node_field_power=field_power)
            return result
        m = np.real(np.sum(E.conj()*j, axis=1))/field_power
        constitutive = np.linalg.norm(j-m[:, None]*E)/max(np.linalg.norm(j), 1e-300)
        receiver_residual = np.linalg.norm(model.S@j-y)/max(np.linalg.norm(y), 1e-300)
    legal = bool(np.all((m >= box[0]) & (m <= box[1])))
    result.update(candidate=m, reconstructed_current=j,
                  constitutive_relative_residual=float(constitutive),
                  receiver_relative_residual=float(receiver_residual),
                  node_field_power=field_power,
                  legal=legal,
                  status='CANDIDATE_REQUIRES_INDEPENDENT_FORWARD' if legal and constitutive <= constitutive_rtol
                  else 'REJECTED_CONSTITUTIVE_OR_BOX',
                  validation='NOT_RUN')
    return result


def fixed_active_polynomial(model, target, active, fixed_material=None, box=(.5, 4.)):
    """Build the exact square-free numerator equations on a declared slice.

    Active coordinates use x in [0,1].  Inactive coordinates are fixed by the
    caller before observing outcomes, defaulting to the declared box midpoint.
    """
    active = np.asarray(active, dtype=int)
    m0 = np.full(model.n, .5*(box[0]+box[1])) if fixed_material is None else np.asarray(fixed_material, dtype=float).copy()
    if m0.shape != (model.n,):
        raise ValueError('fixed_material must have full material shape')
    m0[active] = box[0]
    kernel = ReducedKernel(model, m0, active=active)
    poly = kernel.multiaffine_numerators(target, np.full(len(active), box[1]-box[0]))
    poly.update(kernel=kernel, box=np.asarray(box), variable='m_active=box_lower+(box_upper-box_lower)*x',
                fixed_material=m0.copy())
    return poly


def validate_candidate(model, candidate, target, box=(.5, 4.), response_rtol=1e-10,
                       state_rtol=1e-10, sigma_floor=1e-6, label='algebraic_candidate'):
    """Fresh full nonlinear solve; fitted currents are never reused."""
    state = model.forward(candidate, independent=True, label=label)
    response = float(np.linalg.norm(state['y']-target)/max(np.linalg.norm(target), 1e-300))
    legal = bool(np.all((state['m'] >= box[0]-1e-12) & (state['m'] <= box[1]+1e-12)))
    return dict(m=state['m'], response_relative=response,
                state_relative=state['state_relative_residual'],
                resolvent_sigma_min=state['resolvent_sigma_min'], legal=legal,
                qualified=bool(legal and response <= response_rtol
                               and state['state_relative_residual'] <= state_rtol
                               and state['resolvent_sigma_min'] >= sigma_floor),
                independent_validation='FRESH_NUMPY_SOLVE')
