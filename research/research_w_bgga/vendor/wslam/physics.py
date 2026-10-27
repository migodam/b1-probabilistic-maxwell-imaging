"""Campaign-neutral algebra and paid derivative actions; no geometry or truth I/O."""
import numpy as np
from skfem import asm
from .fem_backend import _forcing


def acquisition(outer=False):
    """Declared independent SLAM line-source acquisition (metres/lambda0)."""
    tx = np.array([[-1.5, 0.], [1.5, 0.]], dtype=float)
    if outer:
        tx = np.vstack((tx, [[-2.25, 0.], [2.25, 0.]]))
    rx = np.column_stack((np.linspace(-2., 2., 16), np.zeros(16)))
    return tx, rx


def pack(a, sigma=.01):
    """Frequency, source, receiver C order; all real then all imaginary."""
    if not np.isfinite(sigma) or sigma <= 0:
        raise ValueError('sigma must be finite and positive')
    a = np.asarray(a).ravel()
    if not np.isfinite(a).all():
        raise ValueError('data must be finite')
    return np.r_[a.real, a.imag]/sigma


def packj(a, sigma=.01):
    """Complex (flattened measurement, real parameter) Jacobian whitening."""
    a = np.asarray(a)
    if a.ndim != 2 or not np.isfinite(a).all():
        raise ValueError('Jacobian must be a finite matrix')
    if not np.isfinite(sigma) or sigma <= 0:
        raise ValueError('sigma must be finite and positive')
    return np.vstack((a.real, a.imag))/sigma


def _state(p, s):
    if s._backend_token is not p._backend_token:
        raise ValueError('state belongs to another backend')


def _directions(p, s, derivative):
    _state(p, s)
    d = np.asarray(derivative(s.material_points.copy(), p.freq), dtype=np.complex128)
    if d.ndim != 2 or d.shape[0] != len(s.epsilon) or d.shape[1] < 1 or not np.isfinite(d).all():
        raise ValueError('derivative(points, frequency) must return finite (nquad,nparam)')
    return d


def coefficient_jacobian(p, s, derivative):
    """Exact volume receiver chain rule, paying nr transpose RHS per frequency."""
    if s.receiver_mode != 'volume':
        raise ValueError('receiver-transpose chain rule requires volume receivers')
    d = _directions(p, s, derivative)
    nr = len(p.rx)
    ns = s.field.shape[1]
    p.book.reserve(p.phase, 'fem_general_receiver_transpose', rhs=nr,
                   factorizations=0, freq=p.freq, backend='wslam_fem_pml')
    rhs = np.empty((len(p._free_dofs), nr), dtype=np.complex128, order='F')
    for r in range(nr):
        load = (p.k**2*(s.epsilon-1)*p._receiver_green[r]).reshape(p._material_shape)
        rhs[:, r] = asm(_forcing, p.material_basis, load=load)[p._free_dofs]
    # Bilinear complex transpose here, not the real-adjoint conjugate transpose.
    v = s.factor.solve(rhs, trans='T')
    defects = np.linalg.norm(s._matrix.T@v-rhs, axis=0)/np.maximum(np.linalg.norm(rhs, axis=0), np.finfo(float).tiny)
    residual = float(defects.max())
    p._event('backend_solve', label='general_receiver_transpose', rhs=nr,
             transpose='T', residual=residual, factor_reused=True)
    if not np.isfinite(v).all() or not np.isfinite(residual) or residual > p.residual_tolerance:
        raise RuntimeError('receiver-transpose residual exceeds tolerance')
    dofs = np.zeros((p.ndofs, nr), dtype=np.complex128)
    dofs[p._free_dofs] = v
    recv = p._receiver_green.T + np.column_stack([
        np.asarray(p.material_basis.interpolate(dofs[:, r])).ravel() for r in range(nr)])
    j = np.empty((ns*nr, d.shape[1]), dtype=np.complex128)
    for source in range(ns):
        j[source*nr:(source+1)*nr] = (p.k**2*(s.material_weights*s.field[:, source])[:, None]*recv).T@d
    return j


def material_jvp_jacobian(p, s, derivative):
    """Exact J, choosing min(nparam*ns,nr) paid RHS for volume output.

    Probe extraction uses direct coefficient JVPs. Returns data Jacobian only;
    full fields/factors are never exported by this wrapper.
    """
    d = _directions(p, s, derivative)
    if s.receiver_mode == 'volume' and len(p.rx) <= d.shape[1]*s.field.shape[1]:
        return coefficient_jacobian(p, s, lambda points, freq: d)
    return np.column_stack([p.jvp(s, lambda points, v=d[:, i]: v).ravel()
                            for i in range(d.shape[1])])


def apply_target_mass(p, wall_state, target_state, total_field=None):
    """Assemble k² M_delta applied to a total field, on identical quadrature.

    No PDE action. Both states must be from the same mesh/backend and source set.
    The optional total_field is sampled at material quadrature points.
    """
    _state(p, wall_state)
    _state(p, target_state)
    if wall_state.field.shape != target_state.field.shape or not np.array_equal(wall_state.incident, target_state.incident):
        raise ValueError('paired states require identical illumination')
    field = target_state.field if total_field is None else np.asarray(total_field, dtype=np.complex128)
    if field.shape != target_state.field.shape or not np.isfinite(field).all():
        raise ValueError('total_field must match finite material-quadrature source array')
    delta = target_state.epsilon-wall_state.epsilon
    rhs = np.empty((len(p._free_dofs), field.shape[1]), dtype=np.complex128, order='F')
    for source in range(field.shape[1]):
        load = (p.k**2*delta*field[:, source]).reshape(p._material_shape)
        rhs[:, source] = asm(_forcing, p.material_basis, load=load)[p._free_dofs]
    return rhs


def paired_feedback_residual(p, wall_state, target_state):
    """Check exact wall-factor feedback identity without solving another PDE."""
    rhs = apply_target_mass(p, wall_state, target_state)
    delta_dofs = (target_state.scattered_dofs-wall_state.scattered_dofs)[p._free_dofs]
    defect = wall_state._matrix@delta_dofs-rhs
    relative = np.linalg.norm(defect, axis=0)/np.maximum(np.linalg.norm(rhs, axis=0), np.finfo(float).tiny)
    if not np.isfinite(relative).all():
        raise RuntimeError('nonfinite target feedback residual')
    return dict(relative_residuals=relative, maximum_relative_residual=float(relative.max()),
                passed=bool(relative.max() <= p.residual_tolerance),
                identity='L_wall (u_scattered_full-u_scattered_wall)=k2 M_delta E_full')


def paired_data_decomposition(wall_state, target_state):
    """Exact full paired response difference; no Born approximation or new solve."""
    if wall_state._backend_token is not target_state._backend_token:
        raise ValueError('paired states require identical backend')
    if wall_state.data.shape != target_state.data.shape or not np.array_equal(wall_state.incident, target_state.incident):
        raise ValueError('paired states require identical illumination')
    return dict(wall=wall_state.data.copy(), delta=target_state.data-wall_state.data,
                full=target_state.data.copy())


def target_feedback_terms(p, wall_state, target_state):
    """Split exact target forcing into wall illumination and repeated feedback.

    Algebra only. The returned free-DOF source arrays satisfy full=first+feedback;
    they are not receiver predictions or a standalone four-block approximation.
    """
    full = apply_target_mass(p, wall_state, target_state)
    first = apply_target_mass(p, wall_state, target_state, wall_state.field)
    return dict(full_rhs=full, first_interaction_rhs=first,
                repeated_feedback_rhs=full-first)


def point_receiver_adjoint(p, s):
    """Complex bilinear point-receiver transpose fields, one paid RHS per RX.

    Returns full global-DOF receiver fields and maximum RHS-relative residual.
    The caller may reuse these exact fields for all coefficient/shape directions.
    """
    _state(p, s)
    if s.receiver_mode != 'probe':
        raise ValueError('point receiver adjoint requires probe extraction')
    nr = len(p.rx)
    p.book.reserve(p.phase, 'fem_point_receiver_transpose', rhs=nr,
                   factorizations=0, freq=p.freq, backend='wslam_fem_pml')
    rhs = np.asfortranarray(p._receiver_probe[:, p._free_dofs].T.toarray(),
                           dtype=np.complex128)
    v = s.factor.solve(rhs, trans='T')
    defects = np.linalg.norm(s._matrix.T@v-rhs, axis=0)/np.maximum(
        np.linalg.norm(rhs, axis=0), np.finfo(float).tiny)
    residual = float(defects.max())
    p._event('backend_solve', label='point_receiver_transpose', rhs=nr,
             transpose='T', residual=residual, factor_reused=True)
    if not np.isfinite(v).all() or not np.isfinite(residual) or residual > p.residual_tolerance:
        raise RuntimeError('point receiver-transpose residual exceeds tolerance')
    qnodes = np.zeros((p.ndofs, nr), dtype=np.complex128)
    qnodes[p._free_dofs] = v
    return qnodes, residual


def _receiver_nodes(p, qnodes):
    q = np.asarray(qnodes, dtype=np.complex128)
    if q.shape != (p.ndofs, len(p.rx)) or not np.isfinite(q).all():
        raise ValueError('qnodes must be finite (global_ndof,nreceiver)')
    return q


def point_coefficient_jacobian(p, s, derivative, qnodes=None):
    """Exact point-output coefficient J; optional paid transpose-field reuse."""
    d = _directions(p, s, derivative)
    if s.receiver_mode != 'probe':
        raise ValueError('point coefficient Jacobian requires probe extraction')
    if qnodes is None:
        qnodes, _ = point_receiver_adjoint(p, s)
    qnodes = _receiver_nodes(p, qnodes)
    qquad = np.column_stack([
        np.asarray(p.material_basis.interpolate(qnodes[:, r])).ravel()
        for r in range(len(p.rx))])
    j = np.einsum('n,ns,nr,nj->srj', p.k**2*s.material_weights,
                  s.field, qquad, d, optimize=True)
    return j.reshape(s.field.shape[1]*len(p.rx), d.shape[1])


def surface_shape_jacobian(p, s, qnodes, boundary_points, normals, ds,
                           radius_weights=None, *, contrast):
    """Point-output scalar-TM interface translation/radius shape expression.

    Columns are x translation, y translation, then each supplied radius normal
    velocity. `contrast` is epsilon_inside-epsilon_outside, scalar or one value
    per boundary point. `ds` are positive boundary quadrature weights; normals
    are outward unit normals. `radius_weights` is (npoint,) or (npoint,nradius)
    containing the normal velocity for each radius parameter. No conjugation.
    This formula needs independent shape/finite-difference qualification by the
    caller for the chosen polygon geometry and boundary quadrature.
    """
    from .fem_backend import line_green
    _state(p, s)
    if s.receiver_mode != 'probe' or s.illumination['kind'] != 'line':
        raise ValueError('shape helper requires line-source point outputs')
    qnodes = _receiver_nodes(p, qnodes)
    points = np.asarray(boundary_points, dtype=float)
    normal = np.asarray(normals, dtype=float)
    weights = np.asarray(ds, dtype=float)
    if (points.ndim != 2 or points.shape[1] != 2 or len(points) == 0
            or normal.shape != points.shape or weights.shape != (len(points),)
            or not np.isfinite(points).all() or not np.isfinite(normal).all()
            or not np.isfinite(weights).all() or np.any(weights <= 0)):
        raise ValueError('finite points/normals (n,2) and positive ds (n,) required')
    if not np.allclose(np.linalg.norm(normal, axis=1), 1., rtol=1e-10, atol=1e-12):
        raise ValueError('boundary normals must be unit length')
    jump = np.asarray(contrast, dtype=np.complex128)
    if jump.ndim == 0:
        jump = np.full(len(points), jump, dtype=np.complex128)
    if jump.shape != (len(points),) or not np.isfinite(jump).all():
        raise ValueError('contrast must be finite scalar or (npoint,)')
    velocity = normal.copy()
    if radius_weights is not None:
        radius = np.asarray(radius_weights, dtype=float)
        if radius.ndim == 1:
            radius = radius[:, None]
        if radius.ndim != 2 or radius.shape[0] != len(points) or not np.isfinite(radius).all():
            raise ValueError('radius normal velocity must be finite (npoint,nradius)')
        velocity = np.column_stack((velocity, radius))
    probe = p.basis.probes(points.T)
    # Preserve the state's declared source normalization/amplitudes exactly.
    incident = line_green(points, s.illumination['positions'], p.k)*np.asarray(s.illumination['amplitudes'])[None, :]
    field = incident + probe@s.scattered_dofs
    receiver = probe@qnodes
    j = np.einsum('n,ns,nr,nj->srj', p.k**2*weights*jump,
                  field, receiver, velocity, optimize=True)
    return j.reshape(s.field.shape[1]*len(p.rx), velocity.shape[1])


def weighted_material_gradient(ev, world, raw_residual, gamma, sigma=.01):
    """Ordinary exact-data VJP, two wall + two real pose coordinates.

    The bilinear transpose load is gamma*conj(residual)/sigma². This is
    equivalent to a conjugate-adjoint formulation with its conjugations placed
    consistently. One transpose RHS per illumination, reusing the forward LU.
    It is a standard full-adjoint action, not a feedback/ROM innovation.
    """
    from .world import boundary_quadrature
    from .fem_backend import line_green
    p,s=ev.backend,ev.state
    _state(p,s)
    residual=np.asarray(raw_residual,dtype=np.complex128)
    ns=s.field.shape[1]
    if (s.receiver_mode!='probe' or not ev.metadata['target_included']
            or residual.shape!=s.data.shape or not np.isfinite(residual).all()
            or not np.isfinite(gamma) or sigma<=0):
        raise ValueError('Weighted point VJP requires a finite target-present point-output state')
    weights=gamma*np.conj(residual)/sigma**2
    rhs=np.asfortranarray(p._receiver_probe[:,p._free_dofs].T@weights.T)
    p.book.reserve(p.phase,'weighted_point_material_vjp',rhs=ns,
                   factorizations=0,freq=p.freq,forward_factor_reused=True)
    v=s.factor.solve(rhs,trans='T')
    defects=np.linalg.norm(s._matrix.T@v-rhs,axis=0)/np.maximum(np.linalg.norm(rhs,axis=0),np.finfo(float).tiny)
    backward=float(defects.max())
    p._event('backend_solve',label='weighted_point_material_vjp',rhs=ns,
             transpose='T',residual=backward,factor_reused=True)
    if not np.isfinite(v).all() or backward>p.residual_tolerance:
        raise RuntimeError('Weighted material-adjoint backward residual exceeds frozen tolerance')
    nodes=np.zeros((p.ndofs,ns),dtype=np.complex128);nodes[p._free_dofs]=v
    qquad=np.column_stack([np.asarray(p.material_basis.interpolate(nodes[:,i])).ravel() for i in range(ns)])
    n=int(ev.metadata['circle_polygon_n'])
    d=world.material_directions(s.material_points,p.freq,n,kind='shared')[:,:2]
    material=np.real(np.einsum('n,ns,ns,nj->j',p.k**2*s.material_weights,s.field,qquad,d,optimize=True))
    points,normal,ds,radial=boundary_quadrature(world.pose,world.radius,n)
    probe=p.basis.probes(points.T)
    incident=line_green(points,s.illumination['positions'],p.k)*np.asarray(s.illumination['amplitudes'])[None,:]
    field=incident+probe@s.scattered_dofs
    qboundary=probe@nodes
    jump=world.contrast+1j*world.target_loss/p.freq
    pose=np.real(np.einsum('n,ns,ns,nj->j',p.k**2*ds*jump,field,qboundary,normal,optimize=True))
    return np.r_[material,pose],dict(backward_residual=backward,paid_rhs=ns,
        full_factor_reused=True,coordinate_order=['epsinf','ell','pose_x','pose_y'],
        formula='Re sum gamma conj(raw residual) dF / sigma_real²',
        type='ORDINARY_FULL_ADJOINT_NO_ROM')
