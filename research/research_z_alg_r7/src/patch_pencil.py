"""Coherent sources for a known uniformly loaded set of scalar VIE cells.

Column fields are complex. Response rows multiply feeds without conjugation.
Whitened feeds have Euclidean norm equal to the square root of radiated power.
"""
import numpy as np
from numpy.polynomial import polynomial as poly
from scipy.linalg import null_space
from scipy.optimize import brentq


def whitening(gram):
    eig, vec = np.linalg.eigh(gram)
    if eig[0] <= 0:
        raise ValueError("Choose an explicitly powered positive-definite feed subspace")
    return vec / np.sqrt(eig)


def chain_kernel(left0, left1, degree):
    """Return coefficient vectors of p(t) with (left0+t*left1)p(t)=0."""
    rows, cols = left0.shape
    matrix = np.zeros((rows * (degree + 2), cols * (degree + 1)), complex)
    for j in range(degree + 1):
        matrix[j*rows:(j+1)*rows, j*cols:(j+1)*cols] = left0
        matrix[(j+1)*rows:(j+2)*rows, j*cols:(j+1)*cols] = left1
    if not rows:
        return np.eye(matrix.shape[1]), np.empty(0)
    return null_space(matrix, rcond=2e-12), np.linalg.svd(matrix, compute_uv=False)


def feed_basis(coeff):
    u, singular, _ = np.linalg.svd(coeff, full_matrices=False)
    return u[:, singular > 2e-12 * singular[0]]


def nominal_mode(gain, energy):
    """Maximum squared scalar gain per positive quadratic energy."""
    coef = np.linalg.pinv(energy, hermitian=True, rcond=1e-13) @ gain.conj()
    norm = np.linalg.norm(coef)
    if norm == 0:
        raise ValueError("The calibration receiver sees no mode")
    return coef / norm


def design_pencil(H, whitened_sources, target, center, scale, nominal_receiver):
    """Shortest polynomial focused field, with nominal gain / mean feed power choice."""
    m = H.shape[0]
    E = np.eye(m)[:, target]
    image_perp = null_space(whitened_sources.conj().T).conj().T
    h, r = len(image_perp), E.shape[1]
    field0, field1 = E - center * H @ E, -scale * H @ E
    c, d = image_perp @ field0, image_perp @ field1
    preimage = np.linalg.pinv(whitened_sources)
    ranks = []
    # For r>h, dimension counting guarantees termination by this index.
    # For r<=h, one full-column-rank evaluation rules out any polynomial mode.
    if r <= h and np.linalg.matrix_rank(c) == r:
        return None
    limit = h // (r-h) if r > h else h
    for degree in range(limit + 1):
        basis, singular = chain_kernel(c, d, degree)
        ranks.append(dict(degree=degree, nullity=basis.shape[1], singular=singular))
        if basis.shape[1]:
            break
    else:
        return None
    P = basis.reshape(degree+1, r, -1)
    field = np.zeros((degree+2, m, basis.shape[1]), complex)
    for j, block in enumerate(P):
        field[j] += field0 @ block
        field[j+1] += field1 @ block
    feeds = preimage @ field
    # Uniform average over t in [-1,1] for the energy denominator.
    energy = np.zeros((basis.shape[1], basis.shape[1]), complex)
    for j, first in enumerate(feeds):
        for k, second in enumerate(feeds):
            if (j+k) % 2 == 0:
                energy += first.conj().T @ second / (j+k+1)
    coef = nominal_mode(nominal_receiver @ E @ P[0], energy)
    templates = (P @ coef).T
    feed_coeff = (feeds @ coef).T
    shots = feed_basis(feed_coeff)
    return dict(degree=degree, ranks=ranks, templates=templates, shots=shots,
                coord=shots.conj().T @ feed_coeff,
                relative_synthesis=float(np.linalg.norm(whitened_sources @ feed_coeff
                                    - (field @ coef).T) / np.linalg.norm((field @ coef).T)))


def local_anchor_sources(H, whitened_sources, target, anchors, nominal_rows):
    """Independently choose best locally nuisance-nulled calibration beam at each anchor."""
    E = np.eye(len(H))[:, target]
    off = null_space(whitened_sources.conj().T).conj().T
    preimage = np.linalg.pinv(whitened_sources)
    cols = []
    for q, receiver_row in zip(anchors, nominal_rows):
        field = E - q * H @ E
        possible = null_space(off @ field, rcond=2e-12)
        feed = preimage @ field @ possible
        mode = nominal_mode(receiver_row @ E @ possible, feed.conj().T @ feed)
        a = feed @ mode
        cols.append(a / np.linalg.norm(a))
    return feed_basis(np.column_stack(cols))


def focused_polynomial(data, R_target, design, center, scale):
    prompt = poly.polymul([center, scale], R_target @ design["templates"])
    out = data @ design["coord"]
    out[:len(prompt)] -= prompt
    return out


def _tangent_components(g, coordinates, epsilon):
    """Isolate repeated roots of the inequality for the *supplied binary coefficients*."""
    import sympy as sp
    t = sp.Symbol("t", real=True)
    def magnitude_squared(row):
        re = sum(sp.Rational(float(z.real))*t**j for j, z in enumerate(row))
        im = sum(sp.Rational(float(z.imag))*t**j for j, z in enumerate(row))
        return re**2 + im**2
    noise_term = sum(magnitude_squared(row) for row in coordinates) if epsilon else 0
    f = sp.Poly(magnitude_squared(g) - sp.Rational(float(epsilon))**2*noise_term, t)
    if f.degree() <= 0:
        return [[-1., 1.]] if f.eval(0) <= 0 else []
    boxes = [bounds for bounds, _ in f.intervals(eps=sp.Rational(1, 10**19), inf=-1, sup=1)]
    interior = [pair for pair in boxes if pair not in [(-1,-1), (1,1)]]
    edges = [(sp.Integer(-1),sp.Integer(-1)), *interior, (sp.Integer(1),sp.Integer(1))]
    mid = lambda b: float((b[0]+b[1])/2)
    pieces = [[mid(b), mid(b)] for b in boxes]
    for left, right in zip(edges, edges[1:]):
        if f.eval((left[1]+right[0])/2) < 0:
            pieces.append([mid(left), mid(right)])
    merged = []
    for left, right in sorted(pieces):
        if merged and left <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], right)
        else:
            merged.append([left,right])
    return merged


def noisy_components(g, coordinates, epsilon):
    """All real t in [-1,1] satisfying the focused bounded-error inequality (float)."""
    g = np.asarray(g, complex)
    coordinates = np.asarray(coordinates, complex)
    if epsilon < 0:
        raise ValueError("The noise radius must be nonnegative")
    if epsilon == 0:
        return _tangent_components(g, coordinates, epsilon)
    inequality = poly.polymul(g, g.conj()).real
    for row in coordinates:
        part = poly.polymul(row, row.conj()).real
        inequality = poly.polysub(inequality, epsilon**2 * part)
    def value(t):
        v = poly.polyval(t, g)
        a = poly.polyval(t, coordinates.T)
        return float(abs(v)**2 - epsilon**2 * np.vdot(a, a).real)
    derivative = poly.polyder(inequality)
    critical = [z.real for z in poly.polyroots(derivative) if abs(z.imag) < 2e-8
                and -1 < z.real < 1] if np.any(derivative) else []
    cuts = sorted(set([-1., *critical, 1.]))
    # A stationary contact can consist of one point. A sign test at the
    # floating-point stationary coordinate alone cannot distinguish that case.
    for t in cuts:
        mag = poly.polyval(abs(t), abs(g))
        feed = np.linalg.norm(poly.polyval(abs(t), abs(coordinates.T)))
        evaluated = abs(poly.polyval(t, g)) + epsilon*np.linalg.norm(poly.polyval(t, coordinates.T))
        err = 64*np.finfo(float).eps*(mag+epsilon*feed)
        if abs(value(t)) <= err*(evaluated+err):
            return _tangent_components(g, coordinates, epsilon)
    roots = [t for t in cuts if value(t) == 0]
    for l, u in zip(cuts[:-1], cuts[1:]):
        if (value(l) < 0) != (value(u) < 0) and value(l) and value(u):
            roots.append(brentq(value, l, u, xtol=2e-14))
    edges = sorted(set([-1., *roots, 1.]))
    segments = [[a, b] for a, b in zip(edges[:-1], edges[1:]) if value((a+b)/2) <= 0]
    for t in roots:
        if not any(l <= t <= u for l,u in segments):
            segments.append([t,t])
    merged = []
    for l, u in sorted(segments):
        if merged and l <= merged[-1][1]:
            merged[-1][1] = u
        else:
            merged.append([l, u])
    return merged


def profile_distance(q, y, H, R, measured_fields, target):
    """Exact distance after relaxing all off-target reciprocal currents to arbitrary complex values."""
    T = np.asarray(target)
    N = np.setdiff1d(np.arange(len(H)), T)
    K = np.linalg.solve(np.eye(len(T)) - q * H[np.ix_(T,T)], measured_fields[T])
    prompt = q * R[T] @ K
    nuisance = measured_fields[N] + q * H[np.ix_(N,T)] @ K
    null = null_space(nuisance, rcond=2e-11)
    return np.linalg.norm((y-prompt) @ null), null.shape[1]
