"""Small algebraic VIE diagnostics; NOT a physical Maxwell imaging benchmark.

Checks a new local alternating-sweep consequence, not the already archived
CIE material-readout identity. All dimensions are tiny and no test GT is used
to select an inverse result.
"""
import json
from pathlib import Path
import numpy as np


def main():
    rng = np.random.default_rng(231)
    n, p, nr = 4, 3, 3
    A = rng.normal(size=(n, n)) + 1j*rng.normal(size=(n, n))
    G = 0.025*(A+A.T)
    E = rng.normal(size=(n, p)) + 1j*rng.normal(size=(n, p))
    S = rng.normal(size=(nr, n)) + 1j*rng.normal(size=(nr, n))
    m0 = np.array([0.7, 1.1, 1.6, 2.2])
    J0 = np.linalg.solve(np.eye(n)-m0[:, None]*G, m0[:, None]*E)
    Q0 = E+G@J0
    Y = S@J0
    beta_j = 1.0

    def exact_sweep(m, beta_m):
        w = beta_j/(1+beta_j*m)
        L = np.eye(n)-m[:, None]*G
        B = w[:, None]*L
        rhs = S.conj().T@Y + B.conj().T@(w[:, None]*m[:, None]*E)
        J = np.linalg.solve(S.conj().T@S+B.conj().T@B, rhs)
        Q = E+G@J
        a = np.sum(abs(Q)**2, axis=1)
        b = np.sum(np.real(Q.conj()*J), axis=1)
        c = np.sum(abs(J)**2, axis=1)
        return (b+beta_m*c)/(a+beta_m*b)

    def realify_matrix(Z):
        return np.block([[Z.real, -Z.imag], [Z.imag, Z.real]])

    L = np.eye(n)-m0[:, None]*G
    B = realify_matrix(np.kron(np.eye(p), L))
    D = realify_matrix(np.kron(np.eye(p), S))
    Cc = np.vstack([np.diag(Q0[:, s]) for s in range(p)])
    C = np.vstack([Cc.real, Cc.imag])
    w = beta_j/(1+beta_j*m0)
    W = np.diag(np.tile(w, 2*p))
    Haa = D.T@D+B.T@W@W@B
    Hcc = C.T@W@W@C
    Hca = C.T@W@W@B
    M = np.linalg.solve(Hcc, Hca@np.linalg.solve(Haa, Hca.T))
    sqrtH = np.diag(np.sqrt(np.diag(Hcc)))
    Minmetric = sqrtH@M@np.linalg.inv(sqrtH)
    eps = 1e-5
    beta_values = [0.2, 1.0, 4.0]
    jac_errors = {}
    for bm in beta_values:
        fd = np.column_stack([(exact_sweep(m0+eps*np.eye(n)[i], bm)-exact_sweep(m0-eps*np.eye(n)[i], bm))/(2*eps) for i in range(n)])
        jac_errors[str(bm)] = float(np.linalg.norm(fd-M))
    v = rng.normal(size=n)
    v /= np.linalg.norm(v)
    nonlinear_diff = []
    for t in [0.08, 0.04, 0.02, 0.01]:
        d = np.linalg.norm(exact_sweep(m0+t*v, 0.2)-exact_sweep(m0+t*v, 4.0))
        nonlinear_diff.append(dict(t=t, difference=float(d), difference_over_t2=float(d/t**2)))

    # Exact illustrative failure of sourcewise line-search equivariance.
    H = np.diag([1.0, 4.0])
    grad = np.eye(2)
    U = np.array([[1., 1.], [1., -1.]])/np.sqrt(2.)
    def sourcewise_step(g):
        alpha = np.sum(g*g, axis=0)/np.sum(g*(H@g), axis=0)
        return -g*alpha
    orig = sourcewise_step(grad)
    mixed_unrotated = sourcewise_step(grad@U)@U.T
    out = dict(
        scope="Algebraic reciprocal finite-VIE toy, not a Hankel/dyadic Maxwell benchmark",
        dimensions=dict(material=n, sources=p, receivers=nr),
        exact_root_error=float(np.linalg.norm(exact_sweep(m0, 1.0)-m0)),
        root_linearization_beta_material_errors=jac_errors,
        material_sweep_eigenvalues=np.linalg.eigvalsh(Minmetric).tolist(),
        metric_symmetry_error=float(np.linalg.norm(Minmetric-Minmetric.T)),
        beta_material_off_root_differences=nonlinear_diff,
        sourcewise_step=orig.tolist(),
        rotated_sourcewise_step_transported_back=mixed_unrotated.tolist(),
        sourcewise_non_equivariance=float(np.linalg.norm(orig-mixed_unrotated)),
    )
    path=Path(__file__).with_suffix('.json')
    path.write_text(json.dumps(out, indent=2), encoding='utf-8')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
