#!/usr/bin/env python3
"""Full-rank spectral weak planes need not be integrable: zero Maxwell calls."""
from pathlib import Path
import json
import numpy as np


def jacobian(q, alpha, eps):
    x, y, _ = q
    return np.vstack((np.sqrt(eps) * np.eye(3), [-y, -x, 1.0],
                      [0.0, np.sqrt(alpha), 0.0]))


def normal(q, alpha, eps):
    j = jacobian(q, alpha, eps)
    _, v = np.linalg.eigh(j.T @ j)
    n = v[:, -1]
    return n if n[2] >= 0 else -n


def main():
    records = []
    origin = np.zeros(3)
    for alpha in (0.2, 0.6):
        eps = 0.01
        for h in (1e-3, 1e-4, 1e-5):
            dn = np.column_stack([
                (normal(origin + h * e, alpha, eps)
                 - normal(origin - h * e, alpha, eps)) / (2*h)
                for e in np.eye(3)])
            curl = np.array([dn[2, 1] - dn[1, 2],
                             dn[0, 2] - dn[2, 0],
                             dn[1, 0] - dn[0, 1]])
            value = float(normal(origin, alpha, eps) @ curl)
            expected = -alpha / (1-alpha)
            j = jacobian(origin, alpha, eps)
            analytic_dn = np.array([[0., -1., 0.],
                                    [-1/(1-alpha), 0., 0.], [0., 0., 0.]])
            if abs(value - expected) > 2e-5:
                raise AssertionError((alpha, h, value, expected))
            records.append(dict(alpha=alpha, eps=eps, h=h,
                eigenvalues=np.linalg.eigvalsh(j.T @ j).tolist(),
                rank_J=int(np.linalg.matrix_rank(j)), Dn=dn.tolist(),
                Dn_max_error=float(np.max(np.abs(dn-analytic_dn))),
                n_dot_curl_n=value, analytic_n_dot_curl_n=expected,
                normal_bracket=[0., 0., -value],
                analytic_normal_bracket_z=-expected))
    out = dict(purpose="Non-Maxwell algebraic verification only",
               maxwell_calls=0, records=records)
    path = Path(__file__).resolve().with_name("frobenius_check.json")
    path.write_text(json.dumps(out, indent=2) + "\n")
    print("Verified full-rank nonintegrable weak planes; Maxwell calls = 0.")


if __name__ == "__main__":
    main()
