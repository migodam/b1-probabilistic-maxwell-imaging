"""Toy checks for Research Z-ALG; no Maxwell claim.

Numerical roots are matched along two nonsingular complex loops for
f(z)=z**3-3*z. Analytic statements are separately proved in the report.
"""
import json
import numpy as np
from scipy.optimize import linear_sum_assignment

base = np.array([-np.sqrt(3), 0, np.sqrt(3)], dtype=complex)
result = {}
for center in [2, -2]:
    prev = base.copy()
    worst = 0.0
    minsep = 10.0
    for theta in np.linspace(0, 2 * np.pi, 801)[1:]:
        y = center - center * np.exp(1j * theta)
        roots = np.roots([1, 0, -3, -y])
        _, cols = linear_sum_assignment(abs(prev[:, None] - roots[None, :]))
        nxt = roots[cols]
        worst = max(worst, float(np.max(abs(nxt**3 - 3 * nxt - y))))
        minsep = min(minsep, float(min(abs(nxt[i]-nxt[j]) for i in range(3) for j in range(i))))
        prev = nxt
    result[str(center)] = {
        "permutation_zero_based": [int(np.argmin(abs(base-r))) for r in prev],
        "max_polynomial_residual": worst,
        "minimum_root_separation": minsep,
        "endpoint_distance_to_original_fiber": float(max(min(abs(base-r)) for r in prev)),
        "steps": 800,
    }
print(json.dumps(result, indent=2))
