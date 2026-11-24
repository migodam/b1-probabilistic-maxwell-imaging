"""Check real profiled loss against explicit scalar LS and central differences."""

from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
if not __package__:
    sys.path.insert(0, str(HERE.parents[2]))
sys.path.insert(0, str(HERE.parent / "cross_cie"))

from research.z_alg_xudong_cie.mechanism.design import read_freeze
from research.z_alg_xudong_cie.holdout.data import write_json
from cross_engine import PixelLaw
from profile_law import ProfileLaw


def check():
    seed = read_freeze()["seeds"]["numerical_check"]
    rng = np.random.default_rng(seed)
    def sample(shape):
        return rng.standard_normal(shape) + 1j*rng.standard_normal(shape)
    n, s = 8, 5
    J, B, dJ = (sample((n, s)) for _ in range(3))
    G = sample((n,n)) / np.sqrt(2*n)
    E = B + G @ J
    dE = sample((n,s))
    law, optical = ProfileLaw(B,E,False), PixelLaw(B,E,False)
    gj, ge, value, _ = law.gradient_parts(J,E)
    j, e = J*law.weight, E*law.weight
    b = np.sum(abs(e)**2, axis=1)
    c = np.sum(j.conj()*e, axis=1)
    explicit = np.sum(.5*b/(b+law.delta)*np.sum(abs(j-(c.real/b)[:,None]*e)**2,axis=1))
    phase_gap = np.sum(.5*c.imag**2/(b+law.delta))
    closure = abs(optical.gradient_parts(J,E)[2] - value - phase_gap)
    trials = {}
    for label, dj, de in (("J_only", dJ, np.zeros_like(E)),
                          ("E_only", np.zeros_like(J), dE),
                          ("E_equals_B_plus_GJ", dJ, G@dJ)):
        delta = 1e-5
        plus = law.gradient_parts(J+delta*dj,E+delta*de)[2]
        minus = law.gradient_parts(J-delta*dj,E-delta*de)[2]
        numeric = (plus-minus)/(2*delta)
        analytic = float((np.vdot(gj,dj)+np.vdot(ge,de)).real)
        diff = abs(numeric-analytic)/max(1.,abs(numeric),abs(analytic))
        line = law.line(J,E,dj,de)
        line_error = max(abs(law.line_value(line,t)[0] - law.gradient_parts(J+t*dj,E+t*de)[2])
                         for t in (-.12,0.,.08))
        if diff > 2e-8 or line_error > 2e-12:
            raise ArithmeticError("The independent central-difference or line check failed.")
        trials[label] = dict(analytic=analytic, central_difference=numeric, relative_error=diff,
                             inherited_line_max_error=line_error)
    dark = E.copy()
    dark[0,:] = 0.
    zero_j = J.copy()
    zero_j[1,:] = 0.
    assert all(np.isfinite(x).all() for x in law.gradient_parts(zero_j,dark))
    if abs(explicit-value) > 2e-12 or closure > 2e-12:
        raise ArithmeticError("Explicit profiling or phase decomposition was not reproduced.")
    result = dict(seed=seed, explicit_LS_profile_error=abs(explicit-value), removed_phase_term=phase_gap,
                  optical_minus_profile_closure_error=closure,directions=trials,passed=True)
    write_json(HERE/"gradient_check.json",result)
    return result


if __name__ == "__main__":
    print(check())
