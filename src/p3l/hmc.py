"""Fixed full-dimensional HMC for an explicitly unbounded 66-D target.

The callable is the full log density in the supplied coordinates, including
priors and transformation Jacobians. There is no reflection or adaptation.
``covariance`` is the inverse mass C: K(p)=p.T C p/2 and p~N(0,C^-1).
Invocation counts describe Python calls, not underlying solver costs.
"""
import numpy as np
from scipy.linalg import solve_triangular

from .kernels import ChainState, Transition


class PreconditionedHMC:
    """Standard leapfrog followed by an endpoint flip and Hamiltonian MH.

    Finite Hamiltonian errors are never used as early stopping conditions.
    Nonfinite numerical trajectories return an explicitly unqualified failure;
    that event must be audited before treating a run as an exact reference.
    """

    def __init__(self, target, gradient, step_size, covariance, leapfrog_steps,
                 dimension=66):
        if dimension != 66:
            raise ValueError("HMC requires the full unbounded dimension=66")
        if not np.isfinite(step_size) or step_size <= 0:
            raise ValueError("step_size must be finite and positive")
        if (isinstance(leapfrog_steps, (bool, np.bool_)) or
                not isinstance(leapfrog_steps, (int, np.integer)) or
                leapfrog_steps < 1):
            raise ValueError("leapfrog_steps must be a positive integer")
        C = np.array(covariance, dtype=float, copy=True)
        if (C.shape != (66, 66) or not np.all(np.isfinite(C)) or
                not np.allclose(C, C.T, atol=1e-12, rtol=1e-12)):
            raise ValueError("covariance must be a finite symmetric SPD 66x66 matrix")
        C = .5 * (C + C.T)
        L = np.linalg.cholesky(C)
        C.setflags(write=False)
        L.setflags(write=False)
        self.target, self.gradient = target, gradient
        self.step_size, self.leapfrog_steps = float(step_size), int(leapfrog_steps)
        self.dimension, self.C, self.L = dimension, C, L

    @staticmethod
    def _gradient(value):
        g = np.array(value, dtype=float, copy=True)
        if g.shape != (66,):
            raise ValueError("HMC gradient must have exactly 66 entries")
        return g

    def step(self, current: ChainState, rng):
        nt = ng = completed = attempted = 0
        initial_g = current.gradient
        meta = dict(kernel="PRECONDITIONED_HMC", accepted=False,
                    log_mh_ratio=float("nan"), log_hamiltonian_error=float("nan"),
                    large_energy_error=False, energy_error_diagnostic_threshold=1000.,
                    rejected_proposals=1, repeated_state=True,
                    requested_leapfrog_steps=self.leapfrog_steps,
                    step_size=self.step_size)

        def finish(state, status, divergence=False):
            meta.update(target_evaluations=nt, gradient_evaluations=ng,
                        leapfrog_steps=completed, attempted_leapfrog_steps=attempted,
                        divergence=bool(divergence), divergences=int(divergence),
                        status=status)
            return Transition(state, meta)

        def retained():
            # Keep the useful current gradient even when the proposal rejects.
            if initial_g is not None and np.all(np.isfinite(initial_g)):
                return ChainState(current.q, current.log_density, initial_g)
            return current

        if initial_g is None:
            ng += 1
            initial_g = self._gradient(self.gradient(current.q.copy()))
        else:
            initial_g = self._gradient(initial_g)
        if not np.all(np.isfinite(initial_g)):
            return finish(current, "NUMERICAL_FAILURE_NONFINITE_CURRENT_GRADIENT", True)

        # If C=L L.T, L^-T xi has covariance C^-1.
        p = solve_triangular(self.L.T, rng.normal(size=66), lower=False)
        q = current.q.copy()
        with np.errstate(over="ignore", invalid="ignore"):
            kinetic0 = .5 * (p @ (self.C @ p))
            p += .5 * self.step_size * initial_g
            for i in range(self.leapfrog_steps):
                attempted += 1
                q += self.step_size * (self.C @ p)
                if not np.all(np.isfinite(q)):
                    return finish(retained(), "NUMERICAL_FAILURE_NONFINITE_POSITION", True)
                ng += 1
                g = self._gradient(self.gradient(q.copy()))
                if not np.all(np.isfinite(g)):
                    return finish(retained(), "NUMERICAL_FAILURE_NONFINITE_GRADIENT", True)
                p += (.5 if i == self.leapfrog_steps - 1 else 1.) * self.step_size * g
                completed += 1
                if not np.all(np.isfinite(p)):
                    return finish(retained(), "NUMERICAL_FAILURE_NONFINITE_MOMENTUM", True)
            p = -p  # Makes the deterministic proposal an involution.
            kinetic1 = .5 * (p @ (self.C @ p))
        nt += 1
        lp = float(self.target(q.copy()))
        if np.isnan(lp) or np.isposinf(lp) or not np.isfinite(kinetic0 + kinetic1):
            return finish(retained(), "NUMERICAL_FAILURE_NONFINITE_HAMILTONIAN", True)
        ratio = lp - current.log_density + kinetic0 - kinetic1
        if (not np.isfinite(ratio)) and not np.isneginf(lp):
            return finish(retained(), "NUMERICAL_FAILURE_NONFINITE_ENERGY_ERROR", True)
        meta.update(log_mh_ratio=float(ratio), log_hamiltonian_error=float(-ratio),
                    large_energy_error=bool(np.isfinite(ratio) and abs(ratio) > 1000.),
                    initial_kinetic_energy=float(kinetic0), final_kinetic_energy=float(kinetic1))
        accepted = bool(np.log(max(rng.random(), np.finfo(float).tiny)) < min(0., ratio))
        meta.update(accepted=accepted, rejected_proposals=int(not accepted),
                    repeated_state=not accepted or bool(np.array_equal(q, current.q)))
        state = ChainState(q, lp, g) if accepted else retained()
        status = ("EXACT_INVARIANT_KERNEL_NOT_A_MIXING_CERTIFICATE" if np.isfinite(lp)
                  else "ZERO_DENSITY_ENDPOINT_REJECTED")
        return finish(state, status)
