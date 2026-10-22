# Fixed HMC Gaussian audit

`check_gaussian.py` uses only a full 66-dimensional correlated analytic Gaussian.
It checks the numerical endpoint against an independently assembled leapfrog
matrix, momentum covariance convention, Hamiltonian MH ratio and acceptance
threshold, momentum-flip reversibility, volume preservation by shear factors,
cached-gradient rejection, complete finite high-energy trajectories, explicit
nonfinite numerical failure, and all 66 stationary means and marginal variances
plus a nonzero off-diagonal covariance. The moment tolerances are fixed in the
script; results are written to `gaussian_results.json`.

Run from the P3-L project root:

```
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ../../.venv_nn/bin/python research/p3l_local_posterior_imaging/AUDIT/hmc/check_gaussian.py
```

The fixed matrix supplied as `covariance` is inverse mass C, with kinetic energy
`p.T @ C @ p / 2` and momentum covariance `C^-1`. All target coordinates must
be unbounded. The target callable must include every prior and coordinate
transformation Jacobian. Step size and mass adaptation belong outside the kernel.

Counts are callable invocations; `leapfrog_steps` counts completed position and
gradient/momentum updates, and `attempted_leapfrog_steps` also includes a failed
update. `divergence`/`divergences` currently flag nonfinite numerical failures;
all finite Hamiltonian errors are recorded and never truncate a trajectory.
`large_energy_error` flags a finite absolute Hamiltonian error greater than
1000, purely as a diagnostic, with no effect on integration or acceptance.
That threshold is not a scientific gate.
Numerical failure statuses explicitly withhold an exact-reference certificate.
The analytic Gaussian results establish neither Maxwell correctness nor actual
posterior convergence, calibration, or mixing.
