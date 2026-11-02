# 3D vector backend audit and bounded qualification

The adapter reuses the actual A9/A17 outgoing 3D electric dipole Green, not the
scalar-TM CUDA backend. No historical source was edited, no GPU was used, and
there is no new physical backend claimed here. `src/a9_green.py` contains the
unchanged selected `dyadic_green` and `build_goff` functions from
`Gaussian/A9/code/a9_engine.py`.

## Conventions and physical material mapping

Time dependence is exp(-i omega t), background epsilon=mu=1, k=omega, Z0=1.
Lengths below use lambda0=1. Fields use cell-major/component-minor rows and
illumination columns, `(3N,P)`. Material input is physical relative permittivity
epsilon `(N,3,3)`, real or complex symmetric. Fixed mu=1. Complex symmetry does
not imply passivity; experiments must enforce their own material class.

For density j=p/v and chi=epsilon-I, the exact implemented finite model is

```
s = -1/(3v) + i k^3/(6 pi)
G = v Goff + (-1/3 + i k^3 v/(6 pi)) I
E = Einc + G j
j = chi E
(I - G blockdiag(chi)) E = Einc
Eexc = E - s p
J = -i k j
```

This is the field-native form of Clausius-Mossotti plus radiative reaction
(CM+RR). The equivalent isotropic dipole polarizability is
`alpha=3v chi/[chi+3-i k^3 v chi/(2pi)]`, exactly the A9 mapping.
For tensors, `alpha_CM=3v(epsilon-I)(epsilon+2I)^-1` and
`alpha_RR=(I-i k^3 alpha_CM/(6pi))^-1 alpha_CM`. The implemented field
form avoids contrast inversion and remains well defined at epsilon=I. It
does not confuse exciting Eexc with material internal E. It does not interpret
the radiatively corrected effective contrast as physical epsilon-1.

The static `-I/(3v)` depolarization and imaginary `i k^3/(6pi)` self term are
essential; the diagonal is not zero in the internal-field operator G. This is
the cubic/isotropic CM self prescription, with no finite-cell dynamic
depolarization correction. The continuum target is approximated by dipoles.
No exact cubical voxel volume integral, exact interface field, or continuous
Maxwell error certificate is claimed.

## Radiation, observations, and H

The outgoing `exp(+ikr)/(4pi r)` electric kernel is reused verbatim. Legal
plane waves have unit real direction d and transverse (possibly complex)
polarization a, d dot a=0. Every material state factors one CPU complex128 LU
and solves all illumination columns together; factor and matrix are retained.

External near receivers must lie outside all occupied cubic voxel cells. Both
electric and magnetic full near kernels are evaluated before adding incident
fields. The separate far operator returns amplitudes F such that
`Esc(r n)=exp(ikr) F_E(n)/r+O(r^-2)` and `F_H=n cross F_E`.
It is not a far-radius evaluation of the near receiver operator.

Magnetic field is computed using the analytic curl of the same dipole
radiation kernel: `G_H(d)=k^2 g(d)(1+i/(k r)) [n cross]`. Its symmetric
principal self contribution is zero. This is a dipole/cell-center evaluation,
not a compatible discrete curl complex. Pointwise H and L2 snapshot widths
do not qualify an H(curl) representation claim. Derivative/trace norms and
continuum Hcurl convergence remain **NOT_RUN**. In particular, any magnetic
representation experiment must preserve this qualification boundary.

## One sphere with one refinement

Frozen target: physical epsilon=3, mu=1, radius 0.12 lambda, incident k=z with
x and y polarizations. Cubic center-inclusion grids have h=0.04 and 0.03
lambda (123 and 257 active cells). An independent analytic Riccati-Bessel Mie
calculation uses outgoing h_n^(1), 12 orders, and conventional angular
recurrences. It is independent of DDA discretization, despite sharing standard
NumPy/SciPy libraries. The comparison uses absolute complex far-vector
amplitudes, including phase and scale, and a 24x32 angular quadrature of
scattering cross section.

| Metric | h=.04 | h=.03 |
|---|---:|---:|
| Absolute complex far-vector relative Mie error | 0.1046793 | 0.0365537 |
| Scattering cross-section relative Mie error | 0.2180790 | 0.0598918 |
| Discrete equation relative residual | 1.40e-15 | 2.01e-15 |
| Cell volume / true sphere volume - 1 | +8.76% | -4.13% |

The single refinement reduces both reported errors. Its nonmonotone geometric
volume bias is exposed rather than corrected by rescaling or tuning radius.
The fine grid is a modest approximate sphere control, not an exact sphere.
The material wavelength is lambda/sqrt(3)=.57735 lambda, giving about 19.2
fine-grid cells per internal wavelength; only eight cells span the sphere
diameter. Sphere-boundary curvature/staircase error remains material, as the
independent Mie error quantifies. The free-space outgoing kernel imposes
radiation directly, so no PML is present and no PML-reflection claim is needed.
No arbitrary numerical acceptance threshold was introduced. Near-to-far
amplitude differences at radius 1000 lambda are ~3.18e-4, consistent with a
finite-distance remainder; this is a convention diagnostic, not an independent
physical accuracy test. Far fields are transverse to ~1e-16 by construction.

The run spent 2 LU factorizations and 4 complex RHS. Setup, solve, total wall
and process CPU time are in `results/backend_qualification/QUALIFICATION.json`
and the two solve receipts. No failed numerical attempts occurred. Shell
inspection and Python import/startup are outside the script timer; total
process execution was separately returned by the execution tool (~0.46 s).

## Sources inspected

- `Gaussian/A9/code/a9_engine.py`, specifically original conventions,
  dyadic_green, build_goff, polarizability, and DDAState exciting/current
  mapping; `Gaussian/A9/code/sphere_reversal.py` independently documents the
  Riccati-Bessel/outgoing Mie convention.
- `Gaussian/A17/SHAPES_3D_20261030/run_shapes.py` and
  `Gaussian/A17/SHAPES_3D_BASELINES_20261030/run_comparison.py`: actual 12
  transverse-plane-wave dense vector model and cost conventions. Historical
  same-grid reconstruction results are not reused as qualification evidence.
- `research/research_z_alg_cuda_pipeline/OPTIMIZED_OPERATOR_CATALOG.md` and
  `REPRODUCE.md`: the scalar-TM resident CUDA pipeline is incompatible with
  this physical model; the catalog explicitly routes vector problems to A9.
- Attached `audit/sources/Research_MAX_STRUCT_Sources/Research_MAX_STRUCT_ZH.tex`,
  selectively lines discussing known mu, H, epsilon inverse, and fixed Gauss
  constraints. Its H-space theory does not establish a discrete Hcurl backend.
- [Draine and Goodman, finite-wavelength polarizability beyond CM](https://collaborate.princeton.edu/en/publications/beyond-clausius-mossotti-wave-propagation-on-a-polarizable-point-/),
  and [Mie coefficient algorithm documentation](https://miepython.readthedocs.io/en/3.0.0/07_algorithm.html).

## Reproduction

```
/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python \
  research/research_max_struct_r1/experiments/qualify_backend.py
```

Run from the repository root. The shell executable path contains an apostrophe
and must be appropriately quoted/escaped. The script uses two existing BLAS
threads, NumPy, SciPy and threadpoolctl; no installation is required. It writes
only this new study. Re-running overwrites its qualification artifacts, so
preserve them before a scientifically changed rerun.
