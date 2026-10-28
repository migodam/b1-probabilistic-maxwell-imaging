# Independent material-jet mathematics QA

Status: **PASS_WITH_REGULAR_STRATUM_GATE** after the root's material-effective
GN correction. This audit used the original R4 P1 prior only; no Maxwell module
was imported, no physics solve ran, and no GT was loaded. Detailed reproducible
metrics are in `RESULT.json`; run `validate.py` with the existing Gaussian venv.

## Exact decoder coefficients

The convention `x(h)=x+h v+h² x2+h³ x3` uses Taylor coefficients, not raw
derivatives. The implementation has the correct factorials and cross terms.
For the quintic saturation in its active interval,

`H'(s)=15(1-s²)²/16`,
`H''(s)=15(s³-s)/4`,
`H'''(s)=15(3s²-1)/4`.

Consequently `h2=H' s2+H'' s1²/2` and
`h3=H' s3+H'' s1 s2+H''' s1³/6`, as implemented. The Gaussianized property
coefficients use `phi(x)`, `-x phi(x)`, and `(x²-1)phi(x)` with the same
coefficient convention. The complex real/imaginary material increments multiply
the real fraction coefficients by ordinary complex multiplication.

Independent degree-six interpolation of the original decoder at 17 points
within step 0.015 recovered orders 0,1,2,3 with relative errors
`5.0e-16, 8.6e-14, 3.8e-11, 4.9e-9`. Direct cubic predictions converged at
orders `3.990, 3.995, 3.998` on an active regular stratum and approximately fourth
order on the upper plateau. The lower plateau was exactly constant.

## Complex/real Hilbert handling

`real_material(a)=[Re(a);Im(a)]/24` satisfies
`real_material(a)·real_material(b)=Re(vdot(a,b))/576` exactly in the fixture.
The matrix-action layout residual was `5.7e-16`. Real least squares is correct
because response vectors and Jacobians are already real whitened coordinates.
Its normal residual was `6.2e-14`; the shared-column source-group conflict
identity defect was `1.6e-14`.

The inspected current recurrence is algebraically consistent with
`L=I-diag(chi)G`, `Lj=diag(chi)incident`, and `e=incident+Gj`:

`L j1 = m1 e`,
`L j2 = m2 e + m1 G j1`,
`L j3 = m3 e + m1 G j2 + m2 G j1`.

These coefficient recurrences contain no factorials. With `T=S L^-1`, the
forward action `T@(dchi[:,None]*e)` uses ordinary multiplication; inserting an
extra conjugation there would be wrong. The adjoint construction of T remains
the original adapter's responsibility. This is a formal convention inspection,
not a substitute for the frozen paid JVP/adjoint/jet gates.

## GN finding and verified correction

The original implementation only imposed `null(unit@M)`, which retained
`ker(M)`. This contradicted the frozen material-effective-only correction arm.
A concrete original-prior saturated point has rank 2: the old formula produced
a latent step of 0.5 entirely in the exact decoder kernel, with zero material
change. The reproducible old formula is retained as an audit fixture.

The root changed the correction space to
`Veff @ null_space(unit@M@Veff)`, using the frozen relative `1e-8` and absolute
`1e-10` thresholds. Current code returns zero gauge step for that counterexample.
On a regular point, the signed-progress tangent defect was zero; the linear
material step was 0.01 and latent step 0.0635, within both bounds. Independently
computed quadratic predicted reduction matched the return value exactly.

## Remaining scientific gate

The clipping law is C2, not globally C3. Exactly at `s=1`, inward and outward
cubic behavior differ, so there is no unique third Taylor coefficient. The
module's docstring correctly disclaims a global C3 assertion. The runner must
apply the frozen regular-stratum rule and report seam crossings when assessing
finite-radius predictions. Exact decoder retraction and actual paid physical
prediction tests remain necessary; these zero-physics checks cannot certify
their validity at radius 0.02.
