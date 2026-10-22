# Actual saved reference and conditional failures

**Neither common-parent reference is qualified. The changing-decoder mobility
raised acceptance but did not establish better mixing or a speed benefit.
The Gaussian conditional pilot is severely importance-degenerate and cannot
provide actual-posterior material means or validated strong volumes.**

This is a read-only analysis of completed inference outputs, conditional and
reference code, and existing numerical audits. No GT, online observation,
anchor/geometry input, Maxwell evaluator, new chain, or new conditional samples
were accessed or run. New physics/RHS: **0**. No ledger or parent module was
modified. The parent Codex retains final scientific judgment.

Companion data: `ACTUAL_REFERENCE_AND_CONDITIONAL_FAILURES.json`; reproducible
array-only analysis: `analyze_saved_failures.py`, in this directory. All numbers
below are descriptive of the existing paid runs, not newly qualified posterior
estimates.

## 1. Saved references all fail their recorded gates

The available saved reference groups report:

| Parent / method | Material max Rhat | Material min bulk ESS | Full-grid mean RMS MCSE | Qualified |
|---|---:|---:|---:|---|
| DEV_000 / original MALA | 3.2009 | 4.4545 | 0.055655 | No |
| DEV_028 / original MALA | 5.3921 | 4.1760 | 0.105155 | No |
| DEV_030 / static retracted-init MALA | 2.2900 | 5.0010 | 0.023906 | No |
| DEV_030 / changing-decoder MALA | 4.1857 | 4.2660 | 0.023044 | No |

The recorded targets are Rhat <1.01, material-function bulk/tail ESS ≥1,000,
and full-grid mean RMS MCSE ≤0.001. Each group remains
`UNQUALIFIED_WITHIN_MECHANISM_DIAGNOSTICS`; cross-mechanism reference
verification is also absent. Original-MALA material tail ESS contains
unavailable values in DEV_028 and must remain unavailable, not zero or a pass.

For DEV_030, changing the mobility lowers the reported RMS MCSE by only about
3.6% while making worst material Rhat and minimum bulk/tail ESS worse. Neither
MCSE is reliable evidence of proximity to a true posterior mean in the presence
of this severe nonstationarity/disagreement. The two unqualified pooled material
means differ by **0.0323441 full-grid RMS**; that is a discrepancy diagnostic,
not a measured bias against truth or a trusted posterior.

## 2. Matched DEV_030 starts and actual completed-job costs

All four chain pairs have identical saved weak initialization coordinates,
identical initial native properties, and exactly equal saved initial conditional
Phi. Both methods save 1,024 warmup plus 4,096 production transitions per chain.

Costs use each completed job's `CLOSED.json` solver counters and controlled wall
time. `linear_solve_rhs` was reconciled exactly to full-forward RHS plus
receiver-transfer adjoint RHS plus likelihood-gradient adjoint RHS. Shared
previously paid anchors/backgrounds are reused, so these are incremental job
costs, not cold-start truth/anchor costs. Stage-cumulative totals were **not**
summed. This reviewer did not independently reread the authoritative ledger.

| Method | Four-job RHS | Fresh states / LU | Sum of job wall seconds |
|---|---:|---:|---:|
| Static retracted-init MALA | 329,280 | 20,512 / 20,512 | 472.489 |
| Changing-decoder MALA | 328,896 | 20,512 / 20,512 | 547.624 |

The changing-decoder runs saved **384 RHS (0.117%)**, corresponding to twelve
avoided 32-RHS warmup receiver Jacobians. They consumed **15.9% more summed
job wall time** in these runs. The fresh-state and LU counts did not fall.
Summed job wall time is not node makespan; concurrent load and timing variance
were not controlled here. These measurements establish no timing improvement
or general cost crossover. High acceptance alone cannot supply the missing
material ESS/reference-quality evidence.

## 3. Chain drift and property separation remain visible

The table uses the first and last 1,024 production rows. Native properties
are u1/u2 in [0,1], rather than the physical A/B increments. All ranges below
cover the full production portion. No thinning or rejection deduplication was
performed.

| Method / chain | Acceptance | Phi first→last quarter mean | u1 production range | u2 production range | Quarter-mean Δu1 / Δu2 |
|---|---:|---:|---|---|---|
| Static 0 | 0.4402 | 283.22→281.78 | .0611–.4660 | .2445–.4534 | −.0437 / −.0096 |
| Static 1 | 0.3279 | 291.55→288.70 | .1202–.5318 | .2892–.5525 | +.1866 / +.1196 |
| Static 2 | 0.1045 | 285.09→283.20 | .1172–.5716 | .2125–.5155 | +.1753 / +.1052 |
| Static 3 | 0.2661 | 283.34→282.45 | .6889–.8648 | .5690–.7463 | −.0209 / −.0380 |
| Changing 0 | 0.6721 | 322.93→303.58 | .0802–.2654 | .1584–.2328 | −.0625 / +.0021 |
| Changing 1 | 0.6143 | 432.33→393.70 | .1900–.4278 | .5735–.8159 | −.1531 / −.1661 |
| Changing 2 | 0.5095 | 292.74→290.99 | .5301–.7433 | .3241–.4334 | +.0656 / +.0505 |
| Changing 3 | 0.5898 | 276.40→277.85 | .8941–.9744 | .6655–.7110 | +.0213 / +.0109 |

Static chain 3's u1 range never overlaps the other static chains' ranges.
Changing chains remain even more separated in u1, with chain 0 below .266 and
chain 3 above .894. Changing chain 1's Phi drops about 38.6 between production
quarter means. These are concrete failure diagnostics; quarter means are not
independent replications and no significance test is implied.

The saved arrays passed two inexpensive integrity checks: every rejected
production proposal retained both q and x exactly, and all saved transformed
log densities satisfy
`log_density + Phi + .5*sum(x*x) = 0` to at most about `3.6e-14`.
Thus the observed failure is not explained by dropped rejection repetitions
or an omitted probit Gaussian prior in the saved target values.

## 4. Forty Gaussian conditional batches: strong weight collapse

The B pilot has four development parents, five weak Gaussian dimensions and
two fixed weak contexts, with the native properties explicitly fixed at their
anchors. This is **not** an inferred distribution over the outer (a,u).
Each batch draws a 128-point Gaussian strong reference, then examines the same
8/32/128 prefixes with exact target weights.

| Weak d / strong s | K128 ESS min / median / max | Median Gaussian mean−mode RMS | Median SNIS mean−Gaussian RMS |
|---|---|---:|---:|
| 2 / 62 | 1.000 / 1.000 / 1.002 | .09878 | .11958 |
| 4 / 60 | 1.000 / 1.000 / 1.066 | .09477 | .11861 |
| 8 / 56 | 1.000 / 1.005 / 1.475 | .07965 | .09823 |
| 16 / 48 | 1.002 / 1.267 / 1.919 | .04921 | .05470 |
| 24 / 40 | 1.000 / 1.029 / 5.269 | .02640 | .03574 |

**38/40** batches have ESS <2; **27/40** have a maximum normalized weight >.99.
Even the best batch has only 5.269 effective weights out of 128. These are
empirical weight-concentration diagnostics, not independent posterior samples
or reliable teacher moments.

At K8, median material-mean distance to K128 is .07245 RMS; at K32 it is
.04656 RMS, with maxima .16672 and .13425. The median absolute log-normalizer
change from K8 to K128 is 134.20, and from K32 to K128 is 2.606; the respective
maxima are 4,512.58 and 564.88. A K32 material difference of zero can simply mean
that the same single weight dominates both prefixes. The K128 difference to
itself is identically zero and provides no convergence evidence.

Larger outer dimension reduces the median Gaussian mean−mode discrepancy and
the median SNIS−Gaussian discrepancy in this finite pilot, as the strong block
shrinks. ESS does **not** improve monotonically: its median is lower at d24 than
at d16, and every dimension remains inadequate. Frozen frames/outer points
change with dimension; this is a descriptive trend, not a controlled proof that
dimension alone causes the failures or that d24 is a validated choice.

## 5. Mode convergence and GN volume cannot be treated as exact marginal evidence

**0/40** finite conditional mode candidates meet the saved convergence rule
`strong_gradient_inf < 1e-4`. Zero-context gradients have median .06941,
range .0003735–1.8395; shifted contexts have median 99.741,
range .8774–30,726.17. Some zero contexts are near stationary, so this does not
mean all mode errors have comparable magnitude. The capped shifted searches
are often plainly unfinished.

The implementation uses

\[
H_{GN}=J_b^TJ_b+I,\qquad
\log Z_{GN}=-U(b_0)+\frac s2\log(2\pi)-\frac12\log\det H_{GN}.
\]

At a nonstationary candidate b0, a local quadratic retaining its linear term
would instead have

\[
U(b_0+\delta)\approx U(b_0)+g^T\delta+\tfrac12\delta^TH_{GN}\delta,
\]

\[
\log Z_{quad}=-U(b_0)+\tfrac12g^TH_{GN}^{-1}g
+\frac s2\log(2\pi)-\tfrac12\log\det H_{GN},
\quad E_{quad}b=b_0-H_{GN}^{-1}g.
\]

Thus `conditional.py:270–278` defines a legitimate normalized Gaussian
**proposal/reference** centered at the finite candidate, but its unshifted
`logZ_GN` is not the integral of even that linear-term-retaining local quadratic
unless g is negligible. Exact posterior curvature also includes residual second
derivatives omitted by GN. Without convergence and curvature/tail controls,
these values cannot be certified Laplace marginals. This concern does **not**
invalidate ordinary importance weighting with the explicitly declared r(b).

Twenty matched zero/shift comparisons verify the implementation identity

\[
\Delta\log Z_{GN}=\Delta\log p_{profile}
-\tfrac12\Delta\log\det H_{GN}.
\]

The GN log-volume component varies from **−8.6322 to +2.8468** across those pairs
(median −.10124). It is supported to say that the **approximate GN volume
factor varies**, and that finite conditional optimization can change a slice
height substantially. It is not supported to say how much the **true** strong
conditional volume changes, that profile bias was corrected, or that the true
weak marginal has been recovered. Different d also changes density dimension,
so raw normalizer values cannot simply be ranked across d as comparable point
probabilities.

Finite IS logZ minus GN logZ ranges from large negatives to a positive
**+332.872** at DEV_026/d24/shift, whose strong gradient is 30,726 and ESS is 1.
Neither the sign nor magnitude is a measured GN error against an exact integral:
the finite IS log estimator is itself unstable and generally biased. Ordinary
unnormalized IS estimates remain nonnegative/unbiased in expectation under
full-support reference draws; that theorem does not make this realized
128-particle log estimate or self-normalized material mean reliable.

## 6. Direct Gaussian material means are reference observables only

The scalar conditional field construction, precision-Cholesky orientation,
fixed-property treatment and complex modulus variance are consistent with the
separate analytic controls. `conditional_material_moments` correctly labels its
outputs `GAUSSIAN_CONDITIONAL_REFERENCE` and does not claim spatial covariance
or local support correction.

Those means and variances are exact decoder integrals for the supplied r(b).
They are **not** E[chi | a,u,y] for the actual Maxwell conditional unless that
reference law is established as the true conditional or its error is controlled.
The observed importance concentration and unfinished mode searches prevent that
interpretation here. They are also not the full posterior mean: the outer
distribution over a and u has not been inferred. The conditional Gaussian
mean−mode discrepancy is therefore a real property of this reference, not yet
a measured posterior Jensen improvement or Bayes-risk gain. SNIS corrections
with ESS near one do not resolve this distinction.

## 7. Changing-decoder proposal: no new target algebra defect found

`FrozenMaterialMobility` uses one paid fixed-anchor material-response factor
and the **current decoder derivative** to create a deterministic SPD precision.
The likelihood and source gradient remain full-target computations. The current
position-dependent MALA implementation includes both endpoint Gaussian log
determinants and quadratic terms, retains current gradients on rejection, and
freezes its rule/step size after warmup. Its analytic normalized-density audit
passes; that is an implementation control, not actual posterior mixing evidence.

Changing the mobility with state does not require a diffusion-divergence term
for correctness of this explicitly MH-corrected Gaussian proposal. No fatal
target/balance error was identified that would explain away the observed
reference failures. Poor efficiency and failure to equilibrate remain actual
negative evidence within these allocated runs.

One metadata caution: dynamic checkpoints say `metric_frozen=true` after warmup.
For this method the **rule and anchor factor** are frozen; the precision matrix
continues changing with x. Its final `METRIC.npz` covariance is the initialized C,
not a single actual production covariance for the changing-mobility kernel.
Reports must preserve that distinction.

## Supported scope at this point

The existing evidence supports a negative bounded-pilot conclusion: both actual
reference mechanisms remain unqualified; this changing-decoder configuration
does not establish a cost/mixing advantage; the finite Gaussian conditional/IS
construction cannot yet provide trustworthy posterior moments or marginal
volume comparisons. It does not establish a universal impossibility of Bayesian
Maxwell imaging, conditional integration, or other proposal designs. There is
no GT risk, calibration, teacher, learned-density or Bayes-action improvement
claim available from this review.
