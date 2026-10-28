# Independent pre-locked review

Reviewed the preserved source snapshots in `source_snapshot/`; provenance is
recorded in `REVIEWED_SOURCES.json`. **Initial status: HOLD pending three code
fixes and the final pre-locked freeze.** The root acknowledged the findings and
is implementing the fixes. This review did not modify running sources, invoke
physics, load GT, or launch workers.

## Findings requiring correction before locked generation

1. **Anchor cap bypass.** `direct_anchor` originally counts only `get()` calls
   against the frozen cap 64. The Jacobian restoration and final best-LU
   restoration call `evaluate` directly. An isolated fake-target fixture made
   64 counted calls and a 65th uncounted restoration; the status still reported
   `fresh_state_attempts=64`. All fresh restoration attempts must either fit
   within the 64-state cap or be covered by an explicit pre-locked restoration
   policy, with total fresh calls reported accurately. No Maxwell was used to
   demonstrate this bug.

2. **Finite jet tests omitted from probe verdict.** The original probe records
   the two finite jet prediction trials but does not include their execution
   status or errors in `passed`. Its actual saved trials are excellent
   (`2.22e-9` and `1.39e-8`, both below the frozen 0.1 threshold), so this is a
   latent enforcement gap rather than evidence that the existing probe failed.
   Require both finite trials to execute, produce finite errors, and satisfy
   the frozen threshold before allowing locked work. The separate finite-radius
   two-adjacent-radii/0.02 claim gate remains necessary.

3. **Retraction can exceed the trust subproblem bounds.** The original corrector
   scales the unprojected GN step, then exactly retracts to the material-distance
   shell. It correctly recomputes the predicted decrease using the retracted
   step, but never rechecks that step's latent or material trust bounds.
   The root requested a full decoder-only check of the existing DEV4 output;
   results are in `DEV_RECONSTRUCTED_TRUST_AUDIT.json`.

   All **150 executed DEV correction steps** were reconstructed from saved paid
   endpoint x, using the original decoder and derivative. Every reconstructed
   latent step matched the saved log. Four executed steps exceeded a bound:

   | Parent / arm | Radius / iteration | Exceedance |
   |---|---|---|
   | R6_CDEV_0002 / C_JETS | 0.0025 / 1 | linear material 0.001258808 > 0.00125 (+0.7046%) |
   | R6_CDEV_0002 / C_JETS | 0.005 / 0 | latent 0.503981516 > 0.5 (+0.7963%) |
   | R6_CDEV_0002 / C_JETS | 0.01 / 0 | latent 0.500999501 > 0.5 (+0.1999%) |
   | R6_CDEV_0002 / C_JETS | 0.01 / 1 | latent 0.753775652 > 0.75 (+0.5034%) |

   The other 146 executed steps had no latent or linear material exceedance.
   None of the 150 steps exceeded the finite decoder material displacement
   bound; the single material exceedance concerns the linearized norm actually
   constrained by `gn_step`. Preserve DEV4 as pre-fix evidence and identify the
   affected rows explicitly. Do not silently call the old run trust-compliant.

## Isolation and measurement-only anchors

An AST inspection enumerated every `np.load` path in the reviewed runner.
For locked parents, `load_parent` reads only
`raw/online/<parent>.npz`, requires exactly observation/sigma/parent fields,
checks the original sigma, and calls the fresh direct anchor from x=0. Its
optimization criterion is the original `Phi + ||x||²/2`, and fallback selection
uses the lowest finite observed objective in that same measurement-only history.
No locked online call opens `offline_truth` or reads a material label.

Generation and parent processing are distinct command/process invocations.
`generate` writes truth privately and writes only y/sigma/parent to the online
file. Historical material coordinates and geometry are admitted only for the
explicitly exposed `R6_CDEV_` development branch. The online isolation here is
verified source/call-path isolation, not filesystem permission isolation.

`seal` requires all eight locked `COMPLETE.json` files and records the endpoint
artifacts before scoring. There is no scorer in this runner. The reviewed code
does not itself enforce that final REV1 source/contract sealing and DEV4 review
precede `generate`; this must be completed before locked generation, ideally
with an explicit gate. The available `E2_PREFINAL_REVISION1.json` records DEV-only
changes and is not yet a final post-DEV4 source freeze. No locked-data-selected
hyperparameter change is justified by this review.

## Physics, source groups and accounting

- Original P1 prior, gain 1, sigma `0.00119211798514519`, complex128, and full
  source-major Re32/Im32 whitening are preserved. Gaussianizing Uniform2 with
  normal CDF and using identity Gaussian prior precision is correct.
- Group A is source 0..3, exactly the first 256 real packed entries; group B is
  source 4..7. Group-selected correction still performs the full eight-source
  forward and full 32-receiver Jacobian before slicing the objective. It does
  not claim acquisition savings from discarded channels.
- Common compensation uses the same lower-order material path, common real
  strong-coordinate columns, and separate A/B/joint real least squares. The
  added group-A weak direction is evaluated through the complete acquisition.
- Main arms share material-distance matching, signed progress, original prior,
  the same number of trial-state slots, and objective-based feasible selection.
  Corrector actual/predicted objective comparisons use the matching source
  subset and original prior. Prediction is recomputed after nonlinear distance
  retraction; it is not taken from the obsolete unprojected step.
- One live LU owner is retained by Target; a new evaluation replaces the old
  Point state with factor-free endpoint arrays and clears receiver caches.
  Older base/best Points are used for saved arrays/J only. Fresh solves remain
  reserved through the sole flock worker and the independent 50000-RHS ledger.
- The saved first probe has **200 RHS = 8 base + 24 jet + 8 source adjoint +
  144 central differences + 16 finite jet endpoints**. That is 168 forward RHS
  and 32 derivative RHS. Its receiver T was validated from the same historical
  paid material, so an additional 32-RHS receiver solve was correctly avoided.
- `jet_cost` records the main direction's shared 24-RHS setup; the additional
  group-A weak jet setup is a separate control and remains included in the
  complete parent ledger cost. Reporting must not double-count the shared setup
  once per radius or omit control costs from the parent total.

## Remaining gates and limits

The prior clipping map is only C2 at saturation seams. Keep the regular-stratum
condition, seam-crossing counts, nonzero first-order response and finite-radius
physical prediction checks. The five path nodes cannot certify behavior
between nodes. Neither zero-physics algebra nor current derivative QA certifies
posterior probability, task-scale material nonuniqueness, connected fibers, or
owner acceptance.

This audit itself cost **0 new RHS / 0 new states / 0 new LU**, accessed no GT,
and imported no A30/A24 physical solver modules.
