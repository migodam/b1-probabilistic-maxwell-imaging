# Adaptive continuous-density diagnostic: exposed DEV only

Authorized scope: S1 raw matched CFM, anchor CFM_none and bounded-geometry CFM, each at fixed exposed DEV00000 and00001. One analysis-only Book I,900s total, CPU1,2GiB RSS, zero RHS/optimizer updates. No new point, fallback solver, tolerance relaxation or budget extension is allowed. No locked/confirmation truth is read. Existing fresh inference watcher, NFE32 wrapper, models, normalizers and source files are unchanged. New numerical code is confined to `adaptive_cfm_audit.py` and `audit_r4_adaptive_cfm.py`.

## Equation and signs

In the frozen TRAIN-standardized state s, `ds/dt=v(s,t,c)` and `dl/dt=-div_s v`. Backward DOP853 integrates1→0 from `(s1,l1=0)`. The resulting `l_back=integral_0^1 div(v)dt` is the forward log-Jacobian along that trajectory. Therefore

`logq_raw(x)=logp_base(s0)-l_back+log|det(d s1/d x)|`.

For a linear expansion `v=a*s`, the subtraction is `-66a`, the expected volume-expansion penalty. The outer determinant includes frozen frame/scale: `log|det F|-sum(log x_scale)`; anchor centering is determinant1. These models all use coordinate_none, so there is no GN coordinate determinant. The measure is original Gaussianized dx66; no material-decoder or native-uniform probit Jacobian belongs in this measure.

Each backward RHS evaluates all66 diagonal autograd derivatives, with no stochastic trace or parameter training. Forward cycle integration needs velocity only: integrate0→1 from the recovered s0, invert the original affine map and report maxabs, RMS and relative L2 error in raw x66. Cycle errors are reported without inventing a new pass threshold. The primary numerical threshold was frozen beforehand: two tolerance NLLs must agree within.05nat on each of the six model/object cases.

## Fixed tolerances and precision

DOP853 only, independently at `(atol1e-5,rtol1e-4)` and `(atol1e-6,rtol1e-5)`. Original stored Float32 parameter/buffer values are copied to CPU Float64; original API-normalized Float32 observation condition and computational-base values are captured then promoted. The true Gaussianized target point is mapped by the frozen affine in Float64. This diagnoses the continuous real-valued network with those original stored coefficients; it does not replace the frozen Float32 deployment algorithm.

Fixed32 comparisons separate: original public Float32 RK4; Float64 RK4 at its rounded standardized input; Float64 RK4 at the true standardized input. This separates arithmetic/input rounding from adaptive-integration differences. The earlier fixed16/32/64 audit (3/4 failed) remains intact; it used the older unmatched small raw CFM, not the matched raw model here. It is not silently relabeled as a same-checkpoint comparison.

## Array provenance and gates

Each run writes `FROZEN_AUDIT.json` before GT/density computation, binding3 checkpoint identities, original source identities,2 exact exposed parents, tolerances, threshold and scope. Input NPZs retain raw Gaussianized target, true/rounded standardized target, original normalized condition, base mean/Cholesky, center/frame/scale and an explicit exposed-numerical-audit role. Input JSONs bind original observation, explicit DEV anchor export and isolated DEV truth file byte/mtime identities. Only native `q` is read from those two exposed truth files; no material truth, physical state, decoder label or heldout array enters a solve. Native properties are transformed by ndtri without clipping.

All completed, failed and budget-partial records are retained. A solver failure is not an NLL-agreement failure; unfinished cases remain unestablished. Partial work never authorizes an extra tolerance/object or a changed model. Checkpoints and original source identities are rechecked after successful full completion. Final records include RHS/VJP counts, timings, cycle arrays, agreement status and unchanged scientific limits.

This is an integrator diagnostic. It cannot establish distribution quality, exact finite-RK-map density, exact MH density, a locked NLL result, or posterior qualification. An adaptive solve agreeing on two exposed objects would not alter any frozen fresh-cohort sampling NFE or model-selection rule.
