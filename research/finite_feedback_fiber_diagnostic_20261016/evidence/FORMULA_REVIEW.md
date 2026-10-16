# Independent finite-feedback audit

Read-only source inspection; no physics, tests, imports, or source edits. Inspected `src/finite_feedback/{common,e0,diagnose}.py` and A30 `src/a30/physics.py`. Parent retains scientific judgment.

## Verified formulas and accounting

- A30 uses `L(chi)=I-diag(chi)G`, `E=incident+Gj`, `F=Sj`. The code's first increment is `L0^-1 diag(dx)E0`, and subsequent terms apply `L0^-1 diag(dx)G` on the **left** of the previous increment. This is the correct order. The endpoint-LU rational increment is exactly `L1^-1 diag(dx)E0`; both response and current identities are explicitly audited.
- Source-major data order is per source: 32 real receiver values then 32 imaginary values. `pack`, `real_jac`, and `chart_jac` agree. Material coordinates are `[Re chi, Im chi]/24`; the Jacobian's factor 24 and real noise whitening `sqrt(2)/sigma` are correct.
- Reduced-SVD strong energy is `sum((Vstrong m)^2)/||m||²`; subtracting from 1 correctly includes weak singular directions **and the unrepresented rectangular nullspace**. Relative cutoffs are diagnostic sensitivity conventions, not posterior or finite-fiber certificates. N3/R-base displacement is oracle geometry; TRAIN anchors have hGT=0 and cannot establish recovery.
- Decoder/Maxwell second derivative formula is correct: `T[D2chi E0 + 2 dchi G L0^-1(dchi E0)]`. Decoder second derivative is finite-difference approximated, so the comment “Exact second response derivative” overstates numerical exactness.
- Simpson-4/8 integration uses updated full-state Jv and correct weights; agreement checks the fundamental-theorem identity along the oracle straight path, not equal-data connectivity.
- Worst-case states: 17 shared-case oracle states (including one shared GT), plus 14 for each of two TRAIN anchors = **45**.
- Worst-case complex RHS: full states 360; four base plus four continuation receiver transfers 256; oracle finite series/rational 672; TRAIN finite series/rational 448; oracle path Jv 128; chart curvature 16 = **1880**. Thus headroom is 3 states / 168 RHS under 48 / 2048. Early continuation stops reduce cost. Every solve routes through A30 reservation; no hiding physical calls detected.

## Issues / limits to address before interpretation

1. **Continuation direction transport bug** (`diagnose.py`, updated-J sign test, currently near line 207): compares `D_current v` against `D_initial prev`, even when `prev` belongs to the preceding updated chart. Store the previous accepted physical tangent (or previous D) and align against it. At step 3 this can spuriously reverse/retrace. This does not invalidate endpoint physics, but can invalidate claimed continuation direction/path length meaning.
2. **Curvature anchor consistency** (near lines 186-188): centered differences subtract `base.chi` while outer points are `chart.decode(u0 ± eps*v)`. Require/check `chart.decode(u0)[0]` equals `base.chi`, and preferably center on the decoded value after the check. An anchor mismatch contaminates the fine D2 by `2 mismatch / (5e-4)^2 = 8e6 mismatch`. Saved stencil relative-change alone is not a substitute for this identity.
3. **Wall cap is not a hard elapsed-time kill** (`common.py`, line 37): elapsed time is checked only on reservations. A long SVD or already-admitted solve can cross 1200 seconds and continue until the next reservation; late work can even finish past the cap. State/RHS caps are hard at admission. Report wall enforcement honestly or add an external deadline mechanism if strict 1200 seconds is required.
4. **Redundant matrix work, no hidden RHS**: each legal anchor recomputes its unchanged base operator spectrum; the shared GT endpoint spectrum is recomputed for each oracle-anchor finite analysis. Both are already in action counters, but contribute wall time. First-order feedback solve is recomputed even though T gives the response directly; this pays for current-level audit and is included in 1880.
5. **Scope of correction**: `correct` uses a coordinate-orthogonal pseudoinverse at the accepted pre-predictor point, caps one normal update, then checks the exact endpoint against the original observation ball. It is a one-step candidate corrector, not converged level-set tracing. The output wording mostly preserves this distinction. Noise-ball compatibility permits nonzero data displacement; no equal-data fiber, posterior sampling, or inverse success follows.

No fatal operator-ordering, packing, mass-metric, Simpson-weight, or RHS-accounting defect found by inspection. Execution and resulting evidence remain unreviewed.
