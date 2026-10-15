# Research T: qualification first, compensated topology exploration second

This is an execution specification, not a record of executed physics.
Default authorization: E0 read-only audit and implementation tests that do not call a Maxwell solver. New physical experiments require a separately recorded authorization. Do not reuse any A, B, A30, or P2 remaining allowance.

## Scientific objective
Discover several physically realizable material fields that explain the SAME fixed observation, without truth-provided cavity count, location, support, or topology label. Do not equate lower residual, oracle MAE, a Jacobian identity, or a new optimizer stationary point with discovery.

Primary hypothesis: a finite structural edit, chosen using constrained compensation feasibility rather than residual descent alone, produces more materially and structurally distinct full-Maxwell-validated solutions per paid budget than ordinary topological-gradient proposals or independent multistart.

A17 is an OPTIONAL corrector/screening component. It is not the default generator and receives no preferential budget.

## Frozen evidence to read
- migodam/a17-maxwell-certified-exchange: closure_final/START_HERE.md; manuscript/CORE_THEORY.tex; CONDITIONAL_EXCHANGE_THEORY.tex; PAPER_TAP_FINAL.tex; PAPER_SUPPLEMENT_FINAL.tex; reports/A17_CLOSURE_REPORT_ZH.md.
- migodam/a30-maxwell-observable-material-manifold, historical revision 717d3d26a290248bd2b9aeb55f3d03694b407b84: A30_E0_E1_DECISION.md; src/a30/{plots,hull_audit,topology_supplement,frozen_a26_metrics,physics,compensation,representations}.py; saved result tables and original arrays.
- migodam/b1-probabilistic-maxwell-imaging, Research A revision f1b674eb4893e39126701ed48f3d3dcef3098c3c: research/research_a_20261015/RESEARCH_REPORT_ZH.md; results/E2/RESULT_SUMMARY.json; delegated/history/SUMMARY.md; solver/corrector source and receipts actually used.
- Research B: research/research_b/RESEARCH_REPORT_ZH.md, record the revision actually read.
- Joint P2 originals and A23/A24 actual-execution decisions, not earlier proposed GO decisions.

Known audit findings to REPRODUCE, not silently assume:
1. A30 figure 05 displays abs(chi-background), not an air mask. Its first Shell case is p250065_weak_a1_snr20.
2. Under the historical 0.15/0.25/0.35 truth-relative thresholds, that truth has 7 foreground components and 0 holes. N3 has 1 component and 0 holes. These are raster counts, not continuous topology certificates.
3. A30 image_metrics.feasible checks only Re(chi) in [-0.5,4] and Im(chi) in [0,0.5]; it does not test data compatibility.
4. Latest A E2 reports the same case's base residual 24.198733448253194 against threshold 23.785449004057032. The ring-background case has base residual 27.065435782542366. Independently establish whether these bases are byte/numerically identical to the older displayed N3 before merging records.
5. A30 Sparse-64's 35.69% is achieved offline capacity, not online recovery. Original pairwise distance 8.84 is not an observed-data residual.
6. Latest A: 5/8 two-dimensional conditions have no valid candidates; the sole extra archive material changes one structure-mask cell, is already present without compensation, and is not significant structural discovery.

## E0: zero new full states and zero new full-solver RHS
Read existing arrays, cached predictions/fields, metadata, and logs only. Algebraic residual checks using existing arrays are permitted; do not secretly reconstruct a missing cache through a solve.

Build one joined row per case AND per candidate:
- exact material-array source, background, acquisition, geometry, grid, units, forward implementation and discretization;
- original observation and noise realization, complex packing, sigma meaning, W, fixed tau and its calibration;
- material constraints, nuisance bounds, shared-nuisance scope and profile state;
- saved full prediction and its material/state identity;
- norm W(F_saved-y), noise-only norm, and available generator-versus-inversion prediction discrepancy;
- backward and constitutive residuals, available output error bounds, uncertainty status;
- raw absolute material/air masks, historical increment-support masks, connected components, 2D holes or 3D tunnels/cavities;
- per-case representation fits, optimizer/KKT/stop receipts and full paid costs.

Do not read truth into an online proposal module. Ground truth is permitted only in a separately isolated offline audit/scoring process.

Use four statuses: VALID_NUMERICAL, INVALID, NOT_AUDITED, MODEL_MISMATCH_SUSPECTED. Do not convert missing values to zero or infer a full-physics certificate from a box-feasible flag.

A truth residual above tau alone does not prove model mismatch: distinguish an outlying noise realization, a prior/constraint exclusion, wrong whitening, and an actual clean-forward discrepancy.

Deliver A30_CASE_AUDIT.csv, MASK_SEMANTICS.md, ROOT_CAUSE_MATRIX.md, MISSING_CACHE_MANIFEST.json, and E0_DECISION.md. If no same-data truth/base witness can be reconstructed without new solves, say so.

## E0-Q: separate base qualification, only after explicit authorization
Do not count failure to obtain a feasible base as proof that alternative feasible worlds do not exist.
Compare bounded full-material and multiscale-material correctors under the same noise model, initialization information, regularization policy, and paid limits. Check line search, active bounds, KKT stopping, and TV-prox caps. Keep all failures.
Suggested ceiling: 16 full-state attempts and 4096 complex full-solver RHS per case per qualification arm. These are proposed ceilings, not permission to spend.

Only full-validated bases enter E1/E2. Report conditional exploration results and unconditional base-qualification rate separately.

## Definitions to freeze before E1
- Material field shared across all illuminations; a frequency law is declared rather than independently refitted at every frequency.
- One legal nuisance vector per candidate world across the entire history. Different worlds may have different nuisance values.
- Air/cavity means a declared absolute constitutive range, not a low value of abs(chi-background).
- A protected cavity has a low-contrast core AND a material guard/separating barrier. Protecting core pixels alone does not prevent a vent to the exterior.
- Minimum feature size and contrast margin are in physical units. A two-cell feature on 24x24 is a raster pilot, not mesh-independent evidence.
- A fixed material distance, structural distance, persistence/threshold rule, and minimum significant edit size, independent of truth.
- Archive membership requires material validity, all-source full-forward compatibility, numerical checks and independent structural verification.

## Three algorithms
A. Finite compensated topological nucleation (primary): propose multiscale birth/death/split/merge edits throughout the allowed domain; never provide truth coordinates. Estimate bounded compensation feasibility; retain protected structure; perform a common nonlinear full-Maxwell corrector; append only validated, distinct materials.
B. Multiscale material-field continuation: propose tangent weak directions with explicit novelty forcing, follow a noise-feasible tube with predictor-corrector and pseudo-arclength when useful, refresh the tangent and active set, log observed topology events. Absence of a connecting path is not a disconnected-component certificate.
C. Nonlocal population search: independent starts and structural jumps with the same corrector, retaining multiple feasible endpoints. Optional deflation acts on anchored inverse/KKT problems, not on the uniquely solvable fixed-material linear Maxwell state equation.

## Screening mathematics
For finite edited material chi_e and legal compensation chart c in C_e:
  psi_e(c) = W(F(chi_e + Q_e c, nu_e(c)) - y).
Use the REAL-stacked Jacobian across all illuminations. Compute a finite edited response, not just a differential estimate, when ranking finalists.

Solve the bounded convex local problem
  min_{c in C_e} ||r_hat_e + A_hat_e c||_2.
Use a primal attained upper bound U and, when available, a valid convex dual lower bound L. A KKT iteration cap is not a lower-bound certificate.

With a genuinely justified uniform error epsilon over C_e:
  U + epsilon <= tau: sufficient local-chart feasibility witness;
  L - epsilon > tau: exclusion of this edit and declared chart only;
  otherwise: uncertain; use full correction or retain an explicit UNKNOWN status.
Empirical tangent checks are not a uniform error radius. Final archive acceptance still uses actual full endpoints.

Optional curvature screen for smooth material homotopies:
  J u = 0 is only first order.
  The component of F''[u,u] outside the permitted corrector range is a second-order obstruction.
For linear susceptibility VIE, use the two feedback injections; for DDA include the second derivative of the actual polarizability map. Do not substitute an a^5 cavity law.

## Equal-budget E1 pilot, after separate authorization
Four prequalified bases; three main arms:
1. ordinary topological-gradient proposal ranking;
2. compensated-feasibility proposal ranking;
3. independent multistart.
All use the SAME corrector, legality rules, stopping policy and archive metric. Arms 1 and 2 receive the same truth-free edit inventory. At zero residual the ordinary topological derivative can tie: use a frozen tie-breaking policy, not an adversarial implementation.
Include same-edit/no-compensation as an explicitly cost-accounted diagnostic, not a free extra search arm.
Suggested per-base-per-main-arm ceiling: 24 distinct full-state attempts AND 6144 complex full-solver RHS; stop when either limit is reached. Charge rejected trials, finite secants, tangent/adjoint solves, cache construction and failure attempts. Record setup and wall time separately. No GPU training.

A17 ablation uses the same proposed edits and the same paid budget. Compare its paired corrector space with equal-real-rank material Krylov, multiscale and full/material-adjoint alternatives. Remove A17 if it has no benefit.

## E2 / E3 / E4 (design now; do not auto-launch)
E2: the eight exposed A30 conditions, with separate E0-Q for invalid bases. Search from a genuinely hole-free qualified base to a physically defined cavity without using truth hole count/location. Report negative results and missing cases. Old exposed objects are not a blind confirmation cohort.
E3: several legal material families; compare local continuation, nonlocal topology jumps and finite correction. Proposed ceiling 48 states and 12288 RHS per family per arm; authorize separately. Record feasible paths separately from unvalidated interpolation segments and endpoint-only jumps.
E4: only after a positive E1/E2 mechanism, design a small nonconcentric vector-Maxwell test with multiple illuminations/polarizations, shared material/nuisance, mesh refinement and independent forward checks. Do not run a 3D campaign from this prompt.

## Outcomes and stopping
Primary: count of genuinely distinct, significant, full-validated structural materials and discoveries per paid state/RHS/wall budget; time/cost to the first such witness.
Secondary offline: complex material MAE, Re/Im MAE, Best-of-K with K-valid coverage and missing values, cavity precision/recall with localization matching, support IoU, boundary errors, and persistence under fixed thresholds/resolution.
Always report residual to OBSERVED data and pairwise prediction distance separately.

Mechanism GO: on at least two independent qualified parents, a significant protected edit is invalid without compensation and becomes valid with compensation, while the matched-cost comparison gives a reproducible discovery advantage. This is a pilot gate, not a population significance claim.
NO-GO for the advertised mechanism if gains are only smaller residuals outside the threshold, oracle improvements, sub-resolution holes, looser nuisance, extra corrector work, or an edit already valid without compensation.
NO_QUALIFIED_BASE is distinct from NO_ALTERNATIVE_WORLD_FOUND and from PROVED_EMPTY_LOCAL_CHART.

## Research B / Pro G factorial control
Separate new information from optimization using old/new observation blocks crossed with no-extra/matched-extra continuous correction. Preserve shared nuisance and charge all new measurement and computational work. Candidate omissions must remain possible even when one candidate survives.

## Required deliverables
E0 audit files; source and configuration manifest; append-only NEW ledger; per-trial failure receipts; three algorithm implementations; no-physics unit tests; pilot authorization gate; plots with raw data; candidate package conforming to PRO_G_INTERFACE.md; final decision distinguishing proven, conditional, observed, heuristic and NOT_RUN claims.

Do not modify historical results or use locked truth for online decisions. Do not enlarge a campaign or train a network merely because a local algebraic check passed.
