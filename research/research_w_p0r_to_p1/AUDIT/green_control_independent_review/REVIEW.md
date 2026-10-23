# Green reuse control: independent source-only review

Reviewed 2026-10-24. Scope: `src/wr/green_reuse_control.py`, with read-only interface checks against `../research_w_continuation/src/wq/fem.py`, `src/wr/physical_witness.py`, and `src/wr/native_world.py`. No import, main, Book, PDE, solver, numerical experiment, observations, labels, hash or source modification was performed. This review does not qualify runtime accuracy or make the scientific decision. Line references below refer to the source inspected in this snapshot.

## Findings affecting cost interpretation

### 1. Krylov warm timing receives a current matrix built by the reference Full arm

**Location:** Green lines 360–365, 495–518; cold convention lines 335–340.

`p.solve(...)` assembles and factorizes the actual candidate matrix, then `A1=sf._matrix`. The Krylov timer starts only at line 495; its MVP closes over that already assembled `A1` at line 499. The reference Full arm pays candidate coefficient evaluation/mass assembly/current-matrix preparation, while the Krylov warm row pays neither a candidate assembly nor construction of a matrix-free `A0 + delta*D` operator. Its current matrix is a useful paid artifact of the pooled experiment, not a free artifact of baseline-LU reuse alone. The current `cold_cost_convention.Krylov` mentions baseline/preconditioner and base-LU applications but does not explicitly price this candidate matrix dependency.

**Consequence:** RHS and stage totals remain real and paid; however, `PRECONDITIONED_KRYLOV.warm_wall_seconds` cannot by itself be presented as the complete standalone new-contrast deployment cost against `OPTIMIZED_FULL`. The baseline preconditioner is available legitimately, but the current candidate operator must also be made available. Report its inherited source/cost or separately include the cost of assembling A1 (or a declared matrix-free A0+delta D implementation). All non-Full output extraction also receives `sf.epsilon`, formed in the Full solve, at lines 380–383. That coefficient-sampling dependency is smaller but belongs in the same explicit shared-candidate preparation convention.

### 2. Standalone cold costs are conventions/component records, not complete disjoint numerical totals

**Location:** Green lines 229–252, 256–287, 320–344, 418–442, 530–540.

The code distinguishes actual pooled stage cost (`SUMMARY.cost`) from per-method warm cost and records useful baseline/primal/dual/R0V/exact components. It does **not** emit a fully assembled standalone cold/warm/fallback cost table. In particular:

- `setup.shared_ROM.cost` overlaps its nested primal, dual and R0V component snapshots; it must not be added to those components again.
- `setup.baseline.cost` also includes unit-update/full-RHS assemblies and baseline archive writing, which are not the minimal standalone baseline prerequisites of every method.
- `shared_ROM` includes both V/W projections and current feedback; applying its whole aggregate to FIELD_GOAL_ROM_V would charge enrichment/current setup that this method does not need. Conversely, crediting only its warm row omits its paid setup.
- `engineering_operation` events time operations only. Summing these alone does not reproduce all enclosing component CPU/wall cost, since reservations, residual checks, ledger snapshots, persistence and untimed projection products also execute.
- Geometry qualification (four Full responses, 16 RHS) is paid and saved, but is outside the separately named baseline/ROM/exact setup snapshots. It remains in actual pooled totals.

**Consequence:** Actual pooled expenditure is preserved; numeric standalone amortization/break-even conclusions require explicit nonoverlapping allocation, including qualification and shared candidate preparation. Do not sum overlapping setup records or infer standalone cold cost from `warm_cost` alone. Existing convention strings help disclose this limitation but are not themselves measured standalone totals.

### 3. There is no operational approximate-method-to-Full fallback cost

**Location:** Green lines 389–413, 465–520, 523–541.

ROM rows preserve their actual Full response error and residuals, but no acceptance threshold chooses an approximate deployment result or triggers/charges an additional Full fallback. GMRES nonconvergence yields `GMRES_NOT_CONVERGED`; an exception stops the branch with partial results. Exact resource rejection yields NOT_RUN and continues the other methods. Full reference responses have already been paid for every reached contrast, independently of which approximate method is being compared.

**Consequence:** These are engineering accuracy/cost rows, not measured accepted-deployment policies with fallback. A poor ROM/nonconverged Krylov row must not count as a successful cheap response, and the already pooled Full reference must not silently count as a free standalone fallback. Any later accepted-cost comparison needs an explicit accuracy rule and a separately allocated incremental Full fallback. No fallback is required to run this descriptive benchmark; the limitation is in interpreting deployment claims.

## Numerical formula/interface checks

### Finite update uses the total field, and the output includes the new-epsilon cross-term

Green lines 234–246 define the fixed target mask, `D=-k² M(mask)` and `bunit=k² M(mask E0_total)`. FEM `solve`, lines 579–611, confirms `s0.field` is incident plus interpolated scattered field, while `s0.scattered_dofs` contains only scattered DOFs. Consequently, with fixed geometry and `epsilon1=epsilon0+delta*mask`, subtraction of the discrete scattered equations gives

`(A0+delta*D) du = delta*k² M(mask E0_total)`.

Using `s0.incident` here would have omitted the baseline scattered-field contribution; the code correctly uses `s0.field`. `binc` is used only for the independent full-RHS audit and correctly uses incident field.

Green `_response`, lines 110–118, forms `total=s0.field+interpolate(du)` and `current=(epsilon1-1)*total`, then uses the original analytic free-Green volume receiver and direct term. Thus it retains both the direct finite material change and the new-epsilon × field-increment cross-term. It does not replace the volume output with FE point receivers. The full backend uses the same total-field/current/volume convention (`fem.py` lines 577–611 and 631–632). The physical probe mask is `<RADIUS` in both Green line 235 and `physical_witness.py` line 36, with conductivity fixed .02. The first-candidate independent volume identity audit is correctly outside Full deployment timing (Green lines 402–413); it is a runtime check, not evidence already obtained.

No sign, total-field or missing cross-term error was found by source inspection. The recorded matrix identity error at line 376 is not itself a stopping threshold; runtime values and the independent residual/volume checks still need review.

### Exact selector/Woodbury formula

Green lines 292–301 explicitly require every nonzero column of D to be among its active rows. With that check, `D=U Cunit U^T` uses all structural support, without rank truncation. Lines 450–458 implement

`du = R0 rhs - R0U (I + C U^T R0U)^(-1) C U^T R0 rhs`,

using the already paid `primal=R0 bunit` by scalar linearity. No inverse of a possibly singular C is required. Resource-scale rejection is explicit and does not relabel a truncation exact. The code permits the dense core/solve cost and memory to remain substantial; these must remain in reported CPU/memory totals.

### Same-current-space A17 and Galerkin really are the same finite projected equation

Green lines 266–282 construct QR-orthonormal V, RV=R0 V, `feedback=V^H D RV`, `current_rhs=V^H bunit`, and `current_gram=V^H V`. The current equation is `z + delta D R0 z = rhs`, with trial/test V and `du=RV*c`. Current Galerkin at lines 473–475 solves `(V^H V + delta feedback)c=delta current_rhs`; A17-labelled current feedback at lines 478–481 replaces `V^H V` with I. They are identical up to finite QR orthogonality and dense-roundoff error. Lines 482–492 explicitly check the field/response agreement and disallow an independent algorithmic gain claim.

The different FIELD_GOAL_ROM_V equation uses field trial/test V (`V^H A1 V`); it is correctly not declared the same trial space as current-load Galerkin. A numerical mismatch larger than the identity tolerance stops the branch, rather than become evidence of a new method.

### Ordinary ROM receives the same paid R0V enrichment

Green lines 268–278 pay rank(V) base solves for RV, construct `W=QR[V,RV]`, and project the ordinary field equation onto W. Lines 465–471 evaluate both FIELD_GOAL_ROM_V and FIELD_GOAL_ROM_ENRICHED_QR_V_R0V. A17 has no exclusive access to these paid enrichment fields. This is a fair enrichment-rights control, although the field/current test equations remain different.

The 12 duals at lines 260–265 are for the **baseline** volume functional: conjugating the assembled original volume-goal load and applying the Hermitian baseline solve is consistent with that functional. Candidate epsilon changes that functional; the enriched basis is not a freshly goal-adapted basis for every contrast. This is disclosed in setup metadata (`baseline original volume functional`) and does not corrupt candidate output extraction, which uses actual new epsilon. Accuracy must come from measured candidate errors, not from the word “goal”.

### Warm means new contrast with baseline reuse, not repeated response lookup

The 16 `linspace(2.94,3.06,16)` values are distinct and exclude the baseline contrast 3. Each candidate invokes a real `p.solve` and the FEM cache only reuses a factor when epsilon arrays are exactly equal (`fem.py` lines 456–464). Every reused method solves/reconstructs a new delta-specific field and extracts its response; no result cache is consulted. Baseline four-RHS primal reuse by scalar linearity in Woodbury is valid for this fixed-mask, affine single-contrast experiment. It is not a five-frequency posterior deployment claim.

## Review boundary

The main numerical formulas and the same-space/enrichment controls pass source inspection. The actionable issues above concern standalone operator preparation, nonoverlapping setup allocation and absent deployment/fallback policy. No real accuracy, speedup, memory feasibility, budget completion or scientific success has been established here. Source, budgets, ledger, jobs and physical arrays were left unchanged; only this review artifact was written.
