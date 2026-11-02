# Independent output, material-class and cost review

This review reads the final saved F arrays/receipts, rank-class analysis, implementation source, backend qualification, D diagnostics, M qualification and T source qualification. It launches **no new forward solve or experiment**. The metadata corrections are reflected in the reviewed `results/families/SUMMARY.json`; no older DOF or defect label is silently treated as authoritative. The core scientific conclusions remain the parent orchestrator's responsibility.

## Verdict

**PASS for the bounded finite-model diagnostic claims.** The saved results support legal tensor ambiguity constructions and the stated source/receiver controls, without an imaging or GT recovery claim. **M primary NOT_RUN; T finite-channel certificate INCONCLUSIVE; ROM NOT_RUN; D2 full-history dynamical classification NOT_RUN.** No higher-level experiment gate is supplied by the DDA qualification or local D1 spectrum.

Material concerns are scope limitations, rather than a missing successful image: generic complex-passive fibers are not real-symmetric/diagonal/block-constant fibers; re-nulling each refined discretization is not continuum convergence of one fixed material pair; source separation is a candidate-pair prediction under an ideal coefficient-power metric; and no actual active measurement, blind inversion, true-target-containing library, or GT recovery is demonstrated.

## Saved arrays and material legality

`audit/F_SAVED_ARRAY_AUDIT.json` independently checks the four stored primary `.npz` artifacts, using only saved values. Every saved field/material array is finite. `branch−anchor` agrees with saved `V` to at most2.22e-16 entrywise. The real analytic branches are exactly transpose symmetric with real minimum eigenvalue1 and loss0. The passive-complex branch reciprocity defects are2.96e-18/2.04e-17; minimum real-part eigenvalues are2.16858/2.16852 and minimum loss eigenvalues .21951/.19901 on the two grids. Thus these saved endpoints remain strictly passive reciprocal tensors; the real analytic class is SPD without strict loss.

The real branches have RMS Frobenius material distances .48849/.45260; their old-source internal-E and port gaps are stored exact zeros. Passive-complex distances are .03820/.02646; their old-source whitened port gaps are9.02e-14/1.46e-13 and the annihilator defects are2.77e-18/5.52e-18. These are floating finite-DDA algebraic checks, not independent continuum error bounds.

The corrected `constitutive_defect=0` refers to constructing `j=(epsilon−I)E`; it is distinct from `same_field_annihilator_defect=||V E_old||` and `unit_ray_excitation=||V E||`. The old noise-ray value .01524 is nonzero weak excitation, not a constitutive-equation violation. A small linear-system residual around1e-15 is also not the discrepancy of an inverse image against unknown data. State DOF metadata now correctly gives36N real electric-field unknowns for six complex3-vector columns;12N applies to the two old columns.

## Rank classes: preserve real/complex and spatial intersections

`results/families/RANK_CLASSES.json` examines one old source, two old sources and those two plus an independent source on the **same saved anchor fields**. Its relative threshold is a roundoff diagnostic, not a tuned weak-mode cutoff.

| Anchor / source subset | Local complex rank | Real span rank | Consequence for symmetric annihilators |
|---|---:|---:|---|
| Analytic real, one old |1|1|3 real symmetric coefficients/cell; two diagonal axes unseen|
| Analytic real, two old |2|2|1 real/complex symmetric coefficient/cell; z-diagonal direction unseen|
| Analytic real, old+independent |3|3|No local same-field symmetric fiber|
| Passive complex, one old |1|2|3 complex symmetric coefficients/cell;1 real symmetric coefficient/cell|
| Passive complex, two old |2|3|1 complex symmetric coefficient/cell; no real-symmetric fiber|
| Passive complex, old+independent |3|3|No local same-field symmetric fiber|

For the analytic anchor the stacked all-cell old fields have block-constant rank2, leaving a uniform z-axis tensor direction. For the generic passive anchor the stacked rank is already3, including the one-source case: the local complex direction varies spatially, so **no block-constant tensor** annihilates the whole domain. Its axis-diagonal null DOF is0. Counts labelled complex DOF are complex dimensions; they cannot be substituted for real dimensions. Full local rank3 only removes this same-field construction; it does not prove full inverse uniqueness or stable reconstruction.

## A / B / C and the noise contract

A adds receivers while retaining the same old illumination span. Its passive gaps remain about1e-13; analytic gaps remain0. B solves actual unitary-mixed old incident columns on retained LU factors, with4 extra RHS per primary case. The2×2 source matrix is unitary, so the isotropic bounded-noise metric is preserved. B does not add an independent excitation direction and its gaps remain0/roundoff.

C optimizes a coherent combination of the four independent library columns, under the **declared identity coefficient-power metric**. This is the leading singular vector of the complete candidate transfer difference, not a GT-ranked source. The real pair's predicted whitened gaps25.94/30.33 exceed the single-source two-ball threshold12. The passive pair's1.471/1.070 gaps do not exceed12. Consequently the generic passive example shows predicted data separation but **not noise-qualified pair resolution at the frozen noise level**. An ideal plane-wave coefficient norm is not a calibrated finite-antenna power budget. No selected new measurement is generated and no actual unknown target is acquired; C remains prediction-only.

The finite-noise weak ray is separate from the exact fiber examples. Its local complex rank is3 with nonzero smallest returned singular value .01148; nevertheless one finite ray has small nonzero excitation. The implemented current-resolvent identity and Neumann denominator retain both propagation factors. Its full chosen interval has `T*a≤.4`; the measured pair gap .01750 is below the floating bound .16613. The actual common observation is a disclosed anchor output plus a seeded direction scaled to exactly half the deterministic radius. Both actual residuals5.19615/5.19611 are below radius10.39230. This supports a **constructive finite-model ambiguity lower bound**, with K uniform only on the stated ray. It is not a Gaussian coverage probability, a recovery upper bound, or a uniform material-class certificate. The arithmetic norms were not outward rounded.

## Refinement and backend qualification

The independently evaluated Mie sphere controls reduce electric far-vector relative error10.47%→3.66% and scattering-cross-section error21.81%→5.99% between the two staircase grids. This is useful bounded DDA evidence, without an arbitrary acceptance threshold or a claim of continuum certification. The staircase volume bias changes sign; two points do not establish an asymptotic rate. The large-radius near/far check and transverse far amplitudes are internal consistency checks, while analytic Mie supplies the distinct reference.

The exact passive fibers separately rebuilt on the coarse and fine fields are different discretized directions. The transported coarse `q` field, multiplied by the same analytic bump to construct `q q^T`, gives a fine-grid whitened gap .001878; it is not exact on the fine grid. Its explicit `INCONCLUSIVE` status is appropriate. The interpolation/extrapolation rule is disclosed; extrapolated cell centers lie outside the bump support. No continuum noise qualification follows from this one refinement.

The backend's cell-center magnetic field and analytic curl radiation kernel do not supply compatible edge/face spaces, physical mass matrices, weak Maxwell residual dual norms, or boundary traces. M's stop therefore follows its stated qualification requirement. There is no substitution of an L2 tail for the requested H(curl) metric, and no demonstrated claim that H is easier to represent. T likewise has not obtained a usable material-uniform exception bound `k<p` or complete discretization/noise/compression δ, so the certificate is **INCONCLUSIVE**, rather than a fitted empirical certificate. ROM's required selfadjoint reciprocal colocated temporal observations are absent.

## D diagnostics and historical interpretation

D1's common beta_material root derivative and local spectral radius .836636 apply to the fixed **full-current, exact-current-block, real-material, interior, no-prior** map on one disclosed4×4 root. The finite-step, prior, active-face and one-PR controls delimit that claim. In particular the single-PR directional map has a nonvanishing additivity defect; it does not inherit that Fréchet spectrum. The active-face example is nonsmooth but does not refute the tested first-order beta coincidence.

D3 preserves the objective under full source-metric covariance, while sourcewise iterations still differ; the difference already appears before PR memory is accumulated. Full metrics plus global search/PR agree to roundoff. This is not a newly superior imaging method. D2's historical scalar logs lack the full current/gradient/direction/line-step recurrence; their wall cuts cannot be classified as basin fixed points, cycles or chaos. Historical R2 8/12 and R3-E10/12, including nine specific successful digits and failed ring/curved_bar, remain historical rather than fresh blind results.

## Receipts, packaging and cost audit

`audit/receipts/INDEX.json` indexes18 standardized receipts: six F records, two backend states, D1, D3, shared D setup, D2/M/T/ROM stop records, and three rendering/metadata records. Every required schema field is retained; each null has a dotted-path reason. Raw numerical files are not overwritten by consolidation. The source/config hashes remain null under the user's instruction. `audit/D_INPUT_PROVENANCE.json` records exact frozen Git source paths and the packaged dependencies.

The D source is now standalone within this study: it reads `data_public/dynamics_operators_4x4.npz` containing **only G/E0/S/grid_shape** and loads the included frozen R3/diversity adapters. The R3 snapshot differs only in dependency import routing; numerical code is unchanged. A subsequent import-only check resolves these dependencies without calling `main`, creating a state, or rerunning an experiment. The original executed run's results/costs remain authoritative.

`audit/FINAL_COST_TOTALS.json` adds each independent outer process exactly once. Backend and F setup/solver/per-case counters are nested, not additional paid processes. F's1.30215s process includes all its task execution and array I/O; its per-case1.26207s sum is not added again. D's .04976s process includes its one shared root and active-face root. Totals are **15 original forward LU,72 forward complex RHS,208 additional retained-LU triangular RHS**, plus546 small exact-current/JVP Cholesky RHS. These categories are not interchangeable. Comprehensive Green-product and small-SVD totals remain null because LAPACK norm/least-squares/legality subcalls were not exhaustively counted.

The three numerical process timers total **1.436489s wall /1.689674s CPU**. Metadata scoring adds .01602s. The two actual F rendering passes are both retained, with2.264885s combined wall; the backend render adds .382026s. Thus all measured process times sum **4.099421s wall /4.243236s CPU**. These are process-time sums, not the elapsed research duration, speed comparisons, or complete execution overhead. The confirmed earlier render is included through `previous_plot_runs` once.

Setup suboperations are retained, but exact setup/online/validation/scoring partitions were not instrumented across the whole run. No online blind inversion executed. Source preparation, interpreter imports/startup failures, shell/tool orchestration, literature access, independent read-only inspection and report composition are incompletely timed and remain null/UNMEASURED; local monetary cost is unknown. No fake zero is assigned to unknown overhead. The three pre-numerical D startup failures have known zero forward LU/RHS and unmeasured wall/CPU. These limits must survive parent consolidation.
