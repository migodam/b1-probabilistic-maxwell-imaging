# R12 research loop — completed experiments and publication audit

## 2026-11-11: inherited Pro reproduction

The unchanged four-case prototype ran in a new output directory. Ring DR relative material L2 .0870014, shell .1777494 and disk .1286312 agree closely with source main trajectories. Bridge differs (.608695 versus source .6049); after the common200-step physical corrector it becomes .728807 while the data residual falls to .0343822. Numerical trajectories are not asserted bitwise portable across SciPy/BLAS. Both original and replay remain.

## Mature comparisons close a missing control

The baseline adapter ran two explicit Fourier schedules on the exact same100-cell finite model. It retained mature sourcewise/global current and material rules, adapted box/grid rather than reconstructing Green. Ring still defeated all five baselines; shell was recovered more accurately by mature methods than by the original DR candidate. Therefore the Pro result has conditional rather than universal value. Details, cost and schedule qualifications are in delegated/baseline/RESULTS_ZH.md.

## Exact R11 bad endpoint intervention

Three saved bad materials were located, rather than regenerated. Full physical initial fields are recomputed; only original Phi including inherited prior selects candidates. Four arms are direct continuation, AP, DR and predetermined smooth perturbation plus direct multistart. Remote replication uses own pinned i9 P-core children. This is exposed mechanism development, not confirmation.

Initial local ring observation: direct continuation remains at Phi .009385376; AP incumbent likewise. DR first raises physical Phi (audit50 .0145825), then falls (audit100 .0020018, audit400 .0000650); the first observed noise-feasible audit occurs at400. The accepted physical material and its correctness are separate questions; joint scoring waits completed controls. This is a finite escape-path observation, not a global convergence proof.

## Cost refinement from algebraic reuse

Identical per-source receiver maps can share a grouped port operation and summed eigenspace residual powers in the same scalar multiplier equation. Optional grouped v2 leaves the projection problem unchanged; one toy fixed-projection benchmark measured2.40x CPU improvement. This is engineering evidence for that workload only. Version1 remains available and the already-started basin/vector runs keep it.

## Independent next questions

- Does the R11 escape survive a same-cost multistart control and the other two actual endpoints?
- Does it transfer beyond the exposed ring to independent mixed/high materials?
- Does the full physical NN interface preserve useful starting structure without a new prior?
- Are matched vector improvements retained under common longer correction and cross-grid mismatch?

These questions drive the ongoing runs; no universal-success threshold is imposed.

## 2026-11-12: completed evidence and root acceptance

The questions above have been investigated, not used as stopping gates. All three saved-endpoint interventions and their controls were scored; a postseal saved-pool oracle (0 new forward calls) narrows the positive claim to the actual audited finite pools. Independent confirmation comprises24 materials×8 methods with all192 terminal predictions sealed before official scoring. DR mean NRMSE .184631 is close to Dual .185456 and above G19 .143891; this does not establish equivalence or universal superiority.

The4-object frozenNN interface gives conditional complementarity, not independent NN-DR confirmation. Both moderate and high vector groups were actually inverted; common long correction removes moderate differences, and high matched rawDR loses accuracy and noise feasibility after the same soft-objective corrector. A pre-registered feasible-candidate minR selector was evaluated separately, with unavailability retained. Fresnel and Austria were executed; full nativeDual continuation and matched32 diagnosis distinguish budget and discretization effects without rescuing the DR claim.

Remaining nearest-method comparisons and finer3D are disclosed as paper limitations, not fabricated finished experiments. The current scientific work is complete at this scope. Final work is full raw-evidence archiving, publication, and public readback. An external-storage I/O failure interrupted publication on11November; restored files are checked against receipts on12November. No new science is run just to recover lost text.
