# R3-DIGITS blind protocol

Geometry is predeclared in `configs/geometry.json` before any ground truth. This is a nonlinear scalar TM volume-integral scattering model with 576 real material pixels, 12 line sources and 48 receivers. All methods receive the same shared operators and complex observations. The primary holdout and corrected DEV are noiseless; the unused original 1% noisy DEV draft is retained with a correction receipt and never scored.

1. Generate two DEV targets with local fonts distinct from holdout. Freeze first DEV predictions before releasing DEV truth. DEV metrics may guide subsequent shared settings.
2. Freeze solver source, configurations, objective, regularizer, initialization and budget; parent explicitly releases holdout generation.
3. Generate one target per digit 1–9 with private entropy, random anonymous IDs and random order. Observation files contain no class, font, mask, material, seed, or mapping. Solver reads only shared operators and public observation files.
4. Run every method on every observation, retaining failures and all K predictions. Top1 uses only the frozen physical objective and shared regularizer; oracle best@K is separate.
5. Freeze all holdout predictions and receipts before scorer opens holdout truth. Then expose seedless material, mask and mapping in `truth_after_freeze` and publish scores and figures. Generator entropy remains private.

Blindness is a worker input discipline, not OS access-control isolation. The experiment worker generates and scores truth; the solver worker must not read private folders or generator truth. The reconstruction prior is generic material bounds/background/smoothness; no digit classifier, template, label or trained image prior may enter it. Nine synthetic font digits do not establish population generalization.

Independent endpoint validation solves `(I-G diag(m))E=E0`, forms current `j=mE`, and predicts `Hj`; optimization uses the separate current equation `(I-diag(m)G)j=diag(m)E0`. This independently checks the same discrete physics, with no claim of cross-model or continuous-PDE certification.

Endpoint stability uses the total-field LU already required by validation: LAPACK ZGECON estimates reciprocal condition in the matrix 1-norm without a full inverse or SVD. Report field norms, normwise backward error, estimated condition, estimator status and timing. Explicit source work is 12 RHS per endpoint; ZGECON performs additional opaque estimator work whose internal triangular-solve RHS are not exposed and are never called zero cost. Validate physical top1 and, if different, the offline oracle best@K endpoint; do not redundantly re-solve every retained iterate.

The archive contains optimizer iterates, not certified solution branches or independent samples. Report pairwise material-RMS range and median, plus greedy physical-rank representatives separated by RMS≥1 as a descriptive count. The threshold is predeclared before holdout generation and is not an acceptance gate. Oracle best@K curves use ascending frozen physical-objective rank and are explicitly GT-selected retrospective summaries.
