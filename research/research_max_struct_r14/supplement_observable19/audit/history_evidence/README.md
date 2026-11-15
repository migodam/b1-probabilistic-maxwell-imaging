# Saved-history evidence, observable19 supplement

This is a mechanical extraction of existing scalar results and actual source. No GT arrays, numerical checkpoints, scorer, optimizer, forward operator or remote service were opened/run. Scientific interpretation remains with the parent. `extract_evidence.py` reproduces the CSVs from the named local inputs. `EXTRACTION_RECEIPT.json` records the latest extraction CPU/wall and `EXTRACTION_EXECUTIONS.json` retains both completed extraction invocations (total .104770 CPU s); reused solver cost is not added to extraction cost. Earlier interactive inspection time was not instrumented and is explicitly unknown.

## Local closed-loop diagnostic: quantities and spaces

The source is `research/research_max_struct_r14/src/local_closed_loop.py`. This is the declared side8, 6-source, 16-receiver, zero-noise/zero-prior anchor diagnostic, not the576-cell benchmark. The public scalar summary records48 independent current spatial directions among64 native Fourier columns (relative current reconstruction gap7.24e-15). `SPACE_LABELS.csv` distinguishes material real64, ambient current complex64×6, independent current48/source, and data complex16×6. All norms copied here are Euclidean/Frobenius source norms; a current-space norm cannot be compared numerically with a material-space norm as if units/spaces matched.

The table does not contain a generic propagation spectral-bound field. It contains the spectral radius of the constructed **material closed-loop map T** and a separate induced norm expressed in the common **B00** material Gram:

|pair|spectral_radius(T)|same_B00_norm|
|---|---:|---:|
|00|.9999945440301764|.999994544030176|
|01|1.0027124491177037|1.010669973434584|
|10|.9999865689800003|1.2802093683600224|
|11|.9999886575226579|1.2413793480142075|

Actual construction: lines79–95 form independent-current normal solves;99–105 form the material Gram/readout and T;101,108–112 use **B['0'] for every pair's displayed norm**, separately from eig(T). These are distinct stored quantities; this extraction supplies no theorem identifying them with a different propagation bound or establishing contraction.

At epsilon=.001 in the predeclared h1 direction:

|pair|current_error_norm (complex64×6)|current_error_material_norm (real64)|material_error_norm (real64)|idealgap_norm (real64)|
|---|---:|---:|---:|---:|
|00|.4984487044|.8129338693|1.38470e-8|.8129338684|
|01|.4984487044|.8316134566|.0019727151|.8318032386|
|10|.2165550373|.3615479273|1.20632e-9|.3615479277|
|11|.2165550373|.3560411092|.0029254845|.3554796885|

Exact definitions at144–161: `ec=dJ−ideal_dJ`; the source `readout(ec)` transports current error into material space; `em=dm−readout(dJ)` is the residual relative to that root-linearized readout; `gap=dm−T h`. Thus `gap=readout(ec)+em` by the recorded algebra. The scalar **norms are not additive**. Do not rename em a certified optimization-error bound: it is the defined finite-response residual. All36 directional rows, including sweep/fallback/box activity, are preserved in `LOCAL_DIRECTIONAL.csv`.

At epsilon=.001, finite-response relative additivity residuals (material/current) are00 .01065976/.01240257,01 .01045653/.01240257,10 .14353222/.15447828,11 .14378724/.15447828. These are hsum−h1−h2 central-response discrepancies with denominator sum of h1/h2 response norms (source169–176); all12 rows retained. The saved notes additionally report zero arc sweeps at exact root versus1–2 perturbed sweeps. These finite-program observations and exact residual closure alone do not establish a Frechet derivative or identify which theoretical hypothesis fails.

## What actual G19 selects

Actual executable chain: R13 `src/full_current.py:17,137–146` imports/calls R11 `prepare_full`; R11 `src/material_readouts.py:97–101` calls R10 `prepare_green_strong` and aliases `matrix_material_readout`; R10 `src/lowrank_denominator.py:113–128` selects `green_strong`.

R10 lines19–29 construct strong **right singular vectors of G** at threshold.1 of largest singular value, separately strong vectors of S, and a twofold alternative. The selected G19 basis is **Vg**, not Vo/twofold. Lines38–50 implement `K=P_G G P_G` through U and B; lines71–79 retain the denominator's full material derivative. Lines86–110 perform bounded full-material L-BFGS (default8 iterations) with scalar-seed and acceptance fallback. Reading the code does not support calling this a19-dimensional material parameterization or a simple change to observation coordinates.

Naming provenance: R10 `METHODS_AND_ABLATIONS.md:12–15,31` calls it Green19 and distinguishes19/30 small denominator ranks from576 material variables. R11 same-named document:3–19 explicitly calls G19 a project experimental name inherited from R10, “full” meaning complete material derivative, not full Green rank. **No R19 study/source is invoked in this chain**;19 is the saved two-dimensional strong-Green rank naming. This is source provenance, not a new basis recommendation. A future proposed material observation/Green parameterization is a different implementation choice not evaluated by these saved results.

## Saved scoring and the two0021 IDs

`INDEPENDENT24_SIX_METHODS_THREE_VIEWS.csv` has432 rows=24 cases×six methods×raw/corrected/primary. Methods retain exact tags: g19,cie,dual_cie,pair01,pair11,dr. Here requested A+ is explicitly **R14 pair01 (A+-Arc)**, not the historical R13 MM tag. Scores are copied from previously authorized `analysis/independent24_v1/ALL_VIEWS.csv`; no score is recomputed. Candidate `CPU_seconds` and `physical_cost_at_candidate` are distinct from `whole_child_*` actual OS CPU/final counters joined from the accepted576-child ledger. All432 cost joins are present. Shared candidate rows repeat a task's whole-child cost; **do not sum view rows as extra solver invocations**.

`R14_NEW_0021_ALL11.csv` contains33 rows, all11 methods and three views. On its primary view, pair11 is best among all11 for each listed metric: NRMSE=.0869648351, IoU=.9251968504, Phi=5.560868703e-5. G19 has NRMSE=.5121844277, IoU=.1910569106, Phi=.03488127453. This saved case does not show a G19 metric win.

`R13_FRESH_0021_SIX_METHODS.csv` separately contains18 saved uniform-initialization rows for **r13_fresh_0021**, six historical methods and three views, from the exact R13 supplement compact table. Its primary NRMSE winner is Dual-CIE (.06660743773); maximum IoU=.7054263566 is shared by Dual-CIE and A+-MM; minimum Phi=5.344605166e-5 is A+-MM. G19's primary NRMSE=.07128711547, IoU=.6791044776, Phi=5.355453597e-5. The IDs are different cohorts, not aliases for one object. `CASE0021_METRIC_WINNERS.csv` retains per-view/metric winners and ties rather than a universal winner claim.

Historical and new budgets/stopping contracts differ; old references are reused evidence, not causal equal-cost comparisons. Corrected-vs-raw effects must also retain method phase semantics (direct has no auxiliary proposal/corrector hierarchy). No threshold or GT-based selection was introduced here.

`SOURCE_MAP.csv` gives exact paths and source line ranges; CSV provenance includes original table and one-based CSV line. Saved scalar scores and costs are reused unchanged; source algorithms and old reports are untouched.
