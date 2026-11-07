# Paired-grid8 V2: generated, inverse not dispatched

Root amended the draft before any new data/result review:8 independent geometries × matched24/cross48to24 ×4 inherited methods =64 scientific tasks. The v1 originalDRAFT files and scripts remain unchanged and are copied under `v1_draft_preserved/`; no v1 observations ever existed. V2 is separate from original48 and3D, not16 independent materials.

Actual generation is COMPLETE with0failed attempts and0inverse runs. External Mac wait4 CPU8.816213s, elapsed9.129347s, peakRSS736,657,408bytes, one BLAS1 lane. Exactly8fine2304LU/96sourceRHS plus8coarse576LU/96sourceRHS,0adjoint; one fine Green setup. Per-object and self timers are nested inside this measured child. No main48GT read, GPU, remote resource or postseal additional truthFW was used. Future model-error scoring reads cached clean fields only.

`GENERATOR_CONFIG_V2.json`, `SOLVER_CONFIG_V2.json`, `SCORING_CONFIG_V2.json` are frozen at2026-11-07T02:30:32.384394Z. Solver has its new actual freeze_utc, parent_solver_freeze_utc, and explicit8×2×4authority. Fixed scoring is576cells,NRMSE denominator11.9, strict support>1, .15/.7 plus original.02..20box. All masks derive from exact fine2×2 blocks: purelow/high/background only if all4 equal2.5/12/.1; everything else is partial-volume boundary.

Data were generated once using new SeedSequence202611071118,8 independent case children and separate geometry/noise substreams per case. Conditions reuse the same standard CN draw but each has1% of its own clean-data RMS. Both use the same fine2×2 area-projected material target. Their noise/results are paired, not new independent replicates. Fine pulse support/scattering differs from coarse, so this is finite disk-model mismatch, not a PDE convergence proof.

Public17NPZ files: shared original24-grid operators plus16 noisy observations under `data_public/crossgrid8_paired/{matched24,cross48to24}/`. No fine materials/phase masks/family/seed are in inverse inputs. Private generation arrays under `sealed_truth/crossgrid8_paired/` were written only by the generator and remain closed until all64public predictions and selections seal. Fine/coarse cleanY is private and never used by inverse.

Tools:

- `generate_paired_grid8.py`: already executed once, refuses existing outputs, receipt/external cost/log preserved.
- `TASKS_V2.json`:64 portable forward-slash public tasks with two conditions,4methods CIE/Dual/Global/G19; same uniform.1 initialization and mature3000aux+max3000corrector/480CPU budget.
- `seal_paired_grid8.py`: implemented, syntax parsed only, NOT EXECUTED. Requires64unique terminal outcomes, stores allfinal+60/120/240/480 public minPhi selections, plus128 common-actual-end-CPU pair definitions (96within-condition method comparisons+32paired-grid same-method comparisons). Every original child remains charged once; views overlap and are nonadditive. CommonB readout does not imply exactly equal actual executionCPU.
- `score_paired_grid8.py`: implemented, syntax parsed only, NOT EXECUTED; requires explicit `--authorize ROOT_EXPLICIT_POST64SEAL_SCORE` after all64seal. Writes64final/prefix/actual-cost metrics, pure-phase/boundary RMSE, case-level grid/method paired differences and8cached model-error rows,0extraFW. It never pools64runs or16conditions as independent materials and makes no small-sample statistical-strength claim.
- `PUBLIC_INVERSE_INPUTS_V2.tgz`: prepared local public-only1,431,070byte package; no upload/remote command by diagnostics. Contains17publicNPZ+TASKS+SOLVER only; no generator/GT/configfamily data. A sourceowner may upload it separately without overwriting any active code. Existing frozen production runner/readout dependencies are reused, not bundled over active source.

Physical area semantics remain as in v1 plan: a=.45h; coarse radius.018,fine.009;4fine disk areas=onecoarse disk area. The target is an equal-area first-order polarization-volume material projection, not exact same scattering support. Both conditions score against that same projected target.

Only original288terminal plus explicit root-release of i9 logicalP0/2/8/10 permits64inverse dispatch, max4BLAS1 workers, noGPU/3D preemption. Diagnostics has not started any inverse queue. All failures/attempt costs must be preserved; no material redraw or scientific retry to rescue outcomes.
