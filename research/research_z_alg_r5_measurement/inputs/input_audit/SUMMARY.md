# Z-ALG R5 input audit summary

The R5 handoff archive is intact: all 25 manifest-listed files match their recorded sizes and SHA-256 values, all JSON parses, and the package contains no unlisted files or unsafe manifest paths.

The recorded historical commit `55724cb1a7bafbd9afa3b95ac817522f8dcaf2d0` exists in the local repository. R2-R4 research trees have no Git content difference from that commit to current HEAD `18620d199a4882c1f189dcd276c6572f7414e353`.

`GENERAL6_00` is fully traceable. Its public observation is `research/research_z_alg_r2/experiments/general6_public/GENERAL6_00.npz`, containing only complex `y` with shape `(1,3)`. Geometry is in `experiments/GENERAL6_CONTRACT.json`.

The Q20 result file is `results/general6_GENERAL6_00_general_quotient20.json`. It resolves 20 quotient paths, physically filters 2 candidates, and both pass fresh full-VIE validation. The 64-start baseline file is `results/general6_GENERAL6_00_ReducedFullJ64.json`; it finds 1 qualified root.

The two Q20 roots reproduce material RMS separation `1.1587903026947952`. Both residuals are approximately `1.9e-14` and `1.7e-14`. The status is `NOT_RIGOROUSLY_CERTIFIED`; the extra root is ambiguity evidence, not additional GT recovery.

The historical counts are machine-traceable: R2 gives Q20 versus 64-start root counts `2/1/1` versus `1/1/1`; the 192 starts and 24 evaluation-cap hits are directly countable from the three raw baseline files.

R3 count provenance is exact: 864 Q6 calls, 538 new legal archive roots (`118+0+183+237`), 462 material-distant roots (`103+0+172+187`), and zero direct updates of the global best trace. The primary source is `research_z_alg_r3/audit/MECHANISM_STATS.json`.

R4 count provenance is also exact: basin has 29 exposed runs plus 18 fresh predictions with zero fresh joint successes; diversity has 360 predictions from a 288-run CPU set and a 72-run CUDA subset, with no overall GlobalDual superiority.

The reusable physics is the scalar-TM VIE `A(m)=I-diag(m)G`, `J=A^{-1}diag(m)E0`, `Y=SJ`, with exact fixed-active reduction in `ReducedKernel`. Q6 is `src/quotient_solver.py`; Q20 is `src/general_quotient_solver.py`, explicitly limited to 6 active cells and 3 sources.

The critical R5 gap is source realizability. Existing code accepts point-source coordinates and constructs a fixed Hankel `E0`, or accepts a precomputed `E0`. It has no antenna/aperture interface, complex transmitter control model, coherent-mixture set, power/energy normalization, or source Gram matrix. A new source layer is required before calling a `Bplus` result physically realizable.

The checkout is sparse. `research/research_z_alg_asymmetric_validation` exists in the reference tree and Git index but is not materialized, so `physics.model_for()` cannot read its contract in the current worktree until that path is populated.

The report contains all bibliography labels L1-L24. Their identifiers are locally recorded in `bibliography.json`, but no DOI/URL was fetched. Several entries were only abstract-, metadata-, or selected-text-checked in the package; those access limits must remain attached to any later literature claim.

No experiment was run, no historical file was modified, and no web verification is claimed.
