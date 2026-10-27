# Independent G1 offline audit and extra scores

Status: **code written, NOT_RUN**. This worker read the completed replay's
seal/receipt and existing scalar score fields, then wrote only
`tools/g1_offline_extra.py` and this delegated handoff directory. No Python,
tests, NPZ evaluation, truth payload, PDE action or new acquisition was run.
The original `results/g1_replay_001` files were not edited.

## Immediate existing-scalar readback

The global seal records 24 cases x7 methods, and the receipt retains
`raw_interval_bound_comparison_status=VIOLATION_PRESENT`. Readback across the
existing scalar fields gives:

| Signed positive violation | Largest existing scalar | Case | Checkpoint |
| --- | --- | --- | --- |
| atom upper | 2.55351295663786e-15 | parent05 / repeat01 | 2949 queries |
| atom lower | 7.771561172376096e-16 | parent08 / repeat00 | 2950 queries |
| unknown-mass upper | 6.661338147750939e-16 | parent10 / repeat01 | 16 and 64 queries |

Exact source files/lines and limits are in `RAW_VIOLATION_READBACK.json`.
This is **not** an independent recomputation. In particular, full-array
likelihood-bound maxima are **NOT_RUN**, so this handoff does not invent a zero.
The observed scale is compatible with possible floating-point effects, but
scale alone is not a roundoff proof. All positive values and the original
violation status remain. No tolerance, clipping, patched bound or PASS rescue
was added. The parent owns the cause assessment and existing GATES decision.

## Tool order and data access

The new CLI defaults to:

```text
python tools/run_job.py PREFLIGHT <fresh-offline-label> python tools/g1_offline_extra.py
```

Default source is `results/g1_replay_001`; default output is the fresh,
non-overwriting `results/g1_offline_extra`. Optional `--replay` and `--output`
must remain within this campaign's results tree, and output cannot be inside
the sealed replay. The parent should use its established runtime and wrapper.
This command has **not** been run by this worker.

The first file is PLAN.json. SOURCE_ACCESS.json is initialized before source
access, and logs access intent before every NPZ open. The tool first verifies
the global seal, all 24 case seals, all 168 method NPZ/JSON pairs, joint tables
and independent reference artifacts. It does not require scientific PASS.

It then independently audits the immutable M4/reference NPZ arrays and writes
the full raw audit before opening truth. Only after that barrier and audit does
it open the fixed old source
`ROOT.parents[2]/research-w-slam/results/e2e3/finite_bank/offline/FINITE_PRIOR_TRUTHS_AND_NOISE.npz`,
with the single key `truth_index` (12 existing parent atoms). SOURCE_ACCESS.json
records the completed seal attestation and a truth-open intent before the open.
No independent_noise, clean_means, seeds, parent-ID labels, new samples or new
truth labels are read/created. Online GRID/BANK/OBSERVATIONS paths are fixed to
the same old finite-bank root; the tool reads only the declared GRID features
and BANK/OBSERVATIONS `data`.

## Raw audit outputs

- `INTERVAL_CHECKPOINT_AUDIT.csv`: every case/checkpoint's signed minimum,
  maximum, argmax atom and positive count for atom lower/upper and likelihood
  lower/upper; signed unknown-mass lower/upper and evidence-log differences.
- `interval_parent_XX_repeat_YY.npz`: every original source interval, complete
  signed-difference vector, and every positive index/value. Negative and zero
  entries are retained. JSON companions hold scalar rows.
- `INTERVAL_AUDIT_SUMMARY.json`: maxima by quantity, positive checkpoint count,
  descriptive raw status, and explicit `PARENT_REVIEW_REQUIRED_NO_PASS_RESCUE`.

Likelihood comparisons use saved finite-model log likelihoods, preserving
unknown lower=-Infinity. Posterior/unknown comparisons use the sealed reference
weights exactly as saved. There is no renormalization, clipping or roundoff
adjustment in this audit.

## Additional offline metrics

For each 24 cases x7 methods and the sealed reference, the tool computes:

1. Squared error of posterior mean position against the existing parent truth.
2. Four-cell equal-area wall material MSE of posterior mean complex chi,
   `chi_cell = support_cell * (epsinf-1+i*ell)`, against that same existing truth.
   This scores the wall ROI only. Known target intrinsic parameters are not
   presented as a recovered full target material image.
3. Full 10-dimensional mean/covariance in order
   `[4 Re chi, 4 Im chi, pose_x, pose_y]`, with full 10x10 difference from the
   reference and its Frobenius/max-entry diagnostics. The 8x2 wall-target cross
   covariance and its reference difference are also saved. These are raw
   physical coordinates, with no fitted scaling or new scalar gate.
4. Old outer TX2/3 log posterior predictive density,
   `logsumexp(log q_i - Phi_outer_i) - (64/2)*log(2*pi*.01^2)`.
   Both real/imag noise quadratures have sigma=.01; all 64 real Gaussian
   constants are included. Log weights preserve tails that displayed weights
   may underflow. These are **old exposed observations, unassimilated in the
   sealed base-TX replay**. No factor, acquisition, observation or blind-test
   claim is added.

Saved float weights are used as saved for moments; their sum and saved-log
normalizer are reported, without silently repairing the online result.
Undefined methods retain explicit statuses/NaN arrays and null scalar metrics.

`extra_parent_XX_repeat_YY.npz` holds all seven full covariance/cross-covariance
arrays and means, the reference arrays, and the existing truth used for the
offline score. JSON companions and METHOD_EXTRA_SCORES.csv provide the 168
case/method scalar rows. REFERENCE_EXTRA_SCORES.csv contains 24 separate
reference rows. EXTRA_SCORE_SUMMARY.json and final stdout give each method's
minimum, maximum and mean for the requested metrics for parent review.

## Accounting and limits

Each scored method/case and reference/case independently evaluates all 2950
saved outer-response rows, including point methods with zero weights elsewhere.
The tool records those logical offline lookups by owner; with all seven methods
defined, the expected total is 24 x8 x2950 =566400 logical rows. The one disk
BANK load of 2950 rows is separately recorded and is not added twice. No new
FE/RHS occurs. Existing historical bank and G0 qualification paid setup are
copied from the sealed receipt, without changing their attribution.

The tool writes only its new output directory, checks peak process RSS against
1,000,000,000 bytes, and creates no new scientific gate or geometry image.
Syntax, execution, full-array audit maxima, new GT/outer score ranges, resource
measurements and completeness of actual new outputs remain **NOT_RUN** until
the parent executes the metered wrapper.

## Explicit old-release relocation for public reproduction

Both `tools/g1_replay.py` and `tools/g1_offline_extra.py` now accept
`--old-bank-root <old-finite-bank-online-directory>`. For the public old Release
layout, the argument is the extracted directory ending in
`research/research_w_slam/results/e2e3/finite_bank/online`, containing the exact
filenames FINITE_GRID.npz, FINITE_BANK.npz and OBSERVATIONS.npz.

Without the flag, both tools validate and use CACHE_MAP's original absolute
paths. There is no search for a public/local replacement or new main-campaign
file. With the explicit flag, only the root changes; filenames and NPZ key
allowlists remain fixed. GRID/BANK/OBSERVATIONS additionally read the scalar
identity metadata `prior_id` and require FINITE_PRIOR_ENUM4_V1; all original
shape checks remain. Source paths reject locked/inference/posteriors components
before access. SOURCE_ACCESS.json records the explicit relocation and original
audited paths without hashing or claiming a numeric payload-equivalence proof.

The offline-only truth path is derived as the selected online directory's
sibling `offline/FINITE_PRIOR_TRUTHS_AND_NOISE.npz`. Only truth_index is read,
after the complete seal check and raw audit. The online replay never opens the
offline truth directory. The main thread's G0_DECISION launch guard in
g1_replay.py was retained unchanged.

Two additional test methods were written for explicit versus unauthorized root
replacement, forbidden paths, prior identity, correct sibling truth derivation
and rejection before a truth NPZ can open without the global seal. They remain
NOT_RUN. The completed `results/g1_replay_001` was neither rerun nor edited.
