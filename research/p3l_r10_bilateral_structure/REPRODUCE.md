# Reproduce the R10 evidence without losing its history

This is the executed local evidence tree. The commands below are documentation for a future replay; no new physics was run to write this file. For immediate verification, read the sealed arrays, canonical statistics and independent QA artifacts. Re-running physics is unnecessary for checking the reported results.

## Runtime and fixed paths

Use the declared interpreter as written; do not resolve its symlink to another environment:

```sh
R10_PY="/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python"
R10_WORKTREE="/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r10-bilateral-structure"
R10_STUDY="$R10_WORKTREE/research/p3l_r10_bilateral_structure"
cd "$R10_WORKTREE"
```

Readback at documentation time: Python3.13.13, NumPy2.5.3, SciPy1.18.1, threadpoolctl3.7.0, Matplotlib3.11.1. The scripts also use Python standard-library modules. CPU, float64/complex128, one physical worker and one BLAS thread are the executed configuration. Entry scripts set OPENBLAS/OMP/MKL/VECLIB/NUMEXPR thread limits before numerical imports; candidate physics additionally uses threadpoolctl. No CUDA, MPS or new service is required.

Historical assets stay read-only at their original sibling paths:

| Dependency | Actual path under `Gaussian/B1/` | Purpose |
|---|---|---|
| P1 decoder | `b1-p3l-r4-conditional-transport/src/p1_material/prior.py` | Gaussian64, Uniform2, quintic occupancy and exact derivatives |
| Target | `b1-p3l-r4-conditional-transport/research/p3l_r4_conditional_transport/configs/TARGET.json` | Original material/noise/acquisition contract |
| Background | `b1-p3l-r4-conditional-transport/research/p3l_r4_conditional_transport/inputs/background/weak.npz` |576 points and fixed original acquisition |
| Exact ROI backend | `b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/src/possibility_volume.py` | Exact nonlinear256 ROI forward and full66 real Jacobian |
| Host cache | `b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/results/stable_roi/SETUP.npz` | C, ER, O, Fh and reconstruction arrays;264 historical RHS576,1 historical LU576 |
| DEV observations and anchors | R9 `results/possibility_volume/inputs/online/` and `parents/R9V_FRESH_000..007/anchor/ANCHOR.npz` | Exposed DEV8; shared already-paid anchors/J |
| DEV truth, scorer only | R9 `private/possibility_volume/truth/` | Never opened by online candidate kernels |

Independent readback of the frozen R9 dependencies is in `audit/independent_qa/DEPENDENCY_FREEZE_READBACK.json`; every listed dependency matched its R9 freeze. Source commit and R9 release provenance are in `EXPERIMENT_CONTRACT.json`. That commit predates some R10 generated files; do not imply the entire completed artifact existed at that commit.

## Evidence layout and archive status

The preserved raw tree is the canonical local evidence. The historical dependency archive now exists at `packaging/R10_dependency_inputs.tar.zst`, with its extracted sibling-layout tree at `packaging/dependencies/`. Manifest and restore notes are inside that directory; literal extraction/readback PASS is recorded in `audit/DEPENDENCY_BUNDLE_READBACK.json`. R10 DEV/fresh evidence archives are packaged separately by the root thread; this document does not invent their filenames or publication state.

```text
config/                         executed frozen protocols, enables, statistics plan
src/                            actual entrypoints and retained failed transfer source
results/transfer/               original attempt, V1, revision2, nodes, arrays, cumulative ledger
results/candidates/             DEV384 records, exact attempts/states, prediction seal, scores
results/candidates_fresh10/     generation/online ledgers,480 records, seals, scores, canonical statistics copies
results/reduced_response/       exposed DEV rank32/64/96/128 compression probe
results/reduced_response_fresh10/ heldout-shape320-node fixed64 kernel probe
results/reduced_response_probability_workload/ DEV836-node actual conditional workload replay
private/candidate_fresh10/      separate original-prior generated truth and private generation states
figures/                       rendered scientific figures; visual-QA receipts under audit/
audit/independent_qa/           independent arrays/ledger/gate/statistical audits and producers
audit/candidate_fresh_prephysics_bookkeeping_fix/ original source/config and empty failed ledger
audit/candidate_fresh_online_adapter_fix/ original source/config and charged interrupted-run receipt
theory/                        verified literature records and pre-run physical review
```

The initial fresh protocol draft is preserved at `audit/candidate_fresh10_protocol_pre_review_DRAFT.json`; it is not an executed scientific protocol. `audit/FRESH_GEOMETRY_PROVENANCE_CORRECTION.json` corrects imported legacy text: fresh geometry diagnostics use newly paid fresh anchor J, not an old R9 anchor J. Source and original notices remain intact.

Absolute paths occur in sealed endpoint records and protocols. For read-only restoration, restore the declared paths or use an explicit, audited path-prefix map in a copy of the metadata. Preserve original metadata. New runs in the staged sibling tree below naturally generate new absolute endpoint paths. Do not silently let relocated scoring read a different run through old absolute paths.

## Read-only verification and canonical statistics

Start with `STATISTICAL_ANALYSIS.md`, `INDEPENDENT_QA.md`, and the JSON/CSV artifacts; no physics is needed. Canonical independent bootstrap and frozen plan are copied byte-for-byte to:

- `results/candidates_fresh10/CANONICAL_PAIRED_BOOTSTRAP.json`
- `results/candidates_fresh10/CANONICAL_STATISTICAL_PLAN.json`
- `results/candidates_fresh10/CANONICAL_STATISTICS_PROVENANCE.json`

The statistics-only producer reads independently rescored parent rows, uses the frozen20,000 shared paired draws/seed and writes a separate replay by default. It neither solves physics nor opens GT:

```sh
"$R10_PY" "$R10_STUDY/audit/independent_qa/recompute_bootstrap.py" --output /tmp/r10-bootstrap-replay.json
```

`audit/independent_qa/recompute.py` is also offline, but it opens truth only after checking the full prediction seal and GT_SCORING_OPENED, then rewrites QA output artifacts. Preserve canonical QA outputs before choosing to run that producer. The present documentation task did not rerun it or recompute an alternative bootstrap.

## Stage future physical replay in a new sibling tree

Most entrypoints infer the study directory from `__file__`; they do not accept a fictitious `--config` flag. Keep the same `Gaussian/B1/<worktree>/research/<study>` nesting so sibling R4/R9 imports continue to resolve. Do not erase original results to satisfy no-replay guards.

```sh
R10_REPLAY_WORKTREE="/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r10-bilateral-structure-replay"
R10_REPLAY="$R10_REPLAY_WORKTREE/research/p3l_r10_bilateral_structure"
mkdir -p "$R10_REPLAY"
cp -R "$R10_STUDY/src" "$R10_REPLAY/"
cp -R "$R10_STUDY/config" "$R10_REPLAY/"
```

Before any replay, clone the protocol `output` and `private_output` path fields to the new replay location, preserving seeds, algorithms, caps, source identity and the original files. For fresh protocol, `output` must be the new `results/candidates_fresh10` because the observation-only loader uses that field. No source edits are needed. An output-path change is a new replay record, not permission to raise budgets. Keep copied protocol approvals/enables and source identities consistent; a separate run receipt should distinguish replay results from this canonical run.

### Candidate DEV and fresh generation/online/scoring

Executed parameters are in `config/candidate_protocol.json` and `candidate_fresh10_protocol.json`, with matching enable JSON. DEV uses8 already exposed R9 parents, K16 and three new full66 arms. Fresh uses10 newly generated parents, six arms, K8. New arms retain12F8J; old-form arms retain6F4J, four radii prefix and latentcap4. They share K, not optimization budget or admissible correction geometry.

These are the actual action names. Run sequentially in the new tree, after cloned output paths and empty output directories are established:

```sh
"$R10_PY" "$R10_REPLAY/src/candidate_generation.py" run
"$R10_PY" "$R10_REPLAY/src/candidate_statistics.py"
"$R10_PY" "$R10_REPLAY/src/candidate_fresh.py" generate
"$R10_PY" "$R10_REPLAY/src/candidate_fresh.py" online
"$R10_PY" "$R10_REPLAY/src/candidate_fresh.py" score
"$R10_PY" "$R10_REPLAY/src/candidate_fresh_figures.py"
```

`freeze` actions exist for the two candidate engines, but deliberately refuse existing configs. Do not call freeze on an already frozen evidence tree. Generator and online are separate processes. Online reads only observation/sigma/parent plus its own anchors; private fresh truth is only opened by `score`, after the complete480-row prediction seal. Fresh role seeds are[20261030,101009,role,index]; DEV target namespace is[20261030,1000701]. The root statistical plan was frozen before fresh generation and does not select arms using DEV truth.

The corrected fresh source supports an explicit first resume phase: it loads an already paid anchor, skips completed RECORD.json files and starts distinct endpoint keys. A second unplanned replay is rejected by existing phase/output guards. Do not use this as a silent retry loop.

### Transfer V1 and revision2

`src/transfer_campaign.py run` implements full256 fixed-shape LU/complex-Schur comparisons, mode/Hessian work and initial quadrature on eight DEV shapes, with conditional work on the first three. Its config is `config/transfer_frozen.json`. The original failed source is preserved separately as `transfer_campaign_failed_v1.py`; do not reenact the failure just to reproduce its already preserved cost.

`src/transfer_revision2.py run` consumes V1 setup/conditional checkpoints and `config/transfer_revision2_frozen_v2.json`. It is a checkpoint continuation, not an independent clean-from-zero driver. It requires the V1 cumulative ledger prefix ending40608RHS256 and the preserved first-attempt receipt. The original rejected `transfer_revision2_frozen.json` is audit history, not the executed config. A new clean V1 run without the historical failure does not have that cumulative cost and must not be falsely relabeled40608.

To replay only revision2, restore the V1 checkpoint files from the raw archive into a separate replay `results/transfer/`, retain a copy of the complete original ledger, and restore its exact line prefix through the V1 `END` event at40608RHS256. Preserve first-attempt receipt and all copied revision2 outputs in an adjacent history directory; move copied revision2 output names away so the no-replay guard does not overwrite them. Then use the actual entrypoints:

```sh
"$R10_PY" "$R10_REPLAY/src/transfer_revision2.py" run
"$R10_PY" "$R10_REPLAY/src/transfer_report.py"
```

That replay reuses paid V1 shape/profile/grid construction explicitly. For a scientifically new clean-from-zero transfer study, freeze a new cumulative-cost contract rather than altering the original continuation counters. Original executed V1 conditional references remain UNQUALIFIED; revision2 has only numerical fixed-shape task qualification, no full-support/tail/rounding certificate and no Full66 posterior.

### DEV compression, heldout-shape compression, actual probability workload

For DEV compression, restore only the R9V_FRESH setup and timing-node NPZ dependencies from `results/transfer/` into the staged transfer directory, with their provenance. `reduced_response.py` imports these already formed B,L,R arrays, then pays for new full Schur, energy bounds, a basis built to128 and models32/64/96/128. Its rank64 result is DEV selection, not the independent cold-to64 validation.

For heldout-shape compression, restore or generate the fresh anchors and observation-only inputs under staged `results/candidates_fresh10/`. The source rebuilds B,L,R and a basis only to64 from those anchors; no GT or future candidate processing occurs. For the actual probability workload, restore the original `R9V_FRESH_000_integration_nodes_revision2.csv` dependency into staged `results/transfer/`. Its836level1 nodes/weights and historical profile/event-grid construction are fixed common history; no refit is performed.

```sh
"$R10_PY" "$R10_REPLAY/src/reduced_response.py" run
"$R10_PY" "$R10_REPLAY/src/reduced_response_fresh.py"
"$R10_PY" "$R10_REPLAY/src/reduced_response_fresh_report.py"
"$R10_PY" "$R10_REPLAY/src/reduced_response_probability_workload.py"
```

Corresponding configs: `reduced_response_frozen.json`, `reduced_response_fresh10_frozen.json`, `reduced_response_probability_workload_frozen.json`. Their outputs use distinct directories listed above. Fixed rank64 was selected from DEV only. The heldout study uses10 fresh anchor shapes and32new Uniform2 amplitudes per shape, seed[20261030,101012,index]. Only M32 is an actually executed complete query workload; other amortization counts are measured-mean projections. The probability replay uses one exposed DEV shape and the same836 saved quadrature nodes, not a new reference or independent material test.

Both reduced mechanisms gate on full256 residual, numerical coercivity and the observed residual/potential bound. Actual exact-reference error is reported after queries, never used to accept/tune. Exact fallback is separately reserved/charged; zero fallback occurred in these runs. The heldout source supports all-host analytic response and rank-breakdown fallback, but neither branch was exercised; no additional miniature physics test was run.

## Frozen budgets and actual costs

Hard caps are stops, not targets. Do not increase them to rescue a conclusion. Global original envelope:350000 conservative RHS,10000LU256,5400s and one worker; branch policies executed0new576 solves. Initial allocations and the root's pre-generation reallocation are preserved in `EXPERIMENT_CONTRACT.json` and `config/BUDGET_REALLOCATION_1.json`.

| Branch | Frozen relevant cap | Actual RHS256 | Actual lower-order RHS | LU256 |
|---|---|---:|---:|---:|
| Transfer including failed/V1/revision2 | cumulative70000RHS256; revision2 extra29392 |69784|0|270|
| DEV candidates |250000RHS256/7500LU/900s initial worker allocation |103824|0|4370|
| Fresh candidates |155000RHS256/5000LU/1200s reallocated branch |118760|0|4313|
| DEV compression |10000cross-order RHS/0LU/180s |1920|7680|0|
| Heldout compression |8000cross-order RHS/400LU/300s |5120|2560|320|
| DEV probability workload |20000cross-order RHS/0LU/120s |6688|6688|0|
| Actual total | original global cap unchanged |306096|16928|9273|

Conservative combined RHS=323024. Lower-order RHS remain labeled by dimension; the combined count is a conservative cap convention, not equal256 computational work. New RHS576/LU576=0. Historical shared setup264RHS576/1LU576 is separate and disclosed for cold/standalone views. Kernel speed tables include required cold shape algebra and gates, not full inversion/candidate generation. Historical profile/grid work is shared history with no separately measured pure-grid time; no end-to-end time is invented.

Fresh candidate cost reconciles100624proposal RHS+18016shared-anchor RHS+80generation RHS+40wasted-adapter RHS=118760, and3806+496+10+1=4313LU. See `CANDIDATE_COST_AND_PROVENANCE.md` and `results/candidates_fresh10/COST_RECONCILIATION.json`. Do not sum cumulative fresh JOB_END wall counters as independent phase times.

## Repair chronology and source identities

1. A DEV launch had a shell quote parse failure before Python;0physics. Corrected launch completed384 predictions, then the separate scorer opened already exposed DEV truth.
2. Fresh pre-review draft counter bookkeeping was corrected before any fresh physics; its config remains a draft archive.
3. First fresh generator launch failed at its initial ledger event because `phase` was passed twice; empty ledger,0RHS,0LU,0truth states. Original source/config and receipt are preserved under `candidate_fresh_prephysics_bookkeeping_fix/`.
4. Fresh generator then completed10private truths/observation-only files. Initial online process omitted `ph.y`, required internally by the historical corrector. It stopped at cumulative2288RHS/86LU, including generation, after parent000anchor, three completed new-arm K0 endpoints and one old raw/J pair.
5. Minimal repair assigned `ph.y`, reused the paid anchor and completed records, used new phase ONLINE_RESUME1 and replayed the interrupted old candidate. The wasted original raw/J remains an extra40RHS/1LU in the global ledger. Scientific methods, seeds, statistics and all six arms were unchanged. No raw-state cache injection or hidden extra corrector freedom was used.
6. Imported fresh geometry notices said old R9J; a companion provenance correction identifies the actual freshly paid anchor J without changing scientific source or matrices.
7. Transfer attempt1 failed on nonfinite JSON output after9552RHS/30LU; its source, receipt and arrays remain. V1 and revision2 cumulative accounting includes that failure once. Rejected revision2 preflight configuration is preserved; only the approved v2 configuration executed.

`failed_actions=0` in physical receipts means no failed physical solve. It does not erase software failures above. Canonical source identities are BLAKE2b readbacks below; no SHA check is part of this workflow. Candidate and heldout/workload configs carry their own executed source identities. Transfer and root DEV-compression rows below are documentation-time readbacks unless a separate pre-run identity is present; do not retroactively call them pre-run seals.

| Source | BLAKE2b-256 current readback | Executed config match |
|---|---|---|
|`src/candidate_generation.py`|`fa45480518960c95369ff64866673517e91a6d839f35522124be062bb6266a9f`|True|
|`src/candidate_fresh.py`|`dc47221f0ab28e21503e8cfd69e989195c38ee61508a10e1ffb81e7d202b13a5`|True|
|`src/candidate_statistics.py`|`c47105be77b69e2055d3fd4566a1b47d354635b9e8755a1d3f9dfc76d29b98b0`|documentation-time readback only|
|`src/transfer_campaign.py`|`11678a51727fb9f8a6478062bb1c3869122368b9d426e59a568414d7732e0305`|documentation-time readback only|
|`src/transfer_campaign_failed_v1.py`|`f818608cb998abf4d44535549228764025254c846c0fe5640b5028c25580c4d1`|documentation-time readback only|
|`src/transfer_revision2.py`|`b33978d02c8d9d80114dcbf96be5931dbaebf1519143bb50e370069298583eae`|documentation-time readback only|
|`src/reduced_response.py`|`eb2b96de4ca8e0b2943af1f8dc760a8fdf5367dfa6160053946efd75783f70d5`|documentation-time readback only|
|`src/reduced_response_fresh.py`|`3e27c674cd644f1b23415c5488fc8f014282e4f4b124a4ea02ae4ea9bc87cab2`|True|
|`src/reduced_response_probability_workload.py`|`bab1e4b8dd286a61d84d6c8afdc92d7f71c0c99accbfb992103a6104613011db`|True|

Preserved prior adapter source/config identities and refreeze times are in the two repair RECEIPT.json files. Original scientific parameters remained frozen; text/bookkeeping repairs after data generation are disclosed rather than backdated.

## Portable dependency package pointers

Canonical explanatory report for the836-node workload is `PROBABILITY_WORKLOAD_REPLAY.md`. Its common historical profile/Hessian/event-grid construction was reused and was not retimed; the measured saving is cold response kernel plus the fixed quadrature workload, not a new inversion or new reference.

The real portable archive is `packaging/R10_dependency_inputs.tar.zst` (192886243 bytes). Its directory `packaging/dependencies/` contains423 exact source files,194545889 payload bytes, including380 non-null historical selected endpoints. `DEPENDENCY_MANIFEST.json` gives original absolute paths, portable relative paths, sizes and BLAKE2b-256 identities; `NOTES.md` gives restoration instructions. Source-copy and full archive-extraction literal comparisons both passed; see `audit/DEPENDENCY_BUNDLE_READBACK.json`. The authoritative historical freeze readback remains `audit/independent_qa/DEPENDENCY_FREEZE_READBACK.json`. The actual package mirrors:

```text
historical/Gaussian/B1/
  b1-p3l-r4-conditional-transport/
    src/p1_material/prior.py
    research/p3l_r4_conditional_transport/configs/TARGET.json
    research/p3l_r4_conditional_transport/inputs/background/weak.npz
  b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/
    src/{possibility_volume.py,possibility_volume_generate.py,stable_roi_operator.py,possibility_statistics.py}
    config/{possibility_volume_protocol.json,possibility_statistics_frozen.json}
    results/stable_roi/SETUP.npz
    results/possibility_volume/inputs/online/R9V_FRESH_000..007.npz
    results/possibility_volume/parents/R9V_FRESH_000..007/{anchor/ANCHOR.npz,COMPLETE.json,...}
    private/possibility_volume/truth/R9V_FRESH_000..007.npz
```

The ellipsis is a concrete reference closure, not the whole R9 tree: include every non-null `selected_path` referenced by those eight R9 COMPLETE.json files for the historical DEV scoring comparison. Keep each referenced file at its mirrored relative path. The historical truth files are needed only to reproduce DEV scoring, not online physics. The R9 anchor/observation files support DEV candidate runs and the DEV probability workload; fresh generation and heldout kernels only require the original decoder/target/background/cache plus R10-generated anchors/observations. R7 config files mentioned in the seed-audit receipt are audit history, not runtime physics dependencies.

The final evidence manifest can point to this archive and its existing literal readback receipt. No public URL is asserted by local archive creation. Restore with `zstd -d -c packaging/R10_dependency_inputs.tar.zst | tar -xf - -C <restore-directory>`; the restore directory must already exist. Restoring the mirrored historical prefix requires either the original absolute root or an audited path-prefix mapping/sibling layout in a replay copy. Original frozen inputs were not rewritten. Packaging is Git-ignored; the original historical and R10 files remain in place.
