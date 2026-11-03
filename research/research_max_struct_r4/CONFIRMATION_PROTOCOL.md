# Six-object R4 confirmation — preserved preparation protocol

Execution note, 2026-11-03: the preparation text below is historical. The finalized configuration is `configs/confirmation_frozen.json`; all 18 tasks completed and were sealed before scoring. Results are in `results/CONFIRMATION_SCORES.md`, with execution/GT boundaries in `audit/confirmation/CONFIRMATION_COMPARISON_AUDIT.md`. Each method achieved 4/6. Later development diagnostics do not modify this denominator. The original draft configuration and source snapshots are retained.

Status: `PREPARED_NOT_GENERATED_NOT_LAUNCHED`. The immutable draft is
`configs/confirmation_draft_frozen.json`; the anonymous 18-slot task template is
`configs/confirmation_manifest_draft.json`. No fresh observation, truth, task ID,
remote package or numerical result has been created by this preparation.

Parent must first finish the development review and freeze the selected onehop
algorithm. Copy the draft to `configs/confirmation_frozen.json`, retain
the draft, set `algorithm_finalized=true`, `generation_authorized=true`, and
`status=FROZEN_BEFORE_GENERATION`. Freeze the actual finalized onehop source with
the configuration before generating. Do not alter scientific settings after
seeing fresh data. The local generation audit preserves the full frozen config;
Before any realization, generation snapshots the exact runtime sources plus
the two generator sources into `audit/confirmation/executed_source_snapshot`.
Packaging refuses config changes and runtime-source byte differences relative
to that pre-generation snapshot. There are no hashes.

## Independent observations and isolation

`experiments/generate_confirmation.py --execute-after-freeze` is the later,
explicit generator entry point. It imports the unchanged R3 `shape_alpha` only.
There are two newly sampled objects from each of `ring`, `curved_bar`, and
`two_islands`, not repeated noise trials on one material. Each has an independent
128-bit entropy seed and independent geometry/noise draw; UUID4 observation IDs
are independent of seeds/families. Family order is randomly permuted; public
records are sorted by opaque ID. No rejection/regeneration is allowed to improve
a material or an outcome. Write-once directories preserve partial failures too.

Reuse the original 24×24 disk-pulse scalar-TM G/E0/S, 12 sources, 48 receivers,
complex128. No operator reassembly, new backend, CUDA, installation, driver or
context change. Material is `.1 + 12*alpha`, box `[.02,20]`; alpha is the existing
192×192-to-24×24 area average. Noise is independent circular CN with
`sigma=.001*||Yclean||/sqrt(576)` and real/imag variance `sigma²/2`; each public
file retains sigma and the 99% chi-square whitened radius.

Public output goes to `data_public/confirmation/`: one shared `operators.npz`,
six observation NPZ files, and `manifest.json`. Private material, alpha, support,
clean Y, realized noise, entropy seeds, and the ID-to-family mapping go only to
`sealed_truth/confirmation/`. The generator writes
`audit/confirmation/GENERATION_RECEIPT.json`, including all six physical LU/72
source RHS, setup/serialization times and any failed attempt. Truth remains
closed until all 18 predictions or explicit failed/timeout dispositions are
sealed and the parent explicitly authorizes independent scoring. No scoring
entry point is supplied in this work package.

## Frozen methods and physical selection

`src/confirmation_solver.py` exposes exactly:

```python
solve(public, method, config, backend='cpu') -> dict
```

Default confirmation methods are `onehop`, `global`, and `native_extended`.
The `block4_cached` compatibility method is also accepted by the adapter but is
not included in the default 18-task manifest. `public` is
the merged exact R2-style allowlisted operator/observation payload; unexpected
keys are refused. No files, generator, labels, truth, or scoring are read by the
adapter. All begin at public host .1 and use the unchanged real material prior
(smooth/TV weights 1e-4, width19.5, epsilon .01) and box [.02,20].

- `onehop`: unchanged `cie_onehop.run_onehop_cie`, diagonal self absorption
  followed by one fixed offdiagonal four-cell residual action; no small LU or
  block inverse. Same mf3/6/10 three-by-500 schedule and two-beta policy as below;
  exact repeated material value/gradient caching enabled.
- `block4_cached` (compatibility only): existing cached 2×2 spatial-cell block absorption; actual
  Green submatrices including self, same fixed alignment, exact-material/beta
  inverse caching and material value/gradient caching. mf3/6/10 each 500 updates;
  first beta1, subsequent beta_current1/beta_material.2.
- `global`: unchanged R2 `cie_global.run_cie`, global PR, the same 1500 schedule.
- `native_extended`: unchanged mature sourcewise native CIE, the same first
  three stages plus mf12 for 6500 updates, beta_current1/beta_material.2;
  8000 updates total. This deliberately has more auxiliary capacity/compute;
  no equal-compute claim is made.

DPCS rtol .001, coefficient_tol0, miniter50, material inner maxiter60,
ftol1e-12/gtol1e-8 remain fixed. All auxiliary candidates are common host plus
true stage endpoints; no internal 500-update checkpoint enters selection.
Re-evaluate only by original normalized physical Phi including the original
prior, then apply the unchanged `cooperative_paths.refine` once to the selected
seed (maxiter1000, maxfun4000, original Phi scaled1e6). Retain the minimum Phi of
host, selected seed, and refined candidate. There is no scientific wall limit,
new stopping rule, GT-driven root selection, derivative suite, or parameter
search.

The return dictionary keeps compact receipts, stage-end history, candidate
scores, cost breakdowns and `trajectory_arrays`. The remote worker removes that
last key and writes its arrays to NPZ: final `m`, `Y_prediction`, `j`, `Q`, all
candidate materials/Phi, pre-corrector seed, columnar complete numeric/string
auxiliary histories and corrector history. JSON does not repeat full trajectories.
`result.status=SEALED_PREDICTION` is the scientific success disposition.
`result.status=FAILED` retains charged costs/error inside the worker's outer
`status=OK` transport receipt; outer OK alone is not a sealed prediction. A
process kill has the existing outer `HARD_TIMEOUT` and incomplete costs.

## Execution, cost and memory

The existing unchanged Windows R4-diversity `remote_runner.py` finite FIFO queue
will use four CPU workers, one thread on one logical sibling per physical P
core. All 18 tasks are CPU, so it creates no CUDA workers or CUDA context.
Its legacy topology guard still requires six physical P cores and its summary
contains fixed `cuda_workers=2` metadata even though zero CUDA tasks/workers are
launched; use `worker_affinity.json`, task backends and the public manifest as
the actual execution record. No global affinity or permanent service is changed.
Planned isolated root: `D:/AI/MAX_STRUCT_R4_20261103`.

Each job has an emergency process cap of 900 s only. A cap hit is a failed,
incomplete disposition, not scientific convergence or a reason to expand the
budget automatically. The existing runner also has its original 120 s worker
startup/assignment-stall guard. No other session or GPU allocation is touched.

Report setup/model cache construction, initial physical evaluation, auxiliary,
physical reranking, corrector and final state costs. Every task pays its receiver
SVD/DPCS/source-metric and Fourier setup; no free cross-task source cache is
claimed. Existing auxiliary counters record cache hits/misses and separate
4×4 factorization/identity-inverse RHS, small block matrix applies, and reduced
products. Onehop instead charges its actual scalar divisions and forward/adjoint small
block matrix RHS, with zero small LU; these are not full physical source RHS. Physical LU/source_RHS/
adjoint_RHS/actual_complex_RHS remain separate. Public loads/NPZ output are outside adapter timing and included in remote
worker totals. Solver/dependency import wall/CPU is separately recorded and
charged once per worker. Earlier remote bootstrap imports are not measured by
the unchanged runner and must be reported as an unavailable setup component;
do not claim complete process startup timing from the worker task timer. Historical operator assembly
is not rerun or priced as new work. Process peak working set/RSS is reported
when available, with worker-lifetime scope; it is not a per-job isolated maximum.
Native exceptions may lack a completed auxiliary operation counter; mark that
cost incomplete, while physical counters remain charged. Failed attempts are
retained. Nothing is silently retried.

## Public-only packaging and later parent handoff

After authorized generation, prepare a new local output directory with:

```sh
python experiments/package_confirmation.py --output /ABSOLUTE/NEW/PUBLIC_PACKAGE
```

The script uses explicit source paths listed in
`CONFIRMATION_FILELIST.md`, checks exact public NPZ keys, copies the existing
remote runner plus same-folder threadpoolctl vendor/license metadata, and writes
`PACKAGE_FILELIST.json` with actual files and byte sizes. It never copies truth,
scoring, generator, old observations/results, credentials, environments, or
whole repository directories. It binds six opaque IDs to three methods only
at this stage; none is invented in the draft. It derives per-task configs with
no private family/seed/label values. Execution remains
`execution_authorized=false` for the parent's final manifest review. Parent
will upload/review/authorize/launch separately; preparation performs none of
those actions. Launch later from the isolated package root using the genuine
existing remote Python environment; no installation step is present.

After all frozen dispositions are sealed, independent parent scoring retains
RMSE/12, reconstruction-only Otsu256 IoU, whitened data radius, original gates,
case-level denominators and failure statuses. Six materials are six independent
objects; 18 methods/tasks are not 18 independent materials. Claim judgment stays
with the parent Codex thread. No research result or novelty is claimed here.
