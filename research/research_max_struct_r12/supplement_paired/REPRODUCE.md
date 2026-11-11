# Reproduce the paired supplement

Source freeze: commit1166fe2256fa0dc9eec5ef7dae2c890dbdb8e9b0. Original numerical backend is the existing R11 shared checkout with R8/R10/R12 compatible sources; exact deployed bundle retained. New projection module runtime/field_current.py is an unchanged copy of R12 v2. Scientific settings and the 44 tasks are frozen under configs/.

Public input preparation (source historical worktrees must exist):

```sh
python prepare_inputs.py
```

The published inputs/ NPZ and manifest provide the resulting public-only arrays directly; regenerating them is optional. Preparation reads saved public predictions and original observation packages, not GT, and does not run Maxwell.

Actual Windows remote execution, with no new software installation:

```text
D:/python/python.exe -u D:/AI/MAX_STRUCT_R12_PAIRED_20261112/src/finite_queue.py --root D:/AI/MAX_STRUCT_R12_PAIRED_20261112 --legacy-repo D:/AI/MAX_STRUCT_R11_DEV_20261107
```

The finite queue uses four existing P-core IDs0,2,4,6, NORMAL priority and one BLAS thread per child. It creates queues/main and RESULTS once; existing results are not overwritten and automatic retries are disabled. For rerunning use a new copied output directory. The retained SSH process stays alive until the finite queue completes; a detached SSH launch did not survive and is recorded under audit/LAUNCH_ATTEMPTS.json.

Single task uses the same root, public task manifest and whole-process budget:

```text
D:/python/python.exe src/run_paired.py --root <new-supplement-copy> --legacy-repo <shared-compatible-repository> --task-id A__r10_twod_0001__dr --cpu 0
```

Record child CPU from process creation, not an offset after imports. Source m0 remains a physical incumbent. New G19 uses the existing legal least-squares projected-current warm start in the original affine DPCS/Fourier family; it resets PR and retains the original material denominator, not physical_complement. Diagnostic multiplication observes the existing projection and is charged.

Whole-phase/atomic overrun: already completed physical states are eligible only before their phase deadline and whole-child cap. Full external CPU includes output serialization after the last usable state; those export costs never purchase an earlier candidate. There is no forced final forward after a cap. Shared upstream G19 and NN preprocessing costs are disclosed separately from incremental rescue cost.

Scoring must use the root-written all44 seal. It is offline and reads GT only after all outputs are sealed. The existing24 cohort is exposed; no new blind confirmation is claimed. New candidates/prefixes are selected from public Phi/R/noise criteria before scoring.

## Offline reproduction from published files

Use a fresh external-disk scratch directory. Copy `configs/`, `inputs/`, `runtime/`, `src/`, and the postseal `offline_truth/` into it. Restore/extract `paired_results_full.zip` there; it contains all original RESULTS and queues but no GT or old scoring outputs. This leaves no pre-existing local acceptance/seal. Then run:

```sh
OPENBLAS_NUM_THREADS=1 python src/verify_and_seal.py --authorize-scoring
OPENBLAS_NUM_THREADS=1 python src/score_and_plot.py --seal ALL44_PREDICTIONS_SEALED.json --controller queues/main/RECEIPT.json --output-dir RESULTS/offline_reproduced
```

The newly generated seal stamps the freshly restored local files; do not pretend Git/ZIP extraction preserves nanosecond mtimes from the original author's disk. `offline_truth/` is explicitly for authorized offline scoring and is never read by run_paired.py. Its provenance retains original GT paths and predeclared amplitudes/support definitions. To reproduce additional historical-context plots, retain the parent R12 study with its existing public decomposition files; absence of that context does not change the new44-task scores.

Actual scorer environment was the existing Gaussian/.venv_nn Python3.13 environment on Apple Silicon; actual numerical child runs used Windows Python3.10.9/NumPy2.2.6/SciPy1.15.3, recorded in audit/REMOTE_ENVIRONMENT.txt. No scientific solver was rerun to correct the v1 ledger: v2 prioritizes completed external parent CPU, preserves the same materials/prefixes and re-renders labeled figures.
