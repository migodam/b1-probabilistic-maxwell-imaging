# Frozen confirmation inference benchmark protocol

State: **PREPARED; TIMINGS NOT RUN**. Execute only after the top-level orchestrator verifies that physical work has ended. The harness requires `--physics-ended` and holds the existing physical-worker lock throughout measurement. It never imports a physical solver, trains a model, or materializes new truth.

Canonical environment is CPU with one PyTorch intra-op and one inter-op thread and one numerical-library thread. MPS/CUDA availability is recorded only; unavailable MPS is `NOT_AVAILABLE`. CPU remains canonical.

Four arms use frozen seed 2026102601: StrongCNN, optimized original K32, deterministic matched Q64, and distilled width128/n32768. Other two training seeds have functional replay checks already; they are not timing replicates. Batch sizes are 1 and 64. Identical observation, residual, q0, chi0, J, and parent IDs are used for every arm. Default inputs are the first64 old VAL parents, read without materializing truth. `--parents-manifest configs/parents_A.json` selects the first64 new online contexts through a strict field allowlist and verifies raw/context observation equality. It never opens `raw/offline_only/`.

Each batch/arm receives at least 30 independent Action constructors, timing checkpoint load and the first full inference separately, plus their sum. These are fresh Action loads in a warm process with ordinary OS file caches; there is no claim of purged-disk or process cold start. Arm order rotates between rounds. Steady full-predict and provider-only timing use at least 30 repeats after an explicit warmup; provider/full measurement order alternates. StrongCNN has no separate moment provider and its provider cost is `NOT_APPLICABLE`, never fabricated as zero.

Every timing series retains raw seconds and p10/p50/p90 per whole batch, with p50 separately divided by batch size for per-parent reporting. Repeats describe runtime variability on the same data, not training-seed or statistical scientific independence. First-inference measurements use 30 fresh model loads; steady measurements reuse loaded actions. The K32 implementation computes the original Gaussian encoder once per batch and preserves all four chunk8 local Gaussian draws per parent. Exact CPU predictions versus the unoptimized historical route already passed on all VAL256 parents.

Full prediction includes input checks, moment generation, frozen online feature construction, learned head, and output packing. Provider-only includes four-channel moments and all K32 draws/decoding. Acquisition, Native96/Probit64 anchor construction, I/O, model imports, and physical solves are outside steady timing. Read-only historical teacher/training costs remain separate. Full physics and end-to-end aggregated parent runtime must be reported separately using the physical ledger and parent receipts.

Memory reports are cumulative process peak RSS, including the raw per-load cumulative values; they are not falsely attributed as isolated arm allocation peaks. Environment, source regime, backend availability, canonical seed, repeat count, and parent IDs accompany results.

Prepared invocation after physical completion (from repository root):

```text
../../.venv_nn/bin/python research/p3l_r5_b_independent_confirmation/src/benchmark_confirmation.py --physics-ended
```

Optional new online context invocation appends `--parents-manifest research/p3l_r5_b_independent_confirmation/configs/parents_A.json`. The output defaults to `audit/timing/CONFIRMATION_TIMING.json` under the independent confirmation root. No timing has been executed while preparing this protocol.
