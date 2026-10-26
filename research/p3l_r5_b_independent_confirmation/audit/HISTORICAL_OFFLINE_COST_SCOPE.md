# Historical offline cost scope

Metadata-only reconstruction; this audit adds **0 RHS, 0 training, 0 truth scoring**. No large ledger scan, checkpoint loading, array payload reading, or new benchmark. Exact source paths and individual records are in [HISTORICAL_OFFLINE_COST_SCOPE.json](HISTORICAL_OFFLINE_COST_SCOPE.json).

| Historical component | Recovered cost | Attribution |
|---|---|---|
| Original joint TRAIN32768 + VAL2048 truth simulation | **278,528 RHS; 34,816 LU; 400.392 s** sum of whole-job controlled wall | Three R2 `A/jobs/data_*/CLOSED.json` plus DATA_MANIFEST; shared dataset later reused by R4. Excludes historical test panels. |
| Strong Native96+Probit64 TRAIN1024 contexts | **3,913,008 RHS; 2,835.898 s** sum of per-parent wall | R4 `results/contexts/train/*.json`; same R2 parents, so do not charge another truth generation. |
| Strong Native96+Probit64 VAL256 contexts | **1,015,096 RHS; 740.835 s** sum of per-parent wall | R4 `results/contexts/validation/*.json`; 4,928,104 RHS combined with TRAIN. |
| Gaussian teacher provider actually used by R5 | **52.423 s** training operation; **53.700 s** enclosing MPS job wall; **35.841 s** operation CPU | R4 `models/density/extended_diag_gaussian_s2026102501/RESULT.json` and matching `I/jobs/.../CLOSED.json`; 682,943 parameters, 6000 updates. These times overlap; do not add operation and whole-job values. |
| Original K32 matched heads, three original seeds | **361.423 s** sum of whole-job controlled CPU wall | R4 `I/jobs/riskfusion_diag_gaussian_s2026102501..3/CLOSED.json`; separate from teacher provider training. |
| Original K32 moment features TRAIN/VAL | **4.074 + 2.049 = 6.123 s** whole-job wall | R4 `I/jobs/fusion_diag_gaussian_features_train/validation/CLOSED.json`; shared feature cache. |
| R5 deterministic matched heads, three seeds | **50.989 s** sum of operation wall | `stage_b_material_risk/COST_LEDGER.json`, `experiments/deterministic_matched_*` only. |
| R5 cached teacher generation n32768 | **29.153 s** CPU wall, **0 RHS**, no truth labels read | `TEACHER_GENERATION_32768.json`; generated once and reused by width64 and width128 three-seed students. Do not charge once per seed. |
| Final width128 student providers, three seeds | **32.383 s** sum training operation time | Three `learned_context/distilled_w128_n32768_s*/RESULT.json`; providers train on32768. |
| Final width128 matched heads, three seeds | **50.373 s** sum training operation time | Same three RESULTs; heads TRAIN1024, VAL256. |
| Final width128 provider+head total | **82.756 s** sum recorded MPS operation wall; **84.771 s** enclosing whole-job wall | `DISTILL_W128_JOB_CLOSED_n32768.json`; totals overlap, not additive. |

Final seed pairs (provider / head seconds) are 2026102601: 11.174 / 16.998; 2026102602: 11.077 / 15.589; 2026102603: 10.132 / 17.785. `GPU_seconds` in these historical records measures wall while using MPS, not measured GPU utilization or energy.

Shared sunk costs above must be shown separately from deployment: a new parent still requires the frozen full-fidelity anchor acquisition, followed by model action inference. Historical `COST_BENCHMARK.json` and `FULL_COST_BENCHMARK.json` exclude paid anchor acquisition and interpreter imports by their own scope. Their online timings cannot alone be reported as end-to-end timing.

`TEACHER_AMORTIZATION.json` reports approximately90,032 parents for **29.153 s teacher-generation-only** amortization using its historical provider-only batch64 saving. It explicitly excludes original Gaussian training, joint data/physics, student/head training, anchor acquisition, power/hardware differences. It is not a total offline or end-to-end break-even. This audit makes no performance judgment.

Unrecovered or non-attributable cost is **not free**: original known-background acquisition; joint repacking/normalization overhead outside recovered records; architecture searches, failed runs and earlier fit versions; standalone imports, packaging, and monetary/energy normalization. Do not allocate the entire R4 campaign to this selected teacher, and do not silently omit shared research cost when reporting a full development total. Historical R5 COST_LEDGER retains whole-campaign460.019 s fit operations and556.488 s recorded operations, plus retained resource failures; these broad totals are not the final width128-only cost.

Source guard review: `src/physics_confirm.py:read_online` now checks exact allowlist,32×8 shape, finite complex observation, exact original P1 target_id, and fixed sigma. All14 source files listed by the original confirmation FREEZE currently match that manifest's bytes and existing BLAKE2B identities; only small source files were checked, no checkpoint/large-array hashes or SHA256. See [POST_GENERATION_SOURCE_GUARD_REVIEW.json](independence/POST_GENERATION_SOURCE_GUARD_REVIEW.json). This is an audit snapshot; no code was changed by this worker.
