# Maxwell implementation references — consult when relevant

See the [experiment handbook](experiment_handbook/OPERATORS_ZH.md) for callable entry points and hardware, plus [baseline errata](experiment_handbook/BASELINES_AND_ERRATA_ZH.md) for reproduced and adapted comparisons.

These are pointers to existing implementations and measured workloads, not mandatory startup reading or universal hardware/algorithm choices. Historical source files and results remain the authority for their exact assumptions. Do not rerun their benchmarks merely to choose an implementation.

| Relevant task | Existing source | Scope to preserve |
|---|---|---|
| Scalar-TM CUDA performance | research/research_z_alg_cuda_pipeline/OPTIMIZED_OPERATOR_CATALOG.md and REPRODUCE.md | Recorded RTX4060/Windows workload, not a universal GPU or GT-recovery claim. |
| Mixed 576-cell TRF/L-BFGS/GlobalDual workload | Same study: manifest_4_hybrid.json, src/hybrid_svd.py, SharedOperators, ObjectiveLane | Compatible SciPy interface; measured 4-lane hybrid SVD is a useful starting point only for this workload. |
| Repeated material changes with known host/support | Original R2 active-block and R9 fixed-host reductions | Coordinate/active-complement and rank assumptions must hold. A dense macro-region parameter is not automatically a rank-one material update. |
| Vector current/material or vector-grid FFT work | A17 compressed endpoint/Cholesky/Woodbury and A18 Route A/B, as indexed in the operator catalog | Different models from scalar-TM; require an appropriate adapter. |
| Fixed-feedback scalar-TM objective | research/research_fbg_r2_adaptive_feedback/MATHEMATICAL_THEORY.md (F02/F03/F09/F12), src/feedback_fit.py | B stays fixed during an inner solve; variable B(m) needs its derivative. Attenuation is not a basin/GT certificate. |
| CIE/Fourier-TSOM comparison for new basin claims | research/research_fbg_r3_experiments/src/cie_tsom_baseline.py, CIE_TSOM_BASELINE_AUDIT.md, audit/E3_BASELINE_PORT.md; src/cie_tsom_cuda.py when applicable | Real-material 2D adaptation, original stopping/RHS conditions, DPCS/APCS and corrector budgets. Similar parameter counts do not imply identical admissible current families. |

For a justified new CUDA task, prefer existing library primitives and compatible operators; reuse retained LU/multiple RHS, shared immutable operators, and analytic derivatives. Choose concurrency from relevant measurements. Synchronize only where dependencies require it; do not infer SM reservations or recovery gains from stream count/utilization. None of this authorizes taking another session's GPU or installing a scheduling service.

The historical hybrid timing (56.7s vs99.5s for 631 evaluations,635LU,26880 complex RHS) is a scoped engineering record. Fixed feedback and CIE studies contain both successes and failures; read those records when making a related claim. Keep CPU jobs on their existing path when CUDA is irrelevant.
