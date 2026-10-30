# Maxwell operator reuse

For Maxwell performance work, consult
`research/research_z_alg_cuda_pipeline/OPTIMIZED_OPERATOR_CATALOG.md` and
`research/research_z_alg_cuda_pipeline/REPRODUCE.md` before reimplementing an
operator or introducing a new execution backend.

- Reuse compatible fixed physics, retained LU/multiple RHS, analytic adjoints,
  exact Schur kernels and mature GPU libraries. Keep the original material and
  field/current conventions and scope assumptions.
- The new scalar-TM CUDA executor provides one shared G/E/S upload, one CUDA
  process, private stream lanes, GPU-resident forward/gradient/Gram/regularizer,
  pinned host returns, continuous FIFO handoff and independent result writing.
  Its measurements concern the stated 576-cell / 12-source / 48-receiver problem.
- For the measured mixed workload, start from `manifest_4_hybrid.json`: four
  streams plus the process-local mature cuSOLVER/SciPy SVD adapter. Small SVDs
  stay on CPU; a large SVD uses one GPU slot only when available, otherwise CPU
  immediately. The same objective/LU/RHS work took 56.7 s versus 99.5 s with
  one lane. This is scoped execution evidence; preserve SciPy compatibility.
- Prefer the measured queue/stream and batch settings to indiscriminately
  launching one CUDA process per case. Do not synchronize the whole device for
  every small operation. Preserve the necessary event dependency before a host
  optimizer reads its result.
- CUDA streams are not SM reservations. CUDA Graphs and forced cuSOLVER are
  optional experiments, not the measured default. Performance depends on the
  matrix size, precision, library, batch and optimization work.
- Reuse published measurements; do not repeat identity/error-floor checks.
  Include setup and failures when a materially different workload needs timing.
  Do not change frozen research files or infer better imaging from faster code.
- Other local vector-Maxwell/FFT assets are catalogued with source paths and
  interface restrictions. Use an appropriate adapter rather than changing the
  physical problem. This guidance does not authorize new dependencies, driver
  changes, taking another session's GPU, or starting a permanent service.
