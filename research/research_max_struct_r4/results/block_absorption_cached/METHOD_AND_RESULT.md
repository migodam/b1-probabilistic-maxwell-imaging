# Cached inverse implementation replay

Exactly one authorized development replay was completed. No truth payloads or material scores were read; no numerical derivative suite, parameter tuning, repeated timing experiment, publication or commit was performed. Previous direct-block and self-only results and executed-source snapshots are preserved.

## Implementation scope

`block_solver='cached_inverse'` is explicit; the default remains `direct`. For a cache miss, the code calls the existing NumPy batched `solve(K,I4)` on each actual 4×4 block. It retains B,K,inverse and a copy of m in an LRU containing at most 2 beta entries, one latest exact material per beta. Hits require equal beta and `np.array_equal(m)`. Residual/adjoint and H/T applications use inverse matrix products, with the adjoint applying inverse^H. No new factorization algorithm, LU kernel, precision, spatial block, or physics backend was introduced.

The optional value/gradient memo stores at most 2 exact material points within each fixed-J/Q readout and resets before the next readout. It only reuses exact repeated calls; gradients are copied when returned from cache. The mathematical derivative function is AST-identical to the original executed source. Seed, prior, bounds, source metrics, tolerances, 1500 outer updates and max 1000 original-Phi corrector remain unchanged. No explicit rounding or quantization is added. Inverse multiplication changes floating-point operation order, so the trajectory and endpoints are retained separately.

## Measured physical outcome and timing

| Solver | Original Phi including prior | Relative physical residual | Auxiliary wall seconds | Process wall seconds |
|---|---:|---:|---:|---:|
| direct | 6.50627883468e-06 | 0.000967014434494 | 43.246437 | 59.423727 |
| cached_inverse | 6.49499025009e-06 | 0.000969536827283 | 35.428370 | 46.374740 |

Observed auxiliary wall reduction was 18.078%; process wall reduction was 21.959%. This is one replay, not a general speed guarantee. In particular, the unchanged physical corrector itself changed from 15.398363s/1026 evaluations to 10.322086s/1019 evaluations. Consequently the entire time reduction cannot be attributed to cached auxiliary kernels. The old process also included its 0.007975s wiring check, which was not repeated.

The cached final full linear-system relative residual is 1.97427331829e-15. The final prediction differs from direct by relative material L2=0.016447170931 and max absolute material difference=0.335748722735; these compare two predictions, not truth. Both final correctors reached 1000 iterations, all three outer stages reached 500, and the cached run had 597 inner solves reach 60 iterations. The physical-fit scale is retained, but no GT-recovery equivalence or superiority is asserted.

## Actual costs

Direct auxiliary work used 18342720 small4×4 factorizations and 261872640 factor-solve RHS. Cached work used 8366544 factorizations/solves with exactly4 RHS per block, totaling 33466176 inverse-solve RHS. Its additional 252054144 matrix-apply RHS are separately counted and are not triangular solves. The inverse cache had 4998 exact hits, 58101 misses, and max resident entries=2. Value/gradient memoization served 3000 exact repeat hits; actual material evaluations=58599, inner L-BFGS iterations=53627. No inner fallback occurred.

The new full run charged 1026 physical Maxwell LU and 24540 physical complex RHS. Process CPU=43.600083s; auxiliary CPU=33.190796s. Startup succeeded; failed_attempts=[]; frozen operator assembly was not repeated. Snapshot and summary postprocessing have separate nonphysics receipts. Prior original costs remain recorded, rather than overwritten by the faster observation.

`RECEIPT.json`, `AUXILIARY_EXECUTION.json`, `AUXILIARY.json`, `CORRECTOR.json`, `FROZEN_CONFIG.json`, `STARTUP.json`, `RUN.log`, three stage NPZs, candidate NPZs, `prediction.npz`, and `executed_source_snapshot/` preserve the run. `COMPARISON.json` includes the configuration/source isolation evidence and all measured counts. Use the scientific Python with CPU threads 1 to run `experiments/block_absorption.py --variant block_absorption_cached` in a fresh output directory; the runner refuses overwrite.
