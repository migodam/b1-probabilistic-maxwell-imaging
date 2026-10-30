# Cost accounting and limitations

`COST_LEDGER.csv` is generated once from completed run manifests and diagnostic receipts. It contains 47 actual solver outputs, including failed and truncated attempts; old frozen R3 work is not billed again as new computation. `COST_TOTALS.json` contains sums and final fresh-arm totals.

Known exposed counts: 54,027 full576 LU; 1,377,960 physical-LU/kernel complex RHS; 16 full576 complex Schur decompositions; 23,328 full576 triangular transfer RHS; 1,547 Full-J builds; 3,770 source Gram factorizations; 673 exact Hessian-vector actions; 36 real material-Jacobian SVDs. Solver numerical wall sums to 729.47 seconds and CPU to 710.00 seconds. These are sums of concurrent jobs, not sprint elapsed time.

There are no unreported zero-cost formula experiments. The diagnostic geometry run, sufficient-certificate conservatism experiment, Hessian attempts, dense-face representation change, factor generation, fresh generation and failed font process are separate ledger rows. Two capped ARPACK attempts retain all paid HVP/RHS. The missing font failed before any physical operation or saved target.

The physical-LU/kernel RHS counter includes the N+ns dense rational kernel setup. The transfer evaluations use a precomputed dense Schur representation; each back substitution still has dimension576. Their RHS are separately listed and their time is included, rather than called free small reduced operations. The max_rhs guard applies to the physical-LU/kernel counter; transfer back substitutions are controlled by actual wall. Equal caps therefore are not an equal-total-flops certificate.

LU and source/adjoint/receiver RHS counts come from the single common frozen physics backend. LAPACK inverse-condition estimation uses opaque internal work, which is not assigned zero RHS. SciPy TRF's internal dense algebra is included in measured wall/CPU but its private factorization calls are not individually exposed. No invented count substitutes for those unknowns.

Solver timers include mathematical preparation, forward, derivatives, Schur, optimization and failed trials. Per-job JSON/NPZ serialization, process startup and some scoring/plotting do not have separate captured CPU timers. Batch start/completion UTC is preserved. Known costs are therefore not a claim of complete operating-system accounting. Model/API monetary usage was not available; explicit model/effort requests and the capacity failure are recorded.

All final methods share a 60-second numerical cap, 2500 LU cap, 60000 exposed physical-LU/kernel RHS cap, and the same inputs, material box and regularizer. Many jobs terminate or hit a physical cap before60s. They are not artificially padded or granted fresh budget after seeing GT. Resource-prefix quality plots select by online physical objective before retrospective GT scoring.

Rebuild ledger without any physical solve:

```bash
python aggregate_costs.py
```
