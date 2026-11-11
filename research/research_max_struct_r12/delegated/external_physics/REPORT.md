# Austria coarse-target and saved DR-shadow diagnostic

Actual reused acquisition metadata is forward150²/inverse32² at400MHz. Earlier authored64²-generation labels were metadata errors; code always reused this actual PUBLIC_DATA. Frozen historical contracts are preserved with explicit correction. Existing R8 grid-adequacy record concerns epsilon3.5 and cannot answer epsilon4. No saved cleanY exists in the inspected public data, so exactly one authorized offlineGT32² forward was performed using the existing backend, with costs below. No inversion, radius change, fine-grid recomputation or noise reconstruction occurred.

The actual coarse32² GT (chi=epsilon−1=3) has residual **0.225256**, compared with the declared.12 ball: **outside**. This target therefore fails the current ball on the coarse model. It does not prove every coarse material infeasible: different fitted Dual materials reach residualbelow.12. This confounds any reading of noise-feasible fit as correct material and prevents assigning the entire failure to an optimizer.

| DR step | Relative A–C gap | State norm / incident | A data residual | C data residual | Physical candidate residual | C propagation defect / incident |
|---|---:|---:|---:|---:|---:|---:|
|50|0.127820|1.172860|0.120000|0.851157|0.817670|0.766424|
|500|0.025018|1.225870|0.120000|0.343287|0.233462|0.274051|
|1000|0.019980|1.228262|0.120000|0.248380|0.231668|0.190535|
|5000|0.013352|1.220175|0.120000|0.227310|0.180285|0.178996|

A shadow satisfies propagation and the data ball, while C shadow satisfies the shared real constitutive constraint. A small physical material residual, the C-shadow data residual and a vanishing A–C shadow gap are different statements. The saved finite timeline determines whether this actual run reached a common intersection; it is not a certificate for all possible paths or unlimited iterations.

New offlineGT cost: {'cached_operator_imports': 0, 'FULL_FORWARD': 1, 'actual_complex_RHS': 20, 'LU': 1, 'source_RHS': 20, 'full_factorization_dimension': 1024}. Four saved-state timeline measurements add8Green/8receiver block products,160source columns each; no extra Maxwell factorization. Whole diagnostic CPU 0.537144s includes operator setup, the sole solve, saved-array measurements and output. GTprediction file is labeled OFFLINE_GT and never enters any algorithm pool.
