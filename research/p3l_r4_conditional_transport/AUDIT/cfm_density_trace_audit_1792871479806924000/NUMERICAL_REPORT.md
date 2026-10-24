# CFM numerical continuous-density audit

The frozen0.05-nat32/64 agreement criterion failed on3of4 preselected exposed DEV objects. All12 full-divergence evaluations were finite. No tolerance was relaxed and no extra NFE/model search was launched.

| Parent | log q16 | log q32 | log q64 | absolute32/64 difference | frozen criterion |
|---|---:|---:|---:|---:|---|
| R4DEV_00000 |-112.216639|-112.130648|-112.245592|0.114944|FAIL|
| R4DEV_00001 |-81.862566|-81.742548|-81.780131|0.037582|within tolerance|
| R4DEV_00002 |-74.950152|-75.024828|-74.965876|0.058952|FAIL|
| R4DEV_00003 |-72.167742|-72.172533|-72.292429|0.119896|FAIL|

Log densities use the same raw Gaussianized dx66 measure with the frozen TRAIN normalizer Jacobian. Full66 divergence is deterministic, with no stochastic trace. These are numerical continuous-ODE densities, not determinants of the finite RK4 sampling map and not exact-MH proposal densities. Agreement at the frozen tolerance is NOT_ESTABLISHED; posterior qualification remains NOT_ESTABLISHED. Four already exposed DEV points do not constitute a fresh holdout or full distributional evaluation.

The selected existing CFM checkpoint remained unchanged. No training, physical state, forward, Jacobian, or RHS occurred. The independent analysis-only Book receipt records about5.07seconds wall,4.85seconds process CPU and236896256bytes peak RSS, below the600-second/2GiB caps. Frozen contract, all raw numbers and convergence/cost records are in FROZEN_AUDIT.json, TRACE_RESULTS.jsonl and RESULT.json.
