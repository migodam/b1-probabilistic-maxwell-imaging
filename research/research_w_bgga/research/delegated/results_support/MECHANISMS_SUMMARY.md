# G0 mechanism figure code handoff

`tools/g0_mechanism_figures.py`: **CODE_READY_NOT_EXECUTED**. Syntax execution, rendering, model runs and physical solves by this worker: **NOT_RUN**.

Default output directory: `figures/g0_mechanisms/`.

- PNG/PDF stems: `OPERATOR_ERROR_VS_QUERY_CPU`, `OPERATOR_SHARED_SETUP`, `STRESS_SCHUR_AND_RETENTION_POSE_ONLY`, `STRESS_SCHUR_AND_RETENTION_POSE_INTRINSIC`, `STRESS_SCHUR_AND_RETENTION_POSE_INTRINSIC_INSTRUMENT`, `FINITE15_WALL_VARIANCE`, `FINITE_MEAN_JOINT_VS_FIXED_J0`, `C4_FUTURE25_AND_INDEPENDENT_ARGMAX`.
- CSVs: `OPERATOR_QUERY_ERROR_CPU.csv`, `OPERATOR_SHARED_SETUP.csv`, `OPERATOR_EQUATIONS.csv`, `CORE_STABILITY.csv`, `STRESS_SPECTRA.csv`, `STRESS_STATUS.csv`, `FINITE15_WALL_VARIANCE.csv`, `FINITE_MEAN_JOINT_AND_FIXED_J0.csv`, `FINITE_PAIR_COMPONENTS.csv`, `STRESS_EQUATION_DUALITY_MESH.csv`, `C4_FUTURE25_DISTANCES.csv`, `C4_FUTURE25_ARGMAX.csv`, `MESH_COUNTS.csv`.
- Runtime summary: `research/delegated/results_support/MECHANISM_FIGURE_SUMMARY.json`. Incomplete input groups yield `NOT_RUN_INCOMPLETE_INPUT`; their CSV missing entries remain blank with recorded status. Plotting/read scope contains the listed diagnostic JSONs and public qualification/background metadata; no truth input or physical backend is imported.

Current JSON snapshot: 35 operator rows (7 cases × 5 methods), 12 stress fixtures, 9 physical worlds. C2/C3 are matched observations of the same three worlds. `results/g0_stress/COMPLETE.json` is now present, as are all three C4 all25 + independent-argmax records.

| Recorded quantity | Minimum | Maximum |
| --- | ---: | ---: |
| Full whitened error vs full | 0 | 0 |
| Woodbury whitened error | 2.1811469588200535e-13 | 2.7724517530658224e-13 |
| Schur-from-paid-ports whitened error | 4.883216157791773e-13 | 3.2664906630703492e-12 |
| Krylov whitened error | 7.884574490477737e-13 | 3.770503469458644e-12 |
| Goal ROM whitened error, both partitions | 4.3763480630568017e-13 | 6.517087414700959 |
| Full query solve CPU / s | 0.3097920000000016 | 0.3876279999999994 |
| Woodbury query solve CPU / s | 0.49776100000000056 | 0.7052240000000012 |
| Schur query solve CPU / s | 0.8082819999999984 | 1.0268120000000014 |
| Krylov query solve CPU / s | 0.2942429999999998 | 0.8963230000000024 |
| Goal ROM query solve CPU / s | 0.009515000000000384 | 0.013475999999997157 |
| Seven separately recorded setup kernels CPU / s | 0.02827799999999847 | 7.639944 |
| finite15 wall variance fraction | 1.4735475451786636e-9 | 0.999159081147629 |
| mean_joint current-data distance | 0.0003952231365232535 | 5.873621293590298 |
| Same-fixed-j0 task distance | 0.040099282823976755 | 6.704126798454728 |
| C4 stored max25 distance | 1.960877127764863 | 3.0681658287576212 |
| C4 independent distance at stored argmax | 1.9610258808204364 | 3.068786044164234 |
| Recorded JVP/adjoint relative error | 3.0083933604967155e-15 | 5.432973531044375e-15 |
| Stress white mesh response difference | 0.00694155478078847 | 0.009768751680817273 |

| Nuisance scope | Stored Schur eigenvalue range | Stored retention eigenvalue range | `UNRESOLVED` fixture count |
| --- | --- | --- | ---: |
| pose_only | [-2.3249776512791963e-17, 4863.8760620152625] | [3.5966642179449857e-33, 1] | 5 |
| pose_intrinsic | [0, 4502.399265627917] | [1.6600694384867478e-33, 0.9890841679460807] | 10 |
| pose_intrinsic_instrument | [0, 4247.606234876779] | [0, 0.9474518353139595] | 9 |

Undefined/missing: C1_0 original JVP/adjoint scalar and original corrector trial history were not preserved before serialization failure. Excluded A-support directions retain undefined retention; current stored arrays contain no fixture-wide `UNDEFINED_EMPTY_A_SUPPORT` status. Core singular values are recorded as `NOT_COMPUTED`; the core condition-product value is an inverse-norm estimate. Rendering and runtime verification of this new script remain **NOT_RUN**. Query timing is separate from shared setup timing. No stage decision is written here.
