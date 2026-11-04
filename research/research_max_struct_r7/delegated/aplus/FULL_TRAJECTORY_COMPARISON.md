# A+ disclosed full-trajectory development

Status: COMPLETE; truth access NONE; recovery scoring NOT_RUN.

All runs use the same 3000 auxiliary iterations and original 3000-iteration corrector allowance. Ratios compare to the same-object original-A control. Original Phi is the online observation/prior objective, not a GT recovery metric.

| Object | Readout | CPU s | Aux CPU s | Phi | Phi / A | CPU / A | LU | RHS | inner loss increases |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| obs_7b8d5df517f1497d | original_a | 43.298 | 31.566 | 5.8843846e-06 | 1.000000 | 1.000 | 1394 | 33408 | 0 |
| obs_b20d5a55677d4eca | original_a | 72.748 | 53.468 | 6.4736232e-06 | 1.000000 | 1.000 | 1410 | 33792 | 0 |
| obs_40ac5b79e8a74ac3 | original_a | 109.314 | 69.500 | 0.0038128817 | 1.000000 | 1.000 | 2918 | 69984 | 0 |
| obs_00f185ff5525406b | original_a | 58.850 | 39.188 | 3.6897772e-06 | 1.000000 | 1.000 | 1245 | 29832 | 0 |
| obs_7b8d5df517f1497d | mobius | 79.494 | 50.726 | 5.8842795e-06 | 0.999982 | 1.836 | 1435 | 34392 | 0 |
| obs_b20d5a55677d4eca | mobius | 82.589 | 53.600 | 6.4716745e-06 | 0.999699 | 1.135 | 1514 | 36288 | 0 |
| obs_40ac5b79e8a74ac3 | mobius | 108.462 | 53.418 | 0.0038128816 | 1.000000 | 0.992 | 2857 | 68520 | 0 |
| obs_00f185ff5525406b | mobius | 76.891 | 52.461 | 3.6937243e-06 | 1.001070 | 1.307 | 1215 | 29112 | 0 |
| obs_7b8d5df517f1497d | projected | 76.383 | 50.421 | 5.8850418e-06 | 1.000112 | 1.764 | 1384 | 33168 | 0 |
| obs_b20d5a55677d4eca | projected | 86.983 | 47.674 | 0.0058324021 | 900.948649 | 1.196 | 2160 | 51792 | 0 |
| obs_40ac5b79e8a74ac3 | projected | 96.135 | 45.107 | 0.0038128826 | 1.000000 | 0.879 | 2906 | 69696 | 0 |
| obs_00f185ff5525406b | projected | 63.525 | 45.469 | 3.6985799e-06 | 1.002386 | 1.079 | 1061 | 25416 | 0 |
| obs_7b8d5df517f1497d | mobius_guarded | 69.907 | 47.223 | 5.8845369e-06 | 1.000026 | 1.615 | 1360 | 32592 | 0 |
| obs_b20d5a55677d4eca | mobius_guarded | 56.082 | 39.866 | 6.4716745e-06 | 0.999699 | 0.771 | 1514 | 36288 | 0 |
| obs_40ac5b79e8a74ac3 | mobius_guarded | 82.069 | 41.279 | 0.0038128812 | 1.000000 | 0.751 | 2905 | 69672 | 0 |
| obs_00f185ff5525406b | mobius_guarded | 59.962 | 41.024 | 3.6917134e-06 | 1.000525 | 1.019 | 1235 | 29592 | 0 |
| obs_7b8d5df517f1497d | mobius_arc_candidates | 45.361 | 26.054 | 5.8894223e-06 | 1.000856 | 1.048 | 1218 | 29184 | 0 |
| obs_b20d5a55677d4eca | mobius_arc_candidates | 46.525 | 26.490 | 6.4738536e-06 | 1.000036 | 0.640 | 1402 | 33600 | 0 |
| obs_40ac5b79e8a74ac3 | mobius_arc_candidates | 62.756 | 25.950 | 0.003812882 | 1.000000 | 0.574 | 2767 | 66360 | 0 |
| obs_00f185ff5525406b | mobius_arc_candidates | 43.140 | 26.450 | 3.6890965e-06 | 0.999816 | 0.733 | 1322 | 31680 | 0 |

Each output directory contains full auxiliary scalar histories, stage joint m/J/Q/A, paid physical stage states and selected/final predictions, corrector evaluation/accepted histories, and process CPU/LU/RHS receipts. Source and configuration are frozen separately.
