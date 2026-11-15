# Leading-real-mode finite32 observations

Declared artificial8x8root only; no inverse/test-GT inputs. One fixed eigenvector/epsilon, no search.

SavedIDEAL01leadingeigenvalue=1.002712449117704; real/imaginarytol1e-10, actual imaginarypart0. Savedrho^32=1.0905485931747956. All8 first-step native/reset array comparisons are bitwise equal; this was recorded after all trajectories and was not an admission gate.

| pair | sign | memory policy | actual material gain at32 | sign-aligned mode coefficient at32 |
|---|---:|---|---:|---:|
| 00 | 1 | native_PR | 0.685240892225 | 0.62525200753 |
| 00 | 1 | reset_before_every_update | 0.622358970361 | 0.569236151883 |
| 00 | -1 | native_PR | 0.685380898395 | 0.625432864627 |
| 00 | -1 | reset_before_every_update | 0.622509905679 | 0.569417308059 |
| 01 | 1 | native_PR | 0.691965439916 | 0.623356432953 |
| 01 | 1 | reset_before_every_update | 0.633511387308 | 0.574444029766 |
| 01 | -1 | native_PR | 0.69214457244 | 0.62357265353 |
| 01 | -1 | reset_before_every_update | 0.63368251609 | 0.574635852484 |
| 10 | 1 | native_PR | 0.473169842407 | 0.449932215117 |
| 10 | 1 | reset_before_every_update | 0.397007574632 | 0.378801522328 |
| 10 | -1 | native_PR | 0.473242167868 | 0.450011960619 |
| 10 | -1 | reset_before_every_update | 0.397100722905 | 0.378896779964 |
| 11 | 1 | native_PR | 0.482624250127 | 0.459583319152 |
| 11 | 1 | reset_before_every_update | 0.404568806793 | 0.38531340363 |
| 11 | -1 | native_PR | 0.482753372459 | 0.459725357391 |
| 11 | -1 | reset_before_every_update | 0.404659961713 | 0.385405856512 |

All16 trajectories/512updates completed;48fresh64-dimensionalLU and288complexsourceRHS. TimedexperimentCPU=1.520583s, wholeprocessCPU=1.554212s. Additional data-only plot fee is PLOT_RECEIPT.json; prior read-only preflight CPU was not separatelyinstrumented and is explicitly not zero-imputed.

These declaredfinite-scale observations differ from repeating the exact-current/exact-material ideal01map. They establish no FrechetJacobian, infiniteLyapunov exponent, globalstability/recovery, or noveltheorem. Original PR/readout guard diagnostics, sourceweights, fullstates and all paidphysics residuals remain in the arrays/tables. Parent owns interpretation.
