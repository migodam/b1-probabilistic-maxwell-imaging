# External anchors: observed readout for parent scientific review

All12 declared runs completed with no failed attempts, extra radius, rerun, or result-dependent replacement. All12 predictions were sealed before first heldout/pixel-GT score. `QUEUE.json` records all children; no queue remains active. Actual external solver-child CPU totals402.666742s including imports/output; per-arm inner clocks sum399.330165s, an overlapping lower bound. Postseal score/report CPU1.247248s; no new physical solves in scoring. Existing Austria observation-generation receipt is separately retained as reused provenance, not new execution.

| Anchor | direct | Dual | AP | DR |
|---|---:|---:|---:|---:|
| Fresnel single heldout predictive residual | .09150 | .09150 | .09150 | .09152 |
| Fresnel twin heldout predictive residual | .09857 | .09856 | .09868 | .09855 |
| Austria final material NRMSE | .98656 | 1.41964 | .98941 | 2.13198 |
| Austria training data residual | .54568 | .10726 | .54567 | .12267 |

The measured predictive fits are similar after the shared physical corrector. They do not establish correct absolute material. Final fitted Fresnel peaks are chi≈9–13 while noise covariance, receive-horn factors and absolute field units are unresolved. The single/twin measured cases have no pixelGT. Raw radius.15||Ytrain|| was fixed before inversion as MODEL_DISCREPANCY_DIAGNOSTIC, not an iid-noise certificate. No heldout gain, radius or illumination was refitted.

The Austria32²/150²-generation control retains a useful failure: AP has nearly the same poor training basin as direct; DR improves physical fit but its material is substantially worse. Frozen Dual has the best training fit, but its NRMSE grows from.72891 at the auxiliary selected stage to1.41964 after the identical physical corrector; DR grows2.03499→2.13198. Dual and DR hit chi20 spikes. No earlier material was selected using these scores. This is evidence that objective decrease and material recovery differ, with final scientific/novelty interpretation reserved to parent Codex.

Bounds0..20, exact chi=epsilon_r−1, source masks, empty-target multipole calibration, inherited prior, auxiliary work and common corrector are fully stated in CONTRACT.json/REPORT.md. Austria epsilon4 means targetchi3, notchi4. It is a reduced32² adaptation, not a full author reproduction. OriginalR8 baseline algorithms and operators were reused; no replacement backend/training was introduced. Projection groupedv2, scale6, max5000/audit50 were fixed for this cohort.

Keep both image ranges: `EXTERNAL_MATERIALS.png` and `AUSTRIA_MATERIALS.png` use display0–4 to reveal structure and thus saturate high spikes; `EXTERNAL_MATERIALS_FULL_RANGE.png` displays the complete0–20 inversion range and retains all raw/corrected failures. Both the full-range image and Austria comparison were visually inspected; labels and grids are readable. `FRESNEL_HELDOUT_PREDICTIONS.png` plots predictions without fitted gain. Full derived arrays preserve all amplitudes and all auxiliary/corrected endpoints; AP/DR complete audit states are also retained.

Primary handoff files: REPORT.md, SOURCE_INVENTORY.json, CONTRACT.json, QUEUE.json, PREDICTION_SEAL.json, SCORES.json, COST_LEDGER.json, run_external.py, score_report.py, RESULTS/{case}/{method}/{derived.npz,RECEIPT.json}, and projection histories/states. Original .exp files remain only at their original sibling R8 paths; none were copied or redistributed.

Metadata correction: actual inherited Austria forward grid is150² (inverse32²). See METADATA_AMENDMENT.json; original frozen metadata remains preserved as history.
