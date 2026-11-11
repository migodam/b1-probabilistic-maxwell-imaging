# Final Austria physics diagnostics: parent handoff

Source correction: the inherited R8 original Austria observation is **150²→32²**, not64²→32². `external/METADATA_AMENDMENT.json` and `external_continued/METADATA_AMENDMENT.json` preserve the original frozen metadata mistake and point to actual PUBLIC_DATA/DATA_RECEIPT values. Reports and Austria figure labels were corrected; original scientific predictions/receipts were not altered.

Exactly one separately authorized offlineGT32² forward used the existing operator:1LU,20sourceRHS, residual qualification1.80e-15. The actual coarseGT has noisy-data residual **.2252565**, outside the unchanged.12 ball. Existing inherited grid-adequacy evidence concerns epsilon3.5, not this epsilon4 case; it was not substituted. No clean fine-grid observation was saved, so no fine solve or reconstructed noise was introduced. `external_physics/RESULTS.json` and REPORT.md retain this evidence and cost.

Saved original DR shadows at50/500/1000/5000 have relative A–C gaps.127820/.025018/.019980/.013352. State norm/incident stays approximately1.17–1.23 on these sampled times. The A shadow satisfies data residual.12 and propagation; the C shadow satisfies shared-real constitutive closure but its propagation defect remains.178996 incident norms at5000 and data residual.227310. The last audited physical material residual is.180285, while the selected raw audit is better and the later shared physical corrector ends.122674. Thus the actual bounded run did not reach a common intersection; affine data fit alone is not successful physical reconstruction. No global impossibility statement follows.

Root conditionally authorized a matched32 control because coarseGT failed the original ball. Matched generation reuses the already paid clean32 prediction above, with fresh exact10%fixednorm noise seed2026111232; no redundant generation forward. GeneratorGT and public observations are separated, all3 predictions sealed before scoring. True matched GT residual **.0995131** lies inside.12. Same contrast targetchi3, box0..20, zero start, prior, scale6, unchanged.12 rule, full native Dual allowances and180aux/60corrector budgets; no retuning.

| Matched32 method | Selected auxiliary NRMSE | Corrected NRMSE | Final residual |
|---|---:|---:|---:|
| direct |.970767|.970767|.566303|
| full native Dual |.606380|1.081429|.093104|
| DR5000 |2.190625|2.253906|.134903|

Direct stopped by projected-gradient criterion; Dual and DR common correctors reached their1000iteration caps. Full Dual used668actualPR steps under inherited coefficient stopping. DR completed5000 steps. Both Dual/DR reachedchi20 spikes. Bad images and all saved fields remain retained. Full-range MATCHED32_MATERIALS_FULL_RANGE.png was visually inspected; its0–4 detail companion explicitly labels its saturation range.

These observations isolate two issues: the original150→32 target/model mismatch makes the.12 ball reject the coarseGT, while the matched control still shows bounded search/selection failure despite an explicitly feasible true target. The matched Dual physical corrector worsens material while improving fit, so data/prior/selection interpretation remains relevant even without discretization mismatch. Scientific claims, novelty and final verdict remain parent Codex decisions; this is OFFLINE_PHYSICS_DIAGNOSTIC after exposed results, not original-paper reproduction or fresh confirmation.

Matched solver external-child CPU146.655343s; fresh noise generation.001657s; scorer/figure CPU.335322s. The sole offlineGT forward is charged once in external_physics (about.537s plus setup/output clock), not again as generation. Matched queue is COMPLETE; no jobs remain active. Primary artifacts: REPORT.md, SCORES.json, COST_LEDGER.json, CONTRACT.json, GENERATION_RECEIPT.json, PREDICTION_SEAL.json, QUEUE.json, RESULTS/, and the two material images.
