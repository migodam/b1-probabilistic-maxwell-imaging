# Native Dual full-allowance budget control

This is a post-v1-score budget diagnostic requested by root, not a fresh preregistered confirmation. Allv1 endpoints remain immutable. Three new Dual arms use full inherited allowances: Fresnel500/500/500/1500 from actual R8 fresnel_invert.py full stages; Austria300perstage from paper_benchmark.py default --iterations300 and actual case mf5/8/10,beta8/2/.5,currentbeta8. Both use the same180proposal+60corrector CPU caps, realbox0..20, original starts/prior, unchanged DPCS and original physical-Phi selection. Austria literalDPCS15 remains the inherited benchmark contract; R9 noise-adaptive selection is a separate previously studied policy and was not silently substituted.

| Case | Train residual | Heldout residual / Austria NRMSE | Raw Austria NRMSE | CPU s |
|---|---:|---:|---:|---:|
|dielTM_dec8f|0.087643|0.091499|NA|84.80|
|twodielTM_8f|0.091377|0.098563|NA|84.82|
|austria_eps4_400MHz|0.093899|1.032829|0.6137737900303937|52.30|

AUSTRIA_OBJECTIVE_MATERIAL_DIAGNOSTIC.json decomposes each saved auxiliary/corrected Phi into its exact inherited prior and data term. A physical corrector can lower the registered Phi and increase material error even from a recognizable Dual auxiliary image. This implicates data/prior/selection adequacy in addition to search behavior; it does not prove a universal optimal replacement prior. There is no GT-based early stopping or alternate endpoint selection.

FRESNEL_CALIBRATION_FORWARD_EVIDENCE.json gives exact inherited source paths and numbers: empty-target validation medians.0251475723/.0240806883; equal-area VIE vs independent circular-interface reference errors at24² are.0243024323/.0189421309, dropping at36² to.0129375430/.0120232929. Those forward differences quantify discretization at nominal geometry, not measured-material accuracy. The initial schematic-axis interpretation failed (>100% measurement error) and is preserved; a later training-only pose diagnostic gave nominal measured train errors.1889438838/.3060280655. Optional known-reference gain transfer gave heldout.1607348202/.2961524435 and was not used by inversions. Source/receive horns and absolute field units remain unseparated; covariance unknown. Therefore similar common Fresnel heldout fit is predictive evidence, not physical material recovery.

All bad images and original v1 scores remain preserved. No delta, box, DPCS or calibration retuning occurred. Full-range material comparisons are in FULL_DUAL_CONTROL.png. New3 predictions were sealed before control scores, with explicit acknowledgement that earlier v1 scores were already known.

Metadata correction: actual inherited Austria forward grid is150² (inverse32²). See METADATA_AMENDMENT.json; original frozen metadata remains preserved as history.
