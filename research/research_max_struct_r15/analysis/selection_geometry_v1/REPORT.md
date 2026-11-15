# Frozen-view selection geometry

Post-seal, already exposed same-cohort offline description. All264 task rows and24 cached coarse reference rows are retained. No candidate is changed or reselected.

Selected primary/raw/corrected/secondary comparison rows: 1056. Fixed-view inspection rows: 1320. Distinct eligible in-ball candidate pairs: 186.

GT_REFERENCE24 contains inverse-grid reference D=.5||Ygt_coarse−Y||²/||Y||², original prior R and Phi=D+R, and residual/delta membership. The reference is reused from old_endpoint_geometry; it is not a recomputed GT forward and does not convert crossgrid mismatch into a noise certificate.

SELECTED_VS_GT_REFERENCE uses the main analysis archived selections unchanged, including explicit secondary→primary fallback. Signed deltas are selected minus cached inverse-grid reference; missing values have reasons.

FEASIBLE_SAVED_VIEW_PAIRS enumerates every distinct candidate-ID pair in the five fixed sealed views only. It uses actual m/Y_prediction, original legality/eligibility and the original public ball radius. Duplicated view aliases/fallbacks do not create extra candidates. No GT is used to choose the maximum-distance pair.

Two distinct archived legal/eligible predictions inside the same finite discrete noise ball only witness finite discretized noise ambiguity at that radius. They do not establish continuous-Maxwell or GT nonidentifiability, or an online candidate selector.

Maximum material-distance and maximum data-distance pairs are OFFLINE_EXISTENCE_DIAGNOSTIC records, never online selectors. Absence of two stored in-ball views is only absence of this particular archived witness, not a uniqueness certificate.

Actual new physical LU=0; actual new physical complexRHS=0; no state/data_value_gradient/forward/SVD/optimizer calls. CPU and wall are recorded in RECEIPT.json. All scientific interpretation remains with the parent.
