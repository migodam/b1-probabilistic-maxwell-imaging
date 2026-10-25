# Pre-inference review repairs

All changes below were completed before generating the finite-prior observations or inspecting their posterior results. Gates, physical worlds and noise are unchanged.

1. Task-weighted triggering now uses the frozen current-reference task-risk threshold. The reference computation is explicitly attributed; this is a reference-dependent diagnostic and cannot claim cheap certification or acceleration.
2. Bank completion, dimensions, calibrated gain, frequency, acquisition, prior positivity/normalization and 72 unique observation IDs are checked before inference.
3. The ordinary iterative point control now performs two wall-then-pose coordinate MAP sweeps with all historical factors. It can revisit wall states. It is finite-prior block MAP, not an implementation of general continuous iterative layer stripping.
4. Partial posterior saves contain a complete mask and terminal boolean; downstream evaluation must reject incomplete arrays. Ordinary exact incremental updates retain log weights without flooring vanished tails.

The source review is `research/delegated/package_audit/E2_CODE_REVIEW.md`. Parent Codex retained responsibility for formulas and interpretation.
