# Reproduce R5

Scientific source/data results are append-only. Never run the historical R2–R4 drivers in their completed output directories.

## Environment actually used

Apple Silicon macOS, Python3.13.13, NumPy2.5.3, SciPy1.18.1, mpmath1.3.0; matplotlib version is saved in final verification. CPU only, BLAS threads1. No CUDA, service, Julia, PHCpack, neural network or orchestration package was installed. Existing interpreter:

    /Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python

Required standard dependencies are numpy, scipy, mpmath, matplotlib. The historical R2 src package in the sibling research_z_alg_r2 directory is a read-only dependency. Its files are hashed in inputs/solver_freeze.json. This R5 folder can be moved with its sibling R2 dependency.

## Commands actually executed, in order

From repository root, with the interpreter above (shell-quote its path):

    python research/research_z_alg_r5_measurement/src/debug.py
    python research/research_z_alg_r5_measurement/src/test_formulae.py
    python research/research_z_alg_r5_measurement/src/freeze.py
    python research/research_z_alg_r5_measurement/src/generate.py
    python research/research_z_alg_r5_measurement/src/run_gate1.py
    python research/research_z_alg_r5_measurement/src/qualify_predictions.py
    python research/research_z_alg_r5_measurement/src/freeze_predictions.py
    python research/research_z_alg_r5_measurement/src/score_v2.py
    python research/research_z_alg_r5_measurement/src/account_and_plot.py

Two engineering failures are preserved: an absent complex-interval report helper during postprocessing; and score.py treating interval real/imag dictionaries as scalar complex numbers. Final qualify_predictions.py and score_v2.py are the corrected entry points. They did not change targets, discovery, selections, retentions or success criteria. The unsafe-to-interpret preliminary qualification output remains marked as before_serialization_audit; final qualification.json checks a safe inward binary64 root box.

## Reproduce the same now-exposed six targets

    python research/research_z_alg_r5_measurement/src/replay.py --output-name research_z_alg_r5_replay_20261118

The new sibling must not exist. It copies stored public observations and private response/GT data, runs the same frozen algorithms, then scores. It is an EXPOSED_REPLAY, never a fresh blind replication. Source choices and floating point estimates should agree up to numerical library/platform differences; timings need not. No target regeneration is performed. The replay entry was syntax/integration audited; the entire scientific batch was not repeated solely to verify a wrapper.

For a new blind batch, do not silently invoke generate.py again and pool it with R5. This contract caps the current study at six targets. A new cohort requires a separately identified study, justified by new evidence; current scoped verdict is CONDITIONAL GO; see the final report and separately labelled Codex amendment.

## Data layout and interpretation

- inputs/solver_freeze.json: discovery source and protocol hashes before target generation.
- inputs/prediction_freeze.json: predictions, measurement requests/observations, retention results and original score source before GT access.
- inputs/scoring_repair_receipt.json: exact decoder-only repair, original preserved.
- results/T*/clean and noisy: all Q endpoints, filters,64 TRF/reciprocal endpoints, certificate attempts, source menu, prediction freeze, requests, observed values, candidate/tube decisions, same-source TRF and source-completion baseline.
- private/: GT, generator seed and simulated measurement bank. They are now deliverable for reproducibility; being published/exposed does not retroactively change the original process isolation record.
- candidate_registry.csv: every accepted material coordinate and post-freeze GT error. Rejected algebraic endpoints remain in Q.json/filter.
- qualification.json: local root certificates, safe box serialization audit, local source-response intervals, global UNKNOWN labels.
- cost_accounting.json, FINAL_COST.json: full/reduced LU/RHS, interval inversions and CPU; do not sum overlapping nested event CPU.

Source noise is a fixed amplitude1e-4 phase perturbation per acquisition. Identical source weights share identical simulated new-noise realizations across methods. Each strategy is a counterfactual acquisition policy over the same synthetic object; totals of all policies are not a single deployment acquisition cost. Deployment comparisons are1 added measurement versus2/3 for completion, each unit calibrated radiation energy and unit duration.

## Verification actually performed

Formula and finite difference tests; independent historical full Maxwell forward; node/pole negative controls; Krawczyk tests; source Gram angular quadrature; source/solver freeze hash checks; no historical R2–R4 modifications; full six-target experiment and independent post-freeze scoring; final data/accounting audit; visual inspection of all five scientific figures. Original local-only coverage is superseded by the Codex central-data certificates. Noise-feasible coverage remains UNKNOWN and mature homotopy remains NOT_RUN.


## Codex amendment replay

The amendment is exposed development on the original data, plus one declared H12 applicability probe. It is not a second blind cohort. All certificate/physics results are at amendment_codex/results. Use a new sibling directory:

    python research/research_z_alg_r5_measurement/amendment_codex/src/replay_amendment.py --output-name research_z_alg_r5_amendment_replay

This copies the frozen inputs and baseline inference needed by the development designs, then runs certification, noise diagnostics, constrained-source acquisition/correction, H12 and presentation. Copied original results must not be counted as newly executed. The wrapper was syntax/path audited, not batch re-executed solely for wrapper validation. Individual actual runs are listed in its steps array. Old prototype v0 and the pre-centered noise attempt are retained in the delivery; the replay executes corrected final versions only.

Arithmetic: mpmath.iv.prec=128 for the amendment. No 512-bit or unbounded precision rescue. Core scripts use append-only save and intentionally fail if an output exists. Launch replays in fresh folders. The standalone archive includes the minimal read-only sibling research_z_alg_r2 dependency and its provenance.

The first report and handoff are preserved as INITIAL_* historical artifacts. Current authoritative interpretation is RESEARCH_REPORT_ZH.md, CODEX_TO_PRO_HANDOFF.md, amendment_codex/THEORY_CERTIFICATES.md and amendment_codex/results/CLAIM_STATUS.json. The original EXPERIMENT_PROTOCOL.md is deliberately unchanged; the user amendment and development scope are recorded separately.
