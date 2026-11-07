# Generalization statistics plotting handoff

Status: **IMPLEMENTED / SCIENTIFIC FIGURES NOT_RUN**. No new TEST score input was opened or synthesized during implementation. No truth array, prediction array, fit or experiment is used by this script. Only the two assigned files were added.

After final TEST methods/predictions are frozen and the R4 scorer has written its actual JSON/CSV:

```sh
python research/research_gn_nn_preservation_r4/experiments/plot_generalization.py
```

Optional model-selection diagnostic, explicitly VAL only:

```sh
python research/research_gn_nn_preservation_r4/experiments/plot_generalization.py --val-scores research/research_gn_nn_preservation_r4/results/VAL_SCORES.json
```

`--scores` accepts final `results/TEST_SCORES.json`, or its companion `.csv` with the R4 JSON provenance present. The script requires the correct R4 schema, TEST scope and TEST rows; VAL selection input independently requires VAL scope. With no final TEST scores it writes only `figures/generalization/STATUS.json` with `NOT_RUN_NO_FINAL_TEST_SCORES`; it does not draw fictional figures.

## Statistical contract

- All actual score-plan rows remain in the exported observations, coverage table and success denominator, including FAIL and NOT_RUN. The script does not manufacture unplanned CIE/baseline runs or 108-parent coverage. Every run/stratum displays its actual planned mother count and NOT_RUN observation count.
- A missing top1 uses its **already-scored common anchor NRMSE** only for a separately labelled material diagnostic. A score from another run for the same observation may supply that identical common-anchor value. This never changes FAIL into success. If no scored anchor is available for an observation, it remains unresolved; the full-stratum diagnostic mean/interval is withheld instead of computing a selected-only mean. Scalar missing-data coverage is exported.
- Clean/noisy observations are averaged **inside each mother first**. Bootstrap samples mothers, fixed seed `20261107`, 500 repetitions. Paired changes compare the diagnostic result against the same mother's common anchor. Rate intervals also use mother resampling. These intervals do not treat clean/noisy or three model seeds as extra independent materials.
- `generalization_seed7_primary.png` separates new-instance, unseen-class and unseen-style populations, showing paired mother NRMSE changes, physical-data success, quality joint success, and joint success within the actual total55 budget. NRMSE panels are explicitly top1/anchor-fill diagnostics. `amplitude_seed7_primary.png` repeats this for A=.1/.2/.4. Unseeded traditional runs are included only to their actual score-plan extent.
- Seed17/27 have a separate confirmation figure. `seed_matched_mothers.png` connects the same mother only when that run family has all three seed records. Complete triplets and actual per-seed planned coverage are labelled. Separate `seed_confirmation_paired.csv` stores per-mother confirmation-minus-seed7 changes; incomplete families are not silently promoted to full confirmation.
- Online cost points use only measured `online_total_wall_seconds`, `online_total_LU`, `online_total_RHS` from the scorer. They show mean measured cost with measured-parent coverage beside each point; missing costs are never filled by the cap or zero. The quality and rate denominators remain the full actual plan even when cost measurement is incomplete. These x-coordinates are therefore explicitly conditional on measured-cost coverage, not an inferred complete-plan cost.
- `one_time_costs.csv` separately invoices actual training and setup metadata when present; it does not amortize them into per-case points or sum ambiguous reused setup scopes. `training_validation_curves.png` reads only TRAIN/VAL history and plots the declared material checkpoint-selection objective. Optional network/loss/alpha comparisons use VAL scores only and carry **VAL MODEL SELECTION ONLY** in filename and title.

## Outputs and limitations

The future output folder includes observation/parent diagnostics; group, amplitude and cost summaries; actual coverage; matched seed pairs; one-time cost metadata; PNG figures; and a provenance/status receipt. No arrays are loaded and no new material error is calculated. The script uses NumPy and Matplotlib already used by this study; no dependency install is required.

Validation at delivery: Python syntax parsed without executing the script or importing experiment modules. Scientific figure rendering and visual inspection are **NOT_RUN** until real final scores exist. No claim about training completion, unseen-class success, CIE completeness or scientific acceptance is made by this handoff.
