# Reproduce cached finite-noise/prior audit

From the existing repository root and its existing compatible Python environment:

```
python research/research_max_struct_r11/experiments/diagnose_bothblocks16_noise_prior.py
```

Requires the already sealed BothBlocks16 predictions and authorized postseal GT package. It checks that every cached Y_clean/Y_prediction exists and raises rather than calling a new forward. The script reproduces all numerical JSON/CSV audit tables. The Chinese factual summary is editorial prose derived from those tables. EXTERNAL_PROCESS_COST.json preserves the measured wait4 cost of the completed audit; rerunning the command is not needed to inspect the saved results.

All80 original choices, successful/failed statuses and thresholds remain unchanged. Outputs are OFFLINE_DIAGNOSTIC and cannot be used retroactively for online candidate selection.
