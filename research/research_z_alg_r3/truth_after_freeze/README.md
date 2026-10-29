# Disclosed material targets for replay

These targets were disclosed only after the corresponding frozen predictions
were scored. They are for reproducibility and are not hidden or new blind
targets. DEV is exploratory; V6 holdouts became development evidence for the
separately frozen V7 revision; V7 used four fresh NEWCONFIRM12 targets.

Every NPZ contains exactly one array, `m`. Target generation seeds remain in
the ignored local `private_truth/` directory and are not included here. Public
optimization seeds in method configs are separate and remain public.

Replay from the R3 directory, writing a new output path:

```sh
python experiments/score_predictions.py --split confirm --freeze-manifest audit/CONFIRM_V7_FINAL_FREEZE.json --truth-root truth_after_freeze --disclosed-truth-replay --output /tmp/r3_v7_replay_scores.json
```

The replay flag permits changed file modification times after a Git clone and
checks the recorded prediction byte lengths. It explicitly labels the result
as disclosed-truth replay; it cannot establish fresh blindness. It does not
change frozen predictions or run a physics solver.
