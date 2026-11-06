# Reproduction

Use an existing NumPy/SciPy/Pillow/Torch environment. Keep BLAS single-threaded. The frozen B1 physics sources are dependencies; FAEMIS linear-scattering measurements are not interchangeable with this dataset.

From this study directory, in a fresh reproduction copy:

```sh
$PYTHON experiments/prepare_data.py
$PYTHON experiments/run_anchors.py --split train --workers 2
$PYTHON experiments/run_anchors.py --split val --split test --workers 2
$PYTHON experiments/run_anchors.py --split val --split test --workers 2 --weakreg 1e-5
$PYTHON experiments/train_models.py --data-ready --output models/pilot256 --seeds 7 --device mps
```

The generator preserves existing data. Optional training expansion uses `prepare_data.py --expand-train 1024` only after the declared validation gate, without replacing the initial cohort or holdout.

The existing Windows transport supports four private P-core workers, without installing a service or changing other sessions:

```sh
python3 src/remote_jobs.py upload
python3 src/remote_jobs.py run --args '--split val --split test --workers 4 --core-mask 85'
python3 src/remote_jobs.py fetch
```

Uploads contain only public observations, operators and source. The transport reuses the local `Gaussian/A9/code/remote.py`; other environments must provide their own SSH configuration and may set `GN_NN_REMOTE_ROOT`.

`src/predict.py --help` describes the checkpoint/public-anchor inference interface. Predictions never read truth. For a completed run:

```sh
$PYTHON experiments/run_evaluation.py evaluate --runs RUN --controls continue_trf weakreg_trf --splits val
$PYTHON experiments/run_evaluation.py seal --runs RUN --controls continue_trf weakreg_trf --splits val --output evaluation/VAL_SEALED.json
$PYTHON experiments/run_evaluation.py score --seal evaluation/VAL_SEALED.json --output evaluation/VAL_SCORES.json --oracle
$PYTHON experiments/run_evaluation.py plots --output evaluation/VAL_SCORES.json
```

Final testing replaces `val` with `test` only after freezing model choices and predictions. Oracle current readout is an offline representation diagnostic, never an online candidate. Training/physics reruns are unnecessary when only redrawing already frozen results; use saved predictions, scores and receipts.
