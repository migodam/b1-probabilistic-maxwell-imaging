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

## Final1024 cohort and one frozen model

The confirmation source is `audit/CONFIRMATION_PLAN_FROZEN.json`. Its primary arms are seed7; seeds17/27 are confirmations, and the two seed7 cascades are separate ablations. All checkpoints are retained, rather than choosing one using TEST truth. In a **fresh reproduction copy**, after all2048TRAIN anchors and128VAL anchors are available:

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
$PYTHON experiments/train_models.py --data-ready --output models/reproduce1024 --seeds 7 17 27 --device auto
$PYTHON experiments/train_models.py --data-ready --output models/reproduce1024 --arms v_physics j_null --seeds 7 --cascade-passes 2 --device auto
$PYTHON experiments/predict_split.py --checkpoint models/reproduce1024/v_physics_seed7/best.pt --run reproduce_v_physics_seed7 --split test --device auto
$PYTHON experiments/run_evaluation.py evaluate --runs reproduce_v_physics_seed7 --controls continue_trf weakreg45_anchor --splits test
$PYTHON experiments/run_evaluation.py seal --runs reproduce_v_physics_seed7 --controls continue_trf weakreg45_anchor --splits test --output evaluation/REPRODUCED_SEALED.json
$PYTHON experiments/run_evaluation.py score --seal evaluation/REPRODUCED_SEALED.json --output evaluation/REPRODUCED_SCORES.json
```

Do not overwrite the archived confirmation names when rerunning. Source fonts can differ between platforms; the saved material/current arrays and public observations are the authoritative frozen dataset. CPU, MPS and CUDA training may not be bit-identical. The original runs used MPS with real-pair complex64 operations; original final physics used float64/complex128.

## A new observation, with no GT

For the same frozen acquisition, provide an NPZ with only complex `Y` of shape48×12 and scalar `noise` (0 or0.001). The entry point composes the same counted GN, trained proposal and full-VIE correction modules. It returns candidates and online top1, not a GT error or a uniqueness certificate:

```sh
$PYTHON experiments/infer_observation.py --observation new_observation.npz --operators data/operators.npz --checkpoint models/confirm1024/v_physics_seed7/best.pt --output new_observation_result --device auto
```

## Redraw / re-score downloaded frozen artifacts

Extract all release parts into one directory, preserving their `research/` paths, then enter this study directory. The frozen score retains the original location, so historical replay explicitly remaps paths. This is a re-score of already disclosed data, not a new blind test:

```sh
$PYTHON experiments/run_evaluation.py score --study "$PWD" --replay-frozen --seal evaluation/CONFIRMATION_SEALED.json --output evaluation/HISTORICAL_REPLAY_SCORES.json
$PYTHON experiments/run_evaluation.py plots --study "$PWD" --replay-frozen --output evaluation/HISTORICAL_REPLAY_SCORES.json
$PYTHON src/summarize_results.py --scores evaluation/HISTORICAL_REPLAY_SCORES.json --output results/historical_replay
```

No formula/error-floor reruns, retraining or new GN solving are needed to inspect the published conclusions. The remote helper is optional and specific to the owner's existing SSH setup; reproduction does not require a new daemon, CUDA installation or scheduler.
